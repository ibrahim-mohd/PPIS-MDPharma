# Written by Mohd Ibrahim
# Technical University of Munich

import os
import re
import json
import time
import math
import itertools
import logging
import tempfile
import argparse
import subprocess
import numpy as np
import networkx as nx
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import warnings
warnings.filterwarnings("ignore")

# --------------------- COLORS ---------------------
NODE_COLORS = {
    "Aromatic": "orange",
    "HydrogenDonor": "magenta",
    "HydrogenAcceptor": "green",
    "Hydrophobic": "cyan",
    "PositiveIon": "red",
    "ExclusionSphere": "gray",
    "NegativeIon": "blue"
}


# --------------------- LOGGING SETUP ---------------------
def setup_logging(output_dir, verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(Path(output_dir) / "generate_and_screen.log")
        ],
    )


# --------------------- ARGUMENT PARSING ---------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate subpharmacophores and immediately screen them with Pharmer."
    )

    # Generation arguments
    parser.add_argument('-j', dest='pharmacophore_file', type=str, required=True,
                        help='Master pharmacophore JSON file')
    parser.add_argument('-min_node', dest='min_node', type=int, default=6,
                        help='Minimum number of nodes to consider')
    parser.add_argument('-max_node', dest='max_node', type=int, default=None,
                        help='Maximum number of nodes to consider')
    parser.add_argument('-top', dest='top_percentage', type=float, default=100,
                        help='Percentage of top-ranked models to keep')
    parser.add_argument('-ntop_limit', dest='ntop_limit', type=int, default=50000,
                        help='If more than this many models exist, use the top ntop_limit pharmacophore models.')
    parser.add_argument('-o', dest='output_dir', type=Path, required=True,
                        help='Output directory for results and logs')

    # Pharmer screening arguments
    parser.add_argument("--path_exe", "-p_exe", type=Path, default="/usr/local/bin/pharmer.static",
                        help="Full path to Pharmer executable")
    parser.add_argument("--database-file", "-df", type=Path,
                        help="Text file containing database paths (one per line)")
    parser.add_argument("--database", "-d", action="append", type=Path,
                        help="Database path (can be used multiple times)")
    parser.add_argument("--workers", "-np", type=int, default=4,
                        help="Number of parallel workers")
    parser.add_argument("--extension", "-e", type=str, choices=["sdf", "mol2"],
                        default="sdf", help="Output file extension")
    parser.add_argument("--max_hits", "-max", type=int, default=10000,
                        help="Stop after total hits exceed this number")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose output")

    return parser.parse_args()


# --------------------- DATABASE HANDLING ---------------------
def load_database_paths(args):
    paths = []
    if args.database:
        paths.extend(args.database)
    if args.database_file and args.database_file.exists():
        with open(args.database_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    paths.append(Path(line))
    # unique and valid
    seen, valid = set(), []
    for p in paths:
        if p.exists() and p not in seen:
            seen.add(p)
            valid.append(p)
        else:
            logging.warning(f"Skipping invalid db path: {p}")
    return valid


# --------------------- GENERATION FUNCTIONS ---------------------
def get_distance(p1, p2): return np.linalg.norm(np.array(p1) - np.array(p2))


def build_graph_from_pharmacophore(pharma_data):
    G = nx.Graph()
    exclusion_spheres = []
    index = 0
    for site in pharma_data["points"]:
        if site["name"] != "ExclusionSphere":
            index += 1
            G.add_node(
                index,
                name=site["name"],
                position=np.array([site["x"], site["y"], site["z"]]),
                radius=site["radius"],
                size=site["size"],
                requirement=site["requirement"],
                enabled=site["enabled"],
                vector=site["vector"],
                color=NODE_COLORS[site["name"]],
                score=site["score"],
                chainID=site["chainID"],
            )
        else:
            exclusion_spheres.append(site)

    for i, j in itertools.combinations(G.nodes, 2):
        distance = get_distance(G.nodes[i]['position'], G.nodes[j]['position'])
        score_sum = G.nodes[i]['score']['normed_score'] + G.nodes[j]['score']['normed_score']
        G.add_edge(i, j, distance=distance, score=score_sum)
    return G, exclusion_spheres


def generate_graphs(G, r, top_percentage=0, ntop_limit=0):
    
    graphs = []
    
    for nodes in itertools.combinations(G.nodes, r):
        newG = G.subgraph(nodes).copy()
        mst_dist = sum(d["distance"] for _, _, d in nx.minimum_spanning_tree(newG).edges(data=True))
        count_A = sum(1 for _, a in newG.nodes(data=True) if a.get("chainID") == "A")
        count_B = sum(1 for _, a in newG.nodes(data=True) if a.get("chainID") == "B")
        sym = min(count_A, count_B) / max(count_A, count_B) if max(count_A, count_B) > 0 else 0
        total_score = sum(a["score"]["normed_score"] for _, a in newG.nodes(data=True))
        graphs.append([newG, dict(min_spanning_tree_dist=mst_dist, score=total_score,
                                  symmetry=sym, combined_score=sym * total_score)])
    sorted_g = sorted(graphs, key=lambda x: x[1]["combined_score"], reverse=True)
    if top_percentage:
        N = min(ntop_limit, int(np.ceil(len(sorted_g) * top_percentage / 100)))
        return sorted_g[:N]
    return sorted_g


def graph_to_pharmacophore_json(graph, info, exclusion_spheres, output_path, filename):
    data = {"extra_info": info, "points": []}
    for _, attr in graph.nodes(data=True):
        node = {}
        for k, v in attr.items():
            if k == "position":
                node.update({"x": float(v[0]), "y": float(v[1]), "z": float(v[2])})
            elif k not in ["color", "value"]:
                node[k] = v
        data["points"].append(node)
    data["points"].extend(exclusion_spheres)
    os.makedirs(output_path, exist_ok=True)
    with open(Path(output_path) / filename, "w") as f:
        json.dump(data, f, indent=4)


# --------------------- PHARMER FUNCTIONS ---------------------
def run_pharmer(input_file, output_file, db_path, exe_path):
    cmd = [str(exe_path), "dbsearch", f"-dbdir={db_path}", f"-in={input_file}", f"-out={output_file}"]
    nhits, pharmer_time, success = 0, 0.0, False
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        hits_match = re.search(r"NumResults:\s*(\d+)", result.stdout)
        time_match = re.search(r"Time:\s*([\d.]+)", result.stdout)
        if hits_match:
            nhits = int(hits_match.group(1))
        if time_match:
            pharmer_time = float(time_match.group(1))
        success = True
    except subprocess.CalledProcessError as e:
        logging.error(f"Pharmer failed: {e}")
        if e.stderr:
            logging.error(e.stderr)
    return nhits, pharmer_time, success


# --------------------- MAIN WORKFLOW ---------------------
def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.output_dir, args.verbose)

    logging.info("Starting generation + screening workflow")

    db_paths = load_database_paths(args)
    if not db_paths:
        logging.error("No valid database paths provided.")
        return 1

    with open(args.pharmacophore_file) as f:
        pharma_data = json.load(f)

    G, exclusion_spheres = build_graph_from_pharmacophore(pharma_data)
    max_node = args.max_node if args.max_node else len(G.nodes)

    summary = {}
    total_hits = 0
    start_time = time.perf_counter()

    # Process from largest graphs down to smallest
    for r in range(max_node, args.min_node - 1, -1):
        
        start_time_graph = time.perf_counter()
        
        subgraphs = generate_graphs(G, r, args.top_percentage, args.ntop_limit)
        logging.info(f"Processing {len(subgraphs)} ranked subgraphs of size {r}")

        n_hits = 0
        total_graphs = len(subgraphs)
        successful_runs = 0
        total_pharmer_time = 0.0
                #################
        output_graph_dir  = args.output_dir / f"n{r}"
        output_json_dir   = output_graph_dir / "json"
        output_sdf_dir    = output_graph_dir / "sdf"
        
        output_graph_dir.mkdir(parents=True, exist_ok=True)
        output_json_dir.mkdir(parents=True, exist_ok=True)
        output_sdf_dir.mkdir(parents=True, exist_ok=True)
                # ---- ranked batch processing ----
        batch_size = args.workers * 2  # tune if needed
        for start_idx in range(0, len(subgraphs), batch_size):
            batch = subgraphs[start_idx:start_idx + batch_size]

            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                futures = [
                    executor.submit(
                        screen_and_save,
                        graph, info, exclusion_spheres,
                        db_paths, args.output_dir,
                        args.path_exe, r, start_idx + idx + 1, args.extension
                    )
                    for idx, (graph, info) in enumerate(batch)
                ]

                for future in as_completed(futures):
                    nhits, pharmer_time, success = future.result()
                    if success:
                        successful_runs += 1
                        total_pharmer_time += pharmer_time
                        n_hits += nhits
                    total_hits += nhits
            
            
            # stop after completing a full ranked batch
            if total_hits >= args.max_hits:
                logging.info(
                    f"Reached MAX_HITS limit ({args.max_hits}) after completing ranked batch."
                )
                break  # stop batching for this r-level


        total_time_graph = time.perf_counter () - start_time_graph
        # record summary for this size
        summary[f"n{r}"] = {
            "total_graphs": total_graphs,
            "successful_runs": successful_runs,
            "number_of_hits": n_hits,
            "total_pharmer_time": total_pharmer_time,
            "wall_time":total_time_graph
        }

        # stop completely after finishing this r-level
        if total_hits >= args.max_hits:
            logging.info(
                "Stopping further screening — top-ranked subgraphs fully processed for current size."
            )
            break
        # remove empty directores
        for dir_path in [output_sdf_dir, output_json_dir, output_graph_dir]:
            try:
                dir_path.rmdir()  # only removes if empty
            except OSError:
                pass

        logging.info(
            f"Completed n{r}: {summary[f'n{r}']['number_of_hits']} total hits "
            f"in {summary[f'n{r}']['total_pharmer_time']:.3f}s"
        )
        

        
    elapsed = time.perf_counter() - start_time
    write_summary(summary, elapsed, args.output_dir, db_paths)
    logging.info(f"Total hits: {total_hits}")
    return 0

 
def screen_and_save(graph, info, exclusion_spheres, db_paths, output_dir, exe_path, r, idx, extension):
    """Screen a generated pharmacophore, save only if hits found."""
    hits_total = 0
    total_time = 0.0
    

     # open temp file in text mode for JSON
    with tempfile.NamedTemporaryFile(suffix=".json", mode="w", encoding="utf-8", delete=False) as tmp:
        tmp_path = Path(tmp.name)
        data = {"extra_info": info, "points": []}
        for _, attr in graph.nodes(data=True):
            node = {}
            for k, v in attr.items():
                if k == "position":
                    node.update({"x": float(v[0]), "y": float(v[1]), "z": float(v[2])})
                elif k not in ["color", "value"]:
                    node[k] = v
            data["points"].append(node)
        data["points"].extend(exclusion_spheres)
        json.dump(data, tmp)

    success_any = False
    
    for db_path in db_paths:
        out_file = Path(output_dir / f"n{r}/sdf") / f"{db_path.stem}_n{r}_{idx}.{extension}"
        nhits, pharmer_time, success = run_pharmer(tmp_path, out_file, db_path, exe_path)
        total_time += pharmer_time
        if nhits == 0 or not success:
            out_file.unlink(missing_ok=True)
        else:
            success_any = True
        hits_total += nhits

    if hits_total > 0:
        graph_to_pharmacophore_json(graph, info, exclusion_spheres,
                                    output_path=Path(output_dir) / f"n{r}/json",
                                    filename=f"n{r}_{idx}.json")
        logging.info(f"n{r}_{idx}: {hits_total} hits found.")
    else:
        logging.debug(f"n{r}_{idx}: no hits found.")
    tmp_path.unlink(missing_ok=True)
    return hits_total, total_time, success_any


# --------------------- SUMMARY WRITER ---------------------
def write_summary(summary, total_time, output_dir, database_paths):
    summary_file = Path(output_dir) / "screening_summary.dat"
    header_fmt = "{:^20} {:^15} {:^15} {:^25} {:^25} {:^25}\n"
    row_fmt = "{:^20} {:^15} {:^15} {:^25.2f} {:^25.2f} {:^25.2f}\n"

    with open(summary_file, "w") as f:
        f.write(f"Pharmer Screening Summary - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 130 + "\n")
        f.write(f"Databases used: {[db.name for db in database_paths]}\n")
        f.write(header_fmt.format("Pharmacophore", "Total graphs", "Total hits",
                                  "Total Wall Time(s)", "Avg. Pharmer time (s)", "%age of Graphs with Hits"))
        f.write("-" * 130 + "\n")

        total_hits = 0
        total_graphs = 0
        total_successful = 0

        for model, data in summary.items():
            success_rate = (data['successful_runs'] / data['total_graphs'] * 100) if data['total_graphs'] > 0 else 0
            avg_pharmer_time = (data['total_pharmer_time'] / data['successful_runs']) if data['successful_runs'] > 0 else 0

            f.write(row_fmt.format(
                model, data['total_graphs'], data['number_of_hits'],
                data ['wall_time'], avg_pharmer_time, success_rate
            ))

            total_hits += data['number_of_hits']
            total_graphs += data['total_graphs']

        f.write("-" * 130 + "\n")
        f.write(f"Total Hits Across All Models: {total_hits}\n")
        f.write(f"Total Execution Time: {total_time:.2f} seconds\n")
        f.write(f"Output Directory: {output_dir}\n")
        f.write(f"Databases Searched: {[str(db) for db in database_paths]}\n")
         
    return summary_file


if __name__ == "__main__":
    import sys
    sys.exit(main())


