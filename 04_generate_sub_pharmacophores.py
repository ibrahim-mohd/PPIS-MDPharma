# Written by Mohd Ibrahim
# Technical University of Munich
import os
import math
import json
import pickle
import itertools
from tqdm import tqdm
import argparse
import numpy as np
import networkx as nx
import warnings

warnings.filterwarnings("ignore")


NODE_COLORS = {
    "Aromatic": "orange",
    "HydrogenDonor": "magenta",
    "HydrogenAcceptor": "green",
    "Hydrophobic": "cyan",
    "PositiveIon": "red",
    "ExclusionSphere": "gray",
    "NegativeIon": "blue"
}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate all possible pharmacophore graphs from a master file")
    parser.add_argument('-j', dest='pharmacophore_file', type=str, default='pharma.json',
                        help='Master pharmacophore JSON file')
    parser.add_argument('-min_node', dest='min_node', type=int, default=5,
                        help='Minimum number of nodes to consider in combinations')
    parser.add_argument('-top', dest='top_percentage', type=int, default=20,
                        help='Percentage of top-ranked models to keep')
    parser.add_argument('-ntop_limit', dest='ntop_limit', type=int, default=20,
                        help='If the total number of pharmacophore is less or equal to this limit, then the top percentage value and instead writes this many sub-pharmacophores')
    parser.add_argument('-o', dest='out_path', type=str, default='./',
                        help='Output path for JSON files')
    return parser.parse_args()


def get_distance(point1: np.ndarray, point2: np.ndarray) -> float:
    return np.linalg.norm(np.array(point1) - np.array(point2))


def create_directory(path: str):
    """Safely create a directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)


def build_graph_from_pharmacophore(pharma_data: dict):
    """Create a NetworkX graph from a pharmacophore JSON structure."""
    G = nx.Graph()
    exclusion_spheres = []
    index = 0

    for site in pharma_data["points"]:
        if site['name'] != 'ExclusionSphere':
            index += 1
            G.add_node(
                index,
                name=site['name'],
                position=np.array([site['x'], site['y'], site['z']]),
                radius=site['radius'],
                size=site['size'],
                requirement=site['requirement'],
                enabled=site['enabled'],
                svector=site['svector'],
                vector=site['vector'],
                minsize=site['minsize'],
                maxsize=site['maxsize'],
                vector_on=site['vector_on'],
                color=NODE_COLORS[site['name']],
                score=site['score'],
                chainID=site['chainID']
            )
        else:
            exclusion_spheres.append(site)

    # Add all edges
    for i, j in itertools.combinations(list(G.nodes), 2):
        score_sum = G.nodes[i]['score']['normed_score'] + G.nodes[j]['score']['normed_score']
        distance = get_distance(G.nodes[i]['position'], G.nodes[j]['position'])
        G.add_edge(i, j, distance=distance, score=score_sum)

    return G, exclusion_spheres


def generate_graphs(G: nx.Graph, r: int, top_percentage: int = 0,ntop_limit: int=0 ):
    """Generate all subgraphs with r nodes and score them."""
    graphs = []

    for nodes in itertools.combinations(list(G.nodes), r):
        new_graph = G.subgraph(nodes).copy()
        min_spanning_tree_dist = sum(data["distance"] for _, _, data in nx.minimum_spanning_tree(new_graph, weight='distance').edges(data=True))
        count_A = sum(1 for _, attr in new_graph.nodes(data=True) if attr.get('chainID') == 'A')
        count_B = sum(1 for _, attr in new_graph.nodes(data=True) if attr.get('chainID') == 'B')
        symmetry = min(count_A, count_B) / max(count_A, count_B) if max(count_A, count_B) > 0 else 0
        total_score = sum(attr['score']['normed_score'] for _, attr in new_graph.nodes(data=True))

        graphs.append([
            new_graph,
            dict(
                min_spanning_tree_dist=min_spanning_tree_dist,
                score=total_score,
                symmetry=symmetry,
                combined_score=symmetry * total_score
            )
        ])

    sorted_graphs = sorted(graphs, key=lambda x: x[1]['combined_score'], reverse=True)

    if top_percentage:
        N = max(ntop_limit, int(np.ceil(len(sorted_graphs) * top_percentage / 100)))
        return sorted_graphs[:N]

    return sorted_graphs


def graph_to_pharmacophore_json(graph: nx.Graph, extra_info: dict, exclusion_spheres: list, output_path: str, filename: str):
    """Convert a NetworkX graph back into a pharmacophore JSON and save it."""
    data = {"extra_info": extra_info, "points": []}

    for _, attr in graph.nodes(data=True):
        node_dict = {}
        for key, value in attr.items():
            if key == 'position':
                node_dict.update({"x": float(value[0]), "y": float(value[1]), "z": float(value[2])})
            elif key in ['color', 'value']:
                continue
            else:
                node_dict[key] = value
        data["points"].append(node_dict)

    data["points"].extend(exclusion_spheres)
    create_directory(output_path)
    with open(os.path.join(output_path, filename), "w") as f:
        json.dump(data, f, indent=4)


def main():
    args = parse_args()
    with open(args.pharmacophore_file, 'r') as f:
        pharma_data = json.load(f)

    G, exclusion_spheres = build_graph_from_pharmacophore(pharma_data)

    all_graphs = {}
    for r in tqdm(range(args.min_node, len(G.nodes) + 1), desc="Generating graphs"):
        subgraphs = generate_graphs(G, r, top_percentage=args.top_percentage, ntop_limit=args.ntop_limit)
        all_graphs[f"n{r}"] = subgraphs

    for key, graphs in tqdm(all_graphs.items(), desc="Writing JSON files"):
        for idx, (graph, info) in enumerate(graphs):
            filename = f"{key}_{idx + 1}.json"
            graph_to_pharmacophore_json(graph, info, exclusion_spheres, output_path=os.path.join(args.out_path, key), filename=filename)


if __name__ == "__main__":
    main()

