# Written by Mohd Ibrahim
# Technical University of Munich

"""
Parallel pharmacophore screening using Pharmer.
Searches multiple pharmacophore models against compound databases.
"""

import subprocess
import re
import logging
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import argparse
import sys
#from typing import Tuple, Dict, List, Optional

   
# --------------------- LOGGING SETUP ---------------------
def setup_logging(output_dir, verbose = False):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(output_dir / "screening.log")
        ]
    )

# --------------------- ARGUMENT PARSING ---------------------
def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parallel pharmacophore screening")

    # pharmer exe path
    
    parser.add_argument("--path_exe","-p_exe", type=Path, default="/usr/local/bin/pharmer.static",
                       help="full path to pharmer executable")

        
    # input output and database paths
    parser.add_argument("--database-file", "-df", type=Path,
                       help="Text file containing database paths (one per line)")
    
    parser.add_argument("--input_dir", "-i", type=Path, required=True, help="Input directory with pharamcophore subdirectories")

    parser.add_argument("--output_dir", "-o", type=Path, required=True, help="Output directory for results and logs (required)")

    ###################################################################################
    
    parser.add_argument("--workers","-np", type=int, default=4, help=f"Number of parallel workers 4")
    
    parser.add_argument("--pattern", type=str, default="n*", help="Pattern to match pharmacophore directories (default: n*)")

    parser.add_argument( "--extension", "-e", type=str, choices=["sdf", "mol2"],  # restrict to valid options
        default="sdf", help="Output file extension for results (default: sdf). Can be 'sdf' or 'mol2'.")

    parser.add_argument("--database", "-d", action="append", type=Path, help="Database path (can be used multiple times)")

    parser.add_argument("--max_hits","-max", type=int, default=10000, 
                        help="Stop screening when total hits reach this number, it  may be much higher than this limit , read comments at line 380 in the code")
 
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Enable verbose output")
    return parser.parse_args()

# --------------------- DATABASE PATH HANDLING ---------------------
def load_database_paths(args):
    """Load database paths from command line arguments and/or file."""
    
    database_paths = []
    
    # Add paths from command line arguments
    if args.database:
        database_paths.extend(args.database)
    
    # Add paths from database file
    if args.database_file:
        if not args.database_file.exists():
            logging.error(f"Database file not found: {args.database_file}")
            sys.exit(1)
        
        try:
            with args.database_file.open('r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if line and not line.startswith('#'):  # Skip empty lines and comments
                        db_path = Path(line)
                        if db_path not in database_paths:  # Avoid duplicates
                            database_paths.append(db_path)
            logging.info(f"Loaded {len(database_paths) - (len(args.database) if args.database else 0)} database paths from {args.database_file}")
        except Exception as e:
            logging.error(f"Error reading database file {args.database_file}: {e}")
            sys.exit(1)
    
    # Remove duplicates while preserving order
    seen = set()
    unique_paths = []
    for path in database_paths:
        if path not in seen:
            seen.add(path)
            unique_paths.append(path)
    
    return unique_paths

# --------------------- VALIDATION FUNCTIONS ---------------------
def validate_environment(input_dir, output_dir, database_paths, pharmer_exe_path):
    """Validate that required executables and directories exist."""
    if not pharmer_exe_path.exists():
        logging.error(f"Pharmer executable not found: {pharmer_exe_path}")
        return False
    
    valid_paths = []
    for db_path in database_paths:
        if db_path.exists():
            valid_paths.append(db_path)
        else:
            logging.warning(f"Database path not found, skipping: {db_path}")
    
    if not valid_paths:
        logging.error("No valid database paths found")
        return False
    
    # Update database_paths to only include valid ones
    database_paths[:] = valid_paths
    
    if not input_dir.exists():
        logging.error(f"Pharmacophore directory not found: {input_dir}")
        return False
    
    output_dir.mkdir(parents=True, exist_ok=True)
    return True

# --------------------- HELPER FUNCTIONS ---------------------
def run_pharmer(input_file, output_file, database_path, pharmer_exe_path):

    # This pharmer command
    # pharmer dbsearch -dbdir=LocalDatabase -in input.json -out ouput.sdf
    
    cmd = [
        str(pharmer_exe_path),
        "dbsearch",
        f"-dbdir={database_path}",
        f"-in={input_file}",
        f"-out={output_file}"
    ]
    
    nhits = 0 # number of ligands satisfying the pharmacophroe model
    pharmer_time = 0.0
    success = False # incase issues for json file,

    try:
        logging.debug(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=None)  # 1 hour timeout
        
        # Parse output
        hits_match = re.search(r"NumResults:\s*(\d+)", result.stdout)
        time_match = re.search(r"Time:\s*([\d.]+)", result.stdout)
        
        if hits_match:
            nhits = int(hits_match.group(1))
        if time_match:
            pharmer_time = float(time_match.group(1))
            
        success = True # everythign ran fine whether found a hit or not
        
        logging.debug(f"Completed {input_file.name}: {nhits} hits in {pharmer_time:.2f}s")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Pharmer failed on {input_file.name}: {e}")
        if e.stderr:
            logging.error(f"Pharmer stderr: {e.stderr.strip()}")
  

    return nhits, pharmer_time, success
    

def cleanup_empty_files(output_file):
    """Remove empty or low-hit output files."""
    try:
        if output_file.exists():
            if output_file.stat().st_size == 0:
                output_file.unlink(missing_ok=True)
                logging.debug(f"Removed empty file: {output_file}")
    except OSError as e:
        logging.warning(f"Could not check/remove {output_file}: {e}")

def find_pharmacophore_dirs(input_dir, pattern = "n*"):
    """Find and sort pharmacophore directories matching pattern."""
    directories = []
    for d in input_dir.iterdir():
        if d.is_dir() and re.fullmatch(pattern.replace("*", ".*"), d.name):
            # Extract numeric part for sorting
            match = re.search(r"(\d+)", d.name)
            if match:
                n_value = int(match.group(1))
                directories.append((n_value, d))
    
    return sorted(directories, key=lambda x: x[0], reverse=True)

# --------------------- MAIN PROCESSING ---------------------
def process_pharmacophore_model(input_dir, output_base, database_paths, pharmer_exe_path,max_workers, extension=".sdf"):
    
    """Process a single pharmacophore model directory against multiple databases."""
    
    pharma_model      = input_dir.name
    model_output_path = output_base / pharma_model
    model_output_path.mkdir(exist_ok=True)
    
    json_files = list(input_dir.glob("*.json"))

    # print info
    logging.info(f"Processing {pharma_model} with {len(json_files)} graphs against {len(database_paths)} databases")
    
    model_data = {
        "successful_runs":0,
        "total_graphs": len(json_files),
        "number_of_hits": 0,
        "wall_time": 0.0,
        "total_pharmer_time": 0.0,
        "databases_used": [db.name for db in database_paths]
    }

    # put a timer to print final time
    model_start = time.perf_counter()
    
    # Process each database sequentially for this model
    for db_path in database_paths:
        
        logging.info(f"Searching against database: {db_path.name}")
        
        # Parallel execution for this database
        futures = {}
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for json_file in json_files:
                out_file = model_output_path / f"{db_path.name}_{json_file.stem}.{extension}"
                
                futures[executor.submit(run_pharmer, json_file, out_file, db_path,pharmer_exe_path)] = (json_file, out_file)
                
            
            for future in as_completed(futures):
                json_file, out_file = futures[future]
                try:
                    # can set time out to a day i.e 86400 seconds
                    # 
                    nhits, pharmer_time, success = future.result(timeout=None)
                    
                    if success:
                        model_data["successful_runs"] += 1
                        model_data["number_of_hits"]  += nhits
                        model_data["total_pharmer_time"] += pharmer_time
                        
                        if nhits > 0:
                            logging.info(f"{json_file.name} against {db_path.name}: {nhits} ligands found")
                        else:
                            logging.debug(f"{json_file.name} against {db_path.name}: no hits")
                        
                        # Cleanup if empty files
                        if nhits ==0:
                            cleanup_empty_files(out_file)
                    else:
                        # remove empty files
                        cleanup_empty_files(out_file)
                        
                except Exception as e:
                    logging.error(f"Failed to process {json_file.name} against {db_path.name}: {e}")
                    json_file, out_file = futures[future]
                    cleanup_empty_files(out_file)
    
    model_end = time.perf_counter()
    model_data["wall_time"] = model_end - model_start
    
    # Remove empty output directory
    try:
        if not any(model_output_path.iterdir()):
            model_output_path.rmdir()
            logging.debug(f"Removed empty directory: {model_output_path}")
    except OSError:
        pass
    
    logging.info(f"Completed {pharma_model}: {model_data['successful_runs']}/"
                f"{model_data['total_graphs']} successful, {model_data['number_of_hits']} total hits, "
                f"in {model_data['wall_time']:.2f}s")
    
    return model_data

def write_summary(summary, total_time, output_dir, database_paths):
    """Write screening summary to file."""
    summary_file = output_dir / "screening_summary.dat"
    
    # Enhanced formatting with more statistics
    #header_fmt = "{:<15} {:<15} {:<12} {:<15} {:<20} {:<20}\n"
    #row_fmt    = "{:<15} {:<15} {:<12} {:<15.2f} {:<20.2f} {:<20.2f}\n"
    header_fmt = "{:^20} {:^15} {:^15} {:^25} {:^25} {:^25}\n"
    row_fmt    = "{:^20} {:^15} {:^15} {:^25.2f} {:^25.2f} {:^25.2f}\n"

    with summary_file.open("w") as f:
        f.write(f"Pharmer Screening Summary - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 130 + "\n")
        f.write(f"Databases used: {[db.name for db in database_paths]}\n")
        f.write(header_fmt.format("Pharmacophore", "Total graphs", "Total hits", 
                                "Total Wall Time(s)", "Avg. Pharmer time (s)", "%age of Graphs Screened"))
        f.write("-" * 130 + "\n")
        
        total_hits = 0
        total_graphs = 0
        total_successful = 0
        
        for model, data in summary.items():
            success_rate = (data['successful_runs'] / data['total_graphs'] * 100) if data['total_graphs'] > 0 else 0
            avg_pharmer_time = (data['total_pharmer_time'] / data['successful_runs']) if data['successful_runs'] > 0 else 0
            
            f.write(row_fmt.format(
                model, data['total_graphs'], data['number_of_hits'], 
                data['wall_time'], avg_pharmer_time,success_rate
            ))
            
            total_hits  += data['number_of_hits']
            total_graphs += data['total_graphs']
         
        f.write("-" * 130 + "\n")
        f.write(f"Total Hits Across All Models: {total_hits}\n")
        f.write(f"Total Execution Time: {total_time:.2f} seconds\n")
        f.write(f"Output Directory: {output_dir}\n")
        f.write(f"Databases Searched: {[str(db) for db in database_paths]}\n")
    
    return summary_file 

def main():
    """Main execution function."""
    args = parse_arguments()
    
    output_dir = args.output_dir
    input_dir  = args.input_dir
    MAX_HITS   = args.max_hits
    pharmer_exe_path = args.path_exe
    
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(output_dir, args.verbose)
    
    logging.info("Starting pharmacophore screening")
    
    # Load database paths from various sources
    database_paths = load_database_paths(args)
    
    if not database_paths:
        logging.error("No database paths specified")
        return 1
    
    logging.info(f"Using {len(database_paths)} database(s): {[db.name for db in database_paths]}")
    
    if not validate_environment(input_dir, output_dir, database_paths, pharmer_exe_path):
        return 1
    
    # Find pharmacophore directories
    pharma_subdirs = find_pharmacophore_dirs(input_dir, args.pattern)
    if not pharma_subdirs:
        logging.error(f"No pharmacophore directories found matching pattern '{args.pattern}'")
        return 1

    # start time
    start_time = time.perf_counter()
    summary = {}

    
    total_hits = 0
    
    try:
        for n_value, input_dir in pharma_subdirs:
            model_data  = process_pharmacophore_model (input_dir, output_dir, database_paths,pharmer_exe_path, args.workers,args.extension)
            summary[input_dir.name] = model_data
            
            total_hits += model_data["number_of_hits"]

            # The --max_hits value is a soft global limit checked only *between* pharmacophore models.
            #As a result, the total number of hits may exceed the specified MAX_HITS threshold, 
            #sometimes by a large margin. For example:
            #If MAX_HITS=10000, and pharmacophore n7 produces 3000 hits,
            #screening will continue with n6. If n6 produces 20000 hits,
            #the final total will be 23000 — exceeding the limit by design.             
            
            if total_hits >= MAX_HITS:
                logging.info(f"Reached MAX_HITS limit ({MAX_HITS}). Stopping further screening.")
                break
            
    except KeyboardInterrupt:
        logging.info("Screening interrupted by user")
        return 1
    except Exception as e:
        logging.error(f"Unexpected error during screening: {e}")
        return 1
    
    end_time = time.perf_counter()
    total_time = end_time - start_time
    
    # Write summary
    summary_file = write_summary(summary, total_time, output_dir, database_paths)
    
    logging.info(f"Screening complete in {total_time:.2f} seconds")
    logging.info(f"Summary written to {summary_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

