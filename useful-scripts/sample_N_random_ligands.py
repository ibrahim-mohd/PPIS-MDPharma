import gzip
import random
import argparse
from tqdm import tqdm

def get_number_of_molecules(input_file,extension):
    
    pattern = dict (sdf="$$$$", mol2="@<TRIPOS>MOLECULE")
    
    """Count the number of molecules in a gzipped SDF file."""
    count = 0
    with gzip.open(input_file, 'rt') as f:
 
        for line in tqdm(f, desc="Calculating number of ligand in the parent file"):
            if line.strip() == pattern [extension]:
                count += 1
    print (f"\n Input file has {count} ligands\n")
    return count

def sample_ligands(input_file, sample_size, extension):
    """
    Sample ligands randomly from a gzipped SDF/MOL2 file using reservoir sampling.
    Returns a list of ligands.
    """
    pattern = dict (sdf="$$$$", mol2="@<TRIPOS>MOLECULE")
    total_molecules = get_number_of_molecules(input_file, extension)
    fraction = sample_size / total_molecules

    reservoir = []
    current_ligand = []

    with gzip.open(input_file, 'rt') as f:
        for line in tqdm(f, desc=f"Random sampling of {sample_size} ligands from parent file"):
            current_ligand.append(line)
            
            if line.strip() == pattern [extension]:  # End of ligand
                if random.random() <= fraction:
                    reservoir.append(current_ligand)
                current_ligand = []

    return reservoir

def main():
    """Parse arguments and perform sampling."""
    parser = argparse.ArgumentParser(
        description="Fetch a random subset of compounds from a large gzipped SDF file")
    
    parser.add_argument( "-i", "--input_file", type=str, required=True, help="Input filename of gzipped SDF/MOL2 file")
    parser.add_argument( "-o", "--output_name", type=str, required=True, help="output_file of gzipped file with .gz extension")
    parser.add_argument( "-n", "--sample_size", type=int, default=100_000,help="Number of ligands to sample")
    parser.add_argument( "--extension", "-e", type=str, choices=["sdf", "mol2"], # restrict to valid options
        default="sdf", help="Extensin of input files (default: sdf). Can be 'sdf' or 'mol2'.")
    
    args = parser.parse_args()

    # Sample ligands
    reservoir = sample_ligands(args.input_file, args.sample_size, args.extension)
    print(f"Sampled {len(reservoir)} ligands")

    # Write sampled ligands to output gzipped SDF
    with gzip.open(args.output_name, 'wt') as out:
        for ligand in reservoir:
            out.writelines(ligand)

if __name__ == "__main__":
    main()


