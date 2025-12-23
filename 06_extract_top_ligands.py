from collections import defaultdict
import copy
import argparse
import pickle
import subprocess
import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
import os
######
import warnings
# Suppress specific warnings from MDAnalysis
warnings.filterwarnings("ignore")#, category=UserWarning, module="MDAnalysis.coordinates.PDB")



# List of your individual dictionaries
def sort_and_join_dicts (dict_list,output_path, output_obj="master_dict.pkl",sorting_key="bsa_harm"):
    
    if len (dict_list.strip().split()) == 1: # if just one dictionary
        
        with open(dict_list, 'rb') as f:
            
            dict_ = pickle.load(f)
        
            sorted_merged_dict =   {
            key: dict(sorted(
                value.items(),
                key=lambda item: item[1][sorting_key],
                reverse=True
            ))
            for key, value in dict_.items()  }
    
          
    else:
        
        dicts = []
        for file in dict_list.strip().split ():
            with open(file, 'rb') as f:
                d = pickle.load(f)
                dicts.append(d)
                
        # Initialize a new merged dictionary
        # merged_dict = defaultdict(dict)
        # Initialize a new merged dictionary
        merged_dict = defaultdict(dict)
        
        # Merge all dictionaries
        for d in dicts:
            for top_key, molecules in d.items():
                if top_key not in merged_dict:
                    merged_dict[top_key] = copy.deepcopy(molecules)
                else:
                    merged_dict[top_key].update(molecules)
        
        # Convert back to normal dict if you don't want defaultdict
        merged_dict = dict(merged_dict)
        
        # sort the merged dictionaries
        merged_dict = {
            key: dict(sorted(
                value.items(),
                key=lambda item: item[1][sorting_key],
                reverse=True
            ))
            for key, value in merged_dict.items()
        }
        
        ##
        sorted_merged_dict = {str(k): merged_dict[str(k)] for k in sorted(int(key) for key in merged_dict)}
        
    for key, val in sorted_merged_dict.items ():
        print (f"Total {len(sorted_merged_dict[key])} in graph {key}")
        
    # save the file
    out_file = os.path.join (output_path, output_obj)
    with open (out_file, "wb") as outfile:
        pickle.dump (sorted_merged_dict, outfile)
        
    return sorted_merged_dict

 
def create_universe (name, positions):
    n_atoms = len(positions)
    u_new = mda.Universe.empty(n_atoms=n_atoms,
                             n_residues=n_atoms,
                             atom_resindex=np.arange (n_atoms),
                             residue_segindex=np.arange (n_atoms),
                             n_segments=n_atoms,
                             trajectory=True) # necessary for adding coordinate


    u_new.add_TopologyAttr('name',   name)
    u_new.add_TopologyAttr('elements', name)
    u_new.add_TopologyAttr('resid', n_atoms*["1"])
    u_new.add_TopologyAttr('resname', n_atoms*["LIG"])
    u_new.add_TopologyAttr('segid', n_atoms*["1"])
    u_new.add_TopologyAttr('chainID', n_atoms * ["C"])
    u_new.atoms.positions = positions
    
    return u_new

def parse_single_sdf_file(file_path):
    """
    Parse SDF file, returning a dictionary:
    {ZINC_ID: {"xyz": [...], "atom_radii": [...], "rmsd": ...}}
    """
    zinc_dict = {}
    with open(file_path, "r") as f:
        
        data, atom_name = [], []
        zinc_id = None
        current_rmsd = None
        read_rmsd = False
        
        for line in f:
            line = line.strip()
            if not line:
                continue

            if "ZINC" in line:
                zinc_id = line
                if zinc_id not in zinc_dict:
                    zinc_dict[zinc_id] = {"xyz": [], "name": [], "rmsd":[]}

            elif "$$$$" in line:
                if zinc_id is not None:
                    zinc_dict[zinc_id]["xyz"].append(np.array(data))
                    zinc_dict[zinc_id]["name"] = atom_name
                    zinc_dict[zinc_id]["rmsd"].append(current_rmsd)
                    
                data, atom_name = [], []
                current_rmsd = None
                read_rmsd = False

            elif len(line.split()) == 9: 
                # Thepharmer output sdf has 9 coloumns for coordinate sectino
                parts = line.split()
                try: 
                    xyz = [float(x) for x in parts[:3]]
               
                    data.append(xyz)
                    atom_name.append (parts[3])
                except(ValueError):
                    pass    
                    
            elif "rmsd" in line.lower():
                read_rmsd = True

            elif read_rmsd:
                try:
                    current_rmsd = float(line)
                except ValueError:
                    current_rmsd = None
                read_rmsd = False            
            
                
    return zinc_dict
    

 

def extract_ligands(sorted_merged_dict, out_path, top_N=10, top_Ngraphs=20):
    
    sdf_dict_list = {}
    Saved_keys = []
    # sort according to graph keys, since graphs are ranked by thte keys 

    #

    #
    gro_file_previous =  None # to laod mda universe only if its a differnt gro file
    
    sorted_keys = sorted(sorted_merged_dict.keys(), key=int)
    total_ligands = 0
    for graph_index, key in enumerate(sorted_keys):
        if graph_index >= top_Ngraphs: 
            break
        
        path = os.path.join(out_path, key)
        os.makedirs(path, exist_ok=True)

        n_ligands = min(top_N, len(sorted_merged_dict[key]))
        
        for index, (zinc_id, value) in tqdm (enumerate(sorted_merged_dict[key].items()), desc=f"Processing Graph {key}" ):
            total_ligands += 1
            if zinc_id in Saved_keys:
                continue
            if index >= n_ligands:
                break

            ligand_path = os.path.join(path, f"ligand{index+1}")
            os.makedirs(ligand_path, exist_ok=True)

            sdf_name = os.path.basename(os.path.normpath(value["sdf_name"]))
   
            if sdf_name not in sdf_dict_list:
                sdf_dict_list[sdf_name] = parse_single_sdf_file(value["sdf_name"])

            ligand_index = np.argmin (sdf_dict_list[sdf_name][zinc_id]['rmsd'])
            ligand_xyz = sdf_dict_list[sdf_name][zinc_id]['xyz'][ligand_index]
            name = sdf_dict_list[sdf_name][zinc_id]['name']

            # Create ligand PDB
            u_lig = create_universe(name, ligand_xyz)
            ligand_pdb_path = os.path.join(ligand_path, f"{zinc_id}.pdb")
            u_lig.atoms.write(ligand_pdb_path)

            
            # Load protein and write merged PDB
            gro_file_current =   value["gro"]

            # only if the gro file is differetn from previous one we open a new one
            # loading MDA universe  every time is slow
            if not gro_file_current == gro_file_previous:
                u = mda.Universe(value["tpr"], gro_file_current)
                u.atoms.fragments [0].chainIDs = "A"
                u.atoms.fragments [1].chainIDs = "B"

            
            protein_pdb_path = os.path.join(ligand_path, "protein.pdb")
            u.select_atoms("protein").write(protein_pdb_path, bonds=None)
            gro_file_previous = gro_file_current

            # Run Open Babel to generate mol2
            
            
            ligand_mol2_path = os.path.join(ligand_path, "ligand.mol2")
            command = f"obabel {ligand_pdb_path} -O {ligand_mol2_path} -p 7.0 --partialcharge gasteiger --minimize --steps 50 --sd"
            result = subprocess.run(command, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                print(f"OpenBabel failed for {ligand_pdb_path}:")
                print(result.stderr)
                continue  # skip this ligand

            if not os.path.exists(ligand_mol2_path) or os.path.getsize(ligand_mol2_path) == 0:
                print(f"Mol2 file missing or empty for {ligand_pdb_path}, skipping...")
                continue

            # Load ligand mol2 safely
            # IN some cases obabel wirtes multiples residue names which causes probelme with charge derivation
            ligand_mol2 = mda.Universe(ligand_mol2_path)
            ligand_mol2.atoms.residues.resnames = ['LIG']*len(ligand_mol2.atoms.residues.resnames) 
            ligand_mol2.atoms.write (ligand_mol2_path)
            
            total_charge = int(np.round(sum(ligand_mol2.atoms.charges)))
            
            # Write acpype command
            acpype_command = f"python /mnt/second/acpype.py -i  ligand.mol2 -a gaff2 -o gmx -c bcc -n {total_charge}"
            acpype_script_path = os.path.join(ligand_path, "get_ligand_itp.sh")
            with open(acpype_script_path, "w") as f:
                f.write ("# Here put the correct acpype.py path, or you can use the same command for antechamber\n")
                f.write(acpype_command)

            Saved_keys.append(key)
    print (f"extracted {total_ligands} ligands ")

            
def main():
    parser = argparse.ArgumentParser(description="Merge dictionaries and extract top ligands")
    parser.add_argument('-i', dest='dict_list', nargs="+", type=str, help='List of space-separated pickle dictionaries')
    parser.add_argument('-topN', dest='top_N', type=int, default=10, help='Top N ligands per graph')
    parser.add_argument('-Ngraph', dest='top_Ngraphs', type=int, default=5, help='Number of top graphs to process')
    parser.add_argument('-odir', dest='output_path', type=str, default='./final-results', help='Output directory')
    parser.add_argument('-skey', dest='sorting_key', type=str, default='bsa_harm', help='Sorting key: bsa_harm or bsa_ab')    
    args = parser.parse_args()
    
    dict_list = " ".join(args.dict_list)
    output_path = args.output_path
    top_N = args.top_N
    top_Ngraphs = args.top_Ngraphs
    sorting_key = args.sorting_key
  
    # Ensure output directory exists
    os.makedirs(output_path, exist_ok=True)
    
    # Sort and merge dictionaries
    sorted_merged_dict = sort_and_join_dicts(dict_list, output_path, sorting_key=sorting_key)
    
    # Extract ligands
    extract_ligands(sorted_merged_dict, output_path, top_N=top_N, top_Ngraphs=top_Ngraphs)

if __name__ == "__main__":
    
    main() 
