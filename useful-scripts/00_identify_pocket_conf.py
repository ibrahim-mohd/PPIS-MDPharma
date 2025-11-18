import os
import pickle
import warnings
import numpy as np
from tqdm import tqdm
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import argparse
 
from pathlib import Path
import shutil
import json
import subprocess

 

# Suppress warnings from MDAnalysis
warnings.filterwarnings("ignore")

# MDTraj atomic radii
_ATOMIC_RADII = md.geometry.sasa._ATOMIC_RADII  # still private, but safer than hardcoding

def element_symbol_string(name):
    """Guess element symbol from atom name."""
    element = guess_atom_element(name)
    if len(element) == 2:
        element = element[0].upper() + element[1].lower()
    return element

def get_element_symbols(atom_names):
    return [element_symbol_string(name) for name in atom_names]

def Sasa_calc(positions, atom_radii, probe_radius=0.14, n_sphere_points=960):
    """
    Compute SASA using MDTraj's public API.
    positions: np.array of shape (n_atoms, 3) in Angstroms
    atom_radii: list or np.array of atom radii in Angstroms
    Returns: total SASA in Å^2
    """
    # Convert positions from Å to nm
    xyz_nm = positions / 10.0
    xyz_nm = xyz_nm[np.newaxis, :, :]  # shape (1, n_atoms, 3)

    radii = np.array(atom_radii, dtype=np.float32) + probe_radius
    n_atoms = xyz_nm.shape[1]
    out = np.zeros((1, n_atoms), dtype=np.float32)
    atom_indices = np.arange(n_atoms, dtype=np.int32)
    mask = np.ones(n_atoms, dtype=np.int32)

    # Use MDTraj private _geometry._sasa (no public API for exact per-atom)
    from mdtraj.geometry import _geometry
    _geometry._sasa(xyz_nm, radii, n_sphere_points, atom_indices, mask, out)

    return float(np.sum(out))  # total SASA in Å^2
 
# pockets
def sanitize_pocket_spheres (pocket, cutoff=1e-4): 
    pairs, distance = mda.lib.distances.capped_distance (pocket.positions, pocket.positions, max_cutoff=cutoff)
    index0 = pocket.indices [0]
    set_ = {i  for i, j in pairs if i != j} | {j for i, j in pairs if i != j}

    index_to_remove = " ".join (str (x +index0) for x in set_)
   # new_pocket 
    return index_to_remove



def get_bsa(dict_protein, dict_ligand):
    """Compute buried surface areas for ligand with protein chains."""

    ligand_xyz = dict_ligand['xyz']
    ligand_radii = dict_ligand['radii']

    def concat_coords_radii(coords, radii):
        return np.concatenate((coords, np.array(ligand_xyz, dtype=np.float32)), axis=0), \
               np.concatenate((radii, np.array(ligand_radii, dtype=np.float32)), axis=0)

    complex_coords, complex_radii = concat_coords_radii(dict_protein['xyz']['protein'], dict_protein['radii']['protein'])
    AL_coords, AL_radii = concat_coords_radii(dict_protein['xyz']['chainA'], dict_protein['radii']['chainA'])
    BL_coords, BL_radii = concat_coords_radii(dict_protein['xyz']['chainB'], dict_protein['radii']['chainB'])

    sasa_ligand = Sasa_calc(ligand_xyz, ligand_radii)
    sasa_complex = Sasa_calc(complex_coords, complex_radii)
    sasa_AL = Sasa_calc(AL_coords, AL_radii)
    sasa_BL = Sasa_calc(BL_coords, BL_radii)
    sasa_chain_A = dict_protein['sasa']['chainA']
    sasa_chain_B = dict_protein['sasa']['chainB']
    sasa_protein = dict_protein['sasa']['protein']

    bsa_al = sasa_chain_A + sasa_ligand - sasa_AL
    bsa_bl = sasa_chain_B + sasa_ligand - sasa_BL
    bsa_abl = sasa_protein + sasa_ligand - sasa_complex

    return np.round(bsa_al,3), np.round(bsa_bl,3), np.round(bsa_abl,3)


def load_protein (pdb_file, pocket_resname="STP"):
    
    """Load protein and chain data, compute SASA and atom radii."""
    u  = mda.Universe (pdb_file)
    
    protein = u.select_atoms("protein and not name H*")
    chainA  = u.select_atoms ("chainID A and not name H*")
    chainB  = u.select_atoms ("chainID B and not name H*")
 
    def compute_radii_sasa(atoms):
        symbols = get_element_symbols(atoms.names)
        radii = [ _ATOMIC_RADII.get(sym) for sym in symbols ]
        sasa = Sasa_calc(atoms.positions, radii)
        return radii, sasa

    radii_protein, sasa_protein = compute_radii_sasa(protein)
    radii_chainA, sasa_chainA   = compute_radii_sasa(chainA)
    radii_chainB, sasa_chainB   = compute_radii_sasa(chainB)

    protein_dict =  {
        "sasa": {"protein": sasa_protein, "chainA": sasa_chainA, "chainB": sasa_chainB},
        "radii": {"protein": radii_protein, "chainA": radii_chainA, "chainB": radii_chainB},
        "xyz": {"protein": protein.positions, "chainA": chainA.positions, "chainB": chainB.positions} }
     
    # pockets
    
    interface_pocket_ids = u.select_atoms (f"resname {pocket_resname} and around 5 chainID A and around 5 chainID B").residues.resids

    pocket_bsa = {}
    
    for resid in interface_pocket_ids:
        pocket = u.select_atoms (f"resname {pocket_resname} and resid {resid}") 
        
        remove_index = sanitize_pocket_spheres (pocket, cutoff=1e-4)        
        if remove_index:
            pocket = u.select_atoms (f"resname {pocket_resname} and resid {resid} and not index {remove_index}") 
            
        radius = [ _ATOMIC_RADII["C"] for i in range (len(pocket.names))]
         
        pocket_dict= dict (xyz=pocket.positions, radii=radius)

        bsa_al, bsa_bl, bsa_abl = get_bsa(protein_dict, pocket_dict)
        bsa_harmonic_mean = (2 * bsa_al * bsa_bl) / (bsa_al + bsa_bl + 1e-8)

        pocket_bsa [resid] = dict (al=bsa_al, bl=bsa_bl, bsa_abl=bsa_abl, hm=bsa_harmonic_mean, pocket=resid)
    
    sorted_pocket_bsa = dict(
    sorted(pocket_bsa.items(), key=lambda item: item[1]['hm'], reverse=True))
    max_index = list(sorted_pocket_bsa.keys())[0]
    
    return sorted_pocket_bsa [max_index]






def main():
    
    parser = argparse.ArgumentParser(description="Obtain ")
    parser.add_argument('-f', dest='xtc_file', type=Path, default='mol.xtc', help='input xtc file')
    parser.add_argument('-s', dest='tpr_file', type=Path, default='npt.tpr', help='input tpr file')
    parser.add_argument('-n', dest='n_frames', type=int, default=50, help='number of frames used for Fpocket analysis to find top pocket')
    parser.add_argument('-b', dest='begin_time', type=int, default=0, help='begin time (ps) for analysis')
    parser.add_argument('-e', dest='end_time', type=int, default=None, help='end time (ps) for analysis')
    parser.add_argument('-on', dest='number_of_confs', type=int, default=3, help='Number of protein conformation to be used for pharmacophore search')
    parser.add_argument('-out_path', dest='out_path', type=Path, default=None, help='path where the pharmacophore folder will be saved')
    parser.add_argument('-keep', dest='keep_files', type=int, default=1, help='if set to zero removes the pocket analysis files that are no longer needed.')

    args = parser.parse_args()

    tpr_file = args.tpr_file 
    xtc_file = args.xtc_file
    pocket_path = xtc_file.parent/"pockets"
    pocket_path.mkdir(parents=True, exist_ok=True)
    u    =  mda.Universe (tpr_file, xtc_file)
    
    n_frames    = args.n_frames
    begin_time  = args.begin_time
    end_time   = args.end_time
    number_of_confs = args.number_of_confs
    out_path        =args.out_path
###############################################



    dt          = u.trajectory.dt
    begin_frame = int(begin_time/dt)
    if not end_time:
        skip        = int ((u.trajectory.totaltime - begin_time)/(dt*n_frames))
        end_frame   = None
    else:
        skip        = int ((end_time - begin_time)/(dt*n_frames))
        end_frame   = int(end_time/dt)
    ################
    #Pocket dirtion
    conf_dict = {}
     
    
    for conf, ts in tqdm(enumerate (u.trajectory [begin_frame:end_frame:skip]),
                        total=n_frames,
                        desc=f"Extracting {n_frames} frames & analysing pockets"):
    
        # create a direcotry for the given conf
        conf_path = pocket_path/f"conf{conf+1}"
        conf_path.mkdir(parents=True, exist_ok=True)
        # write the full coordinate file
        u.atoms.write (Path(conf_path/f"conf{conf+1}.gro"))
    
        # just the proteins, iwth proper chain Labels is important
        u.atoms.fragments [0].chainIDs = "A"
        u.atoms.fragments [1].chainIDs = "B"
    
        u.select_atoms ("protein").write ( (Path(conf_path/"protein.pdb")), bonds=None)
    
        # Run Fpocket
        subprocess.run(
            ["fpocket", "-f", "protein.pdb"],
            cwd=conf_path,
            stdout=subprocess.DEVNULL,# sets working directory
            check=True            # raises error if fpocket fails
        )
    
        pocket_file = conf_path / "protein_out/protein_out.pdb"
        dict_ =  load_protein (pocket_file)
     
        dict_["conf"] = conf + 1
        dict_ ["pocket_file"] = pocket_file
        dict_ ["conf_file"]   = Path(conf_path/f"conf{conf+1}.gro")
        conf_dict [conf+1] = dict_
        
    
    sorted_pocket_bsa = dict(sorted(conf_dict.items(), key=lambda item: item[1]['hm'], reverse=True))
    #
     
        
    # create directories
    if not out_path:
        pharma_path = xtc_file.parent/"pharmacophore-search"
        pharma_path.mkdir (parents=True, exist_ok=True)
    else:
        pharma_path = out_path/"pharmacophore-search"
        pharma_path.mkdir (parents=True, exist_ok=True)

        
    with open(pharma_path/"pocket_scores.json", "w") as f:
        json.dump(sorted_pocket_bsa, f, indent=4,
                 default=lambda x: x.item() if isinstance(x, (np.integer, np.floating, np.int64, np.float64)) else str(x))  # indent=4 makes it pretty-printe  
        
    
    count = 0
    for conf, value in sorted_pocket_bsa.items():
        
        count += 1
        if count > number_of_confs: break
        
        conf_path = pharma_path / f"conf{count}"
        conf_path.mkdir(parents=True, exist_ok=True)
        with open (conf_path / f"PocketID_is_{value ['pocket']}", "w+") as f: 
            f.write (f"Pocket ID for analysis: {value ['pocket']}\n bsa_al={value['al']} \n bsa_bl={value['bl']} \n bsa_abl={value['bsa_abl']}\n bsa_hm={value['hm']}")
          # Source file
         
     
        # Copy the file into the destination directory
        shutil.copy(value ["pocket_file"], conf_path)
        shutil.copy(value ["conf_file"], conf_path)
        shutil.copy (tpr_file, conf_path) 

    if not args.keep_files:  shutil.rmtree(pocket_path)

        
if __name__ == "__main__":
    main()
