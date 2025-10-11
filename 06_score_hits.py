
import os
import pickle
import warnings
import numpy as np
from tqdm import tqdm
import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element
import argparse
 

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

def parse_single_sdf_file(file_path):
    """
    Parse SDF file, returning a dictionary:
    {ZINC_ID: {"xyz": [...], "atom_radii": [...], "rmsd": ...}}
    """
    zinc_dict = {}
    with open(file_path, "r") as f:
        data, atom_radii = [], []
        current_rmsd = None
        zinc_id = None
        read_rmsd = False

        for line in f:
            line = line.strip()
            if not line:
                continue

            if "ZINC" in line:
                zinc_id = line
                if zinc_id not in zinc_dict:
                    zinc_dict[zinc_id] = {"xyz": [], "atom_radii": [], "rmsd": []}

            elif "$$$$" in line:
                if zinc_id is not None:
                    zinc_dict[zinc_id]["xyz"].append(np.array(data, dtype=np.float32))
                    zinc_dict[zinc_id]["atom_radii"] = atom_radii
                    zinc_dict[zinc_id]["rmsd"].append(current_rmsd)
                data, atom_radii = [], []
                current_rmsd = None
                read_rmsd = False

            elif len(line.split()) >= 4:  # flexible parsing
                parts = line.split()
                try:
                    xyz = [float(x) for x in parts[:3]]
                    name = parts[3]
                    data.append(xyz)
                    try:
                        atom_radii.append(_ATOMIC_RADII[element_symbol_string(name)])
                    except KeyError:
                        atom_radii.append(1.5)  # fallback radius
                except ValueError:
                    continue

            elif "rmsd" in line.lower():
                read_rmsd = True

            elif read_rmsd:
                try:
                    current_rmsd = float(line)
                except ValueError:
                    current_rmsd = None
                read_rmsd = False

    return zinc_dict

def load_protein(tpr_file, gro_file):
    """Load protein and chain data, compute SASA and atom radii."""
    u = mda.Universe(tpr_file, gro_file)
    protein = u.select_atoms("protein and not name H*")
    chainA  = u.atoms.fragments [0].select_atoms ("all and not name H*")
    chainB  = u.atoms.fragments [1].select_atoms ("all and not name H*")
 
    def compute_radii_sasa(atoms):
        symbols = get_element_symbols(atoms.names)
        radii = [ _ATOMIC_RADII.get(sym) for sym in symbols ]
        sasa = Sasa_calc(atoms.positions, radii)
        return radii, sasa

    radii_protein, sasa_protein = compute_radii_sasa(protein)
    radii_chainA, sasa_chainA   = compute_radii_sasa(chainA)
    radii_chainB, sasa_chainB   = compute_radii_sasa(chainB)

    return {
        "sasa": {"protein": sasa_protein, "chainA": sasa_chainA, "chainB": sasa_chainB},
        "radii": {"protein": radii_protein, "chainA": radii_chainA, "chainB": radii_chainB},
        "xyz": {"protein": protein.positions, "chainA": chainA.positions, "chainB": chainB.positions}
    }

def get_bsa(dict_protein, dict_ligand, ligindex=0):
    """Compute buried surface areas for ligand with protein chains."""
    ligand_xyz = dict_ligand['xyz'][ligindex]
    ligand_radii = dict_ligand['atom_radii']

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

    return bsa_al, bsa_bl, bsa_abl

def assign_properties(rmsd, bsa_al, bsa_bl, bsa_abl, sdf_file_path, tpr_file, gro_file):
    bsa_harmonic_mean = (2 * bsa_al * bsa_bl) / (bsa_al + bsa_bl + 1e-8)
    return {
        "bsa_a": bsa_al,
        "bsa_b": bsa_bl,
        "bsa_ab": bsa_abl,
        "rmsd": rmsd,
        "bsa_harm": bsa_harmonic_mean,
        "sdf_name": sdf_file_path,
        "tpr": tpr_file,
        "gro": gro_file
    }

def get_all_ligand_scores(sdf_path, tpr_file, gro_file):
    
    sdf_files = [f for f in os.listdir(sdf_path) if f.endswith(".sdf")]
    graph_list = {name.rsplit(".",1)[0].split("_")[-1] for name in sdf_files}

    dict_protein = load_protein(tpr_file, gro_file)
    bsa_score_all_ligands = {graph:{} for graph in graph_list}

    for sdf_file in sdf_files:
        sdf_file_path = os.path.join(sdf_path, sdf_file)
        dict_ligand_all = parse_single_sdf_file(sdf_file_path)
        graph_id = sdf_file.rsplit(".",1)[0].split("_")[-1]

        for zinc_id, dict_ligand in tqdm(dict_ligand_all.items()):
            bsa_al, bsa_bl, bsa_abl = get_bsa(dict_protein, dict_ligand)
            kwargs = dict(
                rmsd=dict_ligand['rmsd'][0] if dict_ligand['rmsd'] else None,
                bsa_al=bsa_al,
                bsa_bl=bsa_bl,
                bsa_abl=bsa_abl,
                sdf_file_path=sdf_file_path,
                tpr_file=tpr_file,
                gro_file=gro_file
            )
            bsa_score_all_ligands[graph_id][zinc_id] = assign_properties(**kwargs)

    return bsa_score_all_ligands
 

def main():
    """Main execution function: calculate BSA scores for Pharmer hits."""

    parser = argparse.ArgumentParser(description="Calculate buried surface area (BSA) scores for Pharmer hits")

    parser.add_argument( "--sdf_path", "-sdf", type=str, required=True, help="Full path to SDF output directory")
    parser.add_argument("--tpr_file", "-s", type=str, required=True, help="Full path + filename of TPR file")
    parser.add_argument( "--gro_file", "-c", type=str, required=True, help="Full path + filename of GRO file")
    parser.add_argument( "--out_pkl", "-o", type=str, default=None, help="Output pickle filename")
    
    # Parse the arguments
    args = parser.parse_args()

    # Compute BSA scores
    all_scores = get_all_ligand_scores(args.sdf_path, args.tpr_file, args.gro_file)

    # Save results to pickle
    #
    if not args.out_pkl:
        args.out_pkl =  os.path.basename(os.path.normpath(args.sdf_path)) + "_bsa_score.pkl"
        
    with open(args.out_pkl, "wb") as out_file:
        pickle.dump(all_scores, out_file, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"BSA scores saved to {args.out_pkl}")


if __name__ == "__main__":
    main()

