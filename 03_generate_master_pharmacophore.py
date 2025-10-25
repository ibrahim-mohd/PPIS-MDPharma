import pickle
import json
import argparse
import numpy as np
import MDAnalysis as mda
import warnings

warnings.filterwarnings("ignore")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate pharmacophore JSON file from PKL files"
    )

    # Input structure files
    parser.add_argument('-s', dest='tpr_file', type=str, default='protein_noh.pdb', help='TPR file')
    parser.add_argument('-c', dest='coord_file', type=str, default='protein.gro', help='PDB/GRO file')

    # Input PKL object files
    parser.add_argument('-p', dest='features_pkl', type=str, default='dgsolv_sasa.pkl',
                        help='Pickle file with solvation energy information')

    # Thresholds
    parser.add_argument('-dG_th', dest='G_solv_threshold', type=float, default=0.40,
                        help='Threshold for solvation free energy')
    parser.add_argument('-acceptor_th', dest='acceptor_threshold', type=float, default=40.0,
                        help='Threshold for hydrogen acceptor')
    parser.add_argument('-donor_th', dest='donor_threshold', type=float, default=40.0,
                        help='Threshold for hydrogen donor')
    parser.add_argument('-ion_th', dest='threshold_ion_difference', type=float, default=30,
                        help='Threshold for |N_anion - N_cation|')

    # Output file
    parser.add_argument('-o', dest='out_file', type=str, default='master_pharma.json',
                        help='Output JSON file')

    return parser.parse_args()


def distance(r1: np.ndarray, r2: np.ndarray) -> float:
    """Calculate Euclidean distance between two points."""
    return np.linalg.norm(r1 - r2)


def get_pharmacophore_position(center: np.ndarray, p_geometric_center: np.ndarray, radius: float = 3.0) -> np.ndarray:
    """
    Move a point towards the geometric center by a given radius.
    """
    direction = p_geometric_center - center
    direction /= np.linalg.norm(direction)
    return center + radius * direction


def pharmer_element_template(site_type: str,
                             position: np.ndarray,
                             score: dict = None,
                             chainID: str = "",
                             requirement: str = "required",
                             size: int = 1,
                             radius: float = 1,
                             minsize: str = "",
                             maxsize: str = "",
                             enabled: str = "true",
                             svector: str = "null",
                             vector: str = "null",
                             vector_on: int = 0,
                             gen_score: bool = True) -> dict:
    """Return a dictionary representing a pharmacophore element."""
    valid_types = ["Aromatic", "HydrogenDonor", "HydrogenAcceptor",
                   "Hydrophobic", "PositiveIon", "ExclusionSphere", "NegativeIon"]

    if site_type not in valid_types:
        raise ValueError(f"Invalid site type: {site_type}. Must be one of {valid_types}")

    element = {
        "name": site_type,
        "radius": radius,
        "requirement": requirement,
        "size": size,
        "x": float(position[0]),
        "y": float(position[1]),
        "z": float(position[2]),
        "enabled": enabled,
        "svector": svector,
        "vector": vector,
        "minsize": minsize,
        "maxsize": maxsize,
        "vector_on": vector_on
    }

    if gen_score and score is not None:
        element.update({
            "score": score,
            "chainID": chainID
        })

    return element


def assign_chainID(atom_indices: list, chainA_indices: np.ndarray, chainB_indices: np.ndarray) -> str:
    """Assign chain ID based on atom index membership."""
    idx = int(atom_indices)
    if idx in chainA_indices:
        return "A"
    elif idx in chainB_indices:
        return "B"
    return ""


def main():
    args = parse_args()

    # Load protein structure
    u = mda.Universe(args.tpr_file, args.coord_file)
    chainA, chainB = u.atoms.fragments[0], u.atoms.fragments[1]
    exclusion_sites = u.select_atoms("protein and not name H* 1H* 2H*")

    # Load features
    with open(args.features_pkl, 'rb') as f:
        all_features = pickle.load(f)

    hbond = all_features["hbonds"]
    anion_cation = all_features["ions"]
    dgsol_sasa = all_features["sasa_dG"]

    # Assign chain IDs for ions
    for ion_type in ['anion', 'cation']:
        for key in anion_cation[ion_type].keys():
            anion_cation[ion_type][key]['chainID'] = assign_chainID(key, chainA.indices, chainB.indices)

    # Assign chain IDs for hydrophobic sites
    for key in dgsol_sasa.keys():
        resid = int(key)
        if resid in chainA.residues.resids:
            dgsol_sasa[key]['chainID'] = "A"
        elif resid in chainB.residues.resids:
            dgsol_sasa[key]['chainID'] = "B"

    # Assign chain IDs for H-bonds
    for hb_type in ['acceptor', 'donor']:
        for key in hbond[hb_type].keys():
            hbond[hb_type][key]['chainID'] = assign_chainID(key, chainA.indices, chainB.indices)

    # Thresholds
    G_solv_threshold = args.G_solv_threshold
    acceptor_threshold = args.donor_threshold # intentionally reversed, i.e donor threshold used as acceptor threshold, since protein donor corresponds to an acceptor pharmacophore
    donor_threshold = args.acceptor_threshold # same as above. Also to keep compatible with feature selection file where we show the pharmacophroe information
    threshold_ion_difference = args.threshold_ion_difference

    # Limits
    # if more than these limits we  use the top 10 scores
    nmax_hydrophobic = 10
    nmax_donor = 10
    nmax_acceptor = 10
    nmax_cation = 10
    nmax_anion = 10

    # Identify pharmacophore sites
    # Hydrophobic
    possible_hydrophobic_sites = [key for key, val in dgsol_sasa.items() if val['dG'][0] > G_solv_threshold]
    if len(possible_hydrophobic_sites) > nmax_hydrophobic:
        dG_values = {key: dgsol_sasa[key]['dG'][0] for key in possible_hydrophobic_sites}
        possible_hydrophobic_sites = sorted(dG_values, key=dG_values.get, reverse=True)[:nmax_hydrophobic]

    # Donors & acceptors
    possible_acceptors = [k for k, v in hbond['acceptor'].items() if v['count'] > acceptor_threshold]
    if len(possible_acceptors) > nmax_acceptor:
        possible_acceptors = sorted(possible_acceptors, key=lambda x: hbond['acceptor'][x]['count'], reverse=True)[:nmax_acceptor]

    possible_donors = [k for k, v in hbond['donor'].items() if v['count'] > donor_threshold]
    if len(possible_donors) > nmax_donor:
        possible_donors = sorted(possible_donors, key=lambda x: hbond['donor'][x]['count'], reverse=True)[:nmax_donor]

    # Cation/anion differences
    cation_counts = np.array([v['count'] for v in anion_cation['cation'].values()])
    anion_counts = np.array([v['count'] for v in anion_cation['anion'].values()])
    cation_anion_diff = {k: c - a for k, (c, a) in zip(anion_cation['cation'].keys(), zip(cation_counts, anion_counts))}

    possible_cation_sites = [k for k, v in cation_anion_diff.items() if v >= threshold_ion_difference][:nmax_cation]
    possible_anion_sites = [k for k, v in cation_anion_diff.items() if v <= -threshold_ion_difference][:nmax_anion]

    # Geometric center
    pharmacophore_positions = []
    for site_list in [possible_donors, possible_acceptors, possible_cation_sites, possible_anion_sites]:
        if site_list:
            atoms = u.select_atoms(f"protein and index {' '.join(site_list)}")
            pharmacophore_positions.extend(atoms.positions)

    for resid in possible_hydrophobic_sites:
        atoms = u.select_atoms(f"resid {resid}")
        pharmacophore_positions.append(atoms.center_of_mass())

    X = np.array(pharmacophore_positions)
    p_geometric_center = np.mean(X, axis=0)

    # Generate pharmacophore elements
    data = {"points": []}

    def add_elements(possible_sites, site_type_protein, site_type_pharma, feature_dict, u, norm_type="max"):
        if not possible_sites:
            return
        norm_value = max([feature_dict[x]['count'] for x in possible_sites])
        atoms = u.select_atoms(f"protein and index {' '.join(possible_sites)}")
        for site, pos in zip(possible_sites, atoms.positions):
            pos = get_pharmacophore_position(pos, p_geometric_center)
            score = {'score': feature_dict[site]['count'], 'normed_score': feature_dict[site]['count'] / norm_value}
            element = pharmer_element_template(site_type_pharma, pos, score, chainID=feature_dict[site]['chainID'])
            data["points"].append(element)

    # Acceptors -> HydrogenDonor
    add_elements(possible_acceptors, 'acceptor', 'HydrogenDonor', hbond['acceptor'], u)
    # Donors -> HydrogenAcceptor
    add_elements(possible_donors, 'donor', 'HydrogenAcceptor', hbond['donor'], u)
    # Cation
    add_elements(possible_cation_sites, 'cation', 'PositiveIon', anion_cation['cation'], u)
    # Anion
    add_elements(possible_anion_sites, 'anion', 'NegativeIon', anion_cation['anion'], u)

    # Hydrophobic
    if possible_hydrophobic_sites:
        hydrophobic_norm = max([dgsol_sasa[x]['dG'][0] for x in possible_hydrophobic_sites])
        for resid in possible_hydrophobic_sites:
            atoms = u.select_atoms(f"resid {resid} and not backbone and not name H*")
            pos = get_pharmacophore_position(atoms.center_of_mass(), p_geometric_center)
            distances = [np.linalg.norm(p - pos) for p in atoms.positions]
            radius = np.average(distances)
            score = {'score': dgsol_sasa[resid]['dG'][0], 'normed_score': dgsol_sasa[resid]['dG'][0]/hydrophobic_norm}
            size = 6 if np.unique(atoms.resnames)[0] in ['TRP', 'TYR', 'PHE'] else 1
            element = pharmer_element_template("Hydrophobic", pos, score, dgsol_sasa[resid]['chainID'], size=size, radius=radius)
            data["points"].append(element)

    # Exclusion spheres
    for pos in exclusion_sites.positions:
        element = pharmer_element_template("ExclusionSphere", pos, score={'score':1}, chainID='x')
        data["points"].append(element)

    # Write JSON
    with open(args.out_file, "w") as outfile:
        json.dump(data, outfile, indent=4)


if __name__ == "__main__":
    main()

