# Written by Mohd Ibrahim
# Technical University of Munich
import os
import sys
import json
import pickle
import warnings
import argparse
import numpy as np
import MDAnalysis as mda
from datetime import datetime
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate pharmacophore JSON file from PKL files"
    )

    parser.add_argument('-s', dest='tpr_file', type=str, default='protein_noh.pdb', help='TPR file')
    parser.add_argument('-c', dest='coord_file', type=str, default='protein.gro', help='PDB/GRO file')
    parser.add_argument('-p', dest='features_pkl', type=str, default='./combined_analysis.pkl',
                        help='Pickle file with solvation energy information')
    parser.add_argument('-dG_th', dest='G_solv_threshold', type=float, default=0.25, help='Threshold for solvation free energy')
    parser.add_argument('-acceptor_th', dest='acceptor_threshold', type=float, default=35.0, help='Threshold for hydrogen acceptor')
    parser.add_argument('-donor_th', dest='donor_threshold', type=float, default=35.0, help='Threshold for hydrogen donor')
    parser.add_argument('-ion_th', dest='threshold_ion_difference', type=float, default=25, help='Threshold for |N_anion - N_cation|')
 
    parser.add_argument('-o', dest='out_file', type=str, default='master_pharma.json', help='Output JSON file')

     
    # Ion options
    parser.add_argument('-nname', dest='nname', type=str, default='Cl-', help='anion name')
    parser.add_argument('-pname', dest='pname', type=str, default='Na+', help='cation name')
    #   # Solvent residue and atom names
    parser.add_argument('-hsol_name', dest='hsol_name', type=str, default='HW1 HW2', help='atom name of water hydrogens')
    parser.add_argument('-osol_name', dest='osol_name', type=str, default='OW', help='atomname of water oxygens')
    parser.add_argument('-sol_resname', dest='sol_resname', type=str, default='SOL', help='resname of water molecules')

    parser.add_argument('-hbond_direction', dest='hbond_direction', type=int, default=1, help='whether to include the h-bond directions, use 0 for no direction default is 1')
    parser.add_argument('-ignore_nowater_hbond', dest='ignore_nowater_hbond', type=int, default=0, help='whether to ignore h-bond bond assignment where no water is found in h-bond geometry close to the site, use 1 to ignore default is 0 i.e do not ignore')


    return parser.parse_args()

def distance(r1, r2):
    return np.linalg.norm(r1 - r2)

def merge_identical_sites(data):
    """
    Merge HydrogenDonor and HydrogenAcceptor points that share identical
    x, y, z, and vector coordinates by summing their scores.
    """    
    from copy import deepcopy
    
    data = deepcopy(data)
    points = data.get("points", [])

    seen = {}
    result = []

    for p in points:
        name = p.get("name")
        
        # Only deduplicate H-bond features and ionc istes as only those can have identical coordinates
        if name in {"HydrogenDonor", "HydrogenAcceptor", "NegativeIon", "PositiveIon"}:
            # Get vector value
            vector_val = p.get("vector")
            
            # Create a hashable key - ALWAYS ensure vector is converted to tuple if it's a list
            if isinstance(vector_val, list):
                # Convert list to tuple for hashability
                vector_key = tuple(vector_val)
            elif isinstance(vector_val, tuple):
                # Tuples are already hashable
                vector_key = vector_val
            else:
                # For strings, None, or other types
                vector_key = vector_val
            
            key = (
                name,
                p.get("x"),
                p.get("y"),
                p.get("z"),
                vector_key
            )

            # Try to use the key - if it's not hashable, skip deduplication
            try:
                if key in seen:
                    # Merge scores
                    seen[key]["score"]["score"] += p["score"].get("score", 0)

                    if "normed_score" in p["score"]:
                        seen[key]["score"]["normed_score"] += p["score"].get("normed_score", 0)
                else:
                    seen[key] = p
                    result.append(p)
            except TypeError:
                # Key is not hashable (e.g., vector is still a list somehow)
                # Skip deduplication for this point
                result.append(p)
        else:
            # Keep everything else as it is
            result.append(p)

    data["points"] = result
    return data
    
def print_summary(data, path):
    output_file = os.path.abspath(path) + "/summary_master_pharma.dat"
    lines = []
    total_features = len(data['points'])
    lines.append(f"\nTotal number of features: {total_features}")
    print(lines[-1])
    temp_dict = {}
    for json_element in data['points']:
        key = json_element['name']
        if key not in temp_dict:
            temp_dict[key] = dict(chainA=0, chainB=0, total=0)
        if json_element['chainID'] == 'A':
            temp_dict[key]['chainA'] += 1
        else:
            temp_dict[key]['chainB'] += 1
        temp_dict[key]['total'] += 1
    lines.append("Feature summary:")
    lines.append("-" * 50)
    print(lines[-2])
    print(lines[-1])
    for key, value in temp_dict.items():
        line = f"{key:20s} | chainA={value['chainA']:3d} | chainB={value['chainB']:3d} | Total={value['total']:3d}"
        lines.append(line)
        print(line)
    with open(output_file, "w") as f:
        f.write(f"{datetime.now()}: {' '.join(sys.argv)}\n")
        f.write("\n".join(lines))


def get_pharmacophore_position(center: np.ndarray, p_geometric_center: np.ndarray, radius: float = 4.0) -> list:
    direction = p_geometric_center - center
    direction /= np.linalg.norm(direction)
    pos = center + radius * direction
    return [float(x) for x in pos]


def assign_aromatic_position(atom_group, p_geometric_center, max_attempts=10):
    assert len(atom_group) >= 3, "Atom group must contain at least 3 atoms"
    geo_center = atom_group.center_of_geometry()
    for attempt in range(max_attempts):
        i, j, k = np.random.choice(len(atom_group), size=3, replace=False)
        OA, OB, OC = atom_group.positions[i], atom_group.positions[j], atom_group.positions[k]
        AB = OB - OA
        BC = OC - OB
        AB /= np.linalg.norm(AB)
        BC /= np.linalg.norm(BC)
        normal_vector = np.cross(AB, BC)
        norm = np.linalg.norm(normal_vector)
        if norm > 1e-6:
            normal_vector /= norm
            break
    else:
        raise ValueError("Failed to select non-collinear atoms after multiple attempts")
    position1 = geo_center + 3 * normal_vector
    position2 = geo_center - 3 * normal_vector
    position = position1 if distance(position1, p_geometric_center) < distance(position2, p_geometric_center) else position2
    return [float(x) for x in position], [float(x) for x in normal_vector]

def get_hbond_acceptors(site_index, u, d_cutoff=4, angle_cutoff=150, min_angle=100,
                        osol_name="OW", sol_resname="SOL", hsol_name="HW1 HW2"):
    hbonds_acceptors = HydrogenBondAnalysis(
        universe=u,
        donors_sel=f"name {osol_name} and resname {sol_resname}",
        hydrogens_sel=f"name {hsol_name} and resname {sol_resname}",
        acceptors_sel=f"index {site_index}",
        d_a_cutoff=d_cutoff, d_h_a_angle_cutoff=angle_cutoff,
        update_selections=False)
    hbonds_acceptors.run()
    hb_array = hbonds_acceptors.results["hbonds"]
    if hb_array.size == 0:
        if angle_cutoff > min_angle:
            return get_hbond_acceptors(site_index, u=u, d_cutoff=d_cutoff,
                                       angle_cutoff=angle_cutoff-5, min_angle=min_angle,
                                       osol_name=osol_name, sol_resname=sol_resname, hsol_name=hsol_name)
        else:
            return None, None
    water_donor_oxygen = u.select_atoms(f"index {int(hb_array[0][1])}")
    hbond_acceptor = u.select_atoms(f"index {site_index}")
    direction = hbond_acceptor.positions[0] - water_donor_oxygen.positions[0]
    direction /= np.linalg.norm(direction)
    return [float(x) for x in water_donor_oxygen.positions[0]], [float(x) for x in direction]


def get_hbond_donors(site_index, u, d_cutoff=4, angle_cutoff=150, min_angle=100,
                     osol_name="OW", sol_resname="SOL", hsol_name="HW1 HW2"):
    hbonds_donors = HydrogenBondAnalysis(
        universe=u,
        donors_sel=f"index {site_index}",
        hydrogens_sel=f"name H* and around 1.5 index {site_index}",
        acceptors_sel=f"name {osol_name} and resname {sol_resname}",
        d_a_cutoff=d_cutoff, d_h_a_angle_cutoff=angle_cutoff,
        update_selections=False)
    hbonds_donors.run()
    hb_array = hbonds_donors.results["hbonds"]
    if hb_array.size == 0:
        if angle_cutoff > min_angle:
            return get_hbond_donors(site_index, u=u, d_cutoff=d_cutoff,
                                    angle_cutoff=angle_cutoff-5, min_angle=min_angle,
                                    osol_name=osol_name, sol_resname=sol_resname, hsol_name=hsol_name)
        else:
            return None, None
    water_acceptor_oxygen = u.select_atoms(f"index {int(hb_array[0][3])}")
    hbond_donor = u.select_atoms(f"index {site_index}")
    direction = hbond_donor.positions[0] - water_acceptor_oxygen.positions[0]
    direction /= np.linalg.norm(direction)
    return [float(x) for x in water_acceptor_oxygen.positions[0]], [float(x) for x in direction]

def get_positive_negative_position(u, site, site_type_protein, max_cutoff=5, nname="Cl-",pname="Na+"):
    # For an accurate psotion assignemtn we see if there are any boudn in in the current frame
    # This is very useful becasue sometimes many sites are coordinated to teh same ion in which case
    # this funtion assigns identical positions which are then later merged into single site with sum of scores
    # ion name
    ion_name = pname if site_type_protein == "cation" else nname
    ion_site = u.select_atoms (f"index {site}")
    ion_selection  = u.select_atoms (f"name {ion_name}")
    ion_position = None # if we do not find any posotin we return None
    
    for cutoff in np.arange (2,max_cutoff+0.1,0.1):
        pairs, _ = mda.lib.distances.capped_distance(ion_site, ion_selection, max_cutoff=cutoff, box=u.dimensions)
        if len (pairs) !=0:
            ion_position = list(ion_selection.positions [pairs[0][1]])
            break   
    return  ion_position


def pharmer_element_template(site_type: str, position: list, score: dict = None, chainID: str = "",
                             requirement: str = "required", size: int = 1, radius: float = 1,
                             enabled: str = "true", vector="null", gen_score: bool = True) -> dict:
    valid_types = ["Aromatic", "HydrogenDonor", "HydrogenAcceptor",
                   "Hydrophobic", "PositiveIon", "ExclusionSphere", "NegativeIon"]
    if site_type not in valid_types:
        raise ValueError(f"Invalid site type: {site_type}. Must be one of {valid_types}")
    element = {
        "name": site_type,
        "radius": float(radius),
        "requirement": requirement,
        "size": size,
        "x": float(position[0]),
        "y": float(position[1]),
        "z": float(position[2]),
        "enabled": enabled,
        "vector": vector,
    }
    if gen_score and score is not None:
        element.update({
            "score": score,
            "chainID": chainID
        })
    return element


def assign_chainID(atom_indices: list, chainA_indices: np.ndarray, chainB_indices: np.ndarray) -> str:
    idx = int(atom_indices)
    if idx in chainA_indices:
        return "A"
    elif idx in chainB_indices:
        return "B"
    return ""


def main():
    args = parse_args()
    u = mda.Universe(args.tpr_file, args.coord_file)
    chainA, chainB = u.atoms.fragments[0], u.atoms.fragments[1]
    exclusion_sites = u.select_atoms("protein and not name H* 1H* 2H*")

    with open(args.features_pkl, 'rb') as f:
        all_features = pickle.load(f)

    hbond = all_features["hbonds"]
    anion_cation = all_features["ions"]
    dgsol_sasa = all_features["sasa_dG"]

    # Assign chain IDs
    for ion_type in ['anion', 'cation']:
        for key in anion_cation[ion_type].keys():
            anion_cation[ion_type][key]['chainID'] = assign_chainID(key, chainA.indices, chainB.indices)

    for key in dgsol_sasa.keys():
        resid = int(key)
        if resid in chainA.residues.resids:
            dgsol_sasa[key]['chainID'] = "A"
        elif resid in chainB.residues.resids:
            dgsol_sasa[key]['chainID'] = "B"

    for hb_type in ['acceptor', 'donor']:
        for key in hbond[hb_type].keys():
            hbond[hb_type][key]['chainID'] = assign_chainID(key, chainA.indices, chainB.indices)

    # Thresholds
    G_solv_threshold = args.G_solv_threshold
    acceptor_threshold = args.donor_threshold
    donor_threshold = args.acceptor_threshold
    threshold_ion_difference = args.threshold_ion_difference

    nmax_hydrophobic = 10
    nmax_donor = 10
    nmax_acceptor = 10
    nmax_cation = 10
    nmax_anion = 10

    possible_hydrophobic_sites = [key for key, val in dgsol_sasa.items() if val['dG'][0] > G_solv_threshold]
    if len(possible_hydrophobic_sites) > nmax_hydrophobic:
        dG_values = {key: dgsol_sasa[key]['dG'][0] for key in possible_hydrophobic_sites}
        possible_hydrophobic_sites = sorted(dG_values, key=dG_values.get, reverse=True)[:nmax_hydrophobic]

    possible_acceptors = [k for k, v in hbond['acceptor'].items() if v['count'] > acceptor_threshold]
    if len(possible_acceptors) > nmax_acceptor:
        possible_acceptors = sorted(possible_acceptors, key=lambda x: hbond['acceptor'][x]['count'], reverse=True)[:nmax_acceptor]

    possible_donors = [k for k, v in hbond['donor'].items() if v['count'] > donor_threshold]
    if len(possible_donors) > nmax_donor:
        possible_donors = sorted(possible_donors, key=lambda x: hbond['donor'][x]['count'], reverse=True)[:nmax_donor]

    cation_counts = np.array([v['count'] for v in anion_cation['cation'].values()])
    anion_counts = np.array([v['count'] for v in anion_cation['anion'].values()])
    cation_anion_diff = {k: c - a for k, (c, a) in zip(anion_cation['cation'].keys(), zip(cation_counts, anion_counts))}

    possible_cation_sites = [k for k, v in cation_anion_diff.items() if v >= threshold_ion_difference][:nmax_cation]
    possible_anion_sites = [k for k, v in cation_anion_diff.items() if v <= -threshold_ion_difference][:nmax_anion]

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

    data = {"points": []}

    def add_elements(possible_sites, site_type_protein, site_type_pharma, feature_dict, u, norm_type="max"):
        if not possible_sites:
            return
        norm_value = max([feature_dict[x]['count'] for x in possible_sites])
        atoms = u.select_atoms(f"protein and index {' '.join(possible_sites)}")
        for site, position in zip(possible_sites, atoms.positions):
            pos = get_pharmacophore_position(position, p_geometric_center)
            pos = [float(x) for x in pos]
            vector = "null"

            # we position positive and negative sites first by checking if there are any bound ions in the 
            # current frame. If not we simply used the above translated one
            if site_type_protein in ['cation', 'anion']:
                ion_position = get_positive_negative_position (u, site, site_type_protein, nname=args.nname,pname=args.pname)
                if ion_position is not None: pos = ion_position
                
            if site_type_protein == 'donor':
                pos_hb, direction = get_hbond_donors(site, u)
                if pos_hb is None:
                    if args.ignore_nowater_hbond: 
                        continue
                        # Usually these donors are in very hidden places 
                    else:
                        pos_hb, direction = pos, None
                        print(f"Acceptor:{site} no water found. We resort to a simple position assignment i.e as simple translation towards pocket center by 3 Angstroms. To ignorre such site set the ignore_nowater_hbond flag to 0")
                        
                         
                pos = [float(x) for x in pos_hb]
                if args.hbond_direction and direction is not None:
                    vector = [{k: float(v) for k, v in zip(["x", "y", "z"], direction)}]
                else: vector= "null"
            if site_type_protein == 'acceptor':
                pos_hb, direction = get_hbond_acceptors(site, u,osol_name=args.osol_name, sol_resname=args.sol_resname, hsol_name=args.hsol_name)
                
                if pos_hb is None:
                    if args.ignore_nowater_hbond: 
                        continue
                        # Usually these donors are in very hidden places 
                    else:
                        pos_hb, direction = pos, None
                        print(f"Donor:{site} no water found. We resort to a simple position assignment i.e as simple translation towards pocket center by 3 Angstroms. To ignorre such site set the ignore_nowater_hbond flag to 0")
                pos = [float(x) for x in pos_hb]
                if args.hbond_direction and direction is not None:
                    vector = [{k: float(v) for k, v in zip(["x", "y", "z"], direction)}]
                else: vector = "null"
            score = {'score': feature_dict[site]['count'],
                     'normed_score': feature_dict[site]['count'] / norm_value}
            element = pharmer_element_template(site_type_pharma, pos,
                                               score,
                                               chainID=feature_dict[site]['chainID'],
                                               vector=vector)
            data["points"].append(element)

    add_elements(possible_acceptors, 'acceptor', 'HydrogenDonor', hbond['acceptor'], u)
    add_elements(possible_donors, 'donor', 'HydrogenAcceptor', hbond['donor'], u)
    add_elements(possible_cation_sites, 'cation', 'PositiveIon', anion_cation['cation'], u)
    add_elements(possible_anion_sites, 'anion', 'NegativeIon', anion_cation['anion'], u)

    # Hydrophobic / Aromatic
    if possible_hydrophobic_sites:
        hydrophobic_norm = max([dgsol_sasa[x]['dG'][0] for x in possible_hydrophobic_sites])
        for resid in possible_hydrophobic_sites:
            atoms = u.select_atoms(f"protein and resid {resid} and not backbone and not name H*")
            if not atoms:
                continue
            pos = get_pharmacophore_position(atoms.center_of_mass(), p_geometric_center)
            pos = [float(x) for x in pos]
            distances = [np.linalg.norm(p - pos) for p in atoms.positions]
            radius = float(np.average(distances))
            score = {'score': dgsol_sasa[resid]['dG'][0],
                     'normed_score': dgsol_sasa[resid]['dG'][0] / hydrophobic_norm}
            vector = "null"
            site_type = "Hydrophobic"
            size = 1
            if np.unique(atoms.resnames)[0] in ['TRP', 'TYR', 'PHE']:
                size = 6
                site_type = "Aromatic"
                pos_arom, direction = assign_aromatic_position(atoms, p_geometric_center)
                pos = [float(x) for x in pos_arom]
                vector = [
                    {k: float(v) for k, v in zip(['x','y','z'], direction)},
                    {k: float(v) for k, v in zip(['x','y','z'], -np.array(direction))}
                ]
            element = pharmer_element_template(site_type, pos, score,
                                               dgsol_sasa[resid]['chainID'],
                                               size=size,
                                               radius=radius,
                                               vector=vector)
            data["points"].append(element)

    #if not args.hbond_direction:
        # for cases where same water molecule is detected for hbond site position. so if hbond direction is not included
        # these site becomes identical. We merge scuh site by summing up their scores
    # incase sites are identical then merge
    data = merge_identical_sites (data)
        
    print_summary(data, os.path.dirname(args.out_file))
    for pos in exclusion_sites.positions:
        element = pharmer_element_template("ExclusionSphere", [float(x) for x in pos], score={'score':1}, chainID='x')
        data["points"].append(element)
    with open(args.out_file, "w") as outfile:
        json.dump(data, outfile, indent=4)


if __name__ == "__main__":
    main()

