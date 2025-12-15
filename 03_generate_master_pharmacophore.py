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

        # Solvent residue and atom names
    
    parser.add_argument('-hsol_name', dest='hsol_name', type=str, default='HW1 HW2', help='atom name of water hydrogens')
    parser.add_argument('-osol_name', dest='osol_name', type=str, default='OW', help='atomname of water oxygens')
    parser.add_argument('-sol_resname', dest='sol_resname', type=str, default='SOL', help='resname of water molecules')
    
    return parser.parse_args()

def distance(r1, r2):
    return np.linalg.norm(r1 - r2)


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


def assign_position_donor(site, initial_acceptor_pos, u):
    hydrogen_pos = u.select_atoms(f"protein and name H* and around 1.5 index {site}").positions
    hydrogen_pos = hydrogen_pos.mean(axis=0) if len(hydrogen_pos) > 1 else hydrogen_pos[0]
    donor_position = u.select_atoms(f"protein and index {site}").positions[0]
    DH = hydrogen_pos - donor_position
    HA = hydrogen_pos - initial_acceptor_pos
    DH /= np.linalg.norm(DH)
    HA /= np.linalg.norm(HA)
    rotated_HA = HA.copy()
    def rotate_vector(v, axis, theta_deg):
        theta = np.radians(theta_deg)
        axis = axis / np.linalg.norm(axis)
        return (v * np.cos(theta) +
                np.cross(axis, v) * np.sin(theta) +
                axis * np.dot(axis, v) * (1 - np.cos(theta)))
    angle = np.degrees(np.arccos(np.dot(DH, HA)))
    if angle < 150:
        axis = np.cross(DH, HA)
        if np.linalg.norm(axis) != 0:
            axis /= np.linalg.norm(axis)
            rotate_by = 150 - angle
            HA_plus  = rotate_vector(HA, axis,  rotate_by)
            HA_minus = rotate_vector(HA, axis, -rotate_by)
            angle_plus  = np.degrees(np.arccos(np.dot(DH, HA_plus)))
            angle_minus = np.degrees(np.arccos(np.dot(DH, HA_minus)))
            rotated_HA = HA_plus if angle_plus > angle_minus else HA_minus
    updated_acceptor_position = hydrogen_pos + rotated_HA * 3
    hbond_direction = donor_position - updated_acceptor_position
    hbond_direction /= np.linalg.norm(hbond_direction)
    return [float(x) for x in updated_acceptor_position], [float(x) for x in hbond_direction]


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
            if site_type_protein == 'donor':
                pos_hb, direction = get_hbond_donors(site, u)
                if pos_hb is None:
                    # Usually these donors are in very hidden places 
                    #pos_hb, direction = assign_position_donor(site, position, u)
                    print(f"Donor:{site} no water found")
                    continue 
                pos = [float(x) for x in pos_hb]
                vector = [{k: float(v) for k, v in zip(["x", "y", "z"], direction)}]
            if site_type_protein == 'acceptor':
                pos_hb, direction = get_hbond_acceptors(site, u,osol_name=args.osol_name, sol_resname=args.sol_resname, hsol_name=args.hsol_name)
                if pos_hb is None:
                    print(f"Acceptor:{site} no water found")
                    continue
                    direction = position - p_geometric_center
                    direction /= np.linalg.norm(direction)
                    pos_hb = position
                pos = [float(x) for x in pos_hb]
                vector = [{k: float(v) for k, v in zip(["x", "y", "z"], direction)}]
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

    print_summary(data, os.path.dirname(args.out_file))
    for pos in exclusion_sites.positions:
        element = pharmer_element_template("ExclusionSphere", [float(x) for x in pos], score={'score':1}, chainID='x')
        data["points"].append(element)
    with open(args.out_file, "w") as outfile:
        json.dump(data, outfile, indent=4)


if __name__ == "__main__":
    main()

