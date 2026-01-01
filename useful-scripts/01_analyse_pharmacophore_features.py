# Written by Mohd Ibrahim
# Technical University of Munich

import argparse
from tqdm import tqdm
import numpy as np
import subprocess
import pickle
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis


# ---------------- Utility ---------------- #
def skip_header_np_gentxt(filename):
    file = open(filename, "r+").readlines()
    skip_header = 0
    Resids = []
    for f in file:
        if "@ s" in f:
            key = f.split(" ")[3][1:-2]
        if "@" in f or "#" in f:
            skip_header += 1
        if "Res" in f:
            temp_ = f.split(" ")
            while temp_.count('') > 1:
                temp_.remove('')
            Resids.append(temp_[3].removesuffix("\n")[4:-1])
    return skip_header, Resids

  
# ---------------- Main ---------------- #
def main():
    parser = argparse.ArgumentParser(description="Combined Pocket Analysis Pipeline")

    # General inputs
    parser.add_argument('-fl', dest='coord_file', type=str, default='protein_lig.gro', help='protein-ligand gro or pdb file or FPocket output with protein and pockets')
    parser.add_argument('-f', dest='xtc_file', type=str, default='mol.xtc', help='input xtc file')
    parser.add_argument('-s', dest='tpr_file', type=str, default='npt.tpr', help='input tpr file')
    parser.add_argument("-pocket_id", dest="pocket_id", type=str, default=None, help="Pocket ID (0 = ligand-based pocket)")
    parser.add_argument('-b', dest='begin_time', type=float, default=0, help='begin time in ps')
    parser.add_argument('-e', dest='end_time', type=float, default=-1, help='end time in ps')
    parser.add_argument('-skip', dest='frame_skip', type=int, default=1, help='skip every n frames')

    # Ion options
    parser.add_argument('-nname', dest='nname', type=str, default='Cl-', help='anion name')
    parser.add_argument('-pname', dest='pname', type=str, default='Na+', help='cation name')
    # Solvent residue and atom names
    
    parser.add_argument('-hsol_name', dest='hsol_name', type=str, default='HW1 HW2', help='atom name of water hydrogens')
    parser.add_argument('-osol_name', dest='osol_name', type=str, default='OW', help='atomname of water oxygens')
    parser.add_argument('-sol_resname', dest='sol_resname', type=str, default='SOL', help='resname of water molecules')

    # exclude residues
    parser.add_argument('-res_exclude', dest='res_exclude', type=str, default=None, help='residues to exclude for analysis e.g those that already contributes significantly to stabilziation of the complex')
    parser.add_argument('-pocket_resids', dest='pocket_resids', type=str, default=None, help='pocket residues as space separated resids if known')
    ########
    parser.add_argument('-cutoff', dest='contact_cutoff', type=float, default=5, help='all residues within this distance of the ligand or pocket atoms are considered for analysis. Default is 5 Angstrom')

    # Output
    parser.add_argument('-o', dest='pkl_object', type=str, default='combined_analysis.pkl', help='output pickle object')

    args = parser.parse_args()

    # ---------------- Pocket selection ---------------- #
    coord_file = args.coord_file
    pocket_id = args.pocket_id
    res_exclude = args.res_exclude
    contact_cutoff = args.contact_cutoff
    u = mda.Universe(coord_file)
    # to renumber protien resids from 1 to N in a continous way incase for if the resid is reset for different chains ie. 
    # resid starts from 1 for each chain.  
    protein_ = u.select_atoms("protein")
    protein_.residues.resids = np.arange(1, len(protein_.residues.resids) + 1)
    
    if not pocket_id:
        
        if  args.pocket_resids: pocket = u.select_atoms(f"protein and resid {args.pocket_resids}")
            
        else:
            
            ligand_name = u.select_atoms("all and not protein").resnames[0]
            
            if res_exclude:
                pocket = u.select_atoms(f"byres protein and not resid {res_exclude} and (around {contact_cutoff} resname {ligand_name})")
                
            else:
                pocket = u.select_atoms(f"byres protein and (around {contact_cutoff} resname {ligand_name})")
    else:
        
        if res_exclude:
            pocket = u.select_atoms(f"byres protein and not resid {res_exclude} and (around {contact_cutoff} resname STP and resid {pocket_id})")
        else:
            pocket = u.select_atoms(f"byres protein and (around {contact_cutoff} resname STP and resid {pocket_id})")
    
  

        
    # ---------------- SASA & ΔGsolv ---------------- #
    u = mda.Universe(args.tpr_file, args.xtc_file)
    out_ndx = "index.ndx"
    out_sasa = "sasa.xvg"
    out_gsolv = "dgsolv.xvg"

    each_res = [u.select_atoms(f"resid {resid} and protein") for resid in np.unique(pocket.resids)]
    with mda.selections.gromacs.SelectionWriter(out_ndx, mode='w') as ndx:
        ndx.write(u.select_atoms("all"), name='system')
        ndx.write(u.select_atoms("protein"), name='Protein')
        ndx.write(pocket, name='pocket')
        for selection, resid in zip(each_res, np.unique(pocket.resids)):
            ndx.write(selection, name=f'Res{resid}')

    pocket_residues = " ".join(str(x+3) for x in range(len(each_res)))
    gromacs_command = (
        f"gmx sasa -f {args.xtc_file} -s {args.tpr_file} -n {out_ndx} -o {out_sasa} "
        f"-odg {out_gsolv} -b {args.begin_time} -surface 1 -output {pocket_residues}"
    )
    subprocess.run(gromacs_command, shell=True)

    skip_header, pocket_resids = skip_header_np_gentxt(out_gsolv)
    dG = np.genfromtxt(out_gsolv, skip_header=skip_header)
    skip_header, pocket_resids = skip_header_np_gentxt(out_sasa)
    sasa = np.genfromtxt(out_sasa, skip_header=skip_header)

    sasa_dg_dict = {}
    for index, resid in enumerate(pocket_resids):
        residue = u.select_atoms("protein and resid %s" % (resid))
        dG_ave, std_dg = np.average(dG[:, index+2]), np.std(dG[:, index+2])
        sasa_ave, std_sasa = np.average(sasa[:, index+2]), np.std(sasa[:, index+2])
        sasa_dg_dict[resid] = dict(
            resname=residue.resnames[0],
            dG=(dG_ave, std_dg),
            sasa=(sasa_ave, std_sasa)
        )

    # ---------------- H-Bonds ---------------- #
    
    print ("Calculating hydrogen bond donors and acceptors ...Its going to take a while sit by and relax\n")
    
    #u = mda.Universe(args.tpr_file, args.xtc_file)
        
    #u = mda.Universe(args.tpr_file, args.xtc_file)
    dt = u.trajectory.dt
    begin_frame = int(args.begin_time / dt)
    end_frame = int(np.floor(args.end_time/u.trajectory.dt)) if args.end_time > 0 else None
    
    
    pocket_acceptors = pocket.select_atoms("name O* N*")
    h_bond_acceptors = " ".join(str(x) for x in pocket_acceptors.indices)

    h_bond_donors = ""
    donated_hydrogens = ""
    for index in pocket_acceptors.indices:
        if index+1 < len(u.atoms) and u.atoms.names[index+1][0] == "H":
            h_bond_donors += f" {index}"
            donated_hydrogens += f" {index+1}"

    hbonds_acceptors = HydrogenBondAnalysis(
        universe=u,
        donors_sel=f"name {args.osol_name} and resname {args.sol_resname}",
        hydrogens_sel=f"name {args.hsol_name} and resname {args.sol_resname}",
        acceptors_sel=f"index {h_bond_acceptors}",
        d_a_cutoff=3.0, d_h_a_angle_cutoff=150,
        update_selections=False
    )
    hbonds_donors = HydrogenBondAnalysis(
        universe=u,
        donors_sel=f"index {h_bond_donors}",
        hydrogens_sel=f"index {donated_hydrogens}",
        acceptors_sel=f"name {args.osol_name} and resname {args.sol_resname}",
        d_a_cutoff=3.0, d_h_a_angle_cutoff=150,
        update_selections=False
    )
    
    hbonds_acceptors.run(start=begin_frame, stop=end_frame, step=args.frame_skip,verbose=True)
    hbonds_donors.run(start=begin_frame, stop=end_frame, step=args.frame_skip,verbose=True)
    
    #hbonds_acceptors.run(start=begin_frame,verbose=True)
    #hbonds_donors.run(start=begin_frame,verbose=True)

    if not end_frame: end_frame = u.trajectory.n_frames
    n_frames = int((end_frame - begin_frame)/args.frame_skip )    
    #n_frames = u.trajectory.n_frames
    dict_acceptors = {str(i): dict(name=n, count=0) for i, n in zip(pocket_acceptors.indices, pocket_acceptors.names)}
    pocket_donors = u.select_atoms(f"index {h_bond_donors}") if h_bond_donors else u.atoms[:0]
    dict_donors = {str(i): dict(name=n, count=0) for i, n in zip(pocket_donors.indices, pocket_donors.names)}

    for i in hbonds_acceptors.results['hbonds']:
        dict_acceptors[str(int(i[3]))]['count'] += 1/n_frames
    for i in hbonds_donors.results['hbonds']:
        dict_donors[str(int(i[1]))]['count'] += 1/n_frames

    # Normalize by SASA-per-heavy-atom
    for key in dict_acceptors.keys():
        resid = u.atoms[int(key)].resid
        sasa_norm = sasa_dg_dict.get(str(resid), {}).get('sasa', (1,))[0]
        sasa_norm = sasa_norm / len(u.select_atoms(f"resid {resid} and protein and not name H*"))
        dict_acceptors[key]['count'] /= sasa_norm if sasa_norm > 0 else 1

    for key in dict_donors.keys():
        resid = u.atoms[int(key)].resid
        sasa_norm = sasa_dg_dict.get(str(resid), {}).get('sasa', (1,))[0]
        sasa_norm = sasa_norm / len(u.select_atoms(f"resid {resid} and protein and not name H*"))
        dict_donors[key]['count'] /= sasa_norm if sasa_norm > 0 else 1

    hbond_dict = dict(acceptor=dict_acceptors, donor=dict_donors)

    # ---------------- Ion Binding ---------------- #
    print ("Calculating coloumb interaction sites ...\n")


    NA = u.select_atoms(f"name {args.pname}")
    CL = u.select_atoms(f"name {args.nname}")
    ion_sites = pocket.select_atoms("name O* N*")

    residue_dict_na = {str(i): dict(name=n, count=0) for i, n in zip(ion_sites.indices, ion_sites.names)}
    residue_dict_cl = {str(i): dict(name=n, count=0) for i, n in zip(ion_sites.indices, ion_sites.names)}

    for _ in tqdm(u.trajectory[begin_frame:end_frame:args.frame_skip]):
        pairs, _ = mda.lib.distances.capped_distance(ion_sites, NA, max_cutoff=3, box=u.dimensions)
        for i_res, _ in pairs:
            index = ion_sites.indices[i_res]
            residue_dict_na[str(index)]['count'] += 1

        pairs, _ = mda.lib.distances.capped_distance(ion_sites, CL, max_cutoff=3, box=u.dimensions)
        for i_res, _ in pairs:
            index = ion_sites.indices[i_res]
            residue_dict_cl[str(index)]['count'] += 1

    
    # Normalize the frame to percentage of total frame for easy interpretation

    for key in residue_dict_cl.keys(): 
        residue_dict_cl [key] ['count'] = residue_dict_cl [key] ['count']*100/n_frames  
    
    for key in residue_dict_na.keys(): 
        residue_dict_na [key] ['count'] = residue_dict_na [key] ['count']*100/n_frames
        
    
    ion_dict = dict(cation=residue_dict_na, anion=residue_dict_cl)

    # ---------------- Final Output ---------------- #
    output = dict(
        sasa_dG=sasa_dg_dict,
        hbonds=hbond_dict,
        ions=ion_dict
    )
    with open(args.pkl_object, 'wb') as file:
        pickle.dump(output, file)

if __name__ == "__main__":
    main()
