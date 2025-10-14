import MDAnalysis as mda
import numpy as np
import subprocess
import argparse
import warnings

# Configure warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser (description="Remove molecules breakign across periodic box for gromacs trajectories")

#########################################################

parser.add_argument('-f', dest='xtc_file', type=str, default='npt.xtc',help='input xtc file')
parser.add_argument('-s', dest='tpr_file', type=str, default='npt.tpr',help='input tpr file')
parser.add_argument('-o', dest='out_xtc', type=str, default='mol.xtc',help='output xtc file')
parser.add_argument('-on', dest='out_ndx', type=str, default='mol_ndx.ndx',help='output index file')
 

args                = parser.parse_args()
############################################
xtc_file            = args.xtc_file
tpr_file            = args.tpr_file
out_xtc             = args.out_xtc
out_ndx             = args.out_ndx

###################################################
letters = "ABCDEFGHIJKLMNOP" # just for naming chains
# load trajectory
u = mda.Universe (tpr_file, xtc_file)

protein = u.select_atoms ("protein")
protein_chains = [chain for chain in protein.fragments]

# createa n index file
with mda.selections.gromacs.SelectionWriter(out_ndx, mode='w') as ndx:
    ndx.write(u.select_atoms('all'), name='system')
    ndx.write(protein, name='Protein')
    ndx.write(u.select_atoms ("all and not protein"), name='all_not_protein')
    # write index for seperate chains
    for idx, chain in enumerate (protein_chains):
        ndx.write (chain, name=f'chain{letters[idx]}')


# center the molecule on one of the chains
gmx_cmd = (
    f'echo "chainA system" | '
    f'gmx trjconv -f {xtc_file} -s {tpr_file} -center -pbc mol -ur compact '
    f'-n {out_ndx} -o {out_xtc}'
)

print("Running:", gmx_cmd)
subprocess.run(gmx_cmd, shell=True, check=True)


