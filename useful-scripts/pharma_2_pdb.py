import json
import MDAnalysis as mda
import sys
import numpy as np

if len(sys.argv) < 2:
    print("provide the pharmacophore json file as first argument and the output name as second")
    exit()

# input 
input_name = sys.argv [1]

if len (sys.argv) ==2: output_name = "pharmacophore.pdb"

else: output_name = sys.argv [2]




def create_universe (n_atoms, name, resname, positions, masses, resids):

    u_new = mda.Universe.empty(n_atoms=n_atoms,
                             n_residues=n_atoms,
                             atom_resindex=np.arange (n_atoms),
                             residue_segindex=np.arange (n_atoms),
                             n_segments=n_atoms,
                             trajectory=True) # necessary for adding coordinate

    u_new.add_TopologyAttr('name',   name)
    u_new.add_TopologyAttr('resid', resids)
    u_new.add_TopologyAttr('resname', resname)
    u_new.add_TopologyAttr('masses', masses)
    u_new.atoms.positions = positions
    
    return u_new


name_dict = dict (Aromatic="Ar", HydrogenDonor="Dh", HydrogenAcceptor="Ah", Hydrophobic="Ph", PositiveIon="P", 
                ExclusionSphere="Ex",NegativeIon="N")

# = "/home/ibrahim/Documents/Second-docking/refined-set/6no9/pharmacophore.json"

with open(input_name, 'r') as file:
    data = json.load(file)

name = []
positions = []
masses =  []
resname = []
resids = []

for index, d in enumerate (data ['points']):
    
    name.append (name_dict [d ['name']])
    pos = [float(d ['x']), float(d ['y']), float(d ['z'])]
    positions.append (pos)
    resids.append (index+1)
    resname.append ("pharma")
    masses.append(12)


u =  create_universe (n_atoms=len(positions), name=name, resname=resname, positions=positions,
                      masses=masses, resids=resids)



 
u.select_atoms ("all").write (output_name)


