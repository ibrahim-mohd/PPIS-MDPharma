import subprocess
import parmed as pmd
import shutil
import os
import numpy as np   
import argparse
import MDAnalysis as mda

# Initialize the argument parser
parser = argparse.ArgumentParser(description="Setup simulations.")

# Add arguments with single dash flags and default values
parser.add_argument('-p', '--path0', type=str, help="Base directory path")
parser.add_argument('-f', '--pdb_file', type=str, default="protein.pdb", help="protein pdb file")
parser.add_argument('-ff', '--force_field', type=str, default="ff14SB", help="forcefield, ff19SB or ff14SB")
parser.add_argument('-ligand_itp', '--ligand_itp', type=str, default=None, help="provide the full path and name to the ligand itp file")
parser.add_argument('-ligand_gro', '--ligand_gro', type=str, default=None, help="ligand gro file, provide the ligand gro path with full filename")

parser.add_argument('-opc', '--opc_gro', type=str, default="/mnt/second/opc.gro", help="opc strucutre gro file")

parser.add_argument('-sim_time', '--sim_time', type=float, default=100.0, help="simulation time for production in ns")
parser.add_argument('-pdb4amber', '--pdb4amber', type=int, default=1, help="if 1 runs pdb4amber before tleap if 0 do not run. default is 1")


# water and ion model
# Define the lines to add
opc_ions_nonbonded  = """\
[ atomtypes ] ; For water and ions
Na+           11  22.990000  0.00000000  A     0.21844837      0.7047425
Cl-           17  35.450000  0.00000000  A     0.49177609    0.048791716
OW             8  16.000000  0.00000000  A     0.31665521      0.8903586
HW             1   1.008000  0.00000000  A              0              0
MW             0   0.000000  0.00000000  A              0              0
"""

tip_3p_ions_nonbonded = """\
[ atomtypes ] ; For water and ions
;Tip3p  Water
Na+           11  22.990000  0.00000000  A     0.21844837   0.7047425
Cl-           17  35.450000  0.00000000  A     0.49177609   0.048791716
OW            8   16.00      0.0000      A     3.15061e-01  6.36386e-01
HW            1   1.008      0.0000      A     0.00000e+00  0.00000e+00
"""

opc_ions_bonded = """ \
;; ions and stuff
[ moleculetype ]
; molname       nrexcl
SOL             2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW          1       SOL       OW       1       0        16.00000
  2   HW          1       SOL       HW1      1       0.67914   1.00800
  3   HW          1       SOL       HW2      1       0.67914   1.00800
  4   MW          1        SOL      MW       1      -1.35828   0.00000

#ifndef FLEXIBLE

[ settles ]
; i     funct   doh     dhh
1       1       0.08724 0.13712

#else

[ bonds ]
; i     j       funct   length  force.c.
1       2       1       0.08724 502416.0 0.08724        502416.0
1       3       1       0.08724 502416.0 0.08724        502416.0

[ angles ]
; i     j       k       funct   angle   force.c.
2       1       3       1       103.6   628.02  103.6  628.02

#endif


[ virtual_sites3 ]
; Vsite from                    funct   a               b
4       1       2       3       1       0.147722363     0.147722363


[ exclusions ]
1       2       3       4
2       1       3       4
3       1       2       4
4       1       2       3


; The position of the virtual site is computed as follows:
;
;               O
;             
;               V
;         
;       H               H
;
; Ewald OPC:
; const = distance (OV) / [ cos (angle(VOH))    * distance (OH) ]
;         0.01594 nm     / [ cos (51.8 deg)     * 0.0872433 nm    ]
;       then a = b = 0.5 * const = 0.147722363
;
; Vsite pos x4 = x1 + a*(x2-x1) + b*(x3-x1)


;; ions ans stuff

[ moleculetype ]
; Name            nrexcl
Na+          3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 Na+ rtp Na+ q 1.0
    1        Na+      1    Na+    Na+      1 1.00000000  22.990000   ; qtot 1.000000


[ moleculetype ]
; Name            nrexcl
Cl-          3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 Cl- rtp Cl- q -1.0
    1        Cl-      1    Cl-    Cl-      1 -1.00000000  35.450000   ; qtot -1.000000
    
"""

tip3p_ions_bonded = """\
;; ions and stuff
[ moleculetype ]
; molname       nrexcl
SOL             2

[ atoms ]
; id  at type     res nr  res name  at name  cg nr  charge    mass
  1   OW          1       SOL       OW       1      -0.834    16.00000
  2   HW          1       SOL       HW1      1       0.417     1.00800
  3   HW          1       SOL       HW2      1       0.417     1.00800

#ifndef FLEXIBLE

[ settles ]
; OW    funct   doh     dhh
1       1       0.09572 0.15139

[ exclusions ]
1       2       3
2       1       3
3       1       2

#else

[ bonds ]
; i     j       funct   length  force_constant
1       2       1       0.09572 502416.0   0.09572        502416.0
1       3       1       0.09572 502416.0   0.09572        502416.0


[ angles ]
; i     j       k       funct   angle   force_constant
2       1       3       1       104.52  628.02      104.52  628.02

#endif


;; ions ans stuff

[ moleculetype ]
; Name            nrexcl
Na+          3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 Na+ rtp Na+ q 1.0
    1        Na+      1    Na+    Na+      1 1.00000000  22.990000   ; qtot 1.000000


[ moleculetype ]
; Name            nrexcl
Cl-          3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass  typeB    chargeB      massB
; residue    1 Cl- rtp Cl- q -1.0
    1        Cl-      1    Cl-    Cl-      1 -1.00000000  35.450000   ; qtot -1.000000
"""

## MDP files

em_mdp = """\
 
integrator               = steep
nsteps                   = 1000000

; Restrain all protein positions
constraints              = none

; EM criteria and other stuff
emtol                    = 10
emstep                   = 0.01
niter                    = 20
nbfgscorr                = 10

; Output control
nstxout                  = 1000
nstvout                  = 50000
nstfout                  = 0
nstlog                   = 50000
nstenergy                = 50000
nstxtcout                = 20


; Neighborsearching and short-range nonbonded interactions
nstlist                  = 10
ns_type                  = grid
pbc                      = xyz
rlist                    = 1.2

; Electrostatics
coulombtype              = PME
rcoulomb                 = 1.2

; van der Waals
vdw-type                 = cut-off
rvdw-switch              = 1.2
rvdw                     = 1.2

; Apply long range dispersion corrections for Energy and Pressure
DispCorr                  = EnerPres

; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; EWALD/PME/PPPM parameters
pme_order                = 6
ewald_rtol               = 1e-06
epsilon_surface          = 0
optimize_fft             = no

; Temperature and pressure coupling are off during EM
tcoupl                   = no
pcoupl                   = no

"""
npt_restraint_mdp = """\
define                  = {} ;
integrator              = md
dt                      = 0.002
nsteps                  = 100000
nstlog                  = 10000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 500  ;every 100 ps
nstcalcenergy           = 100000
nstenergy               = 100000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = pme
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Potential-shift
rvdw                    = 1.2
;
tcoupl                  = v-rescale
tc_grps                 = Sol_ions Protein_lig
tau_t                   = 1.0  1.0
ref_t                   = 300 300
;
pcoupl                   = berendsen
pcoupltype               = isotropic
tau_p                    = 5.0
compressibility          = 4.5e-05 0
ref_p                    = 1.0 1.0

;
constraints             = h-bonds
constraint_algorithm    = LINCS
;continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = system
;
refcoord_scaling        = com
gen-vel                 = yes
gen-temp                = 300
gen-seed                = 42
"""

npt_ber_mdp = """\
 
integrator              = md
dt                      = 0.002
nsteps                  = 500000
nstlog                  = 10000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 25000  ;every 100 ps
nstcalcenergy           = 100000
nstenergy               = 100000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = pme
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Potential-shift
rvdw                    = 1.2
;
tcoupl                  = v-rescale
tc_grps                 = Sol_ions Protein_lig
tau_t                   = 1.0  1.0; 1.0
ref_t                   = 300 300
;
pcoupl                   = C-rescale
pcoupltype               = isotropic
tau_p                    = 5.0
compressibility          = 4.5e-05 
ref_p                    = 1.0

;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = no
gen-vel                 = yes
gen-temp                = 300
gen-seed                = -1

;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = system
;
refcoord_scaling        = com

"""

npt_pr_mdp = """\
integrator              = md
dt                      = 0.002
nsteps                  = {}
nstlog                  = 10000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 10000  ;every 100 ps
nstcalcenergy           = 10000
nstenergy               = 10000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = pme
rcoulomb                = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Potential-shift
rvdw                    = 1.2
;
tcoupl                  = v-rescale
tc_grps                 = Sol_ions Protein_lig
tau_t                   = 1.0  1.0; 1.0
ref_t                   = 300 300

;
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic 
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = system
;
refcoord_scaling        = com
"""

def partition_topol_file(input_file = "topol.top", out_path = "./"):
    """
    Partitions a file into multiple files based on '[ moleculetype ]' sections.
    Args:
        input_file (str): Path to the input file.
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    output_files = []
    current_section = []
    section_count = 0
    topol_write = False
    topol_lines = []
    for line in lines:
        # Check for a new '[ moleculetype ]' section
        if line.strip().startswith('[ moleculetype ]'):
            # Save the previous section if it exists
            if current_section:
                if section_count ==0:
                    output_filename = "ff_nonbonded.itp"
                else:
                    output_filename = f"protein{section_count}.itp"
                
                with open(out_path + output_filename, 'w') as out_file:
                    out_file.writelines(current_section)
                output_files.append(output_filename)
                section_count += 1
                current_section = []
                
        if line.strip().startswith('[ system ]') or  topol_write:
            topol_lines.append (line)
            topol_write=True
            
        if not topol_write:
            current_section.append(line)

    # Save the last section if any content exists
    if current_section:
        output_filename  = f"protein{section_count}.itp"
        #output_filename = f"moleculetype_{section_count}.top"
        with open(out_path + output_filename, 'w') as out_file:
            out_file.writelines(current_section)
        output_files.append(output_filename)

    return topol_lines, section_count
     

def cmap_correction (amber):

# This function is adapted from  https://github.com/ParmEd/ParmEd/issues/1292 from the user csy0000.
# I have mannually checked its accuracy by comparing with Charmm-gui generated files and found to be identical
    #amber = pmd.load_file(options.prmtop,options.inpcrd)
    
        # Cmap correction
    CX_dir = { 'ALA' : 'XC0', 'ARG' : 'XC1', 'ASH' : 'XC2', 'ASN' : 'XC3', 'ASP' : 'XC4', 'CYM' : 'XC5', 'CYS' : 'XC6', 
               'CYX' : 'XC5', 'GLH' : 'XC7', 'GLN' : 'XC8', 'GLU' : 'XC9', 'GLY' : 'XC10', 'HID' : 'XC5', 'HIE' : 'XC5',
               'HIP' : 'XC5', 'HYP' : 'XC5', 'ILE' : 'XC11', 'LEU' : 'XC5', 'LYN' : 'XC12', 'LYS' : 'XC12', 'MET' : 'XC5',
               'PHE' : 'XC5', 'PRO' : 'XC13', 'SER' : 'XC14', 'THR' : 'XC15', 'TRP' : 'XC5', 'TYR' : 'XC5', 'VAL' : 'XC11' } 
    
    known_types = {}
    for cmap in amber.cmaps:  # loop through all cmap terms
        atom_type = cmap.atom3.atom_type # get atom type
        new_typename = CX_dir[cmap.atom3.residue.name] # get new CA type name (e.g. XC0)
        if new_typename in known_types: # if new type name in known types, directly assign atom to new type
            cmap.atom3.atom_type = known_types[cmap.type]
            cmap.atom3.type = known_types[cmap.type].name
            continue # skip
    
        # if new type name is not in known types, create new type
        new_type = pmd.AtomType(new_typename, atom_type.number, atom_type.mass, atom_type.atomic_number, new_typename, atom_type.charge)
        new_type.epsilon = atom_type.epsilon # copy over epsilon and rmin
        new_type.rmin = atom_type.rmin
        known_types[cmap.type] = new_type
    
        # assgin atom to new type
        cmap.atom3.atom_type = known_types[cmap.type]
        cmap.atom3.type = known_types[cmap.type].name
    # Replace 'topol.top' with your file's nam

    return amber


# Parse the arguments
args = parser.parse_args()

# input arguments
path0        = args.path0
pdb_file     = args.pdb_file
force_field  = args.force_field
ligand_itp  = args.ligand_itp
ligand_gro   = args.ligand_gro
opc_gro      = args.opc_gro
nsteps       = int (args.sim_time*1e3/0.002) # 
#0.002# create a directory to place everything
os.makedirs(path0 + "/initial/", exist_ok=True)    
#except Exception as e:
shutil.copy(os.path.join(path0, pdb_file), path0 + "/initial/")
os.chdir(path0+ "/initial/")

## 
conda_env = "AmberTools23"
activate_command = f"source $(conda info --base)/etc/profile.d/conda.sh && conda activate {conda_env}"
# activate
## Run pdb for amber
subprocess.run (activate_command, shell=True, executable="/bin/bash")
 
#Step 2
if  args.pdb4amber:
    command = f"pdb4amber -i {pdb_file} -o pdb4amber.pdb"
    subprocess.run(command, shell=True, executable="/bin/bash")
    pdb_file = "pdb4amber.pdb"
 
 
# strip the hydrogens
#mda.Universe ("pdb4amber.pdb").select_atoms("all and not name H* 1H* 2H* 3H*").write ("pdb4amber.pdb")
# Step 3: Run tleap with the input file

tleap_input_file = "tleap.in"

tleap_commands = f"""\
source leaprc.protein.{force_field}
pro = loadpdb {pdb_file}
saveamberparm pro protein.prmtop protein.inpcrd
quit
"""

with open(tleap_input_file, "w") as file: 
    file.write(tleap_commands)

tleap_command = f"tleap -f {tleap_input_file}"
result = subprocess.run(tleap_command, shell=True, executable="/bin/bash")
#####################################################################
### if ff19, save with proper c-map correction
parm = pmd.load_file('protein.prmtop', 'protein.inpcrd')

if force_field == "ff19SB": parm = cmap_correction (parm)
    
# Step 4: Convert to gromacs topology files
parm.save('topol.top', format='gromacs', overwrite=True)
parm.save('protein.gro', overwrite=True)


# Step 5. Create topology and itp files
os.makedirs(path0 + "/itp/", exist_ok=True)  
# if there is ligand
## Add ligand topology
if ligand_itp:
    # read ligand itp
    with open(ligand_itp, 'r') as f: lines = f.readlines()

    # create atom type a sperate itp file
    
    for index, line in enumerate (lines):
        if line.strip().startswith('[ moleculetype ]'): 
            split_index = index 
            break
    
    ligand_atom_type = lines [:split_index]

    with open(path0 + "/itp/ligand.itp", 'w+') as f:
        for line in lines [split_index:]:
            f.write (line)

    #get ligand name

    for index, line in enumerate (lines):
        if line.strip().startswith("[ atoms ]"):
            name_index = index-2
            break
    line = lines [name_index].split (" ")
    while line.count(""): line.remove ("")
    ligand_name = line [0]
    
#print (ligand_name)

## water and ion topology file
with open(path0 + "/itp/water_ions.itp", 'w+') as f:
    
    if force_field == "ff19SB": f.write (opc_ions_bonded)
    
    elif force_field == "ff14SB": f.write (tip3p_ions_bonded)    
    
    else: print ("forcefield can be either of ff14SB or ff19SB")


####################3
topol_lines, n_proteins = partition_topol_file(path0 + "/initial/" + "topol.top",out_path = path0 + "/itp/")

## write topology file
alphabets = "ABCDEFGHIJKL"
with open(path0 + "/topol.top", 'w') as out_file:
    
    out_file.write("#include \"./itp/ff_nonbonded.itp\"\n")
    
    for n in range (1,n_proteins+1):
        
        out_file.write(f"\n#include \"./itp/protein{n}.itp\"\n")
        # write restraint file
        out_file.write (f"#ifdef POSRES_{alphabets[n-1]} \n#include \"./itp/posre_{alphabets[n-1]}.itp\"\n#endif\n")
        
    if ligand_itp:
        out_file.write(f"\n#include \"./itp/ligand.itp\"\n")
        out_file.write (f"#ifdef POSRES_lig \n#include \"./itp/posre_lig.itp\"\n#endif\n")

    # Add water and ions
    out_file.write(f"\n#include \"./itp/water_ions.itp\"\n")

    for line in topol_lines: out_file.write (line)

    if ligand_itp:
        out_file.write (f"{ligand_name} \t 1\n")

# modify the non-bonded file to add solvent, ion and ligand atomtypes
with open(path0 + "/itp/ff_nonbonded.itp", 'r') as f: non_bonded_file = f.readlines()
 
for index, line in enumerate (non_bonded_file):

    # the followign is needed beacause we need to add other non bonded after the protein atom types
    # for ff19SB we have also the cmap section so the atoms types should be placed before it
    # for ff14sb we jsut append the end of file.
    if force_field == "ff14SB":  
        type_index = -1
        break
        
    if force_field == "ff19SB":  
        
        if line.strip().startswith("[ cmaptypes ]"):           
            type_index = index
            break
## add water and ion types
with open(path0 + "/itp/ff_nonbonded.itp", "w+") as output_file:
    
    output_file.writelines (non_bonded_file [:type_index])
    
    if force_field == "ff19SB": output_file.write (opc_ions_nonbonded)
    
    elif force_field == "ff14SB": output_file.write (tip_3p_ions_nonbonded)    
    
    else: print ("forcefield can be either of ff14SB or ff19SB")
    ## add ligand
    if ligand_itp: output_file.writelines (ligand_atom_type)
    ## write back the rest
    output_file.writelines (non_bonded_file [type_index:])     

##########################################################################################################
################# Topology and forcefield creation done #############################################
########################################################################################################

###### From here, simulation setup creation is carried out ##########################################33
os.chdir(path0+ "/initial/")

if ligand_itp: 
    u_protein = mda.Universe ( path0+ "/initial/"+"protein.gro")
    u_ligand =  mda.Universe (ligand_gro)
    mda.Merge (u_protein.atoms, u_ligand.atoms).atoms.write (path0 + "/initial/protein_lig.gro")
    gro_file = "protein_lig.gro"
else: 
    gro_file = "protein.gro"

#1.0 create box

command = f"gmx editconf -f {path0}/initial/{gro_file} -c -d 1.2 -bt dodecahedron -o {path0}/initial/box.gro"
subprocess.run (command, shell=True, executable="/bin/bash")

# solvate 
if force_field == "ff19SB":
    command = f"gmx solvate -cp {path0}/initial/box.gro -cs {opc_gro} -p {path0}/topol.top -o {path0}/initial/sol.gro" 
else:
    command = f"gmx solvate -cp {path0}/initial/box.gro -p {path0}/topol.top -o {path0}/initial/sol.gro" 

subprocess.run (command, shell=True, executable="/bin/bash")

# Add ions
u = mda.Universe (f"{path0}/initial/sol.gro")
Nsol = len(u.select_atoms ("resname SOL"))
if force_field == "ff19SB": Nsol /= 4
else: Nsol /= 3
####
N_ions = int(0.15*Nsol/55.55) #
#####################

# write em mdp
os.makedirs(path0 + "/mdp/", exist_ok=True)  
with open (path0 + "/mdp/em.mdp", 'w') as out_file: out_file.write (em_mdp)

command = f"gmx grompp -f {path0}/mdp/em.mdp -c {path0}/initial/sol.gro -p {path0}/topol.top -o {path0}/initial/ion.tpr -maxwarn 1"
subprocess.run (command, shell=True, executable="/bin/bash")

## add ion
os.makedirs(path0 + "/em/", exist_ok=True)  
command = f"echo SOL | gmx genion -s {path0}/initial/ion.tpr -p {path0}/topol.top -nname Cl- -pname Na+ -nn {N_ions} -neutral -o {path0}/em/ion.gro"
subprocess.run (command, shell=True, executable="/bin/bash")

# energy minimize
command_grompp = f"gmx grompp -f {path0}/mdp/em.mdp -c {path0}/em/ion.gro -p {path0}/topol.top -o {path0}/em/em -maxwarn 1"
command_mdrun = f"gmx mdrun -deffnm {path0}/em/em -v"

subprocess.run (command_grompp, shell=True, executable="/bin/bash")
subprocess.run (command_mdrun, shell=True, executable="/bin/bash")

## create index file and individual gro files for creating restraint files
u  = mda.Universe (path0 + "/em/em.tpr", path0 + "/em/em.gro")
restrain_line = "" # for mdp
with mda.selections.gromacs.SelectionWriter(path0 + "/index.ndx", mode='w') as ndx:

    sol_ions = u.select_atoms ("resname SOL or name Cl- Na+")
    protein_lig = u.select_atoms ("all and not (resname SOL or name Cl- Na+)")
    system = u.select_atoms ("all")

    
    ndx.write(system, name='System')
    # write for each protein chain
    for n in range (0, n_proteins):
        ndx.write(u.atoms.fragments [n], name=f'Chain_{alphabets [n]}')
        u.atoms.fragments [n].write (path0 + f"/initial/chain{alphabets [n]}.gro")
        restrain_line = restrain_line + " " + f"-DPOSRES_{alphabets[n]}"
        
    if ligand_itp:
        ligand = u.select_atoms ("all and not protein and not (resname SOL or name Cl- Na+)")
        ndx.write(ligand, name="ligand")
        ligand.write (path0 + "/initial/ligand.gro")
        restrain_line = restrain_line + " " + "-DPOSRES_lig"
        
    ndx.write(protein_lig, name='Protein_lig')
    ndx.write(sol_ions, name='Sol_ions')

## create directories
os.makedirs(path0 + "/npt1/", exist_ok=True)  
os.makedirs(path0 + "/npt2/", exist_ok=True) 
os.makedirs(path0 + "/npt-restr/", exist_ok=True) 
## create mdp files
#1 restraint mdp file
with open(path0 + "/mdp/md-ber-restr.mdp", 'w') as f:   
    f.write (npt_restraint_mdp.format (restrain_line))

# berendsen equilibration file
with open(path0 + "/mdp/md-ber.mdp", 'w') as f:    
    f.write (npt_ber_mdp)

# parrinello-rahman equilibration file
with open(path0 + "/mdp/md-pr.mdp", 'w') as f:    
    f.write (npt_pr_mdp.format (nsteps))


# write a submission script
with open(path0 + "/sub0.sh", 'w') as f:
    first_lines = """\
count=0
gro_file=em/em.gro
for fc in 1000 100 10 1;do
\t count=$((count+1))
"""
    f.write (first_lines)
    for n in range (0, n_proteins):
        line = f"\t echo Backbone | gmx genrestr -f initial/chain{alphabets[n]}.gro -o itp/posre_{alphabets[n]}.itp -fc $fc\n"
        f.write (line)
    if ligand_itp:
        line = f"\t echo 0 | gmx genrestr -f initial/ligand.gro -o itp/posre_lig.itp -fc $fc\n"
        f.write (line)
    f.write ("\t rm itp/*#\n")
    line = "\t gmx grompp -f ./mdp/md-ber-restr.mdp  -n index.ndx -p topol.top -c $gro_file -r $gro_file  -o npt-restr/npt$count -maxwarn 1\n"
    f.write (line)
    line = "\t gmx mdrun -deffnm npt-restr/npt$count -v \n \t gro_file=npt-restr/npt$count.gro \n\n done"
    f.write (line) 

with open(path0 + "/sub1.sh", 'w') as f:
    # berendsen equilibbration
    line = "gmx grompp -f ./mdp/md-ber.mdp  -n index.ndx -p topol.top -c npt-restr/npt4.gro -o npt1/npt\n"
    f.write (line)
    line = "gmx mdrun -deffnm npt1/npt -v \n"
    f.write (line) 
    # production with parrinello rahman-barostat

    line = "gmx grompp -f ./mdp/md-pr.mdp  -n index.ndx -p topol.top -c npt1/npt.gro -o npt2/npt\n"
    f.write (line)
    line = "gmx mdrun -deffnm npt2/npt -v \n"
    f.write (line) 
