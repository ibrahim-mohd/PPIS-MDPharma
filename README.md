# Introduction
This repository provides codes and scripts to obtain protein-protein interaction stabilizers as described in our paper "". In nut-shell, we generate pharmacophore models from MD simulations and search for corresponding ligands that satisfy the pharmacophores in a local data. The scripts here automates the whole procedure.

# Dependencies
## Software Packages

1. **[Pharmer](https://sourceforge.net/projects/pharmer/files/)** – For pharmacophore search and database creation. Executable available via SourceForge.
2. **[GROMACS](https://manual.gromacs.org/documentation/)** – We are mainly using the ``gmx sasa`` routine for calculating the solvation free energy per residue.
3. **[AmberTools](https://ambermd.org/GetAmber.php)** – For post-processing (e.g., obtaining ligand parameters). The Conda version works fine.
4. **[Fpocket](https://github.com/Discngine/fpocket)** – Use this if you do not know the binding pocket or lack a structure file of the protein–protein complex with a known bound ligand.

## Python Packages

1. **[MDAnalysis](https://www.mdanalysis.org/)** – For all analysis.
2. **[MDTraj](https://www.mdtraj.org/1.9.8.dev0/index.html)** – For buried surface area analysis function.
3. **[NetworkX](https://networkx.org/)** – For visualzing and manipulating pharamcophore graphs

# Usage
1. We assume you have an MD simulation trajectory of a protein-protein complex (apo form)  in Gromacs format (we will soon extend for other packages). Make sure to center the protein i.e no breaking over the box edges due to periodic boundary conditions. Use ``gmx trjconv`` function to do so.
2. We first analyze the trajectory for different pharmacophore features: we calculate the Hydrogen bond donor acceptor, solvation free energies of pocket residues and ion binding sites. To do so we need the informatino about the **binding site pocket**. There are two possible options <br>

   (1) You may have an already known ligand in which case you create a ``PDB`` or ``GRO`` file of the whole complex with the ligand in the bound state. You can even dock a ligand and create such a file, the ligand is not important, it should be in the pocket.  <br>
     (2) You use **Fpocket to detect the pockets** in which case you need the Fpocket output which has typically the whole complex and the detected pockets. In this case you also have to specify the pocket ID i.e the residue ID of the pocket at the interface.
```bash
python 01_analyse_pharmacophore_features.py -fl $PWD/protein_out.pdb -f $PWD/mol.xtc -s $PWD/npt.tpr -pocket_id 1 -b 20000 -skip 10 -o output_features.pkl
```
Where, ``protein_out.pdb`` is the Fpocket output and we are intersted int eh pocket with resid (pocket ID) 1. Incase this file already has the ligand bound, no need to specify pocket ID. The flags are similar to gromacs conventions. Also, the default name of solvent resname is ``SOL`` and the solvent oxygen and hydrogen are `OW`, `HW1 HW2`. If that is not the case use extra flags to add the info. 

3. Plot the above features and apply thresholds for feature selection.
```bash
python 02_feature_selection_plot.py -p $PWD/output_features.pkl -c $PWD/mol.gro -s $PWD/npt.tpr -ion_th 50
```
