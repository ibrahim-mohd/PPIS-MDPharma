# Introduction
This repository provides codes and scripts to obtain protein-protein interaction stabilizers as described in our paper "". In nut-shell, we generate pharmacophore models from MD simulations and search for corresponding ligands that satisfy the pharmacophores in a local data. The scripts here automates the whole procedure.

# Dependencies
1. Pharmer (Executable at: https://sourceforge.net/projects/pharmer/files/)
2. Gromacs (https://manual.gromacs.org/documentation/)
3. AmberTools (For post processing e.g to obtain ligand parameters, conda version is fine: https://ambermd.org/GetAmber.php)
4. Fpocket (https://github.com/Discngine/fpocket) incase you do not know the binding pocket or do  not have a structure file of the protein-protein complex with a known bound ligand.
## Python packages
1. MDAnalysis (https://www.mdanalysis.org/)
2. MDtraj (https://www.mdtraj.org/1.9.8.dev0/index.html)
3. networkx (https://networkx.org/)

# Usage
1. We assume you have an MD simulation trajectory of a protein-protein complex (apo form)  in Gromacs format (we will soon extend for other packages). Make sure to center the protein i.e no breaking over the box edges due to periodic boundary conditions. Use ``gmx trjconv`` function to do so.
