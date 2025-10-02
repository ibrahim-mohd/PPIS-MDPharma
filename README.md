# Introduction
This repository provides codes and scripts to obtain protein-protein interaction stabilizers as described in our paper "". In nut-shell, we generate pharmacophore models from MD simulations and search for corresponding ligands that satisfy the pharmacophores in a local data. The scripts here automates the whole procedure.

# Dependencies
1. Pharmer (Executable at: https://sourceforge.net/projects/pharmer/files/)
2. Gromacs (https://manual.gromacs.org/documentation/)
3. AmberTools (For post processing e.g to obtain ligand parameters, conda version is fine: https://ambermd.org/GetAmber.php)
## Python packages
1. MDAnalysis (https://www.mdanalysis.org/)
2. MDtraj (https://www.mdtraj.org/1.9.8.dev0/index.html)
3. networkx (https://networkx.org/)
