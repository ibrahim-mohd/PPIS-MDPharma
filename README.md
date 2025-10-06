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
 ## Files required

Before proceeding, ensure you have the following files prepared:

1. **Trajectory Files**  
   - The `.xtc` and `.tpr` files from the MD trajectory. The complex should be whole with no jumping across the box boundary before usage e.g you can use ``gmx trjconv`` with ``pbc mol`` option.
2. **Structure File for Pocket Identification**  
   You need a `.pdb` or `.gro` file to identify the binding pocket at the interface. There are two possible approaches:

   **Option 1 — Existing Complex**  
   - If you already have a **protein–protein–ligand complex**, you can use it directly.  
   - Any ligand is acceptable, as long as it is **docked in the pocket of interest**.  
   - Provide the `.pdb` or `.gro` file of the **entire complex** (protein–protein–ligand).

   **Option 2 — Detect Pockets Using Fpocket**  
   - If no ligand-bound structure is available, you can use **[Fpocket](https://github.com/Discngine/fpocket)** to identify potential binding pockets.  
   - Supply the **Fpocket output `.pdb` file**, which includes the protein and all detected pockets.  
   - You will also need the **residue IDs** corresponding to the pocket of interest.

3. **Reference Frame**  
   - Select a representative frame from your trajectory to use for **pharmacophore model generation**.  
   - The pharmacophore hits will correspond to this specific configuration.

## Running the commands

1. **Analyze pharamcophore features**
```bash
python 01_analyse_pharmacophore_features.py -fl $PWD/protein_out.pdb -f $PWD/mol.xtc -s $PWD/npt.tpr -pocket_id 1 -b 20000 -skip 10 -o all_features.pkl
```
Where, ``protein_out.pdb`` is the Fpocket output and we are intersted int eh pocket with resid (pocket ID) 1. Incase this file already has the ligand bound, no need to specify pocket ID. The flags are similar to gromacs conventions. Also, the default name of solvent resname is ``SOL`` and the solvent oxygen and hydrogen are `OW`, `HW1 HW2`. If that is not the case use extra flags to add the info. 

2. **Feature selection**

   Plot the above features and apply thresholds to select features. It is best to keep number of features less than 10.
```bash
python 02_feature_selection_plot.py -p all_features.pkl -c $PWD/mol.gro -s $PWD/npt.tpr -acceptor_th 30 -donor_th 40 -dG_th 0.5 -ion_th 30
```

3. **Master pharmacophore model generation**
   
   Generate pharmacophores using the above obtained thresholds
   ```bash
   python 03_generate_master_pharmacophore.py -p all.pkl -c $PWD/mol.gro -s $PWD/npt.tpr -acceptor_th 30 -donor_th 40 -dG_th 0.5 -ion_th 30 -o master_pharmacophore.json
   ```
