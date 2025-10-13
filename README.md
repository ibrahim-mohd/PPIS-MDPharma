## Introduction
This repository provides codes and scripts to obtain protein-protein interaction stabilizers as described in our paper "". In nut-shell, we generate pharmacophore models from MD simulations and search for corresponding ligands that satisfy the pharmacophores in a local data. The scripts here automates the whole procedure.

## Dependencies
#### Software Packages

1. **[Pharmer](https://sourceforge.net/projects/pharmer/files/)** – For pharmacophore search and database creation. Executable available via SourceForge.
2. **[GROMACS](https://manual.gromacs.org/documentation/)** – We are mainly using the ``gmx sasa`` routine for calculating the solvation free energy per residue.
3. **[AmberTools](https://ambermd.org/GetAmber.php)** – For post-processing (e.g., obtaining ligand parameters). The Conda version works fine.
4. **[Fpocket](https://github.com/Discngine/fpocket)** – Use this if you do not know the binding pocket or lack a structure file of the protein–protein complex with a known bound ligand.
5. **[Openbabel](https://openbabel.org/docs/Installation/install.html#install-binaries)** – This is used near the end to obtain the total ligand charges which is required to simulate the top ligands

#### Python Packages

1. **[MDAnalysis](https://www.mdanalysis.org/)** (v-2.7.0) – For all analysis.
2. **[MDTraj](https://www.mdtraj.org/1.9.8.dev0/index.html)** (v-1.10.3)– For buried surface area analysis function.
3. **[NetworkX](https://networkx.org/)** – For visualzing and manipulating pharamcophore graphs
## Create local database of ligands with Pharmer
Refer to the [Pharmer manual](https://sourceforge.net/p/pharmer/code/HEAD/tree/)  (there seems to be just a README file). But in short: 

1. Obtain/install the Pharmer executable.
2. Download the ligands in ``MOL2`` or ``SDF`` format.
3. Concatenate all the ligands into a single file with the corresponding format i.e ``SDF`` or ``MOL2``. Lets say you named the file as ``all_ligands.sdf``
4. Run the following command. Note that pharmer will show an error if the directory ``LocalDatabase`` or whatever you wanna call it in the command below already exists.  
``` bash
pharmer.static dbcreate -dbdir=LocalDatabase -in all_ligands.sdf
```
You have a Pharmer compatible local database. We also **provide a datasets** with **10 Million Random** compounds, one can simply unzip it and use the above command to create a Pharmer compatabile database. Note that the compressed file is  around **6.5 GB**, upon uncompression, its size is **41 GB** and the resulting Pharmer database becomes **~134 GB** in size. So you need alteast **>190 GB** disk space. However, using the provided scripts, ``./useful-scripts/sample_N_random_ligands.py`` you can just use the subset of the 10 M database. To sample ``N`` compounds"
```bash
  python sample_N_random_ligands.py -i random_10M_ligands.sdf.gz -n 100000 -o 100k.sdf.gz
```

## Usage
 #### Files required

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

#### Running the commands

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
   python 03_generate_master_pharmacophore.py -p all_features.pkl -c $PWD/mol.gro -s $PWD/npt.tpr -acceptor_th 30 -donor_th 40 -dG_th 0.5 -ion_th 30 -o master_pharmacophore.json
   ```

4. **Sub-pharmacophore generation (Subsets of master pharmacophore)**

   Generate subsets of master pharmacophore as: first create an output directory named (let's say) ``pharma``

```bash
mkdir ./pharma
```
   ``` bash
   python 04_generate_sub_pharmacophores.py -j master_pharmacophore.json -min_node 4 -top 30 -ntop_limit 50 -o ./pharma
   ```
Where we generate ``max (30 %, 50)`` number of graphs/models

5. **Perform pharmacophore screening**
   
We have everything we need to perform the screening in a local database. Before proceeding make sure to have a local **pharmer compatible** database ready. Refer to the [Database creation section](#Create-local-database-of-ligands-with-pharmer). 

```bash
python 05_perform_screening.py -d $database_path -i $PWD/pharma -o $PWD/search-output -max 10000 -np 4
```
Where, we provided the path to database, the pharmacophores graphs path (i.e the folder with n* like folder e.g n8 n7 n3 etc), the output folder (it creates if not there) and the number of concurrent workers (parallel threads). Apart from the hits, a screening summary file ``screening_summary.dat`` is also produced.

6. **Scoring hits by calculating buried surface area**
   
 The above ligands are scored by calculating the buried surface area with each protein partner.
```bash
python 06_score_hits.py -sdf $sdf_path -s $tpr_file -c $gro_file -o ligand_scores.pkl
```
Where $sdf_path is the path to pharmer output ``SDF`` files. For different pharmacophore model e.g n7, n6 or n5, one creates seperate scoring files. The code will go through all the ``.sdf`` file in the folder. The ``TPR`` and ``GRO`` files must also be provided with the ``GRO`` file being the MD frame for which the pharmacophores hits were obtained.

7. **Select top hits across all models and prepare for simulation**
```bash
python 07_extract_top_ligands.py -i ligand_scores.pkl -topN 40 -Ngraph 20 -odir $PWD/final-results
```
or if you have multiple score files from different pharmacophore models:

```bash
python 07_extract_top_ligands.py -i ligand_scores1.pkl ligand_scores2.pkl ligand_scores3.pkl -topN 40 -Ngraph 20 -odir $PWD/final-results
```
However, I recommend not mixing up results with different pharmacophore modeles with different number of features
