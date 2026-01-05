## Introduction
This repository provides codes and scripts to obtain protein-protein interaction stabilizers as described in our paper "*Protein-Protein Interaction Stabilizers from MD Simulation-derived Pharmacophores, Mohd Ibrahim and Martin Zacharias*". In nut-shell, we generate pharmacophore models from MD simulations and search for corresponding ligands that satisfy the pharmacophores in a local data. The scripts here automates the whole procedure.

**Note:**  In all the following scripts always use full path for input files. For instance, if the `npt.tpr` file is in the current directory use it as `$PWD/npt.tpr`

## Dependencies
#### Software Packages

1. **[Pharmer](https://sourceforge.net/projects/pharmer/files/)** – For pharmacophore search and database creation. Executable available via SourceForge.
2. **[GROMACS](https://manual.gromacs.org/documentation/)** – We are mainly using the ``gmx sasa`` routine for calculating the solvation free energy per residue.
3. **[AmberTools](https://ambermd.org/GetAmber.php)** – For post-processing (e.g., obtaining ligand parameters) and ff19SB or ff14SB paramters for protein. The Conda version works fine.
4. **[Fpocket](https://github.com/Discngine/fpocket)** – Use this if you do not know the binding pocket or lack a structure file of the protein–protein complex with a known bound ligand.
5. **[Openbabel](https://openbabel.org/docs/Installation/install.html#install-binaries)** – This is used near the end to obtain the total ligand charges which is required to simulate the top ligands

#### Python Packages

1. **[MDAnalysis](https://www.mdanalysis.org/)** (v-2.7.0) – For all analysis.
2. **[MDTraj](https://www.mdtraj.org/1.9.8.dev0/index.html)** (v-1.10.3)– For buried surface area analysis function, older version may **NOT** work since they take different number of arguments for the SASA function
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
You have a Pharmer compatible local database. We also **provide a datasets** with **10 Million Random** compounds, one can simply unzip it and use the above command to create a Pharmer compatabile database. Note that the compressed file is  around **6.5 GB**, upon uncompression, its size is **41 GB** and the resulting Pharmer database becomes **~134 GB** in size. So you need alteast **>190 GB** of disk space. However, using the provided scripts, ``./useful-scripts/sample_N_random_ligands.py`` you can use a random subset of the 10 M database. To sample ``N`` compounds
```bash
  python sample_N_random_ligands.py -i random_10M_ligands.sdf.gz -n 100000 -o 100k.sdf.gz
```
Download our database from [10 Million Random Compounds](https://syncandshare.lrz.de/getlink/fiREpkFfuKx4wt3LPu6u82/random_10M_ligands.sdf.gz). Note that it could take hours to create a 10 Million compound database.

## Usage
 #### Files required

Before proceeding, ensure you have the following files prepared:

1. **Trajectory Files**  
   - The `.xtc` and `.tpr` files from the MD trajectory. The complex should be whole with no jumping across the box boundary before usage e.g you can use ``gmx trjconv`` with ``pbc mol`` option. You can use the `./useful-scripts/remove_jumps_xtc.py` script to achieve that, it just takes the `tpr` and `xtc` file as input and outputs a new trajectory with the periodic jumps corrected.
2. **Pocket Identification**  
   There are three possible approaches to identify the pocket residue for analysis:

   **Option 1 — From known pocket residue IDs**
   - You need the **residue IDs of the pocket residues**.
   - Let’s say the resids are `12 34 56 199 200 234 240`. Simply use the `-pocket_resids` flag and provide the resids as space-separated strings: `-pocket_resids "12 34 56 199 200 234 240"`
     
   - **Important:** The code always renumbers protein residues from `1` to `N`, where `N` is the total number of residues in the complex. Make sure to renumber residues accordingly before specifying pocket IDs.

   **Option 2 — Existing Complex**
   - If you already have a **protein–protein–ligand complex**, you can use it directly.
   - Any ligand is acceptable, as long as it is **docked in the pocket of interest**.
   - Provide the `.pdb` or `.gro` file of the **entire complex** (protein–protein–ligand).
   - **Important:** Here the number of residues in the provided complex should be same as in the simulated protein. If the simulated protein is capped with ACE or NME residues, cap this complex first too and use the capped one. You can cap very quickly usign the script at **[Add NME/ACE Caps](https://github.com/ibrahim-mohd/Add-NME-ACE-residues-to-protein-terminal-residues)**

   **Option 3 — Detect Pockets Using Fpocket**
   - If no ligand-bound structure is available, you can use **[Fpocket](https://github.com/Discngine/fpocket)** to identify potential binding pockets.
   - Supply the **Fpocket output `*_out.pdb` file**, which includes the protein and all detected pockets.
   - You will also need the **residue ID(s)** corresponding to the pocket of interest.
   - The script `00_identify_pocket_conf.py` automates this process and selects the MD frame with the most open interface pocket and also prints the pocket ID that you need  in the next step if this method is used.

     Run:
     ```bash
     python 00_identify_pocket_conf.py -f $xtc_file -s $tpr_file -n 500 -b 20000 -on 5 -keep 0 -out_path $PWD/fpocket-confs
     ```

     Where:
     - `-n 500` → number of frames used for pocket detection
     - `-on 5` → top 5 pocket configurations selected
     - Selection is based on the **harmonic mean of buried surface area** with both protein partners

   

4. **Reference Frame**  
   - Select a representative frame from your trajectory to use for **pharmacophore model generation**.  
   - The pharmacophore hits will correspond to this specific configuration i.e all the hits are docked in this particular conformation.
   - You can also use one of the Fpocket configurations that ``00_identify_pocket_conf.py`` outputs.
   - This frame should correspond to a full MD frame with protein, water and ions etc.

### Running the commands

#### 01. Analyze pharmacophore features   
   - ##### Case-1: Known pocket residue IDs
     
     ```bash
     python 01_analyse_pharmacophore_features.py \
     -f $PWD/mol.xtc \
     -s $PWD/npt.tpr \
     -b 20000 \
     -e 100000 \
     -skip 1 \
     -pocket_resids "20 24 28 36 45 47 89 120 130 134 135 192" \
     -nname Cl- \
     -pname Na+ \
     -hsol_name "HW1 HW2" \
     -osol_name OW \
     -sol_resname SOL \
     -cutoff 5.0 \
     -o combined_analysis.pkl  
 - ##### Case-2: Known Protein-Protein-Ligand complex file (`protein_lig.gro`)

  ```bash
     python 01_analyse_pharmacophore_features.py \
     -fl $PWD/protein_lig.gro \
     -f $PWD/mol.xtc \
     -s $PWD/npt.tpr \
     -b 20000 \
     -e 100000 \
     -skip 1 \
     -nname Cl- \
     -pname Na+ \
     -hsol_name "HW1 HW2" \
     -osol_name OW \
     -sol_resname SOL \
     -cutoff 5.0 \
     -o combined_analysis.pkl
  ```
   - ##### Case-3: Fpocket output file (`protein_out.pdb`):
      Same as before but pocket ID needs to be specified (other than 0). Note that in this case sometimes we end up with too many residues, one can also exclude certain residues using the ``-res_exclude`` flag.
  
     ```bash
     python 01_analyse_pharmacophore_features.py \
     -fl $PWD/protein_out.pdb \
     -pocket_id 16\
     -f $PWD/mol.xtc \
     -s $PWD/npt.tpr \
     -b 20000 \
     -e 100000 \
     -skip 1 \
     -nname Cl- \
     -pname Na+ \
     -hsol_name "HW1 HW2" \
     -osol_name OW \
     -sol_resname SOL \
     -cutoff 5.0 \
     -o combined_analysis.pkl
     ```

#### 02. Feature selection (Optional)

Plot the above features and apply thresholds to select features. It is best to keep total number of features less than 20. Note  this part is optional, it only outputs a figure so that one has an idea of what threshold to apply. If this part is skipped, you can start with default thresholds are 35, 0.2 and 25 for H-bond donor/acceptor, hydrophobic/aromatic and negative/positive features for the next step. 
  ```bash
  python 02_feature_selection_plot.py \
    -p $PWD/combined_analysis.pkl \
    -c $PWD/mol.gro \
    -s $PWD/npt.tpr \
    -acceptor_th 30 \
    -donor_th 30 \
    -dG_th 0.2 \
    -ion_th 25 \
    -figsize 12 12 \
    -o $PWD/features.png
  ```

#### 03. Generate master pharmacophore model
   
 Generate pharmacophores using the above obtained thresholds:

```bash
python 03_generate_master_pharmacophore.py \
  -s $PWD/npt.tpr \
  -c $PWD/mol.gro \
  -p $PWD/combined_analysis.pkl \
  -dG_th 0.20 \
  -acceptor_th 30.0 \
  -donor_th 30.0 \
  -ion_th 25 \
  -nname Cl- \
  -pname Na+ \
  -hsol_name "HW1 HW2" \
  -osol_name OW \
  -sol_resname SOL \
  -hbond_direction 1 \
  -ignore_nowater_hbond 0
  -o $PWD/master_pharma.json \
```
Note that `mol.gro` is our reference frame. If `-hbond_direction` is set to 0, Hbond directions are not taken into account, one obtain more hits that is not super depedenent on the current conformation. `-ignore_nowater_hbond` if set to 1 ignores Hbond assignments for sites for which no water in Hbond geometry is found int he current frame (`mol.gro`), otherwise we use a simple translation towards pocket center for Hbond assignment. 

5. **Perform pharmacophore screening**
   
We have everything we need to perform the screening in a local database. Before proceeding make sure to have a local **pharmer compatible** database ready. Refer to the [Database creation section](#Create-local-database-of-ligands-with-pharmer). 

```bash
python generate_subpharmacophores.py \
  -j master_pharmacophore.json \
  -min_node 6 \
  -max_node 12 \
  -top 100 \
  -ntop_limit 50000 \
  -o /path/to/output_dir \
  -p_exe /usr/local/bin/pharmer.static \
  -df database_paths.txt \
  -d /path/to/database1 \
  -d /path/to/database2 \
  -np 4 \
  -e sdf \
  -max 10000 \
  -v

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
However, I recommend not mixing up results with different pharmacophore modeles with different number of features. To quickly set up simulations, refer to the ./useful/scripts/setup_simulation_protein_ligand_ff19_ff14.py script, which automates the preparation of protein–protein and protein–ligand systems.
## References

If you find this useful please cite:
*Protein-Protein Interaction Stabilizers from MD Simulation-derived Pharmacophores*, Mohd Ibrahim and Martin Zacharias (In preparation)
