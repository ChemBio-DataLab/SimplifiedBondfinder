# SimplifiedBondfinder

<img align="left" width="100" height="100" src="https://github.com/ChemBio-DataLab/SimplifiedBondfinder/blob/main/Logo.png">

**simplifiedBondfinder** is an automated framework for classifying and discovering various bonds and interactions in protein structures.

**Versatile algorithm for protein bond discovery: NOS linkages as a case study**\
*Sophia Bazzi, Sharareh Sayyad*, [arXiv](https://arxiv.org).



## Workflow

Here is the workflow of the `simplifiedBondfinder`
![Workflow](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/blob/main/NOS_overview.png)


 ## Structure of the repository

It comprises data acquisition, feature extraction, and machine-learning analysis. 
The code for data acquisition is stored in the [src_data_acquisition](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/src_data_acquisition) directory.

The machine-learning parts and the scripts for generating figures in the manuscript are stored in [SimplifiedBondfinder_github/Figure_scripts](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Figure_scripts)

The Data used for plotting all figures of our manuscript is in the [Data](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Data) directory. Figures for LYS-NOS-CYS, GLY-NOS-CYS, ARG-N$`_{\eta}`$OS-CYS, and ARG-N$`_{\varepsilon}`$OS-CYS linkages are ploted using data stored in [LYS_CYS](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Data/LYS_CYS), [GLY_CYS](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Data/GLY_CYS), [ARG_N_CYS](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Data/Arg_N_CYS), and [ARG_NE_CYS](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Data/ARG_NE_CYS) folders, respectively. In each of these directories, UMAP and Pair-plot folders are present. Files inside the UMAP folder are two raw data files and a "clusters_with_common_filenames.txt". The latter, which is the output of the UMAP method cells in the "Figure_scripts/Figs_paper_\*.ipynb" scripts, collects the file names for different samples in different clusters of UMAP.
 The data for pairplot figures are stored in the "Pair-plot" directory. Generating these files can be done by uncommenting corresponding cells in "Figure_scripts/Figs_paper_\*.ipynb" scripts.



## Protein Data Acquisition Tool

This tool extracts and analyzes structural atomic and residue features, solvent accessible surface area (SASA), and B-factor values from protein structures, utilizing data from the following sources:

- PDB-REDO: Updated and optimized crystallographic structures.
- BDB: Databank of PDB files with consistent B-factors.
- RCSB PDB: Full Reports, with a focus on sulfur-nitrogen interactions.

The tool specializes in identifying specific atomic interactions/bonds, particularly sulfur-nitrogen interactions, and generates a comprehensive dataset for further analysis in structural biology and bioinformatics research.

Additionally, the tool is highly adaptable, allowing users to focus on any type of atomic interactions present in protein structures and to analyze protein structures from the RCSB PDB.

### Features

- Processes protein structures from CIF files pre-downloaded from PDB-REDO using Biopython.
- Analyzes sulfur-nitrogen interactions.
- Computes distances, angles, and torsion angles between atoms.
- Calculates Solvent Accessible Surface Area (SASA) for atoms and residues.
- Downloads and processes data from Protein Data Bank (RCSB PDB) and BDB.
- Identifies real-space R-value Z-score (RSRZ) outliers and "Too-close contacts".

### Main Functions

1. `process_structure`: Main entry point for processing each protein structure.
2. `process_atom` and `process_neighbor_atom`: Extract relevant information from our two main atoms.
3. `process_CN_bond`: Process the carbon atoms adjacent to the nitrogen atom to identify a CN bond.
4. `process_CS_bond`: Process the carbon atoms adjacent to the sulfur atom to identify a CS bond.
5. `structural_features_neighbour_residues`: Computes geometric features for sulfur and nitrogen atoms and identifies the corresponding adjacent residues:
   - Calculates distances between C-S, C-N, and S-N atoms.
   - Computes C-S-N and C-N-S angles.
   - Calculates C-S-N-C torsion angle.
   - Identifies residues within a 4-angstrom vicinity of the alpha carbon atoms of nitrogen- and sulfur-containing residues
6. `download_full_report`: downloads the Full Reports from the RSCB PDB.
7. `RSRZ`: Identify RSRZ outliers from Full Reports.
8. `sasa`: Computes the Solvent Accessible Surface Area (SASA).
   - Calculating the SASA for each sulfur and nitrogen atoms.
   - Determining the SASA for the residues containing these atoms.
9. `download_BDB`: downloads files from BDB.
10.`process_BDB`: processes the BDB files to extract the B-factors of the sulfur and nitrogen atoms.  
11. `too_close_contacts`: collects the close contacts between the sulfur and nitrogen atoms from Full_reports.
12. `writing`: writes output results to a file.

## Required Libraries

-biopython==1.7
-numpy
-pathlib
-requests
- PyPDF2==2.10.5

## Setting Up the Environment

1. This code has been tested with Python 3.6. Ensure you have Python 3.6 installed on your system. 
2. It's recommended to create a virtual environment to avoid conflicts with other projects:
4. Install the required libraries mentioned above.
5. Verify your Biopython version

## Usage

- By default, CIF files should be placed in the current working directory. Alternatively, one can specify a different directory for the CIF files using: D = Path('your_directory_path'). 
- To initiate data acquisition, one should type, e.g., 

```
python protein_data_acqusition.py
```

 *Summary of Code Performance:*
 The script processed the CIF file 6pgb_final.cif in 4.234 seconds. 
 The test was conducted in a single-threaded manner on a high-performance Linux node with Kernel Version 3.10.0-1160.25.1.el7.x86_64. 
 At the time of execution, nearly 13.33 GB of memory was available.

# Automated Feature Space Exploration

This project performs clustering analysis on protein structure data, focusing on cysteine and lysine residues. It uses automated descriptor set generation and selection to identify patterns in protein structures.

## Analysis Pipeline

Our clustering analysis consists of two main stages:



1. [Automated Descriptor Set Generation](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/src_feature_space/Automated_Feature_Space_Production.py): Generates all possible sets of descriptors for protein structures.
2. [Descriptor Set Selection](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/src_feature_space/Automated_Feature_Space_selection.py): Identifies the most promising descriptor sets based on clustering performance.

Following these automated steps, an expert review is required to determine the optimal, minimal descriptor set for further analysis.

## Installation

1. Clone this repository:

```
git clone https://github.com/ChemBio-DataLab/SimplifiedBondfinder.git
```

2. Install the required packages:
 
pip install pandas umap-learn scikit-learn numpy biopython==1.7 numpy pathlib requests PyPDF2==1.26.0

## Usage

### Stage 1: Automated Descriptor Set Generation

1. Prepare your input data:
- Place a tab-separated file named `CYS_non_N_CYS_LYS_No_back.txt` in the project directory.
- Ensure your reference PDB filenames are listed in the `pdb_filenames` variable in the script.

2. Run the script, e.g., try :
   
```
python Automated_Feature_Space_Production.py
```

3. Output:
- `clusters_with_common_filenames.txt`: Detailed clustering results
- `clusters_summary.txt`: Summary of clustering results
- `CYS_non_N_CYS_LYS_No_back_hashes.csv`: Preprocessed dataset with hash columns

### Stage 2: Descriptor Set Selection

1. Ensure you have generated the `clusters_summary.txt` file from Stage 1.

2. Run the selection script, e.g., try:

```
python Automated_Feature_Space_selection.py
```

3. Output:
- `selected_features.txt`: Contains selected feature sets meeting the criteria

4. An extensive review of `selected_features.txt` should be done to identify the smallest feature space that meets the criteria.

## How It Works

1. Data Preparation and Cleaning:
- Select data containing side chain nitrogen of lysine residue (NZ) and sulfur atom of cysteine residue (SG).
- Remove proteins with None values for BDB B-factors, if included.

2. Descriptor Set Generation:
- Generates all possible sets of 3 to 15 descriptors.
- For each set:
  - Applies UMAP for dimensionality reduction
  - Uses MeanShift for clustering
  - Calculates Silhouette Score and Davies–Bouldin Index
  - Compares clustered structures with reference PDB structures

3. Descriptor Set Selection:
- Identifies sets meeting criteria:
  - Silhouette score ≥ 0.5
  - One cluster with exactly 80 common points with PDB filenames

## Customization

- Adjust the `feature_space` list in the generation script.
- Modify UMAP parameters (`min_dist` and `n_neighbors`).
- Change clustering algorithm or parameters.
- Adjust selection criteria in the selection script.


## Advanced ML Analysis and Visualization

After completing the automated descriptor selection process, you can perform advanced analyses, including machine learning and clustering, and generate visualizations using the final set of descriptors. We provide several Jupyter notebooks for this purpose:

|NOS Linkage | ipython notebook | 
|---------|-----------------|
|LYS-NOS-CYS | [Figs_paper_LYS_CYS.ipynb](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Figure_scripts/Figs_paper_LYS_CYS.ipynb)| 
|GLY-NOS-CYS | [Figs_paper_GLY-CYS.ipynb](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Figure_scripts/Figs_paper_GLY-CYS.ipynb) | 
|ARG-N$`_{\eta}`$OS-CYS | [Figs_paper_ARG_N_CYS.ipynb](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Figure_scripts/Figs_paper_ARG_N_CYS.ipynb) | 
|ARG-N$`_{\varepsilon}`$OS-CYS |[Figs_paper_ARG-NE-CYS.ipynb](https://github.com/ChemBio-DataLab/SimplifiedBondfinder/tree/main/Figure_scripts/Figs_paper_ARG-NE-CYS.ipynb)| 


# Citation

<img align="left" width="100" height="100" src="https://github.com/ChemBio-DataLab/SimplifiedBondfinder/blob/main/Logo.png"> If you find this code useful in your work, please cite our article
``**Versatile algorithm for protein bond discovery: NOS linkages as a case study**''
Sophia Bazzi, Sharareh Sayyad, [arXiv](https://chemrxiv.org/engage/chemrxiv/article-details/66e92a7451558a15ef381a5e).


## License
This project is licensed under the terms of the BSD 3-Clause license.


