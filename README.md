# Multiomic Workflow

This repository contains a pipeline for processing single-cell multiomic data (RNA+ATAC). The workflows are implemented primarily in Jupyter Notebooks and Python scripts.

## Table of Contents
- [Overview](#overview)
- [Getting Started](#getting-started)
    - [Installation](#installation)
    - [Required Libraries](#required-libraries)
- [Usage](#usage)
    - [Running the Master Scripts](#running-the-master-scripts)
    - [Setting up your Data and Configurations](#setting-up-your-data-and-configurations)
    - [QC Cutoff Table and Metadata Table](#qc-cutoff-table-and-metadata-table)
    - [Directory Structure](#directory-structure)
- [Notes](#notes)
- [Troubleshooting](#troubleshooting)
- [Contact](#contact)

---

## Overview

This workflow facilitates the analysis of single-cell multiomic datasets (RNA + ATAC). It includes modules for preprocessing, QC, analysis, and visualization.

Key master scripts:
- [Run RNA QC Workflow on Jupyter](RNA_QC/notebooks/run_rna_workflow.ipynb)
- [Run RNA QC Workflow as script](RNA_QC/notebooks/run_rna_workflow.py)
- [Run ATAC QC Workflow on Jupyter](ATAC_QC/notebooks/ATAC_QC_demo.ipynb)
- [Run ATAC QC Workflow as script](ATAC_QC/notebooks/run_atac_qc_workflow.sh)
- [RUN cell annotation Workflow on Jupyter](cell_annotation/notebooks/cell_annotation_demo.ipynb)
- [RUN multiome data integration as Jupyter](multiome_integration/notebook/multiome_integration_demo.ipynb)
- [RUN multiome data integration as script](multiome_integration/notebook/run_multiome_integration_workflow.sh)
- Please see the Jupyter Notebooks for additional step-by-step analyses.

---

## Getting Started

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/KailiBio/Multiomic-workflow.git
   cd Multiomic-workflow
   ```

2. **Install required libraries:**  
   All required dependencies are listed in [`requirements.txt`](requirements.txt).  
   Install them using:
   ```bash
   pip install -r requirements.txt
   ```

   **Note:** In addition to commonly used libraries, you must install [`harmonypy`](https://github.com/dpeerlab/harmonypy):
   ```bash
   pip install harmonypy
   ```
   If you encounter errors for missing packages (e.g. `upsetplot`, `scikit-image`, `igraph`), ensure these are installed as well.

### Required Libraries

- Python ≥ 3.8
- Jupyter Notebook
- harmonypy
- upsetplot
- igraph
- scikit-image
- scanpy, anndata, pandas (and other packages as listed in `requirements.txt`)

---

## Usage

### Running the Master Scripts

The principal entry point for running RNA workflow is:
- `run_rna_workflow.ipynb` (or `run_rna_workflow.py`)
- Launch a notebook via `jupyter notebook run_rna_workflow.ipynb`

This script coordinates data loading, QC, analysis and results export.

### Setting up your Data and Configurations

1. **Input data folder:**  
   For each sample, create directories structured as follows (files from CellRanger):
   ```
   input_data/
     └── sample_id/
           ├── filtered_feature_bc_matrix.h5
           └── metrics_summary.csv
   ```
   - Place your `filtered_feature_matrix` and `metrics_summary.csv` in each respective `sample_id` folder. Also put your .pct_exonic.tsv file there, if you have it.
   - The `config` file specifies data locations and analysis settings.

2. **Configuration file:**  
   Refer to the provided template or example config for correct formatting.  
   Document all sample IDs and paths accurately to avoid errors.

### QC Cutoff Table and Metadata Table

- The pipeline requires a **QC cutoff table** and **metadata table** for each sample.
- Both must be formatted according to the pipeline expectations (see example tables in the repo).
- Double check for typos and formatting issues—invalid entries will cause errors (error messages may not always specify the location of the typo).

### Directory Structure

- Main workflow scripts and notebooks are in the root directory.
- Results (e.g., summary statistics and figures) are initially saved to `/figures`; note that some outputs may be moved after workflow completion.

---

## Notes

- If you do **not** need the `scrinvex` functionality, you may comment out or skip the `scrinvex` directory.
- Important: If you encounter errors, especially relating to missing Python packages, refer to the installation steps and requirements above.
- The location where `summary_stats` is saved may change during workflow execution; check `/figures` and final output directories.

---

## Troubleshooting

- **Missing Package Errors:**  
  If you see ImportErrors for e.g. `upsetplot`, `igraph`, `scikit-image`, manually install as described above.
- **Configuration/Path Errors:**  
  Most failures are due to typos in config, cutoff, or metadata tables. Double check correctly spelled sample IDs and valid paths.
- **Figure Location:**  
  If you don’t find your summary stats, check both `/figures` and the output directories specified in the config.

---

## Contact

For help, bug reports or feature requests, please [open an issue](https://github.com/KailiBio/Multiomic-workflow/issues) on GitHub.

---

**Contributions** and suggestions to improve documentation or code are welcome!

