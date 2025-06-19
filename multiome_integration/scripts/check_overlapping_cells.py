#!/usr/bin/env python3

"""
Author: Kaili Fan
Check overlapping cells between scRNA and scATAC for all samples in a tissue and generate Venn/bar plots.
"""

import os
import sys
import time
import argparse
import pandas as pd
import anndata as ad

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name, print_elapsed_time
from multiome_integration.integration_plots import plot_venn_and_save, plot_overlap_bar

from glob import glob

def find_h5ad(directory, suffix, sample, mode='rna'):
    """
    Find h5ad file by {suffix}_{sample}*.h5ad pattern in directory.
    """
    pattern = f"{suffix}_{sample}*.h5ad"
    files = glob(os.path.join(directory, pattern))
    if not files:
        raise FileNotFoundError(f"No {mode} h5ad file matching {pattern} in {directory}")
    if len(files) > 1:
        print(f"[Warning] Multiple {mode} h5ad files for {sample} in {directory}, using {files[0]}")
    return files[0]

def process_overlap(sample, rna_dir, atac_dir, suffix, tissue_std, figures_dir):
    try:
        rna_path = find_h5ad(rna_dir, suffix, sample, 'rna')
        atac_path = find_h5ad(atac_dir, suffix, sample, 'atac')
    except FileNotFoundError as e:
        print(f"[Warning] {e} (Skipping sample: {sample})")
        return None

    try:
        rna = ad.read_h5ad(rna_path)
        atac = ad.read_h5ad(atac_path)
    except Exception as e:
        print(f"[Error] Could not read h5ad for sample {sample}: {e}")
        return None

    rna_cells = rna.obs_names
    atac_cells = atac.obs_names
    shared_cells = set(rna_cells) & set(atac_cells)

    rna_only = len(set(rna_cells) - shared_cells)
    overlap = len(shared_cells)
    atac_only = len(set(atac_cells) - shared_cells)
    plot_venn_and_save(rna_cells, atac_cells, sample, figures_dir, tissue_std)
    return [rna_only, overlap, atac_only]

def main(config_path):
    config = load_config(config_path)
    rna_dir = config['paths']['rna_h5ad_dir']
    atac_dir = config['paths']['atac_h5ad_dir']
    figures_dir = config['paths']['output_figures_dir']
    os.makedirs(figures_dir, exist_ok=True)

    df = pd.read_csv(config['paths']['sample_metadata'], sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])

    tissue = config['params']['tissue']
    suffix = config['params']['suffix']
    
    if tissue == "---":
        working_df = df
    else:
        working_df = df[df["tissue"] == tissue]
        
    sample_list = working_df['atacID'].unique()
    print(f"Samples: {list(sample_list)}")
    sample_tissue_dict = dict(zip(working_df['atacID'], working_df['tissue']))
    tissues = sorted(working_df["tissue"].unique())
    print(f"Working tissue: {', '.join(tissues)}")

    cell_counts = {}
    for idx, sample in enumerate(sample_list, 1):
        print(f"[{idx}/{len(sample_list)}] Processing sample: {sample}...")
        res = process_overlap(sample, rna_dir, atac_dir, suffix, tissue_std, figures_dir)
        if res:
            cell_counts[sample] = res

    if cell_counts:
        plot_overlap_bar(cell_counts, figures_dir, tissue_std)
        print(f"Overlap statistics complete. Figures saved in {figures_dir}")
    else:
        print("No overlaps found or no valid samples to plot.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check overlapping cells between scRNA and scATAC for all samples in a tissue and plot Venn/bar graphs."
    )
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    main(args.config)