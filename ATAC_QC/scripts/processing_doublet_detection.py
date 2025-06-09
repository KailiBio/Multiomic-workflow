#!/usr/bin/env python3

"""
Author: Kaili Fan
Script to perform doublet detection for snATAC-seq samples, with summary plots.

Usage:
    python scripts/processing_doublet_detection.py config/ATAC_config.yaml runtag
"""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "64"
os.environ["MKL_NUM_THREADS"] = "64"

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import load_config, standardize_tissue_name
from atac_qc.atac_plots import plot_doublet_score_probability


def main(config_path, runtag):
    config = load_config(config_path)
    workdir = config['paths']['workdir']
    h5ad_dir = config['paths']['output_h5ad_dir']
    out_fig_dir = os.path.join(config['paths']['output_figures_dir'], runtag)
    os.chdir(workdir)
    os.makedirs(out_fig_dir, exist_ok=True)

    tissue = config['params']['tissue']
    tissue2 = standardize_tissue_name(tissue)
    n_threads = config['params'].get('n_threads', 16)

    # Read sample/metadata table
    df = pd.read_csv(config['paths']['sample_metadata'], sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    working_df = df[df["tissue"] == tissue]
    sample_list = working_df['atacID'].unique().tolist()

    for i, fileID in enumerate(sample_list, 1):
        print(f"[{i}/{len(sample_list)}] Processing {fileID}...")

        h5ad_path = os.path.join(h5ad_dir, f'{fileID}.filtered.{runtag}.h5ad')
        if not os.path.exists(h5ad_path):
            print(f"  [SKIP] {fileID}: {h5ad_path} not found.")
            continue

        adata = ad.read_h5ad(h5ad_path)

        # 1. tile matrix
        snap.pp.add_tile_matrix(adata, n_jobs=n_threads, 
                                bin_size=config['params'].get('genomic_bin_size', 500))
        
        # 2. feature selection
        snap.pp.select_features(adata, n_features=config['process'].get('n_features', 250000))

        # 3. doublet detection
        snap.pp.scrublet(adata)
        
        # 4. get midpoint and plot
        doublet_probability = adata.obs['doublet_probability']
        probability_midpoint = ((doublet_probability.max() + doublet_probability.min()) / 2).round(1)
        plot_doublet_score_probability(
            adata.obs['doublet_score'], doublet_probability, probability_midpoint,
            tissue2, fileID, outdir=out_fig_dir, runtag=runtag, show=False
        )
        print(f"Default doublet probability cutoff is {probability_midpoint}")
        
        # 5. save
        proc_h5ad_path = os.path.join(h5ad_dir, f"{fileID}.processed.{runtag}.h5ad")
        adata.write(proc_h5ad_path, compression="gzip")
        print(f"  [OK] wrote {os.path.relpath(proc_h5ad_path)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ATAC doublet detection and visualization.")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("runtag", help="Tag for this run (e.g. round3 or v1)")
    args = parser.parse_args()
    main(args.config, args.runtag)