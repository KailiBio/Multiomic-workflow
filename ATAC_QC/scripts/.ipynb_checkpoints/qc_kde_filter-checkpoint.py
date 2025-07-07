#!/usr/bin/env python3

"""
Author: Kaili Fan
Description:
    Cell filtering and KDE plots for snATAC-seq samples.

Usage:
    python scripts/run_atac_qc_step1_filter_and_kde.py config/ATAC_config.yaml runtag
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import anndata as ad
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import load_config, standardize_tissue_name
from atac_qc.atac_plots import plot_kde_filter


def main(config_path, runtag):
    config = load_config(config_path)
    workdir = config['paths']['workdir']
    os.chdir(workdir)
    figdir = os.path.join(config['paths']['output_figures_dir'], runtag)
    os.makedirs(figdir, exist_ok=True)

    tissue = config['params']['tissue']
    suffix = config['params']['suffix']

    n_threads = config['params'].get('n_threads', 16)

    # Load sample info and QC cutoffs
    df = pd.read_csv(config['paths']['sample_metadata'], sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])

    # Load qc cutoff table
    if config['qc']['atac_qc_cutoff_table'].endswith('.xlsx') or config['qc']['atac_qc_cutoff_table'].endswith('.xls'):
        df_cutoff_all = pd.read_excel(config['qc']['atac_qc_cutoff_table'],
                                      sheet_name=config['qc']['sheet_name'], engine='openpyxl')
    else:
        df_cutoff_all = pd.read_csv(config['qc']['atac_qc_cutoff_table'], sep='\t')

    
    if tissue == "---":
        working_df = df
        df_cutoff = df_cutoff_all 
    else:
        working_df = df[df["tissue"] == tissue]
        df_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue]

    tissues = sorted(working_df["tissue"].unique())
    print(f"Working tissue: {', '.join(tissues)}")
    sample_list = working_df['atacID'].unique().tolist()
    df_cutoff.set_index('atacID', inplace=True)
    sample_tissue_dict = dict(zip(working_df['atacID'], working_df['tissue']))

    before_pdf = os.path.join(figdir, f'ATAC_beforeFilter_kde.{suffix}.{runtag}.pdf')
    after_pdf = os.path.join(figdir, f'ATAC_postFilter_kde.{suffix}.{runtag}.pdf')

    with PdfPages(before_pdf) as pdf_before, PdfPages(after_pdf) as pdf_after:
        for i, fileID in enumerate(sample_list, 1):
            print(f'[{i}/{len(sample_list)}] Filtering {fileID}...')
            tissue_std = standardize_tissue_name(sample_tissue_dict[fileID])
            
            try:
                h5ad_path = os.path.join(workdir, "atac_h5ad", f'{fileID}.raw.h5ad')
                if not os.path.exists(h5ad_path):
                    print(f"  [SKIP] {fileID}: {h5ad_path} not found.")
                    continue
                    
                adata = ad.read_h5ad(h5ad_path)
                
                initial_cell_str = f"Initial cell barcodes: {len(adata.obs_names)}"
                x_cutoff = df_cutoff.loc[fileID, "num_fragment"]
                y_cutoff = df_cutoff.loc[fileID, "TSS_enrichment_score"]
                cutoff_str = f"QC cutoffs: n_fragment > {x_cutoff}, TSS_enrichment > {y_cutoff}"
                passed_mask = snap.pp.filter_cells(
                    adata, min_tsse=y_cutoff, min_counts=x_cutoff,
                    max_counts=100000, inplace=False, n_jobs=n_threads
                )
                passed_cells_str = f"Cells passing QC: {np.sum(passed_mask)}"

                # Pre-filter plot
                snap.pp.filter_cells(adata, min_tsse=3, min_counts=100, max_counts=100000,
                                     inplace=True, n_jobs=n_threads)
                plot_kde_filter(
                    adata, x_cutoff, y_cutoff, f'{tissue_std}: {fileID}', 
                    initial_cell_str, cutoff_str, passed_cells_str, pdf_before, show_cutoff_line=True
                )

                # filter & post-filter plot
                snap.pp.filter_cells(adata, min_tsse=y_cutoff, min_counts=x_cutoff,
                                     max_counts=100000, inplace=True, n_jobs=n_threads)
                plot_kde_filter(
                    adata, x_cutoff, y_cutoff, f'{tissue_std}: {fileID}', 
                    initial_cell_str, cutoff_str, passed_cells_str, pdf_after, show_cutoff_line=False
                )

                out_path = os.path.join(config['paths']['output_h5ad_dir'],
                                       f'{fileID}.filtered.{runtag}.h5ad')
                adata.write(out_path, compression="gzip")
                print(f"  [OK] wrote {os.path.relpath(out_path)} with {adata.n_obs} cells\n")
                
            except Exception as e:
                print(f"  [ERROR] {fileID}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ATAC cell filtering and KDE plotting.")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("runtag", help="Tag for this run (e.g. round3 or v1)")
    args = parser.parse_args()
    main(args.config, args.runtag)