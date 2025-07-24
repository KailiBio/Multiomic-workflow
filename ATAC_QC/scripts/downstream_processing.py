#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Final clustering, doublet removal, and QC figure/statistics for snATAC-seq samples.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import anndata as ad
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import load_config, standardize_tissue_name
from atac_qc.atac_plots import cell_count_post_filter_hist, plot_per_sample_umap_clusters

def summarize_and_plot_cell_counts(sample_list, h5ad_dir, fig_dir, runtag, suffix):
    """
    Summarize cell counts and plot histogram for a set of h5ad files.

    Args:
        sample_list (list): List of sample IDs / fileIDs.
        h5ad_dir (str): Directory containing processed h5ad files.
        fig_dir (str): Output directory for figure and table.
        runtag (str): Tag used in h5ad and output file names.
        suffix (str): Suffix for output naming).
    Returns:
        df_num_cells (pd.DataFrame): DataFrame of sampleID and cell count.
    """
    numCells = []
    for fileID in sample_list:
        h5ad_path = os.path.join(h5ad_dir, f'{fileID}.final.{runtag}.h5ad')

        if not os.path.exists(h5ad_path):
            continue
        adata = ad.read_h5ad(h5ad_path)

        numCells.append((fileID, len(adata.obs_names)))
        
    df_num_cells = pd.DataFrame(numCells, columns=["fileID", "numCells"])
    df_num_cells.to_csv(os.path.join(fig_dir, f"ATAC_NumCell.{runtag}.{suffix}.tsv"), sep="\t", index=False)

    cell_count_post_filter_hist(df_num_cells, suffix, runtag, fig_dir)
    
def main(config_path, runtag):
    config = load_config(config_path)
    
    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    outdir = os.path.join(workdir, f'doublet_filter_processing.{runtag}')
    os.makedirs(outdir, exist_ok=True)
    os.chdir(workdir)

    tissue = config['params']['tissue']
    suffix = config['params']['suffix']

    n_threads = config['params'].get('n_threads', 16)

    # Load sample info
    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])

    # Load qc cutoff table
    qc_cutoff_table = config['qc']['atac_qc_cutoff_table']
    if qc_cutoff_table.endswith('.xlsx') or qc_cutoff_table.endswith('.xls'):
        df_cutoff_all = pd.read_excel(qc_cutoff_table,
                                      sheet_name=config['qc']['sheet_name'], engine='openpyxl')
    else:
        df_cutoff_all = pd.read_csv(qc_cutoff_table, sep='\t')
    
    if tissue == "---":
        working_df = df
        df_cutoff = df_cutoff_all 
    else:
        working_df = df[df["tissue"] == tissue]
        df_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue]
        
    tissues = sorted(working_df["tissue"].unique())
    print(f"Working tissue: {', '.join(tissues)}")
    sample_list = working_df['atacID'].unique()
    print(f"Samples: {list(sample_list)}")

    df_cutoff.set_index('atacID', inplace=True)
    sample_tissue_dict = dict(zip(working_df['atacID'], working_df['tissue']))

    # Remove doublets, embedding, clustering, and save back to disk
    for i, fileID in enumerate(sample_list, 1):
        print(f"\n========== Processing {fileID} ({i}/{len(sample_list)}) ==========")
        
        tissue_std = standardize_tissue_name(sample_tissue_dict[fileID])
          
        h5ad_path = os.path.join(output_h5ad_dir, f'{fileID}.processed.{runtag}.h5ad')
        if not os.path.exists(h5ad_path):
            print(f"[ERROR] No h5ad for {fileID} at {h5ad_path}.")
            continue

        print("[INFO] Loading anndata object...")
        adata = ad.read_h5ad(h5ad_path)

        # Remove doublets
        print("[INFO] Filter doublet...")
        doublet_cutoff = df_cutoff.loc[fileID, 'doublet_cutoff']
        df_doublet = adata.obs[['doublet_score', 'doublet_probability']].copy()
        if str(df_cutoff.loc[fileID, 'use_double_probability_filter']) == 'Yes':
            print(f'  using double probability filter: {doublet_cutoff}')
            df_doublet['doublet_call'] = df_doublet['doublet_probability'].apply(lambda x: 'yes' if x > doublet_cutoff else 'no')
            snap.pp.filter_doublets(adata, n_jobs=n_threads, probability_threshold=doublet_cutoff)
        else:
            print(f'  using double score filter: {doublet_cutoff}')
            df_doublet['doublet_call'] = df_doublet['doublet_score'].apply(lambda x: 'yes' if x > doublet_cutoff else 'no')
            snap.pp.filter_doublets(adata, n_jobs=n_threads, 
                                    score_threshold=doublet_cutoff, probability_threshold=None)
            
        df_doublet.to_csv(f'ATAC_doublet_results.{fileID}.tsv', sep='\t', index=True, header=True, 
                          index_label="cell_barcode")

        
        # Dimension reduction and clustering
        print("[INFO] Dimensional reduction and clustering...")
        snap.tl.spectral(adata)
        snap.tl.umap(adata, random_state=0)
        snap.pp.knn(adata)
        snap.tl.leiden(adata)

        # Save updated AnnData
        h5ad_out_path = os.path.join(output_h5ad_dir, f'{fileID}.final.{runtag}.h5ad')
        adata.write(h5ad_out_path)
        print(f"[DONE] wrote {os.path.relpath(h5ad_out_path)}")

    # Summarize and plot cell stats
    print("[INFO] Generating summary table and histogram for filtered cell counts per sample...")
    summarize_and_plot_cell_counts(sample_list, output_h5ad_dir, outdir, runtag, suffix)

    print("[INFO] Generating per-sample UMAP cluster plots...")
    plot_per_sample_umap_clusters(sample_ids=sample_list, h5ad_dir=output_h5ad_dir, run_tag=runtag,
                                  suffix=suffix, output_dir=outdir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ATAC downstream processing: doublet removal, clustering & statistics.")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("runtag", help="Tag for this run (e.g. round3 or v1)")
    args = parser.parse_args()
    main(args.config, args.runtag)