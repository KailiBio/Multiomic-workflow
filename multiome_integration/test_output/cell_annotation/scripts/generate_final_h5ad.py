#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Update and clean h5ad, generate final h5ad and summary plots/statistics.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from cell_annotation.utils import load_config, standardize_tissue_name, move_figures_to_newdir

def update_annotations(adata, cell_lineage, celltype_broad, celltype_fine):
    if 'leiden_new' in adata.obs_names:
        adata.obs['leiden'] = adata.obs['leiden_new']
    adata.obs['cell_lineage'] = adata.obs[cell_lineage]
    adata.obs['celltype_broad'] = adata.obs[celltype_broad]
    adata.obs['celltype_fine'] = adata.obs[celltype_fine]
    return adata

def plot_cell_counts(adata, tissue_std, figdir):
    donor_col = adata.uns.get('donorID_colors', {})

    df = adata.obs[['donorID', 'leiden', 'cell_lineage', 'celltype_broad']]
    
    # save the cell counts
    cell_count_matrix = df.pivot_table(index='celltype_broad', columns='donorID', aggfunc='size', fill_value=0)
    cell_count_matrix.to_csv(os.path.join(figdir, f'{tissue_std}_final_cell_counts.txt'), sep='\t')

    # number of cells, colored by donorID
    sns.set(style='whitegrid')
    ax = cell_count_matrix.plot(
        kind='barh', stacked=True, color=donor_col, figsize=(8, 5), edgecolor='none'
    )
    ax.set_xlabel('Number of cells')
    ax.set_ylabel('Cell type')
    ax.legend(
        title='DonorID',
        bbox_to_anchor=(1.05, 1),    # Move legend outside
        loc='upper left'
    )
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, f"Number_cell_Bar.colorByDonor.{tissue_std}.png"), 
                dpi=300, bbox_inches='tight')
    plt.close()


def save_final_h5ad(adata, workdir, tissue_std):
    #save the whole .h5ad
    adata.write(os.path.join(workdir, f"{tissue_std}.GEX.all.h5ad"))

    # generate a clean final version of .h5ad
    adata_final = adata.copy()
    # Clean obs/var/uns
    obs_keep = ['tissue', 'sampleID', 'donorID', 'cellbarcode', 
                'total_counts', 'n_genes', 'n_genes_by_counts', 'log10_total_counts', 'log10_n_genes_by_counts',
                'total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 
                'pct_counts_hb', 'pct_exon_reads',  'MALAT1_CPM', 'log10_MALAT1_CPM', 
                'doublet_score', 'predicted_doublet', 'doublet_probabilities', 
                'leiden', 'cell_lineage', 'celltype_broad', 'celltype_fine']
    adata_final.obs = adata_final.obs[obs_keep]
    var_keep = ['gene_ids', 'total_counts', 'n_cells', 'n_cells_by_counts', 'mean_counts',
                'feature_types', 'genome', 'mt', 'ribo', 'hb', 'pct_dropout_by_counts',
                'highly_variable', 'means', 'dispersions', 'dispersions_norm']
    adata_final.var = adata_final.var[var_keep]
    uns_keep = ['dendrogram_leiden', 'donorID_colors', 'hvg', 'leiden', 'leiden_colors',
                'log1p', 'neighbors', 'pca', 'rank_genes_groups', 'umap',]
    adata_final.uns = {k: adata_final.uns[k] for k in uns_keep if k in adata_final.uns}
    # save the final .h5ad
    final_h5ad_path = os.path.join(workdir, f"{tissue_std}.GEX.final.h5ad")
    adata_final.write(final_h5ad_path)
    print(f"[INFO] Final h5ad saved to: {final_h5ad_path}")

def run_per_tissue(workdir, tissue, cell_lineage, celltype_broad, celltype_fine):
    
    tissue_std = standardize_tissue_name(tissue)

    # load h5ad
    print("[INFO] Loading anndata object...")
    adata_path1 = os.path.join(workdir, "cell_annotation_auto", f"{tissue_std}_GEX.filtered.processes.autoAnnotated.updated.h5ad")
    adata_path2 = os.path.join(workdir, "cell_annotation_auto", f"{tissue_std}_GEX.filtered.processes.autoAnnotated.h5ad")
    if os.path.exists(adata_path1):
        print(f"[INFO] loading adata from {adata_path1}")
        adata = sc.read_h5ad(adata_path1)
    elif os.path.exists(adata_path2):
        print(f"[INFO] loading adata from {adata_path2}")
        adata = sc.read_h5ad(adata_path2)
    else:
        print(f"[ERROR] AnnData file not found for {tissue_std}.")
        return

    figdir = os.path.join(workdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    os.chdir(workdir)

    adata = update_annotations(adata, cell_lineage, celltype_broad, celltype_fine)

    # plot final umap
    sc.pl.umap(adata, color='leiden', frameon = False, show=False, save=f'.RNA_final.leiden.{tissue_std}.png')
    sc.pl.umap(adata, color='cell_lineage', frameon = False, show=False, save=f'.RNA_final.cell_lineage.{tissue_std}.png')
    sc.pl.umap(adata, color='celltype_broad', frameon = False, show=False, save=f'.RNA_final.celltype_broad.{tissue_std}.png')
    sc.pl.umap(adata, color='celltype_fine', frameon = False, show=False, save=f'.RNA_final.celltype_fine.{tissue_std}.png')

    plot_cell_counts(adata, tissue_std, figdir)

    print("[INFO] Save files and figures...")
    # Save
    save_final_h5ad(adata, figdir, tissue_std)
    # Organize figures
    move_figures_to_newdir(workdir, old="figures", new="final")

    print(f"[INFO] Completed generating final RNA h5ad for {tissue}.")

def main(config_path, cell_lineage, celltype_broad, celltype_fine):
    config = load_config(config_path)
    
    workdir = config['paths']['workdir']
    
    tissue = config['params']['tissue']

    if tissue == "---":
        sample_metadata = config['paths']['sample_metadata']
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
            names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        
        tissues = sorted(df["tissue"].unique())
        print(f"[INFO] Processing MULTIPLE tissues: {tissues}")
        
        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n========== Final annotation for tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")
            try:
                run_per_tissue(workdir, tissue_name, cell_lineage, celltype_broad, celltype_fine)
            except Exception as e:
                print(f"[ERROR] Final h5ad generation failed for {tissue_name}: {e}")
    else:
        print(f"\n========== Final annotation for tissue: {tissue} ==========")
        try:
            run_per_tissue(workdir, tissue, cell_lineage, celltype_broad, celltype_fine)
        except Exception as e:
            print(f"[ERROR] Final h5ad generation failed for {tissue}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate final cell-annotated h5ad (manual correction) and statistics.")
    parser.add_argument("config", help="YAML config describing tissue, dirs, etc.")
    parser.add_argument("cell_lineage", help="obs name for final cell lineage")
    parser.add_argument("celltype_broad", help="obs name for final broad annotation")
    parser.add_argument("celltype_fine", help="obs name for final fine annotation")
    args = parser.parse_args()
    main(args.config, args.cell_lineage, args.celltype_broad, args.celltype_fine)