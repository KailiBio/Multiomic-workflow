#!/usr/bin/env python3

"""
Author: Kaili Fan
Description:
    Update and clean h5ad, generate final h5ad and summary plots/statistics.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from cell_annotation.utils import load_config, standardize_tissue_name, move_figures_to_newdir

def update_annotations(adata, cell_lineage, broad_annotation, fine_annotation):
    if 'leiden_new' in adata.obs_names:
        adata.obs['leiden'] = adata.obs['leiden_new']
    adata.obs['cell_lineage'] = adata.obs[cell_lineage]
    adata.obs['broad_annotation'] = adata.obs[broad_annotation]
    adata.obs['fine_annotation'] = adata.obs[fine_annotation]
    return adata

def plot_cell_counts(adata, tissue_std, figdir):
    donor_col = adata.uns.get('donorID_colors', {})

    df = adata.obs[['donorID', 'leiden', 'cellannotation', 'cellannotation_color']]
    # save the cell counts
    cell_count_matrix = df.pivot_table(index='cellannotation', columns='donorID', aggfunc='size', fill_value=0)
    cell_count_matrix.to_csv(os.path.join(figdir, f'{tissue_std}_final_cell_counts.txt'), sep='\t')

    # By cell type
    sctype_counts = df['cellannotation'].value_counts()
    color_mapping = df.drop_duplicates(subset=['cellannotation']).set_index('cellannotation')['cellannotation_color']
    colors = [color_mapping[label] for label in sctype_counts.index]
    plt.figure(figsize=(7, 6))
    bars = plt.barh(sctype_counts.index, sctype_counts.values, color=colors)
    for bar in bars:
        plt.text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{int(bar.get_width())}',
                 va='center', ha='left', color='black')
    plt.title('Number of cells in each group')
    plt.xlabel('# of Cells')
    plt.grid(False)
    plt.savefig(os.path.join(figdir, f"stat.{tissue_std}.celltype_cellcount.step6.png"), dpi=300, bbox_inches='tight')
    plt.close()

    # By donor, stacked
    donor_ids = df['donorID'].unique()
    n_cols = min(4, len(donor_ids))
    n_rows = int(np.ceil(len(donor_ids) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows), sharey=True, constrained_layout=True)
    axes = axes.flatten()
    for i, donor_id in enumerate(donor_ids):
        ax = axes[i]
        donor_data = df[df['donorID'] == donor_id]
        sctype_counts = donor_data['broad_annotation'].value_counts()
        color_mapping = donor_data.drop_duplicates(subset=['broad_annotation']).set_index('broad_annotation')['broad_annotation_color']
        colors = [color_mapping[label] for label in sctype_counts.index]
        bars = ax.barh(sctype_counts.index, sctype_counts.values, color=colors)
        for bar in bars:
            ax.text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{int(bar.get_width())}',
                    va='center', ha='left', color='black')
        ax.set_title(f"{donor_id}")
        ax.set_xlabel('# of Cells')
        ax.grid(False)
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])
    fig.suptitle('Number of Cells in Each Group by Donor', fontsize=16)
    plt.savefig(os.path.join(figdir, f"stat.{tissue_std}.celltype_cellcount_byDonor.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Stacked bar, all donors
    donor_sctype_counts = df.groupby(['broad_annotation', 'donorID']).size().unstack(fill_value=0)
    donor_sctype_counts.T.plot(kind='barh', stacked=True, figsize=(10, 6),
                               color=[color_mapping[label] for label in donor_sctype_counts.index])
    plt.title('Number of Cells in Each Group by Donor')
    plt.xlabel('# of Cells')
    plt.legend(title='Cellannotation', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, f"stat.{tissue_std}.donor_cellcouns_byCelltype.png"), dpi=300, bbox_inches='tight')
    plt.close()

def save_final_h5ad(adata, workdir, tissue_std):
    adata.write(os.path.join(workdir, f"{tissue_std}_GEX.all.h5ad"), compression="gzip")
    adata_final = adata.copy()
    # Clean obs/var/uns
    obs_keep = ['tissue', 'donorID', 'cellbarcode', 'total_counts', 'n_genes', 'n_genes_by_counts',
                'log10_total_counts', 'log10_n_genes_by_counts', 'doublet_score', 'predicted_doublet',
                'doublet_probabilities', 'total_counts_mt', 'pct_counts_mt', 'total_counts_ribo',
                'pct_counts_ribo', 'total_counts_hb', 'pct_counts_hb', 'pct_exon_reads',  'MALAT1_CPM',
                'log10_MALAT1_CPM', 'leiden', 'cell_lineage', '', 'fine_annotation']
    adata_final.obs = adata_final.obs[obs_keep]
    var_keep = ['gene_ids', 'total_counts', 'n_cells', 'n_cells_by_counts', 'mean_counts',
                'feature_types', 'genome', 'mt', 'ribo', 'hb', 'pct_dropout_by_counts',
                'highly_variable', 'means', 'dispersions', 'dispersions_norm']
    adata_final.var = adata_final.var[var_keep]
    uns_keep = ['dendrogram_leiden', 'donorID_colors', 'hvg', 'leiden', 'leiden_colors',
                'log1p', 'neighbors', 'pca', 'rank_genes_groups', 'umap',]
    adata_final.uns = {k: adata_final.uns[k] for k in uns_keep if k in adata_final.uns}
    fn = os.path.join(workdir, f"{tissue_std}_GEX.final.h5ad")
    adata_final.write(fn, compression="gzip")
    print(f"[INFO] Final sharing h5ad saved to: {fn}")

def run_per_tissue(workdir, tissue, cell_lineage, broad_annotation, fine_annotation):
    tissue_std = standardize_tissue_name(tissue)

    # load h5ad
    adata_path1 = os.path.join(workdir, tissue_std, f"{tissue_std}_GEX.filtered.processes.autoAnnotated.updated.h5ad")
    adata_path2 = os.path.join(workdir, tissue_std, f"{tissue_std}_GEX.filtered.processes.autoAnnotated.h5ad")
    if os.path.exists(adata_path1):
        print(f"[INFO] loading adata from {adata_path1}")
        adata = sc.read_h5ad(adata_path1)
    elif os.path.exists(adata_path2):
        print(f"[INFO] loading adata from {adata_path2}")
        adata = sc.read_h5ad(adata_path2)
    else:
        print(f"[WARN] AnnData file not found: {adata_file}, skipping {tissue}.")
        return

    outdir = os.path.join(workdir, tissue_std, 'final')
    os.makedirs(outdir, exist_ok=True)
    
    figdir = os.path.join(outdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    os.chdir(os.path.join(outdir))

    adata = update_annotations(adata, cell_lineage, broad_annotation, fine_annotation)

    # plot final umap
    sc.pl.umap(adata, color='leiden', frameon = False, show=False, save=f'.final_umap.leiden.{tissue_std}.png')
    sc.pl.umap(adata, color='cell_lineage', frameon = False, show=False, save=f'.final_umap.cell_lineage.{tissue_std}.png')
    sc.pl.umap(adata, color='broad_annotation', frameon = False, show=False, save=f'.final_umap.broad_annotation.{tissue_std}.png')
    sc.pl.umap(adata, color='fine_annotation', frameon = False, show=False, save=f'.final_umap.fine_annotation.{tissue_std}.png')

    plot_cell_counts(adata, tissue_std, figdir)

    # Save
    save_final_h5ad(adata, workdir, tissue_std)
    # Organize figures
    move_figures_to_newdir(workdir, tissue_std, old="figures", new="final_figures")

def main(config_path, cell_lineage, broad_annotation, fine_annotation):
    config = load_config(config_path)
    
    workdir = config['paths']['workdir']
    sample_metadata = config['paths']['sample_metadata']
    
    tissue = config['params']['tissue']

    if tissue == "---":
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
            names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        
        tissues = sorted(df["tissue"].unique())
        print(f"[INFO] QC reassessment for MULTIPLE tissues: {tissues}")
        
        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n========== Final annotation for tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")
            try:
                run_per_tissue(workdir, tissue_name, cell_lineage, broad_annotation, fine_annotation)
            except Exception as e:
                print(f"[ERROR] Final h5ad generation failed for {tissue_name}: {e}")
    else:
        print(f"\n========== Final annotation for tissue: {tissue} ==========")
        try:
            run_per_tissue(workdir, tissue, cell_lineage, broad_annotation, fine_annotation)
        except Exception as e:
            print(f"[ERROR] Final h5ad generation failed for {tissue}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate final cell-annotated h5ad (manual correction) and statistics.")
    parser.add_argument("config", help="YAML config describing tissue, dirs, etc.")
    parser.add_argument("cell_lineage", help="obs name for final cell lineage")
    parser.add_argument("broad_annotation", help="obs name for final broad annotation")
    parser.add_argument("fine_annotation", help="obs name for final fine annotation")
    args = parser.parse_args()
    main(args.config, args.cell_lineage, args.broad_annotation, args.fine_annotation)