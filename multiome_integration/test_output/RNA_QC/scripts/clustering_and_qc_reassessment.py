#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Cluster and re-assess QC metrics for scRNA-seq pipeline.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from rna_qc.utils import load_config, standardize_tissue_name
from rna_qc.rna_plots import assign_donor_colors, move_figures_to_newdir, plot_umap_by_donor, plot_cellcount_per_cluster_barplot, plot_umap_highlight_by_qc_metrics, plot_qc_metrics_violin_by_cluster

def get_h5ad_path(output_h5ad_dir, tissue_std, batch_corrected, runtag):
    if batch_corrected:
        return os.path.join(output_h5ad_dir, f"{tissue_std}_GEX.filtered.RMbatch.{runtag}.h5ad")
    else:
        return os.path.join(output_h5ad_dir, f"{tissue_std}_GEX.filtered.{runtag}.h5ad")

def run_leiden_multi_res(adata, tissue_std, resolutions=[0.1, 0.5, 1.0], figdir=None):
    for res in resolutions:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:0.1f}", resolution=res, flavor="igraph")
        
    sc.pl.umap(
        adata,
        color=[f"leiden_res_{r:0.1f}" for r in resolutions],
        show=False,
        save=f".check_leiden_res.{tissue_std}.png"
    )

def set_best_leiden(adata, tissue_std, best_res=0.5):
    if 'leiden' in adata.obs: adata.obs.drop(columns='leiden', inplace=True)
    if 'leiden' in adata.uns: del adata.uns['leiden']
    if 'leiden_colors' in adata.uns: del adata.uns['leiden_colors']
    sc.tl.leiden(adata, resolution=best_res, flavor="igraph")
    
    sc.pl.umap(adata, color=["leiden"], title = f'{tissue_std}: leiden {best_res}', save=f'.LeidenCluster.{tissue_std}.png')

def save_stats(adata, out_dir, tissue_std, runtag):
    stat_fp = os.path.join(out_dir, f"{tissue_std}_stat_counts.{runtag}.txt")
    with open(stat_fp, 'w') as f:
        f.write("---------------\n")
        f.write("Total Number of Cells\n")
        cellbarcode_count_per_donor = adata.obs.groupby('donorID').size()
        f.write("Number of cellbarcodes per donorID:\n")
        f.write(cellbarcode_count_per_donor.to_string())
        f.write("\n\n")
        f.write("---------------\n")
        f.write("Number of Genes Detected per Cell\n")
        mean_n_genes_by_counts = adata.obs.groupby('donorID')['n_genes_by_counts'].mean().round().astype(int)
        median_n_genes_by_counts = adata.obs.groupby('donorID')['n_genes_by_counts'].median().astype(int)
        f.write("Mean of n_genes_by_counts per donorID:\n")
        f.write(mean_n_genes_by_counts.to_string())
        f.write("\n")
        f.write("Median of n_genes_by_counts per donorID:\n")
        f.write(median_n_genes_by_counts.to_string())
        f.write("\n\n")
        f.write("---------------\n")
        f.write("Number of UMI Counts per Cell\n")
        mean_total_counts = adata.obs.groupby('donorID')['total_counts'].mean().round(0).astype(int)
        median_total_counts = adata.obs.groupby('donorID')['total_counts'].median().astype(int)
        f.write("Mean of total_counts per donorID:\n")
        f.write(mean_total_counts.to_string())
        f.write("\n")
        f.write("Median of total_counts per donorID:\n")
        f.write(median_total_counts.to_string())
        f.write("\n")
    print(f"[INFO] Summary stats saved to {stat_fp}")

def save_processed_adata(adata, output_h5ad_dir, tissue_std, runtag):
    print("[INFO] Saving h5ad...")
    out_h5ad = os.path.join(output_h5ad_dir,  f"{tissue_std}_GEX.filtered.processed.{runtag}.h5ad")
    adata.write(out_h5ad)

def run_per_tissue(workdir, output_h5ad_dir, qc_cutoff_tissue, tissue, donor_colors, runtag,
                                resolutions=[0.1,0.5,1.0], default_res=0.5):
    
    tissue_std = standardize_tissue_name(tissue)

    # check whether use batch corrected h5ad
    qc_cutoff_dict = qc_cutoff_tissue.set_index('donorID').T.to_dict()
    batch_corrected = qc_cutoff_dict[next(iter(qc_cutoff_dict))]['Whether_batch_correction'] == "Yes"
    adata_path = get_h5ad_path(output_h5ad_dir, tissue_std, batch_corrected, runtag)
    if not os.path.exists(adata_path):
        print(f"[ERROR] No h5ad for tissue {tissue} at {adata_path}.")
        return

    print("[INFO] Loading anndata object...")
    adata = sc.read_h5ad(adata_path)
    
    figdir = os.path.join(workdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    os.chdir(workdir)
    
    all_colors = assign_donor_colors(adata.obs, donor_colors)
    
    # check and decide Leiden resolution
    run_leiden_multi_res(adata, tissue_std, resolutions)
    set_best_leiden(adata, tissue_std, best_res=default_res)

    sc.pl.umap(adata, color=["leiden", "donorID"], wspace=0.3, 
               title=[f"{tissue}: {feature}" for feature in ["leiden", "donorID"]],
               show=False, save=f".LeidenCluster-donorID.{tissue_std}.png")
    plot_umap_by_donor(adata, tissue, tissue_std, figdir, all_colors)
    plot_cellcount_per_cluster_barplot(adata, tissue_std, figdir, all_colors)
    
    # QC reassessment
    sc.pl.highest_expr_genes(adata, n_top=20, show=False, save=f".postFilter.{tissue_std}.png")
    
    qc_metrics = ["leiden", "log10_total_counts", "log10_n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", 
                  "pct_exon_reads", "log10_MALAT1_CPM", "doublet_score", "doublet_probabilities"]
    plot_umap_highlight_by_qc_metrics(adata, tissue, tissue_std, figdir, qc_metrics)
    plot_qc_metrics_violin_by_cluster(adata, tissue, tissue_std, figdir, qc_metrics[1:])
    
    # Save outputs
    save_stats(adata, figdir, tissue_std, runtag)
    save_processed_adata(adata, output_h5ad_dir, tissue_std, runtag)
    
    move_figures_to_newdir(workdir, old="figures", new=f"clustering_and_qc_reassessment.{runtag}")
    
    print(f"[INFO] Finished clustering and QC re-assessment for {tissue}.")

def main(config_path, runtag):
    config = load_config(config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    
    tissue = config['params']['tissue']

    donor_colors = config['color'].get("donor_colors")
    
    # Load qc cutoff table
    qc_cutoff_table = config['qc']['rna_qc_cutoff_table']
    if qc_cutoff_table.endswith('.xlsx') or qc_cutoff_table.endswith('.xls'):
        df_cutoff_all = pd.read_excel(qc_cutoff_table,
                                      sheet_name=config['qc']['sheet_name'], engine='openpyxl')
    else:
        df_cutoff_all = pd.read_csv(qc_cutoff_table, sep='\t')

    if tissue == "---":
        sample_metadata = config['paths']['sample_metadata']
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
                         names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        
        tissues = sorted(df["tissue"].unique())
        print(f"[INFO] QC reassessment for MULTIPLE tissues: {tissues}")
        
        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n========== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")
            
            try:
                QC_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue_name]
                run_per_tissue(workdir, output_h5ad_dir, QC_cutoff, tissue_name, donor_colors, runtag)
            except Exception as e:
                print(f"[ERROR] QC re-assessment failed for {tissue_name}: {e}")
    else:
        print(f"\n========== Processing tissue: {tissue} ==========")
        
        try:
            QC_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue]
            run_per_tissue(workdir, output_h5ad_dir, QC_cutoff, tissue, donor_colors, runtag)
        except Exception as e:
            print(f"[ERROR] QC re-assessment failed for {tissue}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster and re-assess QC metrics for scRNA-seq h5ad.")
    parser.add_argument("config", help="YAML config file describing tissue, paths, etc.")
    parser.add_argument("runtag", help="Tag for this run (e.g. round3 or v1)")
    args = parser.parse_args()
    main(args.config, args.runtag)