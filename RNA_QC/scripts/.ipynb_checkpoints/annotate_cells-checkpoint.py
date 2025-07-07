#!/usr/bin/env python3

"""
Author: Kaili Fan
Description:
    Identify DEG and generates dot plots, run automatic annotation with scType.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import urllib.request
import shutil

from rna_qc.utils import load_config, standardize_tissue_name
from rna_qc.rna_plots import move_figures_to_newdir

def filter_marker_genes(marker_gene_list, var_names):
    valid_marker_genes = {}
    for cell_type, marker_genes in marker_gene_list.items():
        valid_genes = [gene for gene in marker_genes if gene in var_names]
        if valid_genes:
            valid_marker_genes[cell_type] = valid_genes
    return valid_marker_genes
    
def load_marker_genes(config):
    if config['annotation']['marker_gene_file'].endswith('.xlsx') or config['annotation']['marker_gene_file'].endswith('.xls'):
        marker_gene_file = pd.read_excel(config['annotation']['marker_gene_file'],
                                      sheet_name=config['annotation']['sheet_name'], engine='openpyxl')
    else:
        marker_gene_file = pd.read_csv(config['annotation']['marker_gene_file'], sep='\t')
        
    marker_gene_file["CellType"] = marker_gene_file["CellType"].str.rstrip()
    marker_gene_list = marker_gene_file.groupby("CellType")["Gene"].apply(list).to_dict()
    return marker_gene_list

def Process_DEG(adata, tissue, tissue_std, valid_marker_genes, outdir):
    # rank and pick DEGs among clusters
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="t-test")
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False, save=f'.top25DEG.{tissue_std}.png')
    sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden", standard_scale="var", n_genes=5,
                                    show=False, save=f'.top5DEG.{tissue_std}.png')

    # save top DEGs
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    output_txt = os.path.join(outdir, "top_genes_per_cluster.txt")
    with open(output_txt, "w") as f:
        for group in groups:
            f.write(f"Top DEGs for cluster {group}:\n")
            genes = result['names'][group][:20]
            f.write(", ".join(genes) + "\n\n")
    print(f"[Info] Saved top DEGs to {output_txt}")

    # UMAPs for top14 DEGs of each cluster
    clusters = adata.obs['leiden'].unique()
    for cluster in clusters:
        dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group=cluster).head(14)["names"]
        titles = [f"{tissue}: leiden-{cluster}\n {feature}" for feature in ["leiden", *dc_cluster_genes]]
        sc.pl.umap(
            adata,
            color=["leiden", *dc_cluster_genes],
            legend_loc="on data",
            ncols=5,
            hspace=0.5,
            title=titles,
            show=False,
            save=f'.topDEG_leiden{cluster}.{tissue_std}.png'
        )

def load_sctype():
    url = "https://raw.githubusercontent.com/kris-nader/sc-type-py/main/sctype_py.py"
    response = urllib.request.urlopen(url)
    script = response.read().decode()
    exec(script, globals())

def run_sctype_annotation(adata, marker_gene_list):
    # Build gs2 (empty/non-immune as in ScType)
    marker_gene_list2 = {cell_type: [''] for cell_type in marker_gene_list}
    
    # Extract scaled data for ScType: cells as columns, genes as rows
    if hasattr(adata.X, "todense"):
        scaled_data = pd.DataFrame(adata.X.todense(), index=adata.obs_names, columns=adata.var_names).T
    else:
        scaled_data = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names).T
        
    es_max = sctype_score(scRNAseqData=scaled_data, scaled=True, gs=marker_gene_list, gs2=marker_gene_list2)
    clusters = adata.obs['leiden'].unique()
    cL_results = pd.concat([process_cluster(cluster, adata, es_max, "leiden") for cluster in clusters])
    sctype_scores = cL_results.groupby('cluster').apply(lambda x: x.nlargest(1, 'scores')).reset_index(drop=True)
    sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'
    # Assign annotation
    adata.obs['sctype_annotation'] = ""
    for cluster in sctype_scores['cluster'].unique():
        cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
        cl_type_value = cl_type['type'].iloc[0]
        adata.obs.loc[adata.obs['leiden'] == cluster, 'sctype_annotation'] = cl_type_value.rstrip()

def plot_final_marker_dotplots(adata, tissue_std, valid_marker_genes):
    # Dotplot for top HVGs (by sctype_annotation)
    sc.pl.rank_genes_groups_dotplot(adata, groupby="sctype_annotation", standard_scale="var", n_genes=5,
                                    show=False, save=f'.DEG2.{tissue_std}.png')
    # Dotplot for provided markers by sctype_annotation
    sc.pl.dotplot(adata, valid_marker_genes, groupby="sctype_annotation", standard_scale="var",
                  show=False, save=f'.markerGenes2.{tissue_std}.png')

def save_results(adata, output_h5ad_dir, out_dir, tissue_std):
    # Save AnnData
    out_h5ad = os.path.join(output_h5ad_dir, tissue_std, f'{tissue_std}_GEX.filtered.processes.autoAnnotated.h5ad')
    adata.write(out_h5ad, compression="gzip")
    
    # Save Main annotation table
    auto_cellannotation_df = adata.obs[['leiden', 'rough_cellcluster_annotation', 'sctype_annotation']].drop_duplicates()
    auto_cellannotation_df.to_csv(os.path.join(out_dir, f'{tissue_std}_GEX.autoAnnotation.tsv'), sep='\t', index=False)

def run_per_tissue(workdir, output_h5ad_dir, output_figures_dir, tissue, marker_gene_list):
    tissue_std = standardize_tissue_name(tissue)

    # Load data
    adata_path = os.path.join(output_h5ad_dir, tissue_std, f"{tissue_std}_GEX.filtered.processed.h5ad")
    if not os.path.exists(adata_path):
        print(f"[WARN] No h5ad for tissue {tissue} at {adata_path}. Skipping.")
        return
    adata = sc.read_h5ad(adata_path)
    
    outdir = os.path.join(workdir, tissue_std)
    os.makedirs(outdir, exist_ok=True)
    
    figdir = os.path.join(output_figures_dir, tissue_std, 'figures')
    os.makedirs(figdir, exist_ok=True)
    os.chdir(os.path.join(output_figures_dir, tissue_std))

    # check marker genes
    valid_marker_genes = filter_marker_genes(marker_gene_list, adata.var_names)
    print(f"[INFO] {tissue}: {adata.shape[0]} cells, {adata.shape[1]} genes; {len(valid_marker_genes)} valid marker sets")

    # 1. DEG and marker plots
    Process_DEG(adata, tissue, tissue_std, valid_marker_genes, outdir)

    # 2. Dotplot with input markers
    sc.pl.dotplot(adata, valid_marker_genes, groupby="leiden", standard_scale="var",
                  show=False, save=f'.markerGenes.{tissue_std}.png')

    # 3. ScType annotation
    load_sctype()
    run_sctype_annotation(adata, marker_gene_list)
    
    sc.pl.umap(adata, color='sctype_annotation', frameon=False, show=False, save=f'.ScType_Annotation.{tissue_std}.png')
    plot_final_marker_dotplots(adata, tissue_std, valid_marker_genes)

    # 4. Save h5ad and annotation table
    save_results(adata, output_h5ad_dir, outdir, tissue_std)

    # 5. Move figures
    move_figures_to_newdir(output_figures_dir, tissue_std, old="figures", new="cell_annotation")

    print(f"[INFO] Completed annotation step for {tissue}.")

def main(config_path):
    config = load_config(config_path)
    
    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    sample_metadata = config['paths']['sample_metadata']
    output_figures_dir = config['paths']['output_figures_dir']

    marker_gene_file = config['annotation']['marker_gene_file']

    # Load marker gene file
    marker_gene_list = load_marker_genes(config)
        
    tissue = config['params']['tissue']
    
    if tissue == "---":
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
            names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        tissues = sorted(df['tissue'].unique())
        print(f"[INFO] Annotating MULTIPLE tissues: {tissues}")
        
        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n----- Annotating tissue: {tissue_name} ({idx}/{len(tissues)}) -----")
            try:
                run_per_tissue(workdir, output_h5ad_dir, output_figures_dir, tissue_name, marker_gene_list)
            except Exception as e:
                print(f"[ERROR] Cell annotation failed for {tissue_name}: {e}")
    else:
        print(f"[INFO] Annotating tissue: {tissue}")
        try:
            run_per_tissue(workdir, output_h5ad_dir, output_figures_dir, tissue, marker_gene_list)
        except Exception as e:
            print(f"[ERROR] Cell annotation failed for {tissue}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cell annotation workflow for scRNA-seq h5ad.")
    parser.add_argument("config", help="YAML config with tissue/path/marker genes spec.")
    args = parser.parse_args()
    main(args.config)