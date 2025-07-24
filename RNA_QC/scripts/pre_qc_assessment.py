#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Computes QC metrics and runs doublet detection, generates plots for further QC filtering.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse
import logging
from sklearn.mixture import BayesianGaussianMixture

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from rna_qc.utils import load_config, standardize_tissue_name
from rna_qc.rna_plots import assign_donor_colors, move_figures_to_newdir, plot_qc_violin, plot_qc_jointplot, plot_qc_cumulative_distribution, plot_doublet_hist, clustering_umap

def calculate_qc_metrics(adata):
    # calculate mt, ribo, hb
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=False)

    # convert to log10
    adata.obs['log10_n_genes_by_counts'] = np.log10(adata.obs['n_genes_by_counts'] + 0.1)
    adata.obs['log10_total_counts'] = np.log10(adata.obs['total_counts'] + 0.1)
    
    # normalize raw counts to CPM and log10
    if not scipy.sparse.issparse(adata.X):
        raw_counts = adata.X
    else:
        raw_counts = adata.X.toarray()
    row_sums = np.sum(raw_counts, axis=1)[:, np.newaxis]
    normalized_counts = raw_counts / row_sums * 1e6
    adata.layers['CPM'] = normalized_counts
    
    # MALAT1
    CPM_df = pd.DataFrame(normalized_counts, index=adata.obs_names, columns=adata.var_names)
    if "MALAT1" in CPM_df.columns:
        adata.obs['MALAT1_CPM'] = CPM_df['MALAT1']
        adata.obs['log10_MALAT1_CPM'] = np.log10(CPM_df['MALAT1']+0.1)
    else:
        adata.obs['MALAT1_CPM'] = np.nan
        adata.obs['log10_MALAT1_CPM'] = np.nan

def get_doublet_probability(
    doublet_scores_sim: np.ndarray,
    doublet_scores: np.ndarray,
    random_state: int = 0,
    verbose: bool = False,
):
    from sklearn.mixture import BayesianGaussianMixture

    X = doublet_scores_sim.reshape((-1, 1))
    gmm = BayesianGaussianMixture(
        n_components=2, n_init=10, max_iter=1000, random_state=random_state
    ).fit(X)

    if verbose:
        logging.info("GMM means: {}".format(gmm.means_))

    i = np.argmax(gmm.means_)
    return gmm.predict_proba(doublet_scores.reshape((-1, 1)))[:,i]

def run_doublet_detection(adata, donor_col):
    
    sc.pp.scrublet(adata, batch_key="donorID")
    
    dfs = []
    for donorID in adata.obs['donorID'].unique():
        try:
            print(f"[INFO] Doublet detection for {donorID}...")
            doublet_scores_sim = adata.uns['scrublet']['batches'][donorID]['doublet_scores_sim']
            doublet_scores = adata[adata.obs['donorID'] == donorID].obs['doublet_score'].to_numpy()
            
            probabilities = get_doublet_probability(
                doublet_scores_sim=doublet_scores_sim,
                doublet_scores=doublet_scores,
                random_state=0,
                verbose=True
            )
            
            obs_names = adata[adata.obs['donorID'] == donorID].obs_names
            df = pd.DataFrame({'obs_names': obs_names, 'probabilities': probabilities})
            dfs.append(df)
            
        except Exception as e:
            print(f"[WARNING] Could not run doublet GMM for {donorID}: {e}")
    
    if dfs:
        all_prob_df = pd.concat(dfs).set_index('obs_names')
        adata.obs['doublet_probabilities'] = all_prob_df['probabilities']
    else:
        adata.obs['doublet_probabilities'] = np.nan

def compress_and_save(adata, output_h5ad_dir, tissue_std):
    # Compress to sparse
    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    if 'CPM' in adata.layers and not scipy.sparse.issparse(adata.layers['CPM']):
        adata.layers['CPM'] = scipy.sparse.csr_matrix(adata.layers['CPM'])
    print("[INFO] Saving h5ad...")
    adata.write(os.path.join(output_h5ad_dir, f'{tissue_std}_GEX.withQC.h5ad'))

def run_per_tissue(working_df, tissue, output_h5ad_dir, outdir, donor_colors, scrinvex_dir):
    tissue_std = standardize_tissue_name(tissue)
        
    adata_path = os.path.join(output_h5ad_dir, f"{tissue_std}_GEX.raw.h5ad")
    if not os.path.exists(adata_path):
        print(f"[ERROR] No h5ad found for tissue {tissue_std} at {adata_path}.")
        return
    
    print(f"[INFO] Loading raw anndata obejct for tissue: {tissue_std}...")
    adata = ad.read_h5ad(adata_path)

    figdir = os.path.join(outdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    os.chdir(os.path.join(outdir))

    all_colors = assign_donor_colors(working_df, donor_colors)

    calculate_qc_metrics(adata)

    # Highest expressed genes
    sc.pl.highest_expr_genes(adata, n_top=20, save=f'.{tissue_std}.png', show=False)
    
    # Violin plots for cell QC
    QC_metrics = ['log10_n_genes_by_counts', 'log10_total_counts', 'pct_counts_in_top_50_genes',
                  'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_exon_reads', 'log10_MALAT1_CPM']

    # Cutoff for percent exon reads
    if scrinvex_dir != None:
        Q75 = np.nanpercentile(adata.obs['pct_exon_reads'], 75)
        Q25 = np.nanpercentile(adata.obs['pct_exon_reads'], 25)
        percent_exon_cutoff = Q75 + 1.5 * (Q75 - Q25)
    else:
        percent_exon_cutoff = None

    # default cutoffs
    metrics_with_cutoffs = {
        'log10_n_genes_by_counts': [np.log10(200), np.log10(5000)],
        'log10_total_counts': [np.log10(25000)],
        'pct_counts_mt': [10],
        'pct_counts_ribo': [10],
        'pct_exon_reads': [percent_exon_cutoff],
        'log10_MALAT1_CPM': [np.log10(10)]
    }

    for metric in QC_metrics:
        plot_qc_violin(adata, metric, tissue, tissue_std, all_colors, metrics_with_cutoffs, figdir)

    # Joint scatter gene/cell counts
    plot_qc_jointplot(
        adata=adata, 
        x_metric='log10_total_counts', 
        y_metric='log10_n_genes_by_counts', 
        title=f"{tissue_std}",
        save_path=os.path.join(figdir, f"RNA_QC.{tissue_std}.cellCounts_geneCounts_scatter.png"),
    )
    
    # Joint per-donor
    for donorID in adata.obs['donorID'].unique():
        donor_data = adata[adata.obs['donorID'] == donorID]
        plot_qc_jointplot(
            adata=donor_data,
            x_metric='log10_total_counts',
            y_metric='log10_n_genes_by_counts',
            title=f"{tissue_std}: {donorID}",
            save_path=os.path.join(figdir, f"RNA_QC.{tissue_std}-{donorID}.cellCounts_geneCounts_scatter.png"),
        )

    # Cumulative plots
    plot_qc_cumulative_distribution(
        adata=adata,
        metrics=QC_metrics,
        tissue=tissue,
        tissue_std=tissue_std,
        all_colors=all_colors,
        figdir=figdir
    )

    # Doublet detection
    run_doublet_detection(adata, all_colors)
    plot_doublet_hist(adata, all_colors, tissue, tissue_std, figdir)

    # Clustering, UMAP, etc.
    clustering_umap(adata, tissue, tissue_std, figdir)

    # Save output h5ad
    compress_and_save(adata, output_h5ad_dir, tissue_std)
    # Move figures to new directory
    move_figures_to_newdir(outdir, old="figures", new="pre_qc_assessment")
    print(f"[INFO] RNA QC pre-assessment complete for {tissue}.")

def main(config_path):
    config = load_config(config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    scrinvex_dir = config['paths'].get('scrinvex_dir', None)
    
    tissue = config['params']['tissue']

    donor_colors = config['color'].get("donor_colors")

    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    
    if tissue == "---":
        tissues = sorted(df["tissue"].unique())
        print(f"[INFO] Running analysis for MULTIPLE tissues: {tissues}")
        
        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n============== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==============")

            try:
                working_df = df[df["tissue"] == tissue_name]
                run_per_tissue(working_df, tissue_name, output_h5ad_dir, workdir, donor_colors, scrinvex_dir)
            except Exception as e:
                print(f"[ERROR] Encountered error for tissue {tissue_name}: {e}")
    else:
        print(f"\n========== Processing tissue: {tissue} ==========")
        
        try:
            working_df = df[df["tissue"] == tissue]
            run_per_tissue(working_df, tissue, output_h5ad_dir, workdir, donor_colors, scrinvex_dir)
        except Exception as e:
            print(f"[ERROR] Encountered error for tissue {tissue}: {e}")
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pre-QC metric/assessment for scRNA-seq h5ad.")
    parser.add_argument("config", help="YAML config file describing tissue, paths, etc.")
    args = parser.parse_args()
    main(args.config)