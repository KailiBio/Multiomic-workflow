#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Computes QC metrics and runs doublet detection, generates plots for further QC filtering.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

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
from rna_qc.rna_plots import assign_colors, assign_donor_colors, move_figures_to_newdir, plot_qc_violin, plot_qc_jointplot, plot_qc_cumulative_distribution, plot_doublet_hist, clustering_umap

def calculate_qc_metrics(adata):
    # calculate mt, ribo, hb
    macaque_mito_gene_list = ['COX1', 'KEG06_p08', 'KEG06_p13', 
                              'KEG06_p05', 'COX3', 'ND1', 'KEG06_p10', 
                              'KEG06_p02', 'ND6', 'KEG06_p07', 'ND4L', 
                              'ND3', 'KEG06_p12', 'KEG06_p04', 'COX2', 
                              'KEG06_p09', 'KEG06_p01', 'ND5', 'KEG06_p06', 
                              'ATP8', 'CYTB', 'ND2', 'KEG06_p11', 'KEG06_p03']
    # listing here for completeness but follows human convention
    #marmoset_mito_gene_list = ['MT-NAD3', 'MT-COX1', 'MT-COX3', 'MT-COB', 'MT-NAD2', 'MT-COX2', 'MT-NAD1', 'MT-NAD4L', 'MT-NAD6', 'MT-ATP8']
    adata.var["mt"] = (
        adata.var_names.str.startswith("MT-") | 
        adata.var_names.isin(macaque_mito_gene_list)
    )
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

    X = doublet_scores_sim.reshape((-1, 1))
    with warnings.catch_warnings():
        # Ignore RuntimeWarnings (divide by zero, overflow, etc.)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        gmm = BayesianGaussianMixture(
            n_components=2, n_init=10, max_iter=1000, random_state=random_state
        ).fit(X)

    if verbose:
        logging.info("GMM means: {}".format(gmm.means_))

    i = np.argmax(gmm.means_)
    return gmm.predict_proba(doublet_scores.reshape((-1, 1)))[:,i]

def run_doublet_detection(adata, donor_col, key):
    
    sc.pp.scrublet(adata, batch_key=key)
    
    dfs = []
    for ID in adata.obs[key].unique():
        try:
            print(f"[INFO] Doublet detection for {ID}...")
            doublet_scores_sim = adata.uns['scrublet']['batches'][ID]['doublet_scores_sim']
            doublet_scores = adata[adata.obs[key] == ID].obs['doublet_score'].to_numpy()
            
            probabilities = get_doublet_probability(
                doublet_scores_sim=doublet_scores_sim,
                doublet_scores=doublet_scores,
                random_state=0,
                verbose=True
            )
            
            obs_names = adata[adata.obs[key] == ID].obs_names
            df = pd.DataFrame({'obs_names': obs_names, 'probabilities': probabilities})
            dfs.append(df)
            
        except Exception as e:
            print(f"[WARNING] Could not run doublet GMM for {ID}: {e}")
    
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

def run_per_tissue(working_df, tissue, output_h5ad_dir, outdir, my_color_palette, scrinvex_dir, nmads, key = "sampleID", qc_cutoff_df=None, mad_scope='per-sample'):
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

    print(f'[INFO] all figure plots by {key}')
    if key == 'donorID':
        all_colors = assign_donor_colors(working_df, my_color_palette, key = key)
        print(f"use colors: {all_colors}")
    elif key == 'sampleID':
        sample_colors={}
        #all_colors = assign_donor_colors(working_df, sample_colors, key = 'rnaID')
        all_colors = assign_colors(sorted(set(working_df['rnaID'])), palette=my_color_palette)
        print(f"use colors: {all_colors}")
    else:
        print("[WARNING] need to edit for colors")

    calculate_qc_metrics(adata)

    # Highest expressed genes
    sc.pl.highest_expr_genes(adata, n_top=20, save=f'.{tissue_std}.png', show=False)
    
    # Violin plots for cell QC
    QC_metrics = ['log10_n_genes_by_counts', 'log10_total_counts', 'pct_counts_in_top_50_genes',
                  'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_exon_reads',
                  'log10_MALAT1_CPM']
    if 'pct_exon_reads' not in adata.obs:
        QC_metrics.remove('pct_exon_reads')
    if ('log10_MALAT1_CPM' not in adata.obs or
        adata.obs['log10_MALAT1_CPM'].isnull().all() or
        ((adata.obs['log10_MALAT1_CPM'] == np.inf) | (adata.obs['log10_MALAT1_CPM'] == -np.inf)).all()):
        QC_metrics.remove('log10_MALAT1_CPM')
        print("[INFO] Skipping MALAT1 plots — gene absent from reference or all values NaN/inf.")

    # Build cutoff bounds from Excel QC table if available, otherwise use defaults.
    # Each metric maps to {'lower': v_or_None, 'upper': v_or_None}; MAD lines are
    # drawn only on sides where a cutoff exists.
    if qc_cutoff_df is not None and not qc_cutoff_df.empty:
        row = qc_cutoff_df.iloc[0]  # cutoffs are per-sample but use first row for tissue-wide lines
        def _safe(val, transform=None):
            if val == '---' or pd.isna(val):
                return None
            v = float(val)
            return transform(v) if transform else v

        metrics_with_cutoffs = {
            'log10_n_genes_by_counts': {
                'lower': _safe(row.get('Min_genes_in_cells'), lambda v: np.log10(max(v, 1))),
                'upper': _safe(row.get('Max_genes_in_cells'), np.log10),
            },
            'log10_total_counts': {'upper': _safe(row.get('Max_counts_in_cells'), np.log10)},
            'pct_counts_mt': {'upper': _safe(row.get('Max_percent_mt_in_cells'))},
            'pct_counts_ribo': {'upper': _safe(row.get('Max_percent_ribo_in_cells'))},
            'pct_exon_reads': {'upper': _safe(row.get('Exon_ratio_cutoffs'))},
            'log10_MALAT1_CPM': {'lower': _safe(row.get('MALAT1_CPM_cutoffs'), lambda v: np.log10(max(v, 1)))},
        }
        print(f"[INFO] Pre-QC cutoff lines from Excel: {metrics_with_cutoffs}")
    else:
        # Fallback: default cutoffs
        if 'pct_exon_reads' in adata.obs and adata.obs['pct_exon_reads'].notna().any():
            Q75 = np.nanpercentile(adata.obs['pct_exon_reads'].dropna(), 75)
            Q25 = np.nanpercentile(adata.obs['pct_exon_reads'].dropna(), 25)
            percent_exon_cutoff = Q75 + 1.5 * (Q75 - Q25)
        else:
            percent_exon_cutoff = None

        metrics_with_cutoffs = {
            'log10_n_genes_by_counts': {'lower': np.log10(200), 'upper': np.log10(5000)},
            'log10_total_counts': {'upper': np.log10(25000)},
            'pct_counts_mt': {'upper': 10},
            'pct_counts_ribo': {'upper': 10},
            'pct_exon_reads': {'upper': percent_exon_cutoff},
            'log10_MALAT1_CPM': {'lower': np.log10(10)},
        }
        print("[INFO] Pre-QC cutoff lines using defaults (no Excel cutoffs provided)")
    for metric in QC_metrics:
        plot_qc_violin(adata, metric, tissue, tissue_std, all_colors, metrics_with_cutoffs, figdir, key, nmads=nmads, add_mad_lines=True, mad_scope=mad_scope)

    # Joint scatter gene/cell counts
    plot_qc_jointplot(
        adata=adata, 
        x_metric='log10_total_counts', 
        y_metric='log10_n_genes_by_counts', 
        title=f"{tissue_std}",
        save_path=os.path.join(figdir, f"RNA_QC.{tissue_std}.cellCounts_geneCounts_scatter.png"),
    )
    
    # Joint per-donor
    for ID in adata.obs[key].unique():
        per_data = adata[adata.obs[key] == ID]
        plot_qc_jointplot(
            adata=per_data,
            x_metric='log10_total_counts',
            y_metric='log10_n_genes_by_counts',
            title=f"{tissue_std}: {ID}",
            save_path=os.path.join(figdir, f"RNA_QC.{tissue_std}-{ID}.cellCounts_geneCounts_scatter.png"),
        )

    # Cumulative plots
    plot_qc_cumulative_distribution(adata, QC_metrics, tissue, tissue_std, all_colors, figdir, key)

    # Doublet detection
    run_doublet_detection(adata, all_colors, key)
    plot_doublet_hist(adata, tissue, tissue_std, figdir, key)

    # Clustering, UMAP, etc.
    clustering_umap(adata, tissue, tissue_std, figdir, key)

    # Save output h5ad
    compress_and_save(adata, output_h5ad_dir, tissue_std)
    # Move figures to new directory
    move_figures_to_newdir(outdir, old="figures", new="pre_qc_assessment")
    print(f"[INFO] RNA QC pre-assessment complete for {tissue}.")

def main(config_path, nmads=5.0, mad_scope='per-sample'):
    config = load_config(config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    scrinvex_dir = config['paths'].get('scrinvex_dir', None)
    
    tissue = config['params']['tissue']

    donor_colors = config['color'].get("donor_colors")
    my_color_palette = config["my_color_palette"]

    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None, index_col=False,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])

    # Load QC cutoff table (optional — used for red lines in pre-QC plots)
    df_cutoff_all = None
    qc_config = config.get('qc', {})
    qc_cutoff_table = qc_config.get('rna_qc_cutoff_table')
    if qc_cutoff_table and os.path.exists(qc_cutoff_table):
        print(f"[INFO] Loading QC cutoff table for pre-QC plots: {qc_cutoff_table}")
        if qc_cutoff_table.endswith('.xlsx') or qc_cutoff_table.endswith('.xls'):
            df_cutoff_all = pd.read_excel(qc_cutoff_table,
                                          sheet_name=qc_config.get('sheet_name', 0), engine='openpyxl')
        else:
            df_cutoff_all = pd.read_csv(qc_cutoff_table, sep='\t')

    if tissue == "---":
        tissues = sorted(df["tissue"].unique())
        print(f"[INFO] Running analysis for MULTIPLE tissues: {tissues}")

        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n============== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==============")

            try:
                working_df = df[df["tissue"] == tissue_name]
                qc_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue_name] if df_cutoff_all is not None else None
                run_per_tissue(working_df, tissue_name, output_h5ad_dir, workdir, my_color_palette,
                               scrinvex_dir, nmads, key="sampleID", qc_cutoff_df=qc_cutoff, mad_scope=mad_scope)
            except Exception as e:
                print(f"[ERROR] Encountered error for tissue {tissue_name}: {e}")
    else:
        print(f"\n========== Processing tissue: {tissue} ==========")

        try:
            working_df = df[df["tissue"] == tissue]
            qc_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue] if df_cutoff_all is not None else None
            run_per_tissue(working_df, tissue, output_h5ad_dir, workdir, my_color_palette,
                           scrinvex_dir, nmads, key="sampleID", qc_cutoff_df=qc_cutoff, mad_scope=mad_scope)
        except Exception as e:
            print(f"[ERROR] Encountered error for tissue {tissue}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pre-QC metric/assessment for scRNA-seq h5ad.")
    parser.add_argument("config", help="YAML config file describing tissue, paths, etc.")
    parser.add_argument("--nmads", type=float, default=5.0,
                        help="number of MADs from the median used to define cutoffs (default: 5).")
    parser.add_argument("--mad-scope", choices=["per-sample", "per-tissue"], default="per-sample",
                        help="Show per-sample or per-tissue MAD lines on violin plots (default: per-sample).")

    args = parser.parse_args()
    main(args.config, nmads=args.nmads, mad_scope=args.mad_scope)
