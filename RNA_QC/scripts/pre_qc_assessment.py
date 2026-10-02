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
from rna_qc.utils import load_config, standardize_tissue_name, setup_logging, require_keys
from rna_qc.rna_plots import assign_colors, assign_donor_colors, move_figures_to_newdir, plot_qc_violin, plot_qc_jointplot, plot_qc_cumulative_distribution, plot_doublet_hist, clustering_umap

# Non-"MT-" mitochondrial gene IDs for species whose reference annotation
# doesn't follow the human "MT-" naming convention, keyed by species
# (case-insensitive match against the sample metadata table's 'species'
# column). Add an entry here for any other non-"MT-" reference you use.
_SPECIES_MITO_GENES = {
    "macaque": ['COX1', 'KEG06_p08', 'KEG06_p13',
                'KEG06_p05', 'COX3', 'ND1', 'KEG06_p10',
                'KEG06_p02', 'ND6', 'KEG06_p07', 'ND4L',
                'ND3', 'KEG06_p12', 'KEG06_p04', 'COX2',
                'KEG06_p09', 'KEG06_p01', 'ND5', 'KEG06_p06',
                'ATP8', 'CYTB', 'ND2', 'KEG06_p11', 'KEG06_p03'],
    # marmoset follows human "MT-" convention so needs no extra entries here:
    # ['MT-NAD3', 'MT-COX1', 'MT-COX3', 'MT-COB', 'MT-NAD2', 'MT-COX2', 'MT-NAD1', 'MT-NAD4L', 'MT-NAD6', 'MT-ATP8']
}

def calculate_qc_metrics(adata, species=None):
    # calculate mt, ribo, hb
    extra_mito_genes = _SPECIES_MITO_GENES.get(str(species).strip().lower(), [])
    if extra_mito_genes:
        logging.info(f"Species '{species}' matched — including {len(extra_mito_genes)} extra "
                      f"non-'MT-' mitochondrial gene IDs for pct_counts_mt.")
    adata.var["mt"] = (
        adata.var_names.str.startswith("MT-") |
        adata.var_names.isin(extra_mito_genes)
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
            logging.info(f"Doublet detection for {ID}...")
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
            logging.warning(f"Could not run doublet GMM for {ID}: {e}")

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
    logging.info("Saving h5ad...")
    adata.write(os.path.join(output_h5ad_dir, f'{tissue_std}_GEX.withQC.h5ad'))

def run_per_tissue(working_df, tissue, output_h5ad_dir, outdir, my_color_palette, scrinvex_dir, nmads, key = "sampleID"):
    tissue_std = standardize_tissue_name(tissue)

    adata_path = os.path.join(output_h5ad_dir, f"{tissue_std}_GEX.raw.h5ad")
    if not os.path.exists(adata_path):
        logging.error(f"No h5ad found for tissue {tissue_std} at {adata_path}.")
        return

    logging.info(f"Loading raw anndata obejct for tissue: {tissue_std}...")
    adata = ad.read_h5ad(adata_path)

    figdir = os.path.join(outdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    os.chdir(os.path.join(outdir))

    logging.info(f'all figure plots by {key}')
    if key == 'donorID':
        all_colors = assign_donor_colors(working_df, my_color_palette, key = key)
        logging.info(f"use colors: {all_colors}")
    elif key == 'sampleID':
        all_colors = assign_colors(sorted(set(working_df['rnaID'])), palette=my_color_palette)
        logging.info(f"use colors: {all_colors}")
    else:
        logging.warning("need to edit for colors")

    species_values = working_df['species'].dropna().unique() if 'species' in working_df else []
    if len(species_values) == 1:
        species = species_values[0]
    else:
        if len(species_values) > 1:
            logging.warning(f"Multiple species values found for tissue {tissue}: {list(species_values)}. "
                             f"Skipping species-specific mitochondrial gene list.")
        species = None
    calculate_qc_metrics(adata, species=species)

    # Highest expressed genes
    sc.pl.highest_expr_genes(adata, n_top=20, save=f'.{tissue_std}.png', show=False)

    # Violin plots for cell QC
    QC_metrics = ['log10_n_genes_by_counts', 'log10_total_counts', 'pct_counts_in_top_50_genes',
                  'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_exon_reads',
                  'log10_MALAT1_CPM']
    if 'pct_exon_reads' not in adata.obs:
        QC_metrics.remove('pct_exon_reads')

    # Cutoff for percent exon reads (check data, not config path)
    if 'pct_exon_reads' in adata.obs and adata.obs['pct_exon_reads'].notna().any():
        Q75 = np.nanpercentile(adata.obs['pct_exon_reads'].dropna(), 75)
        Q25 = np.nanpercentile(adata.obs['pct_exon_reads'].dropna(), 25)
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
        plot_qc_violin(adata, metric, tissue, tissue_std, all_colors, metrics_with_cutoffs, figdir, key, nmads=nmads, add_mad_lines=True)

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
    plot_doublet_hist(adata, all_colors, tissue, tissue_std, figdir, key)

    # Clustering, UMAP, etc.
    clustering_umap(adata, tissue, tissue_std, figdir, key)

    # Save output h5ad
    compress_and_save(adata, output_h5ad_dir, tissue_std)
    # Move figures to new directory
    move_figures_to_newdir(outdir, old="figures", new="pre_qc_assessment")
    logging.info(f"RNA QC pre-assessment complete for {tissue}.")

def main(config_path, nmads=5.0):
    setup_logging()
    config = load_config(config_path)
    require_keys(config, [
        "paths.workdir", "paths.output_h5ad_dir", "paths.sample_metadata",
        "params.tissue", "color", "my_color_palette",
    ], context=config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    scrinvex_dir = config['paths'].get('scrinvex_dir', None)

    tissue = config['params']['tissue']

    donor_colors = config['color'].get("donor_colors")
    my_color_palette = config["my_color_palette"]

    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None, index_col=False,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])

    failed_tissues = []
    if tissue == "---":
        tissues = sorted(df["tissue"].unique())
        logging.info(f"Running analysis for MULTIPLE tissues: {tissues}")

        for idx, tissue_name in enumerate(tissues, 1):
            logging.info(f"============== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==============")

            try:
                working_df = df[df["tissue"] == tissue_name]
                run_per_tissue(working_df, tissue_name, output_h5ad_dir, workdir, my_color_palette,
                               scrinvex_dir, nmads, key = "sampleID")
            except Exception as e:
                logging.error(f"Encountered error for tissue {tissue_name}: {e}")
                failed_tissues.append(tissue_name)
    else:
        logging.info(f"========== Processing tissue: {tissue} ==========")

        try:
            working_df = df[df["tissue"] == tissue]
            run_per_tissue(working_df, tissue, output_h5ad_dir, workdir, my_color_palette,
                           scrinvex_dir, nmads, key = "sampleID")
        except Exception as e:
            logging.error(f"Encountered error for tissue {tissue}: {e}")
            failed_tissues.append(tissue)

    if failed_tissues:
        logging.error(f"pre_qc_assessment failed for {len(failed_tissues)} tissue(s): {failed_tissues}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pre-QC metric/assessment for scRNA-seq h5ad.")
    parser.add_argument("config", help="YAML config file describing tissue, paths, etc.")
    parser.add_argument("--nmads", type=float, default=5.0,
                        help="number of MADs from the median used to define cutoffs (default: 5).")

    args = parser.parse_args()
    try:
        main(args.config, nmads=args.nmads)
    except (FileNotFoundError, KeyError, ValueError) as e:
        logging.error(f"{type(e).__name__}: {e}")
        sys.exit(1)
