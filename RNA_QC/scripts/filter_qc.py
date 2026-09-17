#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Filters cells by user-provided (per sample) QC cutoffs, generates Upset plots, processes filtered data, and saves results.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse
from scipy.stats import median_abs_deviation

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from rna_qc.utils import load_config, standardize_tissue_name, setup_logging, require_keys
from rna_qc.rna_plots import assign_donor_colors, move_figures_to_newdir, plot_upset, run_umap_clustering, plot_doublet_hist

def _compute_mad_for_sample(sample_obs, nmads):
    """Compute MAD-based QC thresholds for a single sample's cells.

    Count metrics (n_genes_by_counts, total_counts) use log1p-transformed
    values for MAD computation (standard practice for right-skewed count
    data; see OSCA, Luecken & Theis 2019). Percentage metrics use raw values.

    Returns a dict of threshold values.
    """
    thresholds = {}

    # --- n_genes_by_counts (Min and Max, log-scale MAD) ---
    log_vals = np.log1p(sample_obs['n_genes_by_counts'].dropna())
    med = np.nanmedian(log_vals)
    mad = median_abs_deviation(log_vals, nan_policy='omit')
    thresholds['Min_genes_in_cells'] = int(max(0, np.floor(np.expm1(med - nmads * mad))))
    thresholds['Max_genes_in_cells'] = int(np.ceil(np.expm1(med + nmads * mad)))

    # --- total_counts (Max, log-scale MAD) ---
    log_vals = np.log1p(sample_obs['total_counts'].dropna())
    med = np.nanmedian(log_vals)
    mad = median_abs_deviation(log_vals, nan_policy='omit')
    thresholds['Max_counts_in_cells'] = int(np.ceil(np.expm1(med + nmads * mad)))

    # --- pct_counts_mt (Max, floor of 5%) ---
    med = np.nanmedian(sample_obs['pct_counts_mt'])
    mad = median_abs_deviation(sample_obs['pct_counts_mt'], nan_policy='omit')
    thresholds['Max_percent_mt_in_cells'] = round(min(100, max(5.0, med + nmads * mad)), 2)

    # --- pct_counts_ribo (Max, floor of 5%) ---
    med = np.nanmedian(sample_obs['pct_counts_ribo'])
    mad = median_abs_deviation(sample_obs['pct_counts_ribo'], nan_policy='omit')
    thresholds['Max_percent_ribo_in_cells'] = round(min(100, max(5.0, med + nmads * mad)), 2)

    # --- Exon_ratio_cutoffs (upper bound only, no lower MAD) ---
    if 'pct_exon_reads' in sample_obs.columns and sample_obs['pct_exon_reads'].notna().any():
        vals = sample_obs['pct_exon_reads'].dropna()
        if len(vals) > 0:
            Q75 = np.percentile(vals, 75)
            Q25 = np.percentile(vals, 25)
            thresholds['Exon_ratio_cutoffs'] = round(min(100, Q75 + 1.5 * (Q75 - Q25)), 2)
        else:
            thresholds['Exon_ratio_cutoffs'] = '---'
    else:
        thresholds['Exon_ratio_cutoffs'] = '---'

    # --- MALAT1_CPM (Min lower bound, Max upper bound, log-scale MAD) ---
    if 'MALAT1_CPM' in sample_obs.columns and sample_obs['MALAT1_CPM'].notna().any():
        vals = sample_obs['MALAT1_CPM'].dropna()
        if len(vals) > 0:
            log_vals = np.log1p(vals)
            med = np.median(log_vals)
            mad = median_abs_deviation(log_vals, nan_policy='omit')
            thresholds['MALAT1_CPM_cutoffs'] = int(max(0, np.floor(np.expm1(med - nmads * mad))))
            thresholds['MALAT1_CPM_max_cutoffs'] = int(np.ceil(np.expm1(med + nmads * mad)))
        else:
            thresholds['MALAT1_CPM_cutoffs'] = '---'
            thresholds['MALAT1_CPM_max_cutoffs'] = '---'
    else:
        thresholds['MALAT1_CPM_cutoffs'] = '---'
        thresholds['MALAT1_CPM_max_cutoffs'] = '---'

    return thresholds


def compute_mad_thresholds(adata, QC_cutoff_df, nmads, key='sampleID', scope='per-sample'):
    """Compute MAD-based QC thresholds.

    Parameters
    ----------
    adata : AnnData
        Annotated data with QC metrics in obs.
    QC_cutoff_df : DataFrame
        Per-sample QC cutoff table (from Excel). A copy is returned with
        cell-level thresholds overwritten by MAD-derived values.
    nmads : float
        Number of MADs from the median for threshold computation.
    key : str
        Column in adata.obs identifying samples (default: 'sampleID').
    scope : str
        'per-sample' computes MAD per sample (default).
        'per-tissue' computes MAD across all cells and applies uniform
        thresholds (more robust, standard practice per OSCA/Luecken & Theis).

    Returns
    -------
    DataFrame
        Modified copy of QC_cutoff_df with MAD thresholds.
        Doublet cutoffs are preserved from the Excel table.
    """
    QC_cutoff = QC_cutoff_df.copy()

    if QC_cutoff.empty:
        raise ValueError(
            "QC cutoff table is empty for this tissue. "
            "Check that the tissue name in the config matches the Excel sheet."
        )

    # Doublet cutoffs: keep from Excel (bimodal distribution makes MAD unstable)
    if 'doublet_cutoffs' not in QC_cutoff.columns:
        QC_cutoff['doublet_cutoffs'] = 1.0

    if scope == 'per-tissue':
        # Compute thresholds once from all cells in the tissue
        thresholds = _compute_mad_for_sample(adata.obs, nmads)

        logging.info(f"Computing tissue-wide MAD thresholds (nmads={nmads}, {adata.n_obs} cells):")
        logging.info(f"  genes=[{thresholds['Min_genes_in_cells']}, {thresholds['Max_genes_in_cells']}],"
                     f" counts<{thresholds['Max_counts_in_cells']},"
                     f" mt<{thresholds['Max_percent_mt_in_cells']},"
                     f" ribo<{thresholds['Max_percent_ribo_in_cells']}]")

        # Apply same thresholds to every sample row
        for idx, row in QC_cutoff.iterrows():
            for col, val in thresholds.items():
                QC_cutoff.at[idx, col] = val
            logging.info(f"  {row['rnaID']}: doublet={row['doublet_cutoffs']}")
    else:
        # Compute thresholds per sample
        logging.info(f"Computing per-sample MAD thresholds (nmads={nmads}):")

        for idx, row in QC_cutoff.iterrows():
            sample_id = row['rnaID']
            sample_mask = adata.obs[key] == sample_id if key in adata.obs.columns else adata.obs['rnaID'] == sample_id
            sample_obs = adata.obs[sample_mask]

            if len(sample_obs) == 0:
                logging.warning(f"  No cells found for {sample_id} — skipping MAD")
                continue

            thresholds = _compute_mad_for_sample(sample_obs, nmads)
            for col, val in thresholds.items():
                QC_cutoff.at[idx, col] = val

            logging.info(f"  {sample_id} ({len(sample_obs)} cells):"
                         f" genes=[{thresholds['Min_genes_in_cells']}, {thresholds['Max_genes_in_cells']}],"
                         f" counts<{thresholds['Max_counts_in_cells']},"
                         f" mt<{thresholds['Max_percent_mt_in_cells']},"
                         f" ribo<{thresholds['Max_percent_ribo_in_cells']},"
                         f" doublet={row['doublet_cutoffs']}")

    return QC_cutoff


def resolve_doublet_column(cutoff, default="Yes"):
    """Pick which column to filter doublets on for one sample.

    Reads that sample's own 'Whether_use_double_GMM_method' QC cutoff value:
    'Yes' -> GMM-derived doublet_probabilities, 'No' -> raw doublet_score.
    Doublet detection/GMM conversion is computed per sample, so this is
    resolved per sample rather than once for the whole tissue.
    Case-insensitive; missing/blank/'---' values fall back to `default`.
    """
    raw = str(cutoff.get('Whether_use_double_GMM_method', default)).strip()
    value = raw.lower()
    if value in ("", "nan", "none", "---"):
        value = default.lower()
    elif value not in ("yes", "no"):
        logging.warning(f"Unrecognized Whether_use_double_GMM_method value '{raw}' "
                        f"(expected 'Yes' or 'No'). Defaulting to '{default}'.")
        value = default.lower()
    return 'doublet_probabilities' if value == 'yes' else 'doublet_score'


def prepare_upset_summary_allQC(sample_data, QC_cutoff_dict, global_obs, default_cutoffs):
    sampleID = sample_data['sampleID'].iloc[0].strip()

    if sampleID not in QC_cutoff_dict:
        raise KeyError(
            f"QC cutoffs missing for sampleID '{sampleID}'. Check agreement between tissue and sample"
        )

    # Create masks
    mask_min_genes = sample_data['n_genes_by_counts'] < int(QC_cutoff_dict[sampleID]['Min_genes_in_cells'])
    mask_max_genes = sample_data['n_genes_by_counts'] > int(QC_cutoff_dict[sampleID]['Max_genes_in_cells'])
    mask_max_counts = sample_data['total_counts'] > int(QC_cutoff_dict[sampleID]['Max_counts_in_cells'])
    if QC_cutoff_dict[sampleID]['Max_percent_mt_in_cells'] == '---':
        mask_pct_counts_mt = sample_data['pct_counts_mt'] > default_cutoffs['pct_counts_mt']
    else:
        mask_pct_counts_mt = sample_data['pct_counts_mt'] > float(QC_cutoff_dict[sampleID]['Max_percent_mt_in_cells'])
    if QC_cutoff_dict[sampleID]['Max_percent_ribo_in_cells'] == '---':
        mask_pct_counts_ribo = sample_data['pct_counts_ribo'] > default_cutoffs['pct_counts_ribo']
    else:
        mask_pct_counts_ribo = sample_data['pct_counts_ribo'] > float(QC_cutoff_dict[sampleID]['Max_percent_ribo_in_cells'])
    if 'pct_exon_reads' in sample_data:
        if QC_cutoff_dict[sampleID]['Exon_ratio_cutoffs'] == '---':
            Q75 = np.nanpercentile(global_obs['pct_exon_reads'], 75)
            Q25 = np.nanpercentile(global_obs['pct_exon_reads'], 25)
            cutoff = (Q75 + (Q75 - Q25))
            mask_pct_exon_reads = sample_data['pct_exon_reads'] > cutoff
        else:
            mask_pct_exon_reads = sample_data['pct_exon_reads'] > float(QC_cutoff_dict[sampleID]['Exon_ratio_cutoffs'])
    if QC_cutoff_dict[sampleID]['MALAT1_CPM_cutoffs'] == '---':
        mask_MALAT1_CPM = sample_data['MALAT1_CPM'] < default_cutoffs['MALAT1_CPM']
    else:
        mask_MALAT1_CPM = sample_data['MALAT1_CPM'] < int(QC_cutoff_dict[sampleID]['MALAT1_CPM_cutoffs'])

    if QC_cutoff_dict[sampleID]['MALAT1_CPM_max_cutoffs'] == '---':
        mask_MALAT1_max_CPM = sample_data['MALAT1_max_CPM'] > default_cutoffs['MALAT1_max_CPM']
    else:
        mask_MALAT1_max_CPM = sample_data['MALAT1_max_CPM'] > int(QC_cutoff_dict[sampleID]['MALAT1_CPM_max_cutoffs'])

    mask_doublets = sample_data['doublet_probabilities'] > float(QC_cutoff_dict[sampleID]['doublet_cutoffs'])
    # Indices
    condition_indices = {
        'Min_genes': sample_data.index[mask_min_genes],
        'Max_genes': sample_data.index[mask_max_genes],
        'Max_counts': sample_data.index[mask_max_counts],
        'pct_counts_mt': sample_data.index[mask_pct_counts_mt],
        'pct_counts_ribo': sample_data.index[mask_pct_counts_ribo],
        'MALAT1_CPM': sample_data.index[mask_MALAT1_CPM],
        'MALAT1_max_CPM': sample_data.index[mask_MALAT1_max_CPM],
        'doublets': sample_data.index[mask_doublets]
    }
    if 'pct_exon_reads' in sample_data:
        condition_indices['pct_exon_reads'] = sample_data.index[mask_pct_exon_reads]
    
    upset_data = pd.DataFrame(index=sample_data.index)
    for cond, indices in condition_indices.items():
        upset_data[cond] = 0
        upset_data.loc[indices, cond] = 1
    upset_data = upset_data[(upset_data.sum(axis=1) > 0)]
    return upset_data.groupby(list(condition_indices.keys())).size()

def prepare_upset_summary_filteredQC(sample_df, QC_cutoff_dict):
    sampleID = sample_df['sampleID'].iloc[0]
    cond = QC_cutoff_dict[sampleID]
    
    # Booleans
    mask_min_genes   = sample_df['n_genes_by_counts'] < int(cond['Min_genes_in_cells'])
    mask_max_genes   = sample_df['n_genes_by_counts'] > int(cond['Max_genes_in_cells'])
    mask_max_counts  = sample_df['total_counts'] > int(cond['Max_counts_in_cells'])
    mask_doublets    = sample_df['doublet_probabilities'] > float(cond['doublet_cutoffs'])
    # Optional
    condition_indices = {
        'Min_genes': mask_min_genes,
        'Max_genes': mask_max_genes,
        'Max_counts': mask_max_counts,
        'doublets': mask_doublets,
    }
    optional_fields = [
        ('Max_percent_mt_in_cells', 'pct_counts_mt'),
        ('Max_percent_ribo_in_cells', 'pct_counts_ribo'),
        ('Exon_ratio_cutoffs', 'pct_exon_reads'),
        ('MALAT1_CPM_cutoffs', 'MALAT1_CPM'),
        ('MALAT1_CPM_max_cutoffs', 'MALAT1_max_CPM')
    ]
    for key, col in optional_fields:
        val = cond[key]
        if val != '---':
            if col in ['MALAT1_CPM']:
                mask = sample_df[col] < int(val)
            elif col in ['MALAT1_max_CPM']:
                mask = sample_df[col] > int(val)
            else:
                mask = sample_df[col] > float(val)
            condition_indices[col] = mask
            
    # Construct binary df
    upset_data = pd.DataFrame(index=sample_df.index)
    for condition, mask in condition_indices.items():
        upset_data[condition] = mask.astype(int)
    upset_data = upset_data[(upset_data.sum(axis=1) > 0)]
    return upset_data.groupby(list(condition_indices.keys())).size()

def export_doublet_calls_by_sample(adata, QC_cutoff_dict, tissue, tissue_std, figdir):
    """
    For each sampleID in adata.obs, determines doublet calls using the supplied QC_cutoff_dict,
    and saves the results as TSV files in the specified figdir.
    """

    for sampleID in adata.obs['sampleID'].unique():
        if sampleID not in QC_cutoff_dict:
            logging.error(f"sampleID {sampleID} not found in QC_cutoff_dict.")
            continue

        cutoff = QC_cutoff_dict[sampleID]
        doublet_cutoff = float(cutoff['doublet_cutoffs'])
        doublet_col = resolve_doublet_column(cutoff)

        adata_sel = adata[adata.obs['sampleID'] == sampleID, :].copy()
        df_doublet = adata_sel.obs[['doublet_score', 'doublet_probabilities']].copy()
        df_doublet['doublet_call'] = df_doublet[doublet_col].apply(
            lambda x: 'yes' if x > doublet_cutoff else 'no'
        )

        doublet_outfile = os.path.join(figdir, f'RNA_doublet_results.{sampleID}.tsv')
        df_doublet.to_csv(doublet_outfile, sep='\t', index=True, header=True, index_label="cell_barcode")

        # plot the doublet distribution with filter cutoff
        plot_doublet_hist(
            adata=adata_sel, tissue=tissue, tissue_std=tissue_std,
            figdir=figdir, key='sampleID', probability_cutoff=doublet_cutoff, stage="filterqc",
            metric_col=doublet_col)


def filter_and_process_adata(adata, df, QC_cutoff_dict, tissue, tissue_std, log_file_dir, key):

    adata_filter = adata.copy()
    
    log_lines = []
    log_lines.append(f'total cellbarcodes for all samples: {len(adata.obs_names)}\n\n')
    
    for ID in adata.obs[key].unique():

        logging.info(f"Filtering {ID}...")
        sample_obs = adata_filter.obs[key]
        adata_process = adata_filter[sample_obs == ID,:].copy()
        adata_remain  = adata_filter[sample_obs != ID,:].copy()
        log_lines.append(f'total cellbarcodes in {ID}: {len(adata_process.obs_names)}\n')

        cutoff = QC_cutoff_dict[ID]
        doublet_col = resolve_doublet_column(cutoff)
        n_before = len(adata_process)

        adata_process = adata_process[adata_process.obs['n_genes_by_counts'] > int(cutoff['Min_genes_in_cells']), :]
        logging.info(f"  after Min_genes ({cutoff['Min_genes_in_cells']}): {n_before} -> {len(adata_process)} cells")
        adata_process = adata_process[adata_process.obs['n_genes_by_counts'] < int(cutoff['Max_genes_in_cells']), :]
        logging.info(f"  after Max_genes ({cutoff['Max_genes_in_cells']}): -> {len(adata_process)} cells")
        adata_process = adata_process[adata_process.obs['total_counts'] < int(cutoff['Max_counts_in_cells']), :]
        logging.info(f"  after Max_counts ({cutoff['Max_counts_in_cells']}): -> {len(adata_process)} cells")
        if cutoff['Max_percent_mt_in_cells'] != '---':
            adata_process = adata_process[adata_process.obs['pct_counts_mt'] < float(cutoff['Max_percent_mt_in_cells']), :]
            logging.info(f"  after pct_mt ({cutoff['Max_percent_mt_in_cells']}): -> {len(adata_process)} cells")
        if cutoff['Max_percent_ribo_in_cells'] != '---':
            adata_process = adata_process[adata_process.obs['pct_counts_ribo'] < float(cutoff['Max_percent_ribo_in_cells']), :]
            logging.info(f"  after pct_ribo ({cutoff['Max_percent_ribo_in_cells']}): -> {len(adata_process)} cells")
        if 'pct_exon_reads' in adata_process.obs:
            if cutoff['Exon_ratio_cutoffs'] != '---':
                adata_process = adata_process[adata_process.obs['pct_exon_reads'] < float(cutoff['Exon_ratio_cutoffs']), :]
                logging.info(f"  after Exon_ratio ({cutoff['Exon_ratio_cutoffs']}): -> {len(adata_process)} cells")
        if cutoff['MALAT1_CPM_cutoffs'] != '---':
            if adata_process.obs['MALAT1_CPM'].notna().any():
                adata_process = adata_process[adata_process.obs['MALAT1_CPM'] > int(cutoff['MALAT1_CPM_cutoffs']), :]
                logging.info(f"  after MALAT1_CPM min ({cutoff['MALAT1_CPM_cutoffs']}): -> {len(adata_process)} cells")
            else:
                logging.warning(f"Skipping MALAT1_CPM min filter for {ID}: all values are NaN (gene may not exist in reference)")
        if cutoff['MALAT1_CPM_max_cutoffs'] != '---':
            if adata_process.obs['MALAT1_CPM'].notna().any():
                adata_process = adata_process[adata_process.obs['MALAT1_CPM'] < int(cutoff['MALAT1_CPM_max_cutoffs']), :]
                logging.info(f"  after MALAT1_CPM max ({cutoff['MALAT1_CPM_max_cutoffs']}): -> {len(adata_process)} cells")
            else:
                logging.warning(f"Skipping MALAT1_CPM max filter for {ID}: all values are NaN")
        # Doublet filter: keep cells with NaN values (doublet detection may have failed)
        doublet_mask = adata_process.obs[doublet_col] < float(cutoff['doublet_cutoffs'])
        doublet_mask = doublet_mask | adata_process.obs[doublet_col].isna()
        adata_process = adata_process[doublet_mask, :]
        logging.info(f"  after doublet ({cutoff['doublet_cutoffs']}): -> {len(adata_process)} cells")
        log_lines.append(f'num of cellbarcodes after QC filtering in {ID}: {len(adata_process.obs_names)}\n\n')
        adata_filter = ad.concat([adata_process, adata_remain])
    log_lines.append(f'\nnum of cellbarcodes after QC filtering: {len(adata_filter.obs_names)}')

    # Save log
    os.makedirs(log_file_dir, exist_ok=True)
    logfile = os.path.join(log_file_dir, f"{tissue_std}_filtering_stat_counts.txt")
    with open(logfile, "w") as f:
        f.write("\n".join(log_lines))
        
    return adata_filter

def downstream_process(adata, tissue_std, figdir, all_colors, include_malat1=True):
    adata.layers["rawcounts"] = adata.X.copy()

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat')

    sc.pl.highly_variable_genes(adata, save=f'.highlyVariableGenes.filterqc.{tissue_std}.png', show=False)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save=f'.var_ratio.filterqc.{tissue_std}.png', show=False)

    # PCA colored by donorID/QC/batch
    pca_list = ["sampleID", "pct_counts_mt","pct_counts_ribo", "pct_exon_reads", "log10_MALAT1_CPM" ]
    if "pct_exon_reads" not in adata.obs:
        pca_list.remove('pct_exon_reads')
    if not include_malat1 or "log10_MALAT1_CPM" not in adata.obs:
        pca_list.remove('log10_MALAT1_CPM')
    sc.pl.pca(
        adata,
        color=pca_list,
        wspace=0.5,
        ncols=2,
        palette=all_colors,
        save=f'.pca_scatter.filterqc.{tissue_std}.png', show=False
    )

def compress_and_save_postqc_h5ad(adata, output_h5ad_dir, tissue_std, runtag):
    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    if 'CPM' in adata.layers and not scipy.sparse.issparse(adata.layers['CPM']):
        adata.layers['CPM'] = scipy.sparse.csr_matrix(adata.layers['CPM'])
    if 'rawcounts' in adata.layers and not scipy.sparse.issparse(adata.layers['rawcounts']):
        adata.layers['rawcounts'] = scipy.sparse.csr_matrix(adata.layers['rawcounts'])

    logging.info("Saving h5ad...")
    adata.write(os.path.join(output_h5ad_dir, f'{tissue_std}_GEX.filtered.{runtag}.h5ad'))

def run_per_tissue(working_df, output_h5ad_dir, workdir, QC_cutoff, tissue, my_color_palette, default_cutoffs, runtag,
                  key, use_mad=False, nmads=5.0, mad_scope='per-sample', sample_colors=None, donor_colors=None):
    tissue_std = standardize_tissue_name(tissue)

    adata_path = os.path.join(output_h5ad_dir, f"{tissue_std}_GEX.withQC.h5ad")
    if not os.path.exists(adata_path):
        logging.error(f"No h5ad for tissue {tissue} at {adata_path}.")
        return

    logging.info("Loading anndata object...")
    adata = sc.read_h5ad(adata_path)

    figdir = os.path.join(workdir, 'figures')
    os.makedirs(figdir, exist_ok=True)

    if use_mad:
        QC_cutoff = compute_mad_thresholds(adata, QC_cutoff, nmads, key=key, scope=mad_scope)
        mad_tsv = os.path.join(figdir, f"{tissue_std}_MAD_QC_cutoffs.nmads{nmads}.tsv")
        QC_cutoff.to_csv(mad_tsv, sep='\t', index=False)
        logging.info(f"MAD QC cutoffs written to {mad_tsv}")

    sc.settings.figdir = figdir

    logging.info(f'all figure plots by {key}')
    if key == 'donorID':
        all_colors = assign_donor_colors(working_df, donor_colors or {}, key=key, fallback_palette=my_color_palette)
        logging.info(f"use colors: {all_colors}")
    elif key == 'sampleID':
        all_colors = assign_donor_colors(working_df, sample_colors or {}, key='rnaID', fallback_palette=my_color_palette)
        logging.info(f"use colors: {all_colors}")
    else:
        logging.warning("need to edit for colors")

    QC_cutoff_dict = QC_cutoff.set_index('rnaID').T.to_dict()

    # Whether_use_double_GMM_method (from QC cutoff table) is resolved per sample:
    # "Yes" filters that sample's cells on the GMM-derived doublet_probabilities,
    # "No" filters directly on the raw doublet_score.
    doublet_metric = pd.Series(np.nan, index=adata.obs_names)
    for sample_id in adata.obs["sampleID"].unique():
        if sample_id not in QC_cutoff_dict:
            continue
        doublet_col = resolve_doublet_column(QC_cutoff_dict[sample_id])
        logging.info(f"  {sample_id}: doublet filtering metric = {doublet_col}")
        sample_mask = adata.obs["sampleID"] == sample_id
        doublet_metric.loc[sample_mask] = adata.obs.loc[sample_mask, doublet_col].values

    # For upset plot input summary
    filter_df = pd.DataFrame({
        "sampleID": adata.obs["sampleID"],
        "n_genes_by_counts": adata.obs["n_genes_by_counts"],
        "total_counts": adata.obs["total_counts"],
        "pct_counts_mt": adata.obs["pct_counts_mt"],
        "pct_counts_ribo": adata.obs["pct_counts_ribo"],
        "doublet_probabilities": doublet_metric
    })
    if 'pct_exon_reads' in adata.obs:
        filter_df["pct_exon_reads"] = adata.obs["pct_exon_reads"]
    if 'MALAT1_CPM' in adata.obs:
        filter_df["MALAT1_CPM"] = adata.obs["MALAT1_CPM"]
        filter_df["MALAT1_max_CPM"] = adata.obs["MALAT1_CPM"]

    # Skip MALAT1 plots/downstream coloring if all samples have both cutoffs set to '---'
    all_malat1_skipped = all(
        str(d.get('MALAT1_CPM_cutoffs', '---')) == '---' and str(d.get('MALAT1_CPM_max_cutoffs', '---')) == '---'
        for d in QC_cutoff_dict.values()
    )

    # 1. Plot upset (all metrics and filtered only)
    logging.info("Generating upset plot...")
    for ID in filter_df[key].unique():
        sample_df = filter_df[filter_df[key]==ID]
        logging.info(f"{ID}")
        up1 = prepare_upset_summary_allQC(sample_df, QC_cutoff_dict, adata.obs, default_cutoffs)
        plot_upset(up1, tissue, key, ID, tissue_std, adata,
                   os.path.join(figdir,f"QC_filtering_upset.all.filterqc.{tissue_std}.{ID}.png"))
        up2 = prepare_upset_summary_filteredQC(sample_df, QC_cutoff_dict)
        plot_upset(up2, tissue, key, ID, tissue_std, adata,
                   os.path.join(figdir, f"QC_filtering_upset.filteringOnly.filterqc.{tissue_std}.{ID}.png"))

    # 2. save doublet results
    logging.info("Export doublet information...")
    export_doublet_calls_by_sample(adata, QC_cutoff_dict, tissue, tissue_std, figdir)

    # 3. Filtering + logging
    logging.info("Filtering...")
    adata_filt = filter_and_process_adata(adata, working_df, QC_cutoff_dict, tissue, tissue_std, figdir, key)
    adata_filt.var = adata.var.copy()
    adata_filt.uns = adata.uns.copy()
    adata = adata_filt.copy()

    # 4. Filter genes with < min_cells, normalize, PCA, cluster, UMAP, save
    if adata.n_obs == 0:
        raise ValueError(f"No cells remaining after QC filtering for {tissue}. Check QC cutoffs.")
    sc.pp.filter_genes(adata, min_cells=int(QC_cutoff_dict[adata_filt.obs[key].unique()[0]]['Min_cells_for_genes']))

    # 5. normalization, feature selection, linear dimensional reduction
    logging.info("Post-filter processing...")
    downstream_process(adata, tissue_std, figdir, all_colors, include_malat1=not all_malat1_skipped)
    run_umap_clustering(adata, tissue, tissue_std, figdir)

    # Save output h5ad
    compress_and_save_postqc_h5ad(adata, output_h5ad_dir, tissue_std, runtag)

    # Move figures to new directory
    move_figures_to_newdir(workdir, old="figures", new=f"filter_qc.{runtag}")

    logging.info(f"Finished QC filtering and post-processing for {tissue}.")


def main(config_path, runtag, use_mad=False, nmads=5.0, mad_scope='per-sample'):
    setup_logging()
    config = load_config(config_path)
    require_keys(config, [
        "paths.workdir", "paths.output_h5ad_dir", "paths.sample_metadata",
        "params.tissue", "color", "my_color_palette", "qc.rna_qc_cutoff_table",
    ], context=config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']

    tissue = config['params']['tissue']

    donor_colors = config['color'].get("donor_colors")
    sample_colors = config['color'].get("sample_colors")
    my_color_palette = config["my_color_palette"]

    # Default QC cutoffs for fallback
    default_cutoffs = {
        'log10_n_genes_by_counts':[np.log10(200), np.log10(5000)],
        'log10_total_counts': np.log10(25000),
        'pct_counts_mt': 10,
        'pct_counts_ribo': 10,
        'MALAT1_CPM':10,
        'MALAT1_max_CPM': 1000000 # or remove entirely
    }

    # Load master sample metadata across tissues/donors
    sample_metadata = config['paths']['sample_metadata']
    logging.info(f"Loading sample metadata {sample_metadata}")
    df = pd.read_csv(sample_metadata, sep='\t', header=None, index_col=False,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])

    # Load qc cutoff table
    qc_cutoff_table = config['qc']['rna_qc_cutoff_table']
    logging.info(f"Loading QC cutoff table {qc_cutoff_table}")
    if qc_cutoff_table.endswith('.xlsx') or qc_cutoff_table.endswith('.xls'):
        df_cutoff_all = pd.read_excel(qc_cutoff_table,
                                      sheet_name=config['qc']['sheet_name'], engine='openpyxl')
    else:
        df_cutoff_all = pd.read_csv(qc_cutoff_table, sep='\t')

    failed_tissues = []
    if tissue == "---":
        tissues = sorted(df["tissue"].unique())
        logging.info(f"Running QC filtering for MULTIPLE tissues: {tissues}")

        for idx, tissue_name in enumerate(tissues, 1):
            logging.info(f"============== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==============")
            try:
                working_df = df[df["tissue"] == tissue_name]
                QC_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue_name]
                run_per_tissue(working_df, output_h5ad_dir, workdir, QC_cutoff, tissue_name, my_color_palette, default_cutoffs, runtag, "sampleID", use_mad=use_mad, nmads=nmads, mad_scope=mad_scope, sample_colors=sample_colors, donor_colors=donor_colors)
            except Exception as e:
                logging.error(f"QC filtering failed for {tissue_name}: {e}")
                failed_tissues.append(tissue_name)
    else:
        logging.info(f"============== Processing tissue: {tissue} ==============")
        try:
            working_df = df[df["tissue"] == tissue]
            QC_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue]
            run_per_tissue(working_df, output_h5ad_dir, workdir, QC_cutoff, tissue, my_color_palette, default_cutoffs, runtag, "sampleID", use_mad=use_mad, nmads=nmads, mad_scope=mad_scope, sample_colors=sample_colors, donor_colors=donor_colors)
        except Exception as e:
            logging.error(f"QC filtering failed for {tissue}: {e}")
            failed_tissues.append(tissue)

    if failed_tissues:
        logging.error(f"filter_qc failed for {len(failed_tissues)} tissue(s): {failed_tissues}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run QC filtering, normalization, feature selection for scRNA-seq h5ad.")
    parser.add_argument("config", help="YAML config file describing tissue, paths, etc.")
    parser.add_argument("runtag", help="Tag for this run (e.g. round3 or v1)")
    parser.add_argument("--use-mad", action="store_true",
                        help="Use MAD-based thresholds instead of Excel cutoffs for cell QC metrics")
    parser.add_argument("--nmads", type=float, default=5.0,
                        help="Number of MADs from median for threshold computation (default: 5.0)")
    parser.add_argument("--mad-scope", choices=["per-sample", "per-tissue"], default="per-sample",
                        help="Compute MAD per sample or across all cells in the tissue (default: per-sample)")
    args = parser.parse_args()
    try:
        main(args.config, args.runtag, use_mad=args.use_mad, nmads=args.nmads, mad_scope=args.mad_scope)
    except (FileNotFoundError, KeyError, ValueError) as e:
        logging.error(f"{type(e).__name__}: {e}")
        sys.exit(1)
