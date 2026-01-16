#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Filters cells by user-provided (per sample) QC cutoffs, generates Upset plots, processes filtered data, and saves results.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse
import warnings

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from rna_qc.utils import load_config, standardize_tissue_name
from rna_qc.rna_plots import assign_donor_colors, move_figures_to_newdir, plot_upset, run_umap_clustering

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

def export_doublet_calls_by_sample(adata, QC_cutoff_dict, figdir):
    """
    For each sampleID in adata.obs, determines doublet calls using the supplied QC_cutoff_dict,
    and saves the results as TSV files in the specified figdir.
    """

    for sampleID in adata.obs['sampleID'].unique():
        if sampleID not in QC_cutoff_dict:
            print(f"[ERROR] sampleID {sampleID} not found in QC_cutoff_dict.")
            continue

        cutoff = QC_cutoff_dict[sampleID]
        doublet_cutoff = float(cutoff['doublet_cutoffs'])
        
        adata_sel = adata[adata.obs['sampleID'] == sampleID, :].copy()
        df_doublet = adata_sel.obs[['doublet_score', 'doublet_probabilities']].copy()
        df_doublet['doublet_call'] = df_doublet['doublet_probabilities'].apply(
            lambda x: 'yes' if x > doublet_cutoff else 'no'
        )
        
        doublet_outfile = os.path.join(figdir, f'RNA_doublet_results.{sampleID}.tsv')
        df_doublet.to_csv(doublet_outfile, sep='\t', index=True, header=True, index_label="cell_barcode")


def filter_and_process_adata(adata, df, QC_cutoff_dict, tissue, tissue_std, log_file_dir, key):

    adata_filter = adata.copy()
    
    log_lines = []
    log_lines.append(f'total cellbarcodes for all samples: {len(adata.obs_names)}\n\n')
    
    for ID in adata.obs[key].unique():

        print(f"[INFO] Filtering {ID}...")
        sample_obs = adata_filter.obs[key]
        adata_process = adata_filter[sample_obs == ID,:].copy()
        adata_remain  = adata_filter[sample_obs != ID,:].copy()
        log_lines.append(f'total cellbarcodes in {ID}: {len(adata_process.obs_names)}\n')

        cutoff = QC_cutoff_dict[ID]

        adata_process = adata_process[adata_process.obs['n_genes_by_counts'] > int(cutoff['Min_genes_in_cells']), :]
        adata_process = adata_process[adata_process.obs['n_genes_by_counts'] < int(cutoff['Max_genes_in_cells']), :]
        adata_process = adata_process[adata_process.obs['total_counts'] < int(cutoff['Max_counts_in_cells']), :]
        if cutoff['Max_percent_mt_in_cells'] != '---':
            adata_process = adata_process[adata_process.obs['pct_counts_mt'] < float(cutoff['Max_percent_mt_in_cells']), :]
        if cutoff['Max_percent_ribo_in_cells'] != '---':
            adata_process = adata_process[adata_process.obs['pct_counts_ribo'] < float(cutoff['Max_percent_ribo_in_cells']), :]
        if 'pct_exon_reads' in adata_process.obs:
            if cutoff['Exon_ratio_cutoffs'] != '---':
                adata_process = adata_process[adata_process.obs['pct_exon_reads'] < float(cutoff['Exon_ratio_cutoffs']), :]
        if cutoff['MALAT1_CPM_cutoffs'] != '---':
            adata_process = adata_process[adata_process.obs['MALAT1_CPM'] > int(cutoff['MALAT1_CPM_cutoffs']), :]
        if cutoff['MALAT1_CPM_max_cutoffs'] != '---':
            adata_process = adata_process[adata_process.obs['MALAT1_CPM'] < int(cutoff['MALAT1_CPM_max_cutoffs']), :]
        adata_process = adata_process[adata_process.obs['doublet_probabilities'] < float(cutoff['doublet_cutoffs']), :]
        log_lines.append(f'num of cellbarcodes after QC filtering in {ID}: {len(adata_process.obs_names)}\n\n')
        adata_filter = ad.concat([adata_process, adata_remain])
    log_lines.append(f'\nnum of cellbarcodes after QC filtering: {len(adata_filter.obs_names)}')

    # Save log
    os.makedirs(log_file_dir, exist_ok=True)
    logfile = os.path.join(log_file_dir, f"{tissue_std}_filtering_stat_counts.txt")
    with open(logfile, "w") as f:
        f.write("\n".join(log_lines))
        
    return adata_filter

def downstream_process(adata, tissue_std, figdir, all_colors):
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

    print("[INFO] Saving h5ad...")
    adata.write(os.path.join(output_h5ad_dir, f'{tissue_std}_GEX.filtered.{runtag}.h5ad'))

def run_per_tissue(working_df, output_h5ad_dir, workdir, QC_cutoff, tissue, donor_colors, default_cutoffs, runtag,
                  key):
    tissue_std = standardize_tissue_name(tissue)

    adata_path = os.path.join(output_h5ad_dir, f"{tissue_std}_GEX.withQC.h5ad")
    if not os.path.exists(adata_path):
        print(f"[ERROR] No h5ad for tissue {tissue} at {adata_path}.")
        return

    print("[INFO] Loading anndata object...")
    adata = sc.read_h5ad(adata_path)

    figdir = os.path.join(workdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    os.chdir(workdir)

    print(f'[INFO] all figure plots by {key}')
    if key == 'donorID':
        all_colors = assign_donor_colors(working_df, donor_colors, key = key)
        print(f"use colors: {all_colors}")
    elif key == 'sampleID':
        sample_colors={}
        all_colors = assign_donor_colors(working_df, sample_colors, key = 'rnaID')
        print(f"use colors: {all_colors}")
    else:
        print("[WARNING] need to edit for colors")

    QC_cutoff_dict = QC_cutoff.set_index('rnaID').T.to_dict()

    # For upset plot input summary
    filter_df = pd.DataFrame({
        "sampleID": adata.obs["sampleID"],
        "n_genes_by_counts": adata.obs["n_genes_by_counts"],
        "total_counts": adata.obs["total_counts"],
        "pct_counts_mt": adata.obs["pct_counts_mt"],
        "pct_counts_ribo": adata.obs["pct_counts_ribo"],
        "doublet_probabilities": adata.obs["doublet_probabilities"]
    })
    if 'pct_exon_reads' in adata.obs:
        filter_df["pct_exon_reads"] = adata.obs["pct_exon_reads"]
    if 'MALAT1_CPM' in adata.obs:
        filter_df["MALAT1_CPM"] = adata.obs["MALAT1_CPM"]
        filter_df["MALAT1_max_CPM"] = adata.obs["MALAT1_CPM"]

    plotlist = ["leiden", "log10_total_counts", "log10_n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo",
                "pct_exon_reads", "log10_MALAT1_CPM", "doublet_score", "doublet_probabilities"]
    if 'pct_exon_reads' not in adata.obs:
        plotlist.remove('pct_exon_reads')
        for d in QC_cutoff_dict.values():
            d.pop('pct_exon_reads', None)

    if 'MALAT1_CPM' not in adata.obs:
        plotlist.remove('log10_MALAT1_CPM')
        for d in QC_cutoff_dict.values():
            d.pop('log10_MALAT1_CPM', None)

    # 1. Plot upset (all metrics and filtered only)
    print("[INFO] Generating upset plot...")
    for ID in filter_df[key].unique():
        sample_df = filter_df[filter_df[key]==ID]
        print(f"{ID}")
        up1 = prepare_upset_summary_allQC(sample_df, QC_cutoff_dict, adata.obs, default_cutoffs)
        plot_upset(up1, tissue, key, ID, tissue_std, adata, 
                   os.path.join(figdir,f"QC_filtering_upset.all.filterqc.{tissue_std}.{ID}.png"))
        up2 = prepare_upset_summary_filteredQC(sample_df, QC_cutoff_dict)
        plot_upset(up2, tissue, key, ID, tissue_std, adata, 
                   os.path.join(figdir, f"QC_filtering_upset.filteringOnly.filterqc.{tissue_std}.{ID}.png"))

    # 2. save doublet results
    print("[INFO] Export doublet information...")
    export_doublet_calls_by_sample(adata, QC_cutoff_dict, figdir)
    
    # 3. Filtering + logging
    print("[INFO] Filtering...")
    adata_filt = filter_and_process_adata(adata, working_df, QC_cutoff_dict, tissue, tissue_std, figdir, key)
    adata_filt.var = adata.var.copy()
    adata_filt.uns = adata.uns.copy()
    adata = adata_filt.copy()

    # 4. Filter genes with < min_cells, normalize, PCA, cluster, UMAP, save
    sc.pp.filter_genes(adata, min_cells=int(QC_cutoff_dict[adata_filt.obs[key].unique()[0]]['Min_cells_for_genes']))

    # 5. normalization, feature selection, linear dimensional reduction
    print("[INFO] Post-filter processing...")
    downstream_process(adata, tissue_std, figdir, all_colors)
    run_umap_clustering(adata, tissue, tissue_std, figdir, plotlist)

    # Save output h5ad
    compress_and_save_postqc_h5ad(adata, output_h5ad_dir, tissue_std, runtag)

    # Move figures to new directory
    move_figures_to_newdir(workdir, old="figures", new=f"filter_qc.{runtag}")

    print(f"[INFO] Finished QC filtering and post-processing for {tissue}.")


def main(config_path, runtag):
    config = load_config(config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']

    tissue = config['params']['tissue']

    donor_colors = config['color'].get("donor_colors")

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
    print(f"[INFO] Loading sample metadata {sample_metadata}")
    df = pd.read_csv(sample_metadata, sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
   
    # Load qc cutoff table
    qc_cutoff_table = config['qc']['rna_qc_cutoff_table']
    print(f"[INFO] Loading QC cutoff table {qc_cutoff_table}")
    if qc_cutoff_table.endswith('.xlsx') or qc_cutoff_table.endswith('.xls'):
        df_cutoff_all = pd.read_excel(qc_cutoff_table,
                                      sheet_name=config['qc']['sheet_name'], engine='openpyxl')
    else:
        df_cutoff_all = pd.read_csv(qc_cutoff_table, sep='\t')

    if tissue == "---":
        tissues = sorted(df["tissue"].unique())
        print(f"[INFO] Running QC filtering for MULTIPLE tissues: {tissues}")
        
        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n============== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==============")
            try:
                working_df = df[df["tissue"] == tissue_name]
                QC_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue_name]
                run_per_tissue(working_df, output_h5ad_dir, workdir, QC_cutoff, tissue_name, donor_colors, default_cutoffs, runtag, "sampleID")
            except Exception as e:
                print(f"[ERROR] QC filtering failed for {tissue_name}: {e}")
    else:
        print(f"\n============== Processing tissue: {tissue} ==============")
        try:
            working_df = df[df["tissue"] == tissue]
            QC_cutoff = df_cutoff_all[df_cutoff_all['Tissue'] == tissue]
            run_per_tissue(working_df, output_h5ad_dir, workdir, QC_cutoff, tissue, donor_colors, default_cutoffs, runtag, "sampleID")
        except Exception as e:
            print(f"[ERROR] QC filtering failed for {tissue}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run QC filtering, normalization, feature selection for scRNA-seq h5ad.")
    parser.add_argument("config", help="YAML config file describing tissue, paths, etc.")
    parser.add_argument("runtag", help="Tag for this run (e.g. round3 or v1)")
    args = parser.parse_args()
    main(args.config, args.runtag)
