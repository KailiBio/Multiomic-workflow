#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Identify consensus ATAC peaks (regulatory elements) from scATAC-seq data by cell type, with summary plots.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import argparse
import logging
import math
import snapatac2 as snap
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name, setup_logging, require_keys
from multiome_integration.integration_plots import plot_celltype_property_barplot

def read_chrom_size(file_path):
    logging.info(f"Reading chromosome sizes from: {file_path}")
    chrom_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            chrom_name, size = parts
            chrom_dict[chrom_name] = int(size)
    logging.info(f"Loaded sizes for {len(chrom_dict)} chromosomes.")
    return chrom_dict

def call_peaks_by_cluster(atac, tissue_std, output_h5ad_dir, nthread=16):
    logging.info(f"Calling peaks by cluster for tissue: {tissue_std} (n_jobs={nthread})")
    snap.tl.macs3(atac, groupby='celltype_std', n_jobs=nthread)
    outfile = os.path.join(output_h5ad_dir, f'ATAC_withPeak.{tissue_std}.h5ad')
    atac.write(outfile)
    logging.info(f"Wrote ATAC with peaks: {outfile}")

    raw_peak_dir = './raw_bed'
    os.makedirs(raw_peak_dir, exist_ok=True)
    for celltype in atac.uns['macs3']:
        fname = os.path.join(raw_peak_dir, f'raw_ATAC_peak.{tissue_std}.{celltype}.txt')
        atac.uns['macs3'][celltype].to_csv(fname, sep='\t', index=False)
        logging.info(f"Raw peaks for {celltype} saved: {fname}")
    logging.info(f"Finished peak calling for clusters in {tissue_std}")
    return atac

def call_and_save_consensus_peaks(atac, chrom_size, tissue_std, half_width=150):
    logging.info(f"Calling consensus peaks for tissue: {tissue_std} (half_width={half_width})")
    consensus_peaks = snap.tl.merge_peaks(atac.uns['macs3'], chrom_size, half_width=half_width)
    consensus_peaks_df = pd.DataFrame(consensus_peaks)
    consensus_peaks_df.columns = consensus_peaks.columns
    
    out_csv = f'consensus_peak_metatable.{tissue_std}.txt'
    consensus_peaks_df.to_csv(out_csv, sep='\t', index=False, header=True)
    logging.info(f"Consensus peaks saved: {out_csv} (n_peaks={len(consensus_peaks_df)})")
    return consensus_peaks_df

def plot_active_peaks(consensus_peaks_df, tissue_std, color_palette):
    logging.info(f"Plotting active peaks for tissue: {tissue_std}")
    true_counts = consensus_peaks_df.iloc[:,1:].sum()
    df_active = pd.DataFrame({
        'celltype': true_counts.index,
        'Counts': true_counts.values
    }).sort_values('celltype')
    out_png = f'./figures/consensusPeak_counts_barplot.{tissue_std}.png'
    out_pdf = f'./figures/consensusPeak_counts_barplot.{tissue_std}.pdf'
    
    plot_celltype_property_barplot(
        df_active, 'Counts',
        label='# of active peaks (k)',
        title=f'Active peaks in {tissue_std} (N={len(consensus_peaks_df)})',
        outfile_png=out_png,
        outfile_pdf=out_pdf,
        colororder=color_palette
    )
    logging.info(f"Barplot for active peaks saved: {out_png}")

def compute_and_save_cellbypeak_umap(atac, consensus_peaks_df, tissue_std, output_h5ad_dir):
    logging.info(f"Computing cell-by-peak matrix and UMAP for {tissue_std}")
    peak_mat = snap.pp.make_peak_matrix(atac, use_rep=consensus_peaks_df['Peaks'])
    logging.info(f"Peak matrix shape: {peak_mat.shape}")
    snap.pp.select_features(peak_mat, n_features=250000)
    snap.tl.spectral(peak_mat)
    snap.tl.umap(peak_mat)
    out_h5ad = os.path.join(output_h5ad_dir, f'ATAC_cellByPeak.{tissue_std}.h5ad')
    peak_mat.write(out_h5ad, compression="gzip")
    logging.info(f"Cell-by-peak matrix saved: {out_h5ad}")
    fig_out = f"_ATAC_cellByPeak.celltype.{tissue_std}.png"
    sc.pl.umap(
        peak_mat, color=["celltype_std"],
        save=fig_out
    )
    logging.info(f"UMAP plot saved: ./figures/{fig_out}")
    return peak_mat

def run_per_tissue(tissue, output_h5ad_dir, chrom_size, color_palette, nthread=16, half_width=150):
    
    tissue_std = standardize_tissue_name(tissue)
    
    atac_path = os.path.join(output_h5ad_dir, f'ATAC.GLUE.{tissue_std}.h5ad')
    logging.info(f"Loading ATAC AnnData: {atac_path}")
    atac = ad.read_h5ad(atac_path)

    # standardize cell type names
    atac.obs['celltype_std'] = [standardize_tissue_name(celltype) for celltype in atac.obs['celltype_glue']]
    
    atac = call_peaks_by_cluster(atac, tissue_std, output_h5ad_dir, nthread=nthread)
    consensus_peaks_df = call_and_save_consensus_peaks(atac, chrom_size, tissue_std, half_width=half_width)
    plot_active_peaks(consensus_peaks_df, tissue_std, color_palette)
    compute_and_save_cellbypeak_umap(atac, consensus_peaks_df, tissue_std, output_h5ad_dir)
    logging.info(f"Done with {tissue}\n")

def main(config_path):
    setup_logging()
    config = load_config(config_path)
    require_keys(config, [
        "paths.workdir", "paths.output_h5ad_dir", "params.tissue",
        "my_color_palette", "references.chrom_size",
    ], context=config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    peak_output_dir = os.path.join(workdir, 'atac_peak')
    os.makedirs(peak_output_dir, exist_ok=True)
    os.makedirs(os.path.join(peak_output_dir, 'figures'), exist_ok=True)
    os.chdir(peak_output_dir)
    logging.info(f"Output path set to: {peak_output_dir}")

    tissue = config['params']['tissue']
    
    color_palette = config["my_color_palette"]
    nthread = config['params'].get('threads', 16)
    half_width = config['params'].get('peak_half_width', 150)

    chromsize_path = config['references']['chrom_size']
    chrom_size = read_chrom_size(chromsize_path)

    if tissue == "---":
        sample_metadata = config['paths']['sample_metadata']
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        
        tissues = sorted(df["tissue"].unique())
        logging.info(f"Calling ATAC Peaks for MULTIPLE tissues: {tissues}")
        
        failed_tissues = []
        for idx, tissue_name in enumerate(tissues, 1):
            logging.info(f"========== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")

            try:
                run_per_tissue(tissue_name, output_h5ad_dir, chrom_size, color_palette,
                               nthread=nthread, half_width=half_width)
            except Exception as e:
                logging.error(f"Encountered error for tissue {tissue_name}: {str(e)}")
                failed_tissues.append(tissue_name)
    else:
        logging.info(f"========== Processing tissue: {tissue} ==========")

        failed_tissues = []
        try:
            run_per_tissue(tissue, output_h5ad_dir, chrom_size, color_palette,
                               nthread=nthread, half_width=half_width)
        except Exception as e:
            logging.error(f"Encountered error for tissue {tissue}: {str(e)}")
            failed_tissues.append(tissue)

    if failed_tissues:
        logging.error(f"identify_regulatory_elements failed for {len(failed_tissues)} tissue(s): {failed_tissues}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify consensus ATAC peaks from scATAC-seq data by cluster/celltype (config-driven).")
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    try:
        main(args.config)
    except (FileNotFoundError, KeyError, ValueError) as e:
        logging.error(f"{type(e).__name__}: {e}")
        sys.exit(1)