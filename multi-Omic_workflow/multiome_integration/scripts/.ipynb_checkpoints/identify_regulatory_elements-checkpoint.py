#!/usr/bin/env python3

"""
Author: Kaili Fan
Identify consensus ATAC peaks (regulatory elements) from scATAC-seq data by cell type, with summary plots.
"""

import os
import sys
import argparse
import math
import snapatac2 as snap
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name
from multiome_integration.integration_plots import plot_celltype_property_barplot

def read_chrom_size(file_path):
    print(f"[INFO] Reading chromosome sizes from: {file_path}")
    chrom_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            chrom_name, size = parts
            chrom_dict[chrom_name] = int(size)
    print(f"[INFO] Loaded sizes for {len(chrom_dict)} chromosomes.")
    return chrom_dict

def call_peaks_by_cluster(atac, tissue_std, output_h5ad_dir, nthread=16):
    print(f"[INFO] Calling peaks by cluster for tissue: {tissue_std} (n_jobs={nthread})")
    snap.tl.macs3(atac, groupby='celltype_std', n_jobs=nthread)
    outfile = os.path.join(output_h5ad_dir, f'ATAC_withPeak.{tissue_std}.h5ad')
    atac.write(outfile, compression="gzip")
    print(f"[INFO] Wrote ATAC with peaks: {outfile}")

    raw_peak_dir = './raw_bed'
    os.makedirs(raw_peak_dir, exist_ok=True)
    for celltype in atac.uns['macs3']:
        fname = os.path.join(raw_peak_dir, f'raw_ATAC_peak.{tissue_std}.{celltype}.txt')
        atac.uns['macs3'][celltype].to_csv(fname, sep='\t', index=False)
        print(f"[INFO] Raw peaks for {celltype} saved: {fname}")
    print(f"[INFO] Finished peak calling for clusters in {tissue_std}")
    return atac

def call_and_save_consensus_peaks(atac, chrom_size, tissue_std, half_width=150):
    print(f"[INFO] Calling consensus peaks for tissue: {tissue_std} (half_width={half_width})")
    consensus_peaks = snap.tl.merge_peaks(atac.uns['macs3'], chrom_size, half_width=half_width)
    consensus_peaks_df = pd.DataFrame(consensus_peaks)
    consensus_peaks_df.columns = consensus_peaks.columns
    out_csv = f'consensus_peak_metatable.{tissue_std}.txt'
    consensus_peaks_df.to_csv(out_csv, sep='\t', index=False, header=True)
    print(f"[INFO] Consensus peaks saved: {out_csv} (n_peaks={len(consensus_peaks_df)})")
    return consensus_peaks_df

def plot_active_peaks(consensus_peaks_df, tissue_std, color_palette):
    print(f"[INFO] Plotting active peaks for tissue: {tissue_std}")
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
    print(f"[INFO] Barplot for active peaks saved: {out_png}")

def compute_and_save_cellbypeak_umap(atac, consensus_peaks_df, tissue_std, output_h5ad_dir):
    print(f"[INFO] Computing cell-by-peak matrix and UMAP for {tissue_std}")
    peak_mat = snap.pp.make_peak_matrix(atac, use_rep=consensus_peaks_df['Peaks'])
    print(f"[INFO] Peak matrix shape: {peak_mat.shape}")
    snap.pp.select_features(peak_mat, n_features=250000)
    snap.tl.spectral(peak_mat)
    snap.tl.umap(peak_mat)
    out_h5ad = os.path.join(output_h5ad_dir, f'ATAC_cellByPeak.{tissue_std}.h5ad')
    peak_mat.write(out_h5ad, compression="gzip")
    print(f"[INFO] Cell-by-peak matrix saved: {out_h5ad}")
    fig_out = f"_ATAC_cellByPeak.celltype.{tissue_std}.png"
    sc.pl.umap(
        peak_mat, color=["celltype_std"],
        save=fig_out
    )
    print(f"[INFO] UMAP plot saved: ./figures/{fig_out}")
    return peak_mat

def run_per_tissue(
    tissue, output_h5ad_dir, chrom_size, color_palette, nthread=16, half_width=150):
    tissue_std = standardize_tissue_name(tissue)
    print(f"\n============== [ATAC] Processing tissue: {tissue} ({tissue_std}) ==============")
    atac_path = os.path.join(output_h5ad_dir, f'ATAC.GLUE.{tissue_std}.h5ad')
    print(f"[INFO] Loading ATAC AnnData: {atac_path}")
    atac = ad.read_h5ad(atac_path)
    atac.obs['celltype_std'] = [standardize_tissue_name(celltype) for celltype in atac.obs['celltype_glue']]
    atac = call_peaks_by_cluster(atac, tissue_std, output_h5ad_dir, nthread=nthread)
    consensus_peaks_df = call_and_save_consensus_peaks(atac, chrom_size, tissue_std, half_width=half_width)
    plot_active_peaks(consensus_peaks_df, tissue_std, color_palette)
    compute_and_save_cellbypeak_umap(atac, consensus_peaks_df, tissue_std, output_h5ad_dir)
    print(f"[INFO] Done with {tissue}\n")

def main(config_path):
    config = load_config(config_path)

    output_h5ad_dir = config['paths']['output_h5ad_dir']
    peak_output_dir = os.path.join(config['paths']['workdir'], 'atac_peak')
    os.makedirs(peak_output_dir, exist_ok=True)
    os.makedirs(os.path.join(peak_output_dir, 'figures'), exist_ok=True)
    os.chdir(peak_output_dir)
    print(f"[INFO] Output path set to: {peak_output_dir}")

    tissue = config['params']['tissue']
    color_palette = config["my_color_palette"]
    nthread = config['params'].get('threads', 16)
    half_width = config['params'].get('peak_half_width', 150)

    chromsize_path = config['references']['chrom_size']
    chrom_size = read_chrom_size(chromsize_path)

    if tissue == "---":
        working_df = pd.read_csv(config['paths']['sample_metadata'], sep='\t', header=None,
                                 names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        tissues = sorted(working_df["tissue"].unique())
        for idx, working_tissue in enumerate(tissues, 1):
            print(f"\n============== [Main] Processing tissue: {working_tissue} ({idx}/{len(tissues)}) ==============")
            try:
                run_per_tissue(working_tissue, output_h5ad_dir, chrom_size, color_palette, 
                               nthread=nthread, half_width=half_width)
            except Exception as e:
                print(f"[ERROR] Encountered error for tissue {working_tissue}: {str(e)}")
    else:
        print(f"\n============== [Main] Processing tissue: {tissue} ==============")
        try:
            run_per_tissue(tissue, output_h5ad_dir, chrom_size, color_palette, 
                               nthread=nthread, half_width=half_width)
        except Exception as e:
            print(f"[ERROR] Encountered error for tissue {tissue}: {str(e)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify consensus ATAC peaks from scATAC-seq data by cluster/celltype (config-driven).")
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    main(args.config)