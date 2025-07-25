#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Summarizes and plots CellRanger QC results.
"""

import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from rna_qc.utils import load_config
from rna_qc.rna_plots import assign_donor_colors

def load_cellranger_summary(sample_dict, datadir):
    """
    Loads and concatenates CellRanger summary tables for all samples.
    """
    all_summary = pd.DataFrame()
    for sampleID, donorID in sample_dict.items():
        summary_file = os.path.join(datadir, sampleID, 'metrics_summary.csv')
        if not os.path.exists(summary_file):
            print(f"[WARNING] {summary_file} not found, skipping.")
            continue
            
        sample_summary = pd.read_csv(summary_file, thousands=',')
        sample_summary.index = [donorID]
        sample_summary['sampleID'] = sampleID
        sample_summary['donorID'] = donorID
        all_summary = pd.concat([all_summary, sample_summary], sort=False)
    return all_summary


def plot_qc_metrics(summary_df, all_colors, tissue, outdir, suffix):
    """
    Plots CellRanger QC metrics.
    """
    columns_to_plot = list(summary_df.columns[0:6]) + list(summary_df.columns[16:19])
    
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 9))
    axes = axes.flatten()
    for i, col in enumerate(columns_to_plot):
        ax = axes[i]
        colors = summary_df['donorID'].map(all_colors)
        values = summary_df[col].copy()

        # plot horizontal bar
        if col == 'Number of Reads':
            values = pd.to_numeric(values, errors='coerce') / 1_000_000
        elif col in ['Valid Barcodes', 'Sequencing Saturation', 'Fraction Reads in Cells']:
            values = values.apply(lambda x: float(str(x).replace('%', '').strip()) if isinstance(x, str) else x).fillna(0)
        ax.barh(summary_df['donorID'], values, color=colors, height=0.5)

        # add text
        for index, value in enumerate(values):
            label = (f'{value:.2f}' if col == 'Number of Reads'
                else f'{value:.2f}%' if col in ['Valid Barcodes', 'Sequencing Saturation', 'Fraction Reads in Cells']
                else f'{value}')
            ax.text(value, index, label, va='center', color='black', fontsize=10)
        xlabel = (col + " (M)" if col == 'Number of Reads'
            else col + " (%)" if col in ['Valid Barcodes', 'Sequencing Saturation', 'Fraction Reads in Cells']
            else col)
        ax.set_xlabel(xlabel)
        
    fig.suptitle(f'CellRanger QC Summary of {tissue}', fontsize=12)
    plt.subplots_adjust(hspace=0.6, wspace=0.6)
    fig_fp = os.path.join(outdir, f'CellRanger_QC_summary.{suffix}.png')
    fig.savefig(fig_fp, dpi=300, bbox_inches='tight')
    plt.close(fig)

def plot_read_mappability(summary_df, all_colors, tissue, outdir, suffix):
    """
    Plots CellRanger read mappability metrics.
    """
    fig, axes = plt.subplots(nrows=7, ncols=1, figsize=(6, 9), sharex=True)
    mappability_cols = list(summary_df.columns[9:16])
    for i, (ax, col) in enumerate(zip(axes, mappability_cols)):
        colors = summary_df['donorID'].map(all_colors)
        values = summary_df[col].apply(lambda x: float(str(x).replace('%', '').strip()) if isinstance(x, str) else x).fillna(0)
        ax.barh(summary_df['donorID'], values, color=colors, height=0.6)
        for index, value in enumerate(values):
            if value > 0:
                ax.text(value, index, f'{value:.2f}%', va='center', color='black', fontsize=10)
        ax.set_xlabel(f'{col} (%)')
    fig.suptitle(f'CellRanger read mappability of {tissue}', fontsize=14)
    plt.tight_layout()
    fig_fp = os.path.join(outdir, f'CellRanger_read_QC_summary.{suffix}.png')
    fig.savefig(fig_fp, dpi=300, bbox_inches='tight')
    plt.close(fig)


def main(config_path):
    config = load_config(config_path)

    workdir = config['paths']['workdir']
    input_dir = config['paths']['input_dir']
    outdir = os.path.join(workdir, "cellranger_qc")
    os.makedirs(outdir, exist_ok=True)

    tissue = config['params']['tissue']
    suffix = config['params']['suffix']

    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    if tissue == "---":
        working_df = df
    else:
        working_df = df[df["tissue"] == tissue]

    tissues = sorted(working_df["tissue"].unique())
    print(f"Working tissue: {', '.join(tissues)}")
    sample_list = working_df['rnaID'].unique()
    print(f"Samples: {list(sample_list)}")

    print(f"\n[INFO] Loading CellRanger summaries for {len(sample_list)} samples...")
    sample_dict = dict(zip(working_df['rnaID'], working_df['donorID']))
    summary_df = load_cellranger_summary(sample_dict, input_dir)
    if summary_df.empty:
        print("[WARNING] No CellRanger summary data could be loaded. Exiting.")
        sys.exit(1)

    # Get donor color map
    donor_colors = config['color'].get("donor_colors")
    all_colors = assign_donor_colors(summary_df, donor_colors)

    print("[INFO] Plotting general QC metrics...")
    plot_qc_metrics(summary_df, all_colors, tissue, outdir, suffix)

    print("[INFO] Plotting read mappability metrics...")
    plot_read_mappability(summary_df, all_colors, tissue, outdir, suffix)

    print(f"\n[INFO] QC summary plots saved in {outdir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarizes and plots CellRanger QC results.")
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    main(args.config)