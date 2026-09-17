#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Summarizes and plots CellRanger QC results.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import os
import sys
import argparse
import logging
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from rna_qc.utils import load_config, setup_logging, require_keys
from rna_qc.rna_plots import assign_colors

def load_cellranger_summary(sample_dict, datadir):
    """
    Loads and concatenates CellRanger summary tables for all samples.
    """
    all_summary = pd.DataFrame()
    for sampleID, donorID in sample_dict.items():
        summary_file = os.path.join(datadir, sampleID, 'metrics_summary.csv')
        if not os.path.exists(summary_file):
            logging.warning(f"{summary_file} not found, skipping.")
            continue
            
        sample_summary = pd.read_csv(summary_file, thousands=',')
        sample_summary.index = [donorID]
        sample_summary['sampleID'] = sampleID
        sample_summary['donorID'] = donorID
        all_summary = pd.concat([all_summary, sample_summary], sort=False)
    return all_summary


def plot_qc_metrics(summary_df, all_colors, tissue, outdir, suffix):
    """
    Plots CellRanger QC metrics in 3x2 pages.
    """
    columns_to_plot = list(summary_df.columns[0:6]) + list(summary_df.columns[16:19])
    pct_cols = ['Valid Barcodes', 'Sequencing Saturation', 'Fraction Reads in Cells']

    # Split into pages of 6
    for page, start in enumerate(range(0, len(columns_to_plot), 6)):
        page_cols = columns_to_plot[start:start+6]
        nrows = (len(page_cols) + 1) // 2
        fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(24, 3 * nrows))
        axes = axes.flatten()

        for i, col in enumerate(page_cols):
            ax = axes[i]
            colors = summary_df['sampleID'].map(all_colors)
            values = summary_df[col].copy()

            if col == 'Number of Reads':
                values = pd.to_numeric(values, errors='coerce') / 1_000_000
            elif col in pct_cols:
                values = values.apply(lambda x: float(str(x).replace('%', '').strip()) if isinstance(x, str) else x).fillna(0)
            ax.barh(summary_df['sampleID'], values, color=colors, height=0.5)

            for index, value in enumerate(values):
                label = (f'{value:.2f}M' if col == 'Number of Reads'
                    else f'{value:.2f}%' if col in pct_cols
                    else f'{value}')
                ax.text(value, index, f' {label}', va='center', color='black', fontsize=10)
            xlabel = (col + " (M)" if col == 'Number of Reads'
                else col + " (%)" if col in pct_cols
                else col)
            ax.set_xlabel(xlabel)

        # Hide empty subplots
        for j in range(len(page_cols), len(axes)):
            axes[j].set_visible(False)

        fig.suptitle(f'CellRanger QC Summary of {tissue}', fontsize=12)
        plt.subplots_adjust(hspace=0.6, wspace=0.6)
        page_suffix = f'.p{page+1}' if page > 0 else ''
        fig_fp = os.path.join(outdir, f'CellRanger_QC_summary{page_suffix}.{suffix}.png')
        fig.savefig(fig_fp, dpi=300, bbox_inches='tight')
        plt.close(fig)

def plot_read_mappability(summary_df, all_colors, tissue, outdir, suffix):
    """
    Plots CellRanger read mappability metrics in 2-column layout.
    """
    mappability_cols = list(summary_df.columns[9:16])
    nrows = (len(mappability_cols) + 1) // 2
    fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(20, 3 * nrows))
    axes = axes.flatten()
    for i, col in enumerate(mappability_cols):
        ax = axes[i]
        colors = summary_df['sampleID'].map(all_colors)
        values = summary_df[col].apply(lambda x: float(str(x).replace('%', '').strip()) if isinstance(x, str) else x).fillna(0)
        ax.barh(summary_df['sampleID'], values, color=colors, height=0.6)
        for index, value in enumerate(values):
            if value > 0:
                ax.text(value, index, f' {value:.2f}%', va='center', color='black', fontsize=10)
        ax.set_xlabel(f'{col} (%)')
    # Hide empty subplots
    for j in range(len(mappability_cols), len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(f'CellRanger read mappability of {tissue}', fontsize=14)
    plt.tight_layout()
    fig_fp = os.path.join(outdir, f'CellRanger_read_QC_summary.{suffix}.png')
    fig.savefig(fig_fp, dpi=300, bbox_inches='tight')
    plt.close(fig)


def main(config_path):
    setup_logging()
    config = load_config(config_path)
    require_keys(config, [
        "paths.workdir", "paths.input_dir", "paths.sample_metadata",
        "params.tissue", "params.suffix", "my_color_palette",
    ], context=config_path)

    workdir = config['paths']['workdir']
    input_dir = config['paths']['input_dir']
    outdir = os.path.join(workdir, "cellranger_qc")
    os.makedirs(outdir, exist_ok=True)

    tissue = config['params']['tissue']
    suffix = config['params']['suffix']

    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None, index_col=False,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    if tissue == "---":
        working_df = df
    else:
        working_df = df[df["tissue"] == tissue]

    tissues = sorted(working_df["tissue"].unique())
    logging.info(f"Working tissue: {', '.join(tissues)}")
    sample_list = working_df['rnaID'].unique()
    logging.info(f"Samples: {list(sample_list)}")

    logging.info(f"Loading CellRanger summaries for {len(sample_list)} samples...")
    sample_dict = dict(zip(working_df['rnaID'], working_df['donorID']))
    summary_df = load_cellranger_summary(sample_dict, input_dir)
    if summary_df.empty:
        logging.warning("No CellRanger summary data could be loaded.")
        logging.warning("Check config file to be sure tissue names are matching. Exiting.")
        sys.exit(1)

    # Get donor color map
    #donor_colors = config['color'].get("donor_colors")
    #all_colors = assign_donor_colors(summary_df, donor_colors)
    #sample_colors = {}
    #all_colors = assign_donor_colors(summary_df, sample_colors, 'sampleID')
    all_colors = assign_colors(sorted(set(summary_df['sampleID'])), palette=config["my_color_palette"])
    

    logging.info("Plotting general QC metrics...")
    plot_qc_metrics(summary_df, all_colors, tissue, outdir, suffix)

    logging.info("Plotting read mappability metrics...")
    plot_read_mappability(summary_df, all_colors, tissue, outdir, suffix)

    logging.info(f"QC summary plots saved in {outdir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarizes and plots CellRanger QC results.")
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    try:
        main(args.config)
    except (FileNotFoundError, KeyError, ValueError) as e:
        logging.error(f"{type(e).__name__}: {e}")
        sys.exit(1)
