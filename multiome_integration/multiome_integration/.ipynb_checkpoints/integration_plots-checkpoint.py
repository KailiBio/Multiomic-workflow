#!/usr/bin/env python3

"""
Author: Kaili Fan
Multiome integration plotting functions.
"""

import os
import numpy as np
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib_venn import venn2
import seaborn as sns
import math

import pandas as pd
import anndata as ad

def plot_venn_and_save(rna_cells, atac_cells, rna_sample, atac_sample, outdir, tissue):
    """
    Make and save a Venn diagram of RNA and ATAC cells for one sample.
    """
    venn2([set(rna_cells), set(atac_cells)], ('RNA cells', 'ATAC cells'))
    plt.title(f"Number of cells in{tissue}\n{rna_sample}\n{atac_sample}")
    sample_new = str(rna_sample).split("-GEX")[0]
    outfile = os.path.join(outdir, f"shared_cells_venn.{sample_new}.png")
    plt.savefig(outfile, dpi=800, bbox_inches='tight')
    plt.close()
    print(f"[plotting] Venn diagram saved: {outfile}")


def plot_overlap_bar(cell_counts, outdir, suffix):
    """
    Make and save bar plot of overlapping and modality-specific cells for all samples in one tissue.

    cell_counts: list
        Each element: [rna_sample, atac_sample, rna_only, overlap, atac_only]
    """

    bar_positions = []
    bar_positions_2 = []
    label_positions = []
    current_position = 0
    counts_1_2, counts_2_2 = [], []
    counts_2_3, counts_3_3 = [], []
    labels = []
    colors = ["#16b900", "#ffa806", "#0e3aa9"]

    for row in cell_counts:
        rna_sample = row['rna_sample']
        atac_sample = row['atac_sample']
        rna_only = row['rna_only']
        overlap = row['overlap']
        atac_only = row['atac_only']

        sample_new = str(rna_sample).split("-GEX")[0]
        bar_positions.append(current_position)
        bar_positions_2.append(current_position + 0.8)
        label_positions.append(current_position + 0.4)
        labels.append(sample_new)

        total_1_2 = rna_only + overlap
        counts_1_2.append((rna_only / total_1_2) * 100 if total_1_2 else 0)
        counts_2_2.append((overlap / total_1_2) * 100 if total_1_2 else 0)

        total_2_3 = overlap + atac_only
        counts_2_3.append((overlap / total_2_3) * 100 if total_2_3 else 0)
        counts_3_3.append((atac_only / total_2_3) * 100 if total_2_3 else 0)

        current_position += 2

    # Only use the actual number of successfully processed bars!
    N = len(labels)
    if N == 0:
        print("[plotting] No samples to plot after filtering.")
        return

    fig_height = max(0.5 * N, 3)
    fig, ax = plt.subplots(figsize=(12, fig_height))

    bar1 = ax.barh(bar_positions, counts_2_2, height=0.5, color=colors[1], label='RNA cells overlap with ATAC')
    bar2 = ax.barh(bar_positions, counts_1_2, height=0.5, left=counts_2_2, color=colors[0], label='RNA-only cells')
    bar3 = ax.barh(bar_positions_2, counts_2_3, height=0.5, color=colors[1], label='ATAC cells overlap with RNA')
    bar4 = ax.barh(bar_positions_2, counts_3_3, height=0.5, left=counts_2_3, color=colors[2], label='ATAC-only cells')

    ax.set_xlabel('Percentage of cells (%)')
    ax.set_yticks(label_positions)
    ax.set_yticklabels(labels)
    ax.set_title(f'Cell counts among scRNA and scATAC grouped by sample for {suffix}')
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    for i in range(N):
        if counts_2_2[i] > 0:
            ax.text(counts_2_2[i] / 2, bar_positions[i], f'{counts_2_2[i]:.1f}%', ha='center', va='center', color='black', fontsize=9)
        if counts_1_2[i] > 0:
            ax.text(counts_2_2[i] + counts_1_2[i] / 2, bar_positions[i], f'{counts_1_2[i]:.1f}%', ha='center', va='center', color='white', fontsize=9)
        if counts_2_3[i] > 0:
            ax.text(counts_2_3[i] / 2, bar_positions_2[i], f'{counts_2_3[i]:.1f}%', ha='center', va='center', color='black', fontsize=9)
        if counts_3_3[i] > 0:
            ax.text(counts_2_3[i] + counts_3_3[i] / 2, bar_positions_2[i], f'{counts_3_3[i]:.1f}%', ha='center', va='center', color='white', fontsize=9)

    plt.tight_layout()
    outpath = os.path.join(outdir, f"multiome_overlapPercent_vennBar.{suffix}.png")
    plt.savefig(outpath, dpi=800, bbox_inches='tight')
    plt.close()
    print(f"[plotting] Overlap bar plot saved: {outpath}")

def plot_multiome_celltype_count_bar(
    rna,
    atac,
    figures_dir,
    suffix,
    color_map={'RNA': '#16b900', 'ATAC': '#0e39a9'}
):
    """
    Generate a horizontal barplot comparing cell counts per cell type between
    RNA and ATAC modalities, and save as PNG/PDF. Useful for multiome integration QC.

    Parameters
    ----------
    rna : AnnData
        RNA AnnData object with 'celltype_glue' in obs.
    atac : AnnData
        ATAC AnnData object with 'celltype_glue' in obs.
    figures_dir : str
        Directory to save the figures.
    suffix : str
        Suffix for output file names (usually tissue or sample name).
    color_map : dict
        (Optional) Mapping for bar colors by modality.
    """
    rna_counts = rna.obs['celltype_glue'].value_counts().sort_index()
    atac_counts = atac.obs['celltype_glue'].value_counts().sort_index()
    dfstat = pd.DataFrame({'RNA': rna_counts, 'ATAC': atac_counts}).fillna(0).astype(int)
    dfstat.index.name = 'celltype_glue'
    df_melted = dfstat.reset_index().melt(
        id_vars='celltype_glue', var_name='modality', value_name='count'
    )
    plt.figure(figsize=(12, 0.5 * len(dfstat)))
    ax = sns.barplot(
        data=df_melted,
        y='celltype_glue',
        x='count',
        hue='modality',
        orient='h',
        palette=color_map
    )
    plt.xlabel('Number of cells')
    plt.ylabel('')
    plt.tight_layout()
    sns.despine()
    plt.savefig(os.path.join(figures_dir, f'multiome_celltype_counts.{suffix}.png'), dpi=300)
    plt.savefig(os.path.join(figures_dir, f'multiome_celltype_counts.{suffix}.pdf'))
    plt.close()

def plot_celltype_property_barplot(
    df, value_col, label, title,
    outfile_png, outfile_pdf, colororder=None
):
    """
    Generate a horizontal barplot to visualize a property (e.g., peak counts)
    for each cell type as a bar. Used for peak statistics or other celltype-level metrics.

    Parameters
    ----------
    df : DataFrame
        DataFrame with columns 'celltype' and the property to plot.
    value_col : str
        Column in df to show as bar length.
    label : str
        X-axis label (e.g., '# of peaks (k)')
    title : str
        Figure title.
    outfile_png : str
        Output PNG file path.
    outfile_pdf : str
        Output PDF file path.
    colororder : list
        (Optional) List of colors for bars.
    """
    plt.figure(figsize=(8,5))
    if colororder is not None:
        colors = colororder
    else:
        cmap = plt.get_cmap('tab20')
        colors = [cmap(i) for i in range(len(df))]
    plt.barh(df['celltype'], df[value_col]/1000, color=colors)
    mean_count = math.ceil(df[value_col].mean()/1000)
    median_count = math.ceil(df[value_col].median()/1000)
    plt.xlabel(label)
    plt.title(title + f"\n(mean = {mean_count}k    median = {median_count}k)")
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.savefig(outfile_png, dpi=800, bbox_inches='tight')
    plt.savefig(outfile_pdf, dpi=800, bbox_inches='tight')
    plt.close()