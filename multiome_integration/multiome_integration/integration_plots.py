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

import anndata as ad

def plot_venn_and_save(rna_cells, atac_cells, sample, outdir, tissue):
    """
    Make and save a Venn diagram of RNA and ATAC cells for one sample.
    """
    venn2([set(rna_cells), set(atac_cells)], ('RNA cells', 'ATAC cells'))
    plt.title(f"Number of cells in {tissue}: {sample}")
    outfile = os.path.join(outdir, f"{tissue}_shared_cells_{sample}.venn.png")
    plt.savefig(outfile, dpi=800, bbox_inches='tight')
    plt.close()
    print(f"[plotting] Venn diagram saved: {outfile}")

def plot_overlap_bar(cell_counts, outdir, tissue):
    """
    Make and save bar plot of overlapping and modality-specific cells for all samples in one tissue.

    cell_counts: dict
        {sample: [rna_only, overlap, atac_only]}
    """
    bar_positions = []
    bar_positions_2 = []
    label_positions = []
    current_position = 0
    counts_1_2, counts_2_2 = [], []
    counts_2_3, counts_3_3 = [], []
    labels = []
    colors = ["#16b900","#ffa806","#0e3aa9"]  # adjust as needed

    for sample, counts in cell_counts.items():
        bar_positions.append(current_position)
        bar_positions_2.append(current_position + 0.8)
        label_positions.append(current_position + 0.4)
        labels.append(sample)

        total_1_2 = counts[0] + counts[1]
        counts_1_2.append((counts[0] / total_1_2) * 100 if total_1_2 else 0)
        counts_2_2.append((counts[1] / total_1_2) * 100 if total_1_2 else 0)

        total_2_3 = counts[1] + counts[2]
        counts_2_3.append((counts[1] / total_2_3) * 100 if total_2_3 else 0)
        counts_3_3.append((counts[2] / total_2_3) * 100 if total_2_3 else 0)

        current_position += 2

    fig_height = max(0.5 * len(cell_counts), 3)
    fig, ax = plt.subplots(figsize=(12, fig_height))
    bar1 = ax.barh(bar_positions, counts_2_2, height=0.5, color=colors[1], label='RNA cells overlap with ATAC')
    bar2 = ax.barh(bar_positions, counts_1_2, height=0.5, left=counts_2_2, color=colors[0], label='RNA-only cells')
    bar3 = ax.barh(bar_positions_2, counts_2_3, height=0.5, color=colors[1], label='ATAC cells overlap with RNA')
    bar4 = ax.barh(bar_positions_2, counts_3_3, height=0.5, left=counts_2_3, color=colors[2], label='ATAC-only cells')

    ax.set_xlabel('Percentage of cells (%)')
    ax.set_yticks(label_positions)
    ax.set_yticklabels(labels)
    ax.set_title(f'Cell counts among scRNA and scATAC grouped by sample for {tissue}')
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    for i in range(len(bar_positions)):
        if counts_2_2[i] > 0:
            ax.text(counts_2_2[i] / 2, bar_positions[i], f'{counts_2_2[i]:.1f}%', ha='center', va='center', color='black', fontsize=9)
        if counts_1_2[i] > 0:
            ax.text(counts_2_2[i] + counts_1_2[i] / 2, bar_positions[i], f'{counts_1_2[i]:.1f}%', ha='center', va='center', color='white', fontsize=9)
        if counts_2_3[i] > 0:
            ax.text(counts_2_3[i] / 2, bar_positions_2[i], f'{counts_2_3[i]:.1f}%', ha='center', va='center', color='black', fontsize=9)
        if counts_3_3[i] > 0:
            ax.text(counts_2_3[i] + counts_3_3[i] / 2, bar_positions_2[i], f'{counts_3_3[i]:.1f}%', ha='center', va='center', color='white', fontsize=9)

    plt.tight_layout()
    outpath = os.path.join(outdir, f"{tissue}_overlapPercent_vennBar.png")
    plt.savefig(outpath, dpi=800, bbox_inches='tight')
    plt.close()
    print(f"[plotting] Overlap bar plot saved: {outpath}")

