#!/usr/bin/env python3

"""
Author: Kaili Fan
ATAC_QC plotting functions.
"""

import os
import numpy as np
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
import seaborn as sns

import anndata as ad

def plot_kde_filter(
    adata,
    x_cutoff,
    y_cutoff,
    sample_info,
    initial_cell_str,
    cutoff_str,
    passed_cells_str,
    pdf=None,
    show=False,
    show_cutoff_line=True,
    ylim_post=(5, 30),
    ylim_pre=(0, 30)
):
    """
    Generalized KDE plot for QC filtering (before or after).
    Args:
        adata             : AnnData object (must have 'n_fragment' and 'tsse' in .obs)
        x_cutoff, y_cutoff: QC cutoff values (fragment, TSS enrichment)
        sample_info       : Sample string (e.g. "{tissue2}: {fileID}")
        initial_cell_str  : e.g. "Initial cell barcodes: 23432"
        cutoff_str        : e.g. "QC cutoffs: n_fragment > 1000, TSS_enrichment > 7"
        passed_cells_str  : e.g. "Cells passing QC: 19324"
        pdf               : PdfPages object to save (optional)
        show              : If True, also plt.show() (optional)
        show_cutoff_line  : If True, plot cutoff crosshairs (usually for before filter)
        ylim_post         : Y axis post-filter range (default (5, 30))
        ylim_pre          : Y axis pre-filter range (default (0, 30))
    Returns:
        fig               : The JointGrid matplotlib figure (so user can extend if needed)
    """

    if 'n_fragment' not in adata.obs or 'tsse' not in adata.obs:
        raise ValueError("AnnData must have 'n_fragment' and 'tsse' in .obs")

    plt.style.use("default")
    fig = sns.JointGrid()

    sns.kdeplot(
        x=np.log10(adata.obs['n_fragment']),
        y=adata.obs['tsse'],
        cmap='Blues',
        fill=True,
        ax=fig.ax_joint
    )

    fig.ax_joint.set_xlim(2, 5)
    fig.ax_joint.set_ylim(ylim_pre if show_cutoff_line else ylim_post)
    fig.ax_joint.set_xlabel("Number of unique fragments (log10)")
    fig.ax_joint.set_ylabel("TSS enrichment score")

    if show_cutoff_line:
        fig.refline(x=np.log10(x_cutoff), y=y_cutoff, color='r', linestyle='--')

    sns.histplot(x=np.log10(adata.obs['n_fragment']), kde=True, ax=fig.ax_marg_x)
    sns.histplot(y=adata.obs['tsse'], kde=True, ax=fig.ax_marg_y)

    main_title_line = sample_info + ('\n' if show_cutoff_line else '\npost-filter\n')
    fig.fig.suptitle(
        f"{main_title_line}\n{initial_cell_str}\n{cutoff_str}\n{passed_cells_str}",
        fontsize=14
    )
    plt.subplots_adjust(top=0.85)
    plt.tight_layout()
    if pdf is not None:
        pdf.savefig(fig.fig)
    if show:
        plt.show()
    plt.close(fig.fig)
    return fig  # in case user wants to use outside


def plot_doublet_score_probability(
    doublet_scores,
    doublet_probabilities,
    probability_midpoint,
    tissue,
    fileID,
    outdir=None,
    runtag=None,
    show=False
):
    """
    Generate histograms of doublet scores/probabilities, colored by doublet classification.
    If show=True, displays with plt.show(); if outdir and runtag are given, saves png.
    """
    doublet_mask = doublet_probabilities > probability_midpoint
    bins = np.histogram_bin_edges(doublet_scores, bins=50)
    hist_non_doublets, _ = np.histogram(doublet_scores[~doublet_mask], bins=bins)
    hist_doublets, _ = np.histogram(doublet_scores[doublet_mask], bins=bins)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    # Doublet scores
    axes[0].bar(bins[:-1], hist_non_doublets, width=np.diff(bins),
                edgecolor='black', color='grey', label='Non-doublets', alpha=0.7)
    axes[0].bar(bins[:-1], hist_doublets, width=np.diff(bins),
                edgecolor='black', color='red', bottom=hist_non_doublets, label='Doublets', alpha=0.7)
    axes[0].set_title(f'Doublets N={str(doublet_mask.sum())}')
    axes[0].set_xlabel('Doublet Score')
    axes[0].set_ylabel('Frequency')
    axes[0].legend()
    axes[0].grid(False)

    # Doublet probabilities
    prob_bins = np.histogram_bin_edges(doublet_probabilities, bins=30)
    prob_hist, _ = np.histogram(doublet_probabilities, bins=prob_bins)
    axes[1].bar(prob_bins[:-1], prob_hist, width=np.diff(prob_bins),
                edgecolor='black', color='grey', alpha=0.7)
    axes[1].axvline(probability_midpoint, color='red', linestyle='--', label='Doublet Threshold')
    axes[1].set_title(f'cutoff = {str(probability_midpoint)}')
    axes[1].set_xlabel('Doublet Probability')
    axes[1].set_ylabel('Frequency')
    axes[1].set_xlim(0, 1)
    axes[1].grid(False)

    fig.suptitle(f'Distribution of doublet scores and probability on {tissue}, {fileID}', fontsize=14)
    plt.tight_layout()

    # Save figure only if outdir and runtag are given
    if outdir is not None and runtag is not None:
        outpath = os.path.join(outdir, f"ATAC_QC_doubletHist.{runtag}.{tissue}-{fileID}.png")
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)


def cell_count_post_filter_hist(
    df_num_cells,
    tissue2,
    runtag,
    out_dir
):
    """
    Generate and save histogram of cell counts per sample after filtering, including total cell count in title.

    Args:
        df_num_cells: DataFrame with columns ["fileID", "numCells"]
        tissue2: Standardized tissue name.
        runtag: QC/processing tag.
        out_dir: Output dir for the plot and stats file.
    """
    tsv_path = os.path.join(out_dir, f"ATAC_NumCell.{runtag}.{tissue2}.tsv")
    png_path = os.path.join(out_dir, f'ATAC_NumCell.{tissue2}.png')
    df_num_cells.to_csv(tsv_path, sep="\t", index=False)

    total_cells = df_num_cells["numCells"].sum()

    fig, ax = plt.subplots()
    sns.histplot(df_num_cells['numCells'], binwidth=100, ax=ax)
    ax.set_xlabel("Number of Cells past QC")
    ax.set_title(f"{tissue2} ({total_cells:,} cells)\n{runtag}")
    fig.savefig(png_path)
    plt.close(fig)

def plot_per_sample_umap_clusters(
    sample_ids,
    h5ad_dir,
    run_tag,
    tissue_name,
    output_dir,
    doublet_rate_key="doublet_rate"
):
    """
    Plot UMAP clusters for each sample and save all plots to a single PDF file.

    Args:
        sample_ids (list of str): List of sample IDs (fileID)
        h5ad_dir (str): Directory containing processed h5ad files.
        run_tag (str): QC/processing tag (e.g. round1, round3, etc.).
        tissue_name (str): Standardized tissue name (for output file naming).
        output_dir (str): Directory to save the output PDF.
        doublet_rate_key (str): Key in .uns for doublet rate (default: "doublet_rate")
    """
    pdf_path = os.path.join(output_dir, f'ATAC_post_filter_clusters.{tissue_name}.pdf')
    with PdfPages(pdf_path) as pdf:
        for sample_id in sample_ids:
            h5ad_path = os.path.join(h5ad_dir, f'{sample_id}.final.{run_tag}.h5ad')
            if not os.path.exists(h5ad_path):
                print(f"  [SKIP] Missing {h5ad_path}")
                continue
            adata = ad.read_h5ad(h5ad_path)

            # Cluster assignments
            clusters = [int(x) for x in adata.obs['leiden']]
            unique_clusters = np.unique(clusters)
            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_clusters)))

            fig, ax = plt.subplots()
            for i, cluster in enumerate(unique_clusters):
                idx = np.where(np.array(clusters) == cluster)[0]
                ax.scatter(
                    adata.obsm['X_umap'][idx, 0],
                    adata.obsm['X_umap'][idx, 1],
                    c=[colors[i]],
                    label=cluster,
                    s=1
                )
            ax.legend(
                bbox_to_anchor=(1.18, 1),
                loc="upper right",
                title="Leiden Cluster",
                frameon=False
            )
            ax.set_xlabel('UMAP-1')
            ax.set_ylabel('UMAP-2')
            median_tss = np.median(adata.obs["tsse"])
            median_frags = np.median(adata.obs["n_fragment"])
            if doublet_rate_key in adata.uns:
                dbl_rate = adata.uns[doublet_rate_key]
                dbl_line = f"Doublet rate: {dbl_rate*100:.1f}%"
            else:
                dbl_line = ""
            ax.set_title(
                f'{tissue_name}: {sample_id}\n'
                f'Median TSS enrichment: {median_tss:.2f}\n'
                f'Median # fragments: {median_frags:.0f}\n'
                f'{dbl_line}',
                loc='center',
                fontsize=10
            )
            plt.subplots_adjust(top=0.7, right=0.7)
            pdf.savefig(fig)
            plt.close(fig)


def plot_umap_by_sample(adata, tissue2, runtag, outdir, sample_colors=None, suffix="bySample"):
    """
    Plot UMAP for all cells, colored by sample.
    """
    samples = adata.obs['sample']
    samplelist = sorted(set(samples))
    fig, ax = plt.subplots(figsize=(8, 5))
    plt.style.use("default")
    for sample in samplelist:
        idx = np.where(np.array(samples) == sample)[0]
        color = sample_colors[sample] if (sample_colors and sample in sample_colors) else None
        ax.scatter(adata.obsm['X_umap'][idx,0], adata.obsm['X_umap'][idx,1],
                   color=color, label=sample, s=0.8, alpha=0.3)  # <--- ONLY THIS CHANGED
    ax.legend(bbox_to_anchor=(1.8, 1), loc="upper right", title="Sample", frameon=False)
    ax.set_xlabel('UMAP-1', fontsize=12)
    ax.set_ylabel('UMAP-2', fontsize=12)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_title(f"dGTEx {tissue2}\n", fontsize=14)
    plt.subplots_adjust(right=0.7)
    outpng = os.path.join(outdir, f'dGTEx_ATAC_UMAP.{tissue2}_{suffix}.{runtag}.png')
    plt.savefig(outpng, bbox_inches='tight', dpi=800)
    plt.close(fig)

def plot_umap_by_donor(adata, tissue2, runtag, outdir, donor_colors=None, suffix="byDonor"):
    """
    Plot UMAP for all cells, colored by donor.
    """
    donors = adata.obs['donorID']
    donorlist = sorted(set(donors))
    fig, ax = plt.subplots(figsize=(6, 5))
    plt.style.use("default")
    for donor in donorlist:
        idx = np.where(np.array(donors) == donor)[0]
        color = donor_colors[donor] if (donor_colors and donor in donor_colors) else None
        ax.scatter(adata.obsm['X_umap'][idx,0], adata.obsm['X_umap'][idx,1],
                   color=color, label=donor, s=0.8, alpha=0.3)  # <-- changed here
    ax.legend(bbox_to_anchor=(1.8, 1), loc="upper right", title="Donor", frameon=False)
    ax.set_xlabel('UMAP-1', fontsize=12)
    ax.set_ylabel('UMAP-2', fontsize=12)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_title(f"dGTEx {tissue2}\n", fontsize=14)
    plt.subplots_adjust(right=0.7)
    outpng = os.path.join(outdir, f'dGTEx_ATAC_UMAP.{tissue2}_{suffix}.{runtag}.png')
    plt.savefig(outpng, bbox_inches='tight', dpi=800)
    plt.close(fig)

def plot_umap_by_samplepanel(adata, tissue2, runtag, outdir, sample_colors=None, ncol=4):
    """
    Plot one UMAP panel per sample (multi-panel plot), output as a single PNG.
    """
    samples = adata.obs['sample']
    samplelist = sorted(set(samples))
    import math
    nrow = math.ceil(len(samplelist)/ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(8*ncol, 9*nrow))
    axes = axes.flatten()
    plt.style.use("ggplot")
    for ax in axes:
        ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1], color="white", s=0.8, alpha=0.3)
    for i, sample in enumerate(samplelist):
        idx = np.where(np.array(samples) == sample)[0]
        color = sample_colors[sample] if (sample_colors and sample in sample_colors) else None
        axes[i].scatter(adata.obsm['X_umap'][idx,0], adata.obsm['X_umap'][idx,1],
                        color=color, label=sample, s=0.8, alpha=0.3)        # <-- changed here
        # legend
        legend_entries = [plt.Line2D([0], [0], marker='o', color=color, markerfacecolor=color, markersize=8, label=sample)]
        axes[i].legend(handles=legend_entries, loc="lower left", frameon=False)
        axes[i].set_xlabel('UMAP-1', fontsize=10)
        axes[i].set_ylabel('UMAP-2', fontsize=10)
        axes[i].set_title(f"{tissue2}\n{sample}\nN={len(idx)}\n", fontsize=12)
    outpng = os.path.join(outdir, f'dGTEx_ATAC_UMAP.{tissue2}_sampleBySide.{runtag}.png')
    plt.savefig(outpng, bbox_inches='tight', dpi=800)
    plt.close(fig)

def assign_colors(keys, palette="tab10"):
    """
    Assign colors from a matplotlib palette to an iterable of keys (sample/donor IDs).
    Returns a dict: {key: color}
    """
    import matplotlib.pyplot as plt
    keys = sorted(keys)
    if isinstance(palette, str):
        colors = plt.get_cmap(palette).colors
    else:
        colors = palette
    # Repeat palette if too few colors
    out = {k: colors[i % len(colors)] for i, k in enumerate(keys)}
    return out