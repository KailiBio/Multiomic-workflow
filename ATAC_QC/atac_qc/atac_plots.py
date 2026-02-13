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
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
import seaborn as sns
import math

import anndata as ad

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import standardize_tissue_name

def assign_donor_colors(df, donor_col, key='donorID'):
    """
    Returns a {donor: color} dict for only donors in df[key].unique().
    For donors not in donor_col, assigns extra colors from a colormap.
    """
    donors_in_data = list(df[key].unique())
    color_map = dict(donor_col) 
    unknown_donors = [d for d in donors_in_data if d not in donor_col]

    if unknown_donors:
        colormap = plt.cm.get_cmap('tab20', len(unknown_donors))
        hex_colormap = [mcolors.rgb2hex(colormap(i)) for i in range(colormap.N)]
        for i, donor in enumerate(unknown_donors):
            color_map[donor] = hex_colormap[i]

    # Only keep donors present in df
    return {donor: color_map[donor] for donor in donors_in_data}

def assign_colors(keys, palette="tab10"):
    """
    Assign hex color codes from a matplotlib palette to an iterable of keys (sample/donor IDs).
    Returns a dict: {key: hex_color}
    """
    keys = sorted(keys)
    if isinstance(palette, str):
        colors = plt.get_cmap(palette).colors
    else:
        colors = palette

    def as_hex(c):
        # If already a hex string, just return
        if isinstance(c, str):
            if c.startswith("#") and (len(c) == 7 or len(c) == 9): 
                return c
        return mcolors.to_hex(c)

    hex_colors = [as_hex(c) for c in colors]

    out = {k: hex_colors[i % len(hex_colors)] for i, k in enumerate(keys)}
    return out
    
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
        sample_info       : Sample string (e.g. "{suffix}: {fileID}")
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
        outpath = os.path.join(outdir, f"ATAC_QC_doubletHist.{tissue}.{fileID}.{runtag}.png")
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)


def cell_count_post_filter_hist(
    df_num_cells,
    suffix,
    runtag,
    out_dir
):
    """
    Generate and save histogram of cell counts per sample after filtering, including total cell count in title.

    Args:
        df_num_cells: DataFrame with columns ["fileID", "numCells"]
        suffix: Suffix for naming.
        runtag: QC/processing tag.
        out_dir: Output dir for the plot and stats file.
    """
    tsv_path = os.path.join(out_dir, f"ATAC_NumCell.{runtag}.{suffix}.tsv")
    png_path = os.path.join(out_dir, f'ATAC_NumCell.{suffix}.png')
    df_num_cells.to_csv(tsv_path, sep="\t", index=False)

    total_cells = df_num_cells["numCells"].sum()

    fig, ax = plt.subplots()
    sns.histplot(df_num_cells['numCells'], binwidth=100, ax=ax)
    ax.set_xlabel("Number of Cells past QC")
    ax.set_title(f"{suffix} ({total_cells:,} cells)\n{runtag}")
    fig.savefig(png_path)
    plt.close(fig)

def plot_per_sample_umap_clusters(sample_ids, h5ad_dir, run_tag, suffix, 
                                  output_dir, doublet_rate_key="doublet_rate"):
    """
    Plot UMAP clusters for each sample and save all plots to a single PDF file.

    Args:
        sample_ids (list of str): List of sample IDs (fileID)
        h5ad_dir (str): Directory containing processed h5ad files.
        run_tag (str): QC/processing tag (e.g. round1, round3, etc.).
        suffix (str): Suffix for output file naming.
        output_dir (str): Directory to save the output PDF.
        doublet_rate_key (str): Key in .uns for doublet rate (default: "doublet_rate")
    """
    pdf_path = os.path.join(output_dir, f'ATAC_post_filter_clusters.{suffix}.pdf')
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
            # Shrink plot to make space for the legend
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.70, box.height])

            for i, cluster in enumerate(unique_clusters):
                idx = np.where(np.array(clusters) == cluster)[0]
                ax.scatter(adata.obsm['X_umap'][idx, 0], adata.obsm['X_umap'][idx, 1],
                    c=[colors[i]], label=cluster, s=1)
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper right", title="Leiden Cluster", frameon=False,
                      borderaxespad=0.)
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
                f'{suffix}: {sample_id}\n'
                f'Median TSS enrichment: {median_tss:.2f}\n'
                f'Median # fragments: {median_frags:.0f}\n'
                f'{dbl_line}',
                loc='center', fontsize=10
            )
            plt.subplots_adjust(top=0.7, right=0.7)
            plt.tight_layout(rect=[0, 0, 0.8, 1])
            pdf.savefig(fig)
            plt.close(fig)


def plot_umap_by_sample(adata, suffix, runtag, outdir, sample_colors=None, key="sampleID"):
    """
    Plot UMAP for all cells, colored by sample.
    """
    samples = adata.obs[key]
    samplelist = sorted(set(samples))
    fig, ax = plt.subplots(figsize=(7, 4))
    plt.style.use("default")
    for sample in samplelist:
        idx = np.where(np.array(samples) == sample)[0]
        color = sample_colors[sample] if (sample_colors and sample in sample_colors) else None
        ax.scatter(adata.obsm['X_umap'][idx,0], adata.obsm['X_umap'][idx,1],
                   color=color, label=sample, s=2, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.8, 1), loc="upper right", title=key, frameon=False)
    ax.set_xlabel('UMAP-1', fontsize=12)
    ax.set_ylabel('UMAP-2', fontsize=12)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_title(f"{suffix}\n", fontsize=14)
    plt.subplots_adjust(right=0.7)
    outpng = os.path.join(outdir, f'ATAC_UMAP_bySample.{suffix}.{runtag}.png')
    plt.savefig(outpng, bbox_inches='tight', dpi=800)
    #outpdf = os.path.join(outdir, f'ATAC_UMAP_bySample.{suffix}.{runtag}.pdf')
    #plt.savefig(outpdf, bbox_inches='tight')
    plt.close(fig)

def plot_umap_by_donor(adata, suffix, runtag, outdir, donor_colors=None, key='donorID'):
    """
    Plot UMAP for all cells, colored by donor.
    """
    donors = adata.obs[key]
    donorlist = sorted(set(donors))
    fig, ax = plt.subplots(figsize=(7, 4))
    plt.style.use("default")
    
    for donor in donorlist:
        idx = np.where(np.array(donors) == donor)[0]
        color = donor_colors[donor] if (donor_colors and donor in donor_colors) else None
        ax.scatter(adata.obsm['X_umap'][idx,0], adata.obsm['X_umap'][idx,1],
                   color=color, label=donor, s=2, alpha=0.3) 
    ax.legend(bbox_to_anchor=(1.8, 1), loc="upper right", title=key, frameon=False)
    ax.set_xlabel('UMAP-1', fontsize=12)
    ax.set_ylabel('UMAP-2', fontsize=12)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_title(f"{suffix}\n", fontsize=14)
    plt.subplots_adjust(right=0.7)
    outpng = os.path.join(outdir, f'ATAC_UMAP_byDonor.{suffix}.{runtag}.png')
    plt.savefig(outpng, bbox_inches='tight', dpi=800)
    #outpdf = os.path.join(outdir, f'ATAC_UMAP_byDonor.{suffix}.{runtag}.pdf')
    #plt.savefig(outpdf, bbox_inches='tight')
    plt.close(fig)

    
def plot_umap_single_tissue_sample_by_side(adata, tissue, suffix, runtag, outdir, sample_colors=None, 
                                           tissue_key='tissue', sample_key="sampleID", ncol=4, show=False):
    """
    For a given tissue, plot a multi-panel UMAP: one panel per sample (of that tissue),
    all cells as white background, sample's cells in color.
    """

    # Get samples for this tissue
    sel = adata.obs[tissue_key] == tissue
    samples = adata.obs.loc[sel, sample_key]
    samplelist = sorted(samples.unique())
    nrow = math.ceil(len(samplelist) / ncol)

    plt.style.use("ggplot")
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 5 * nrow))
    axes = axes.flatten()
    

    # Set up colors
    if sample_colors is None and 'sample_colors' in adata.uns:
        sample_colors = adata.uns['sample_colors']
        
    for i, sample in enumerate(samplelist):
        ax = axes[i]
        
        # All cells in white
        
        ax.scatter(
            adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],
            color="white", s=0.8, alpha=0.3, rasterized=True, zorder=1)

        # Highlight cells of this tissue in each sample
        idx = np.where((adata.obs[tissue_key] == tissue) & (adata.obs[sample_key] == sample))[0]
        color = sample_colors[sample] if (sample_colors and sample in sample_colors) else "#1f77b4"
        ax.scatter(
            adata.obsm['X_umap'][idx, 0], adata.obsm['X_umap'][idx, 1],
            color=color, label=sample, s=0.8, alpha=0.8, rasterized=True, zorder=2
        )
        
        ax.legend(
            handles=[plt.Line2D([0], [0], marker='o', color=color, markerfacecolor=color, 
                                markersize=8, label=sample)],
            loc="lower left", frameon=False)
        ax.set_xlabel('UMAP-1', fontsize=10)
        ax.set_ylabel('UMAP-2', fontsize=10)
        ax.set_title(f"{sample}\nN={len(idx)}", fontsize=12)
        
    # Hide unused axes
    for j in range(len(samplelist), len(axes)):
        axes[j].axis('off')
    plt.suptitle(f"{suffix}:{tissue}", fontsize=16)
    plt.tight_layout(rect=(0,0,1,0.97))

    outpng = os.path.join(outdir, f'ATAC_UMAP_sampleBySide.{suffix}.{standardize_tissue_name(tissue)}.{runtag}.png')
    plt.savefig(outpng, bbox_inches='tight', dpi=800)
    if show:
        plt.show()
    plt.close(fig)
    print(f"Wrote: {outpng}")

        
def plot_umap_per_tissue_by_sample_all(adata, suffix, runtag, outdir, sample_colors=None,
                                       tissue_key="tissue", sample_key="sampleID", ncol=4):
    """
    For each tissue in AnnData, call plot_umap_single_tissue_sample_by_side.
    """
    tissues = sorted(adata.obs[tissue_key].unique())
    for tissue in tissues:
        plot_umap_single_tissue_sample_by_side(adata, tissue, suffix, runtag, outdir, 
                                               sample_colors=sample_colors, 
                                               tissue_key=tissue_key, sample_key=sample_key, ncol=ncol)
        
        
def plot_umap_with_QC(adata, suffix, runtag, outdir, sample_colors=None,
                               tss_bounds=[7,10,12,15,20,30,100],
                               frag_bounds=[1000,2000,3000,4000,5000,10000,100000]):
    """
    For each sample in adata.obs['sample'], plot a 1x3 panel UMAP:
    - sample highlight on gray, 
    - TSS enrichment score coloring,
    - n_fragment coloring.
    Save each as a PNG in outdir.
    """
    
    samples = adata.obs['sample']
    samplelist = sorted(set(samples))
    
    for fileID in samplelist:
        print(f"Plotting for sample: {fileID}")
        plt.style.use("ggplot")
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
        
        # Draw all other samples as white/gray background
        for file in samplelist:
            if file != fileID:
                idx = np.where(np.array(samples) == file)[0]
                ax1.scatter(adata.obsm['X_umap'][idx, 0], adata.obsm['X_umap'][idx, 1], 
                            color="white", s=1, alpha=0.3)
                ax2.scatter(adata.obsm['X_umap'][idx, 0], adata.obsm['X_umap'][idx, 1], 
                            color="white", s=1, alpha=0.3)
                ax3.scatter(adata.obsm['X_umap'][idx, 0], adata.obsm['X_umap'][idx, 1], 
                            color="white", s=1, alpha=0.3)
        idx = np.where(np.array(samples) == fileID)[0]
        
        # ax1: main highlight color
        color = sample_colors[fileID] if (sample_colors and fileID in sample_colors) else "#1f77b4"
        ax1.scatter(adata.obsm['X_umap'][idx, 0], adata.obsm['X_umap'][idx, 1],
                    color=color, label=fileID, s=1, alpha=0.7)
        legend_entries = [plt.Line2D([0], [0], marker='o', color=color, markerfacecolor=color, markersize=8, label=fileID)]
        ax1.legend(handles=legend_entries, loc="lower left", frameon=False)
        ax1.set_xlabel('UMAP-1', fontsize=10)
        ax1.set_ylabel('UMAP-2', fontsize=10)
        ax1.set_title(f"{suffix}\n{fileID}\nN={len(idx)}", fontsize=12)
        
        # ax2: TSS enrichment color
        cmap2 = mpl.cm.YlGnBu
        norm2 = mpl.colors.BoundaryNorm(tss_bounds, cmap2.N, extend='both')
        values2 = adata.obs['tsse'].iloc[idx]
        ax2.scatter(adata.obsm['X_umap'][idx, 0], adata.obsm['X_umap'][idx, 1],
                    c=values2, cmap=cmap2, norm=norm2, s=1)
        cbar2 = plt.colorbar(ScalarMappable(norm=norm2, cmap=cmap2), ax=ax2, orientation='vertical', 
                             label='TSS enrichment score')
        ax2.set_xlabel('UMAP-1', fontsize=10)
        ax2.set_ylabel('UMAP-2', fontsize=10)
        ax2.set_title("Color by TSS enrichment score", fontsize=10)
        
        # ax3: n_fragments color
        cmap3 = mpl.cm.YlGnBu
        norm3 = mpl.colors.BoundaryNorm(frag_bounds, cmap3.N, extend='both')
        values3 = adata.obs['n_fragment'].iloc[idx]
        ax3.scatter(adata.obsm['X_umap'][idx, 0], adata.obsm['X_umap'][idx, 1],
                    c=values3, cmap=cmap3, norm=norm3, s=1)
        cbar3 = plt.colorbar(ScalarMappable(norm=norm3, cmap=cmap3), ax=ax3, orientation='vertical', 
                             label='# fragment')
        ax3.set_xlabel('UMAP-1', fontsize=10)
        ax3.set_ylabel('UMAP-2', fontsize=10)
        ax3.set_title("Color by # fragment", fontsize=10)
        
        # Save/show
        outpng = os.path.join(outdir, f'ATAC_UMAP_withQC.{suffix}.{fileID}.{runtag}.png')
        plt.savefig(outpng, bbox_inches='tight', dpi=800)
        plt.close(fig)

def plot_cells_per_tissue_by_donor(adata, suffix, runtag, outdir, tissue_key="tissue", donor_key="donorID"):
    df = adata.obs.groupby([tissue_key, donor_key], observed=True).size().reset_index(name='n_cells')
    df_pivot = df.pivot(index=tissue_key, columns=donor_key, values='n_cells').fillna(0)
    # Order tissues by total cell count
    df_pivot = df_pivot.loc[df_pivot.sum(axis=1).sort_values(ascending=False).index]

    #n_tissue = df_pivot.shape[0]
    #n_donor = df_pivot.shape[1]
    #fig_height = max(3, min(0.5*n_tissue, 20))
    #fig_width  = max(8, min(1.2 + 0.8*n_donor, 20))
    fig_height = 5
    fig_width = 8
    
    # Get color list for donors, matching the donorID column order in the pivot
    donor_ids = df_pivot.columns.tolist()
    donor_colors = adata.uns[donor_key + '_colors']
    if isinstance(donor_colors, dict):
        donor_colors = [donor_colors[did] for did in donor_ids]

    ax = df_pivot.plot(
        kind='barh',
        stacked=True,
        figsize=(fig_width, fig_height),
        color=donor_colors,
        edgecolor='none'
    )
    plt.xlabel('Number of Cells')
    plt.ylabel('')
    plt.title('Number of Cells per Tissue (colored by Donor)')
    plt.legend(
        title="Donor",
        bbox_to_anchor=(-.5, -0.15), #(1.05, 1),
        loc='bottom center', #'upper left',
        labels=donor_ids
    )
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    outpng = os.path.join(outdir, f'ATAC_cellCount_perTissueByDonor.{suffix}.{runtag}.png')
    plt.savefig(outpng, dpi=300)
    #outpdf = os.path.join(outdir, f'ATAC_cellCount_perTissueByDonor.{suffix}.{runtag}.pdf')
    #plt.savefig(outpdf)
    plt.close()
    print(f"Figure saved as:\n  {outpng}")

def plot_cells_per_donor_per_tissue(adata, suffix, runtag, outdir, tissue_key="tissue", donor_key="donorID"):
    df = adata.obs.groupby([tissue_key, donor_key], observed=True).size().reset_index(name='n_cells')
    # Pivot so index is tissue, columns are donor, values are counts
    df_pivot = df.pivot(index=tissue_key, columns=donor_key, values='n_cells').fillna(0)
    df_pivot = df_pivot.loc[df_pivot.sum(axis=1).sort_values(ascending=False).index]

    #n_tissue = df_pivot.shape[0]
    #n_donor = df_pivot.shape[1]
    #fig_height = max(3, min(0.5*n_tissue, 20))   # At least 3, at most 20
    #fig_width  = max(8, min(1.2 + 0.8*n_donor, 20))  # At least 6, at most 20
    fig_height = 5
    fig_width = 8
    
    
    donor_ids = df_pivot.columns.tolist()
    donor_colors = adata.uns[donor_key + '_colors']
    if isinstance(donor_colors, dict):
        donor_colors = [donor_colors[did] for did in donor_ids]

    ax = df_pivot.plot(
        kind='barh', 
        figsize=(fig_width, fig_height),
        color=donor_colors,
        edgecolor='none'
    )
    plt.xlabel('Number of Cells')
    plt.ylabel('')
    plt.title('Number of Cells per Donor in Each Tissue')
    plt.legend(
        title=donor_key,
        bbox_to_anchor=(-.5, -0.15), #(1.05, 1),
        loc='bottom center', #'upper left',
        labels=donor_ids
    )
    
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    outpng = os.path.join(outdir, f'ATAC_cellCount_perDonorPerTissue.{suffix}.{runtag}.png')
    plt.savefig(outpng, dpi=300)
    #outpdf = os.path.join(outdir, f'ATAC_cellCount_perDonorPerTissue.{suffix}.{runtag}.pdf')
    #plt.savefig(outpdf)
    plt.close()
    print(f"Figure saved as:\n  {outpng}")