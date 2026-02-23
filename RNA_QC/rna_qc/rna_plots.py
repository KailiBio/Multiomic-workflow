#!/usr/bin/env python3

"""
Author: Kaili Fan
ATAC_QC plotting functions.
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import upsetplot
import shutil
import warnings
import scanpy as sc
from scipy.stats import median_abs_deviation
import matplotlib.patches as mpatches

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
    
def move_figures_to_newdir(output_figures_dir, old, new):
    old_path = os.path.join(output_figures_dir, old)
    new_path = os.path.join(output_figures_dir, new)
    if os.path.exists(new_path):
        shutil.rmtree(new_path)
    if os.path.exists(old_path):
        os.rename(old_path, new_path)

def plot_qc_violin(adata, metric, tissue, tissue_std, all_colors, metrics_with_cutoffs, figdir, key, nmads, add_mad_lines=True):
    if metric in adata.obs:
        print(f"[INFO] Plotting violin for: {metric}...")
        all_data = pd.to_numeric(adata.obs[metric], errors='coerce').copy()
        df = adata.obs[[metric, key]].copy()
        df[metric] = pd.to_numeric(df[metric], errors='coerce')
        donor_order = list(all_colors.keys())
        fig, axes = plt.subplots(1, 2, figsize=(2 * len(donor_order), 4), gridspec_kw={'width_ratios': [1, 3]}, sharey=True)

        # --- compute MAD-based thresholds once (drop NaN/None) ---
        clean_data = all_data.dropna()
        med = np.median(clean_data)
        mad = median_abs_deviation(clean_data)
        lower = med - nmads * mad
        upper = med + nmads * mad

        # catch warnings specifically for the plotting block to silence divide by zero/overflow
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            # All data
            sns.violinplot(y=all_data, ax=axes[0], color="gray", inner='box')
            axes[0].set_title(tissue_std)
            axes[0].set_xlabel('')
            axes[0].set_xticks([])
            axes[0].set_xticklabels([])
            axes[0].set_ylabel(metric)
            axes[0].grid(False)
            
            # By donor
            sns.violinplot(x=key, y=metric, data=df, ax=axes[1],
                           palette=all_colors, hue=key, legend=False,
                           order=donor_order, inner='box')

            axes[1].set_title(f'by {key}')
            axes[1].set_xlabel('')
            axes[1].set_xticks([])
            axes[1].set_xticklabels([])
            axes[1].tick_params(axis='x', rotation=0)
            axes[1].grid(False)

        # --- custom legend for sampleID colors ---
        legend_handles = [
            mpatches.Patch(color=color, label=sample)
            for sample, color in all_colors.items()
        ]
        axes[1].legend(
            handles=legend_handles,
            title=key,                 
            bbox_to_anchor=(1.05, 1),
            loc='upper left',
            borderaxespad=0.
        )
        
        # Add cutoffs
        if metric in metrics_with_cutoffs:
            for v in metrics_with_cutoffs[metric]:
                if v is not None and np.isfinite(v):
                    axes[0].axhline(y=v, color='red', linestyle='--', linewidth=1)
                    axes[1].axhline(y=v, color='red', linestyle='--', linewidth=1)

        # --- add MAD-based lower/upper lines ---
        if add_mad_lines:
            for ax in axes:
                ax.axhline(y=lower, color='blue', linestyle='--', linewidth=1, label='MAD lower')
                ax.axhline(y=upper, color='blue', linestyle='--', linewidth=1, label='MAD upper')
        
        fig.suptitle(f'{metric} on {tissue}', fontsize=14)
        plt.tight_layout()
        fig.savefig(os.path.join(figdir, f'RNA_QC_violin.by{key}.{tissue_std}.{metric}.png'), dpi=300, bbox_inches='tight')
        plt.close(fig)

def plot_qc_jointplot(adata, x_metric, y_metric, title, save_path):
    if (x_metric in adata.obs) and (y_metric in adata.obs):
        print(f"[INFO] Plotting jointplot for: {x_metric} - {y_metric}...")
        g = sns.jointplot(x=adata.obs[x_metric], y=adata.obs[y_metric], kind="hex", data=adata.obs)
        g.fig.set_size_inches(6, 5)
        g.ax_joint.set_title(title)
        g.ax_joint.set_xlabel(x_metric)
        g.ax_joint.set_ylabel(y_metric)
        g.ax_joint.grid(False)
        plt.tight_layout()
        g.figure.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(g.fig)

def plot_qc_cumulative_distribution(adata, metrics, tissue, tissue_std, all_colors, figdir, key):
    for metric in metrics:
        if metric in adata.obs:
            print(f"[INFO] Plotting cumulative distribution for: {metric}...")
            plt.figure(figsize=(5, 4))
            sns.ecdfplot(data=adata.obs, x=metric, hue=key, palette=all_colors)
            ax = plt.gca()
            ax.set_xlabel(metric)
            ax.set_ylabel('Proportion of Cells')
            ax.set_title(f'Cumulative Distribution in {tissue}')
            ax.grid(False)
            plt.tight_layout()
            plt.savefig(os.path.join(figdir, f"RNA_QC_cumulative.by{key}.{tissue_std}.{metric}.png"), dpi=300, bbox_inches='tight')
            plt.close()

def plot_doublet_hist(adata, tissue, tissue_std, figdir, key, probability_cutoff=None, stage=None):
    for ID in adata.obs[key].unique():
        print(f"[INFO] Plotting doublet score & probability distribution for: {ID}...")
        doublet_scores = adata[adata.obs[key] == ID].obs['doublet_score']
        doublet_probabilities = adata[adata.obs[key] == ID].obs['doublet_probabilities']

        if doublet_probabilities.isnull().all():
            print(f"[WARNING] No doublet probabilities for {ID}")
            continue
        # check whether have given doublet cutoffs. if no, use midpoint
        if probability_cutoff is None:
            print(f"[INFO] using midpoint as doublet cutoff...")
            cutoff = np.round((doublet_probabilities.max()+doublet_probabilities.min())/2,1)
        else:
            print(f"[INFO] using given doublet cutoff...")
            cutoff = float(probability_cutoff)
        print(f"[INFO] current doublet cutoff is {cutoff}.")
        doublet_mask = doublet_probabilities > cutoff

        bins = np.histogram_bin_edges(doublet_scores, bins=50)
        hist_non_doublets, _ = np.histogram(doublet_scores[~doublet_mask], bins=bins)
        hist_doublets, _ = np.histogram(doublet_scores[doublet_mask], bins=bins)
        fig, axes = plt.subplots(1, 2, figsize=(9, 4))
        axes[0].bar(bins[:-1], hist_non_doublets, width=np.diff(bins), edgecolor='black', color='grey', label='Non-doublets', alpha=0.7)
        axes[0].bar(bins[:-1], hist_doublets, width=np.diff(bins), edgecolor='black', color='red', bottom=hist_non_doublets, label='Doublets', alpha=0.7)
        axes[0].set_title(f'Doublets N={str(doublet_mask.sum())}')
        axes[0].set_xlabel('Doublet Score')
        axes[0].set_ylabel('Frequency')
        axes[0].legend()
        axes[0].grid(False)
        prob_bins = np.histogram_bin_edges(doublet_probabilities, bins=30)
        prob_hist, _ = np.histogram(doublet_probabilities, bins=prob_bins)
        axes[1].bar(prob_bins[:-1], prob_hist, width=np.diff(prob_bins), edgecolor='black', color='grey', alpha=0.7)
        axes[1].axvline(cutoff, color='red', linestyle='--', label='Doublet Threshold')
        axes[1].set_title(f'cutoff = {str(cutoff)}')
        axes[1].set_xlabel('Doublet Probability')
        axes[1].set_ylabel('Frequency')
        axes[1].set_xlim(0, 1)
        axes[1].grid(False)
        fig.suptitle(f'distribution of doublet scores and probability on {tissue}, {ID}', fontsize=14)
        plt.tight_layout()
        stage_part = f".{stage}" if stage else ""
        plt.savefig(os.path.join(figdir, f"RNA_QC_doublet_hist{stage_part}.{tissue_std}-{ID}.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

def clustering_umap(adata, tissue, tissue_std, figdir, key):
    print("[INFO] Clustering and UMAP...")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
        # sc.pp.neighbors(adata) # why is this called twice?

        sc.tl.leiden(adata, flavor="igraph")
        sc.tl.umap(adata)
    sc.pl.umap(
        adata, color=["leiden"], title=f'{tissue}: before filtering',
        save=f'.{tissue_std}.LeidenClusterBeforeFiltering.png', show=False)
    
    plotlist = ["leiden", "log10_total_counts", "log10_n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo",
                "pct_exon_reads", "log10_MALAT1_CPM", "doublet_score", "doublet_probabilities"]
    if "pct_exon_reads" not in adata.obs:
        plotlist.remove("pct_exon_reads")
    if "log10_MALAT1_CPM" not in adata.obs:
        plotlist.remove("log10_MALAT1_CPM")
        
    # All together
    sc.pl.umap(
        adata, color=plotlist, wspace=0.3, ncols=3,
        title=[f"{tissue} (before filtering): {feature}" for feature in plotlist],
        save=f'.QCmetrics_LeidenClusterBeforeFiltering.{tissue_std}.png', show=False)
    
    # By sample/donor
    IDlist = adata.obs[key].unique()
    for ID in IDlist:
        adata_plot = adata[adata.obs[key] == ID,:]
        sc.pl.umap(
            adata_plot, color=plotlist, wspace=0.3, ncols=3,
            title=[f"{tissue} - {ID} (before filtering): {feature}" for feature in plotlist],
            save=f'.QCmetrics_LeidenClusterBeforeFiltering.{tissue_std}-{ID}.png', show=False
        )

def plot_upset(upset_data_summary, tissue, key, ID, tissue_std, adata, savepath):
    """
    Plot an upset diagram for the given donor's data summary and save to file.
    """
    plt.figure(figsize=(8, 4))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        upsetplot.plot(upset_data_summary, show_counts=True, sort_by="cardinality")
    plt.suptitle(f'{tissue} - {ID}\ntotal N={len(adata[adata.obs[key] == ID, :].obs_names)}', fontsize=14)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
    plt.close()

def run_umap_clustering(adata, tissue, tissue_std, figdir, plotlist):
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor="igraph")
    
    sc.pl.umap(adata, color=["leiden"], title = f'{tissue}: leiden', save=f'.filterqc.LeidenCluster.{tissue_std}.png', show=False)

    sc.pl.umap(
        adata,
        color=["leiden", "sampleID"],
        wspace=0.5,
        title = [f"{tissue}: {feature}" for feature in plotlist],
        save=f'.filterqc.LeidenCluster-donorID.{tissue_std}.png',
        show=False
    )

def plot_umap_by_ID(adata, tissue, tissue_std, figdir, all_colors, key):
    IDs = adata.obs[key].unique()
    fig, axes = plt.subplots(1, len(IDs), figsize=(3*len(IDs)+2, 3), squeeze=False)
    for i, ID in enumerate(IDs):
        adata_plot = adata[adata.obs[key] == ID, :]
        sc.pl.umap(adata_plot, color=['leiden'], ax=axes[0, i], show=False, save=False)
        axes[0, i].set_title(f"{ID}\nN={len(adata_plot.obs_names)}")
    fig.suptitle(f"{tissue}")
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, f"UMAP_LeidenCluster_by{key}.{tissue_std}.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

def plot_cellcount_per_cluster_barplot(adata, tissue_std, figdir, all_colors, key):
    df_ = adata.obs[[key, 'leiden']]
    df_counts = df_.groupby([key, 'leiden'], observed=False).size().unstack(fill_value=0).T
    
    colors = [all_colors[c] for c in df_counts.columns]
    
    fig, ax = plt.subplots(figsize=(9, 4))
    df_counts.plot(kind='bar', stacked=False, ax=ax, color=colors)
    ax.set_xlabel('Leiden Cluster')
    ax.set_ylabel('Number of Cells')
    ax.set_title(f"{tissue_std}")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    ax.legend(title=key, bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, f"Number_cellbarcode_cluster_by{key}.{tissue_std}.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)

def plot_umap_highlight_by_qc_metrics(adata, tissue, tissue_std, figdir, qc_metrics, key):
    sc.pl.umap(
        adata, color=qc_metrics, wspace=0.3, ncols=3,
        title=[f"{tissue}: {feature}" for feature in qc_metrics],
        show=False,
        save=f".LeidenCluster-QCmetrics.{tissue_std}.png"
    )
    
    for ID in adata.obs[key].unique():
        adata_plot = adata[adata.obs[key] == ID, :]
        sc.pl.umap(
            adata_plot, color=qc_metrics, wspace=0.3, ncols=3,
            title=[f"{tissue} - {ID}\n: {feature}" for feature in qc_metrics],
            show=False,
            save=f".LeidenCluster-QCmetrics.{tissue_std}-{ID}.png"
        )

def plot_qc_metrics_violin_by_cluster(adata, tissue, tissue_std, figdir, features, key, 
                                      nmads, add_mad_lines=True):
    dfqc = adata.obs[[key, "leiden"] + features]

    # === All clusters together ===
    for feature in features:
        fig, ax = plt.subplots(figsize=(12, 4))

        # Compute MAD thresholds from all cells for this feature
        values = dfqc[feature].dropna()
        med = np.median(values)
        mad = median_abs_deviation(values)
        lower = med - nmads * mad
        upper = med + nmads * mad
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sns.violinplot(
               x="leiden", y=feature, data=dfqc, palette=adata.uns.get('leiden_colors', None),
               hue="leiden", legend=False, inner="box", ax=ax
            )
        ax.set_title(f"{tissue}: {feature}")
        ax.set_xlabel("leiden cluster")
        ax.tick_params(axis="x", rotation=0)
        ax.grid(False)

        # Add MAD lines
        if add_mad_lines:
            ax.axhline(y=lower, color="blue", linestyle="--", linewidth=1)
            ax.axhline(y=upper, color="blue", linestyle="--", linewidth=1)
            
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, f"QCmetric_bycluster.{tissue_std}.{feature}.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)
        
    # === Per donor ===
    if len(adata.obs[key].unique()) > 1:
        for feature in features:
            fig, axes = plt.subplots(len(adata.obs[key].unique()), 1, figsize=(12, 3*len(adata.obs[key].unique())), sharey=True)
            y_min, y_max = dfqc[feature].min(), dfqc[feature].max()

            # Compute global MAD thresholds for this feature
            values = dfqc[feature].dropna()
            med = np.median(values)
            mad = median_abs_deviation(values)
            lower = med - nmads * mad
            upper = med + nmads * mad
            
            for i, ID in enumerate(adata.obs[key].unique()):
                df_plot = dfqc[dfqc[key] == ID]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    sns.violinplot(
                    x="leiden", y=feature, data=df_plot,
                    palette=adata.uns.get('leiden_colors', None), 
                    hue="leiden", legend=False, inner="box", 
                    ax=axes[i])
                axes[i].set_ylim(y_min, y_max)
                axes[i].set_title(f"{tissue} - {ID}: {feature}")
                axes[i].set_xlabel("leiden cluster")
                axes[i].tick_params(axis="x", rotation=0)
                axes[i].grid(False)

                # Add MAD lines
                if add_mad_lines:
                    axes[i].axhline(y=lower, color="blue", linestyle="--", linewidth=1)
                    axes[i].axhline(y=upper, color="blue", linestyle="--", linewidth=1)
                    
            plt.tight_layout()
            plt.savefig(os.path.join(figdir, f"QCmetric_bycluster.{tissue_std}-{key}.{feature}.png"), 
                        dpi=300, bbox_inches='tight')
            plt.close(fig)

