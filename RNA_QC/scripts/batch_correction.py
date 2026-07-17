#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Runs Harmony batch correction across donors.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from rna_qc.utils import load_config, standardize_tissue_name
from rna_qc.rna_plots import move_figures_to_newdir

def run_harmony_batch_correction(output_h5ad_dir, workdir, tissue, runtag, key="sampleID"):
    
    tissue_std = standardize_tissue_name(tissue)

    adata_path = os.path.join(output_h5ad_dir, f"{tissue_std}_GEX.filtered.{runtag}.h5ad")
    if not os.path.exists(adata_path):
        print(f"[ERROR] No h5ad for tissue {tissue} at {adata_path}.")
        return
    
    adata = sc.read_h5ad(adata_path)
    print(f"[INFO] Loaded AnnData for tissue {tissue} ({adata.shape[0]} cells, {adata.shape[1]} genes)")
    
    figures_dir = os.path.join(workdir, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    # Save before Harmony embeddings
    adata.obsm['X_umap_before_harmony'] = adata.obsm['X_umap']
    adata.obsm['X_pca_before_harmony'] = adata.obsm['X_pca']
    # Record cluster assignments if present
    if 'leiden' in adata.obs:
        adata.obs['leiden_before_harmony'] = adata.obs['leiden']
    adata_before = adata.copy()

    # Run Harmony batch correction
    print("[INFO] Running Harmony integration...")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        # Direct harmonypy call — sce.pp.harmony_integrate has a shape bug in the
        # current scanpy_external; call harmonypy directly and place the corrected
        # PCA into X_pca_harmony ourselves.
        import harmonypy as hm
        ho = hm.run_harmony(
            adata.obsm['X_pca'],
            adata.obs[[key]],
            key,
            max_iter_harmony=30,
            verbose=True,
        )
        adata.obsm['X_pca_harmony'] = np.asarray(ho.Z_corr).T  # Z_corr is (n_pcs, n_cells) -> (n_cells, n_pcs)

    # Build UMAP from Harmony-corrected PCA
    adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)

    # Plot before and after Harmony UMAP
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ## before Harmony
    sc.pl.umap(adata_before, color=[key], ax=axes[0], show=False, title="Before Harmony", legend_loc=None)
    axes[0].set_title('Before Harmony')
    axes[0].set_aspect('auto')
    ## after Harmony
    sc.pl.umap(adata, color=[key], ax=axes[1], show=False, title="After Harmony", legend_loc=None)
    axes[1].set_title('After Harmony')
    axes[1].set_aspect('auto')
    # Add back in legend
    import matplotlib.patches as mpatches
    
    categories = adata.obs[key].cat.categories
    colors = adata.uns[f'{key}_colors']
    
    # Create handle list for the legend
    handles = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, categories)]
    
    # Place the legend at the bottom center
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False)
    #fig.legend(handles=handles, loc='center right', bbox_to_anchor=(1.02, 0.5), ncol=1, frameon=False)

    # 4. Adjust layout to leave space at the bottom
    fig.suptitle(f"{tissue}: Harmony Batch Correction", fontsize=16)
    # rect=[left, bottom, right, top] -> reserves the bottom 15% of the canvas for the legend
    plt.tight_layout(rect=[0, 0.15, 1, 1])


    plot_path = os.path.join(figures_dir, f"Harmony_batchCorrection_umap.{tissue_std}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[INFO] Saved batch correction UMAP to {plot_path}")

    # Save Harmony-corrected AnnData
    out_adata_path = os.path.join(output_h5ad_dir, f"{tissue_std}_GEX.filtered.RMbatch.{runtag}.h5ad")
    adata.write(out_adata_path)
    print(f"[INFO] Saved Harmony-corrected AnnData: {out_adata_path}")

    # Move figures to new directory
    move_figures_to_newdir(workdir, old="figures", new=f"batch_correction.{runtag}")

    print(f"[INFO] Harmony batch correction done for tissue: {tissue}")

def main(config_path, runtag):
    config = load_config(config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']

    tissue = config['params']['tissue']
    
    if tissue == "---":

        # Load master sample metadata across tissues/donors
        sample_metadata = config['paths']['sample_metadata']
        df = pd.read_csv(sample_metadata, sep='\t', header=None, index_col=False,
                         names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    
        tissue_names = sorted(df['tissue'].unique())
        print(f"[INFO] Running Harmony batch correction for MULTIPLE tissues: {tissue_names}")
        
        failed_tissues = []
        for idx, tissue_name in enumerate(tissue_names, 1):
            print(f"\n============== Processing tissue: {tissue_name} ({idx}/{len(tissue_names)}) ==============")

            try:
                run_harmony_batch_correction(output_h5ad_dir, workdir, tissue_name, runtag, key = "sampleID")
            except Exception as e:
                print(f"[ERROR] Harmony batch correction failed for {tissue_name}: {e}")
                failed_tissues.append(tissue_name)

        if failed_tissues:
            print(f"[ERROR] batch_correction failed for {len(failed_tissues)} tissue(s): {failed_tissues}")
            sys.exit(1)
    else:
        print(f"\n============== Processing tissue: {tissue} ==============")

        try:
            run_harmony_batch_correction(output_h5ad_dir, workdir, tissue, runtag, key = "sampleID")
        except Exception as e:
            print(f"[ERROR] Harmony batch correction failed for {tissue}: {e}")
            sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Harmony batch correction for scRNA-seq h5ad.")
    parser.add_argument("config", help="YAML config file describing tissue, paths, etc.")
    parser.add_argument("runtag", help="Tag for this run (e.g. round3 or v1)")
    args = parser.parse_args()
    main(args.config, args.runtag)
