#!/usr/bin/env python3

"""
Author: Kaili Fan
Concatenate per-sample AnnData, run joint embedding, clustering, and produce summary UMAPs.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import load_config, standardize_tissue_name
from atac_qc.atac_plots import (
    plot_umap_by_sample, plot_umap_by_donor, plot_umap_by_sample_by_side, plot_umap_with_QC, assign_colors
)

def collect_input_adatas(sample_list, h5ad_dir, runtag):
    """Get AnnData file paths and read them into memory."""
    h5ad_paths = [
        os.path.join(h5ad_dir, f"{sid}.final.{runtag}.h5ad") for sid in sample_list
    ]
    input_adatas = []
    sample_ids = []
    for sid, path in zip(sample_list, h5ad_paths):
        if os.path.exists(path):
            input_adatas.append(snap.read(path))
            sample_ids.append(sid)
        else:
            print(f"  [WARN] Skip missing file: {path}")
    return sample_ids, input_adatas

def create_joint_anndataset(sample_ids, adatas_list, tissue2, runtag, workdir):
    """Create a snap AnnDataSet and return the handle."""
    anndataset = snap.AnnDataSet(
        adatas=[(sid, adata) for sid, adata in zip(sample_ids, adatas_list)],
        filename=os.path.join(workdir, f"dGTEx_ATAC.{tissue2}.{runtag}.h5ads"),
    )
    return anndataset

def update_dataset_metrics(adataset, gencode_gtf, n_threads):
    """Copy per-sample metrics and update TSS enrichment."""
    for metric in ['fragment_paired', 'n_fragment', 'frac_dup', 'frac_mito']:
        if metric in adataset.adatas.obsm:
            adataset.obsm[metric] = adataset.adatas.obsm[metric]
        if metric in adataset.adatas.obs:
            adataset.obs[metric] = adataset.adatas.obs[metric]
    snap.metrics.tsse(adataset, gene_anno=gencode_gtf, n_jobs=n_threads)
    # Generate unique cell ids
    adataset.obs_names = [f"{sa}:{bc}" for sa, bc in zip(adataset.obs['sample'], adataset.obs_names)]

def save_merged_anndata(adataset, h5ad_dir, tissue2, runtag):
    """Convert AnnDataSet to AnnData and save as .h5ad."""
    adata_merged = adataset.to_adata()
    adataset.close()
    # Parse donor IDs for plotting
    new_names = [bc.split(":")[-1] for bc in adata_merged.obs_names]
    adata_merged.obs['donorID'] = [bc.split("_")[0] for bc in new_names]
    adata_merged.obs_names = new_names
    merged_h5ad_path = os.path.join(h5ad_dir, f"{tissue2}.{runtag}.h5ad")
    adata_merged.write(merged_h5ad_path, compression="gzip")
    print(f"Wrote merged h5ad: {merged_h5ad_path}")
    return adata_merged

def main(config_path, runtag):
    config = load_config(config_path)
    workdir = config['paths']['workdir']
    h5ad_dir = config['paths']['output_h5ad_dir']
    fig_dir = os.path.join(config['paths']['output_figures_dir'], runtag)
    os.makedirs(fig_dir, exist_ok=True)
    os.chdir(workdir)

    n_threads = config['params'].get('n_threads', 16)
    tissue = config['params']['tissue']
    tissue2 = standardize_tissue_name(tissue)
    n_features_merge = config.get('merge', {}).get('n_features', 50000)
    gencode_gtf = config['references']['gencode_gtf']

    # Load sample info
    df = pd.read_csv(
        config['paths']['sample_metadata'], sep='\t', header=None,
        names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"]
    )
    working_df = df[df["tissue"] == tissue]
    sample_list = working_df["atacID"].unique().tolist()

    # ---- 1. Collect input AnnData objects
    sample_ids, adatas_list = collect_input_adatas(sample_list, h5ad_dir, runtag)
    if not adatas_list:
        print("No valid AnnData files found. Exiting.")
        return

    # ---- 2. Make joint AnnDataSet, update metrics, and embed
    adataset = create_joint_anndataset(sample_ids, adatas_list, tissue2, runtag, workdir)
    update_dataset_metrics(adataset, gencode_gtf, n_threads)
    print(f'Joint AnnDataSet created: cells={adataset.n_obs} samples={len(adatas_list)}')
    snap.pp.select_features(adataset, n_features=n_features_merge)
    snap.tl.spectral(adataset)
    snap.tl.umap(adataset, random_state=0)

    # ---- 3. Save merged AnnData and plot UMAPs
    adata_merged = save_merged_anndata(adataset, h5ad_dir, tissue2, runtag)

    samples = sorted(set(adata_merged.obs['sample']))
    palette_s = "Set1" if len(samples) < 10 else "tab20"
    sample_colors = assign_colors(samples, palette=palette_s)
    print("Sample colors assigned:", sample_colors)

    donors = sorted(set(adata_merged.obs['donorID']))
    config_donor_colors = config['color'].get("donor_color")
    palette_d = "Set1" if len(donors) < 10 else "tab20"
    palette_donor_colors = assign_colors(donors, palette=palette_d)
    if config_donor_colors:
        donor_colors = {d: config_donor_colors.get(d, palette_donor_colors[d]) for d in donors}
    else:
        donor_colors = palette_donor_colors
    print("Donor colors assigned:", donor_colors)

    plot_umap_by_sample(adata_merged, tissue2, runtag, fig_dir, sample_colors=sample_colors)
    plot_umap_by_donor(adata_merged, tissue2, runtag, fig_dir, donor_colors=donor_colors)
    plot_umap_by_sample_by_side(adata_merged, tissue2, runtag, fig_dir, sample_colors=sample_colors, ncol=min(len(samples), 4))
    plot_umap_with_QC(adata_merged, tissue2, runtag, fig_dir, sample_colors=sample_colors)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate AnnData, embed, and plot joint ATAC UMAPs.")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("runtag", help="Run tag (e.g. round3)")
    args = parser.parse_args()
    main(args.config, args.runtag)