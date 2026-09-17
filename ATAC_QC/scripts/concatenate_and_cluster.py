#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Concatenate per-sample AnnData, run joint embedding, clustering, and produce summary UMAPs.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="Transforming to str index.")
warnings.filterwarnings("ignore", message="n_jobs value .* overridden to 1 by setting random_state. Use no seed for parallelism.")
warnings.filterwarnings("ignore", message=r"^`_import_from_c` is deprecated; use `_import_arrow_from_c` instead\.")


import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import load_config, standardize_tissue_name, setup_logging, require_keys
from atac_qc.atac_plots import assign_donor_colors, assign_colors, plot_umap_by_sample, plot_umap_by_donor, plot_umap_single_tissue_sample_by_side, plot_umap_per_tissue_by_sample_all, plot_umap_with_QC, plot_cells_per_tissue_by_donor, plot_cells_per_donor_per_tissue

def collect_input_adatas(sample_list, h5ad_dir, runtag):
    """Get AnnData file paths and read them into memory."""
    h5ad_paths = [os.path.join(h5ad_dir, f"{sid}.final.{runtag}.h5ad") for sid in sample_list]
    input_adatas = []
    sample_ids = []
    for sid, path in zip(sample_list, h5ad_paths):
        if os.path.exists(path):
            input_adatas.append(snap.read(path))
            sample_ids.append(sid)
        else:
            logging.warning(f"Skip missing file: {path}")
    return sample_ids, input_adatas

def create_joint_anndataset(sample_ids, adatas_list, suffix, runtag, h5ad_dir):
    """Create a snap AnnDataSet and return the handle."""
    anndataset = snap.AnnDataSet(
        adatas=[(sid, adata) for sid, adata in zip(sample_ids, adatas_list)],
        filename=os.path.join(h5ad_dir, f"ATAC.{suffix}.{runtag}.h5ads"),
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

def assign_obs_colors(adata, sample_tissue_dict, config,
                      tissue_key='tissue', donor_key='donorID', sample_key='sampleID'):

    samples = sorted(set(adata.obs[sample_key]))
    sample_colors = assign_colors(samples, palette=config["my_color_palette"])
    logging.info(f"Sample colors assigned: {sample_colors}")

    tissue_colors = assign_donor_colors(adata.obs, config['color'].get("tissue_colors"),
                                        key = tissue_key)
    logging.info(f"tissue colors assigned: {tissue_colors}")

    sample_tissue_colors = {s: tissue_colors.get(sample_tissue_dict[s], tissue_colors[sample_tissue_dict[s]])
                                for s in samples}

    donor_colors = assign_donor_colors(adata.obs, config['color'].get("donor_colors"),
                                        key = donor_key)
    logging.info(f"Donor colors assigned: {donor_colors}")

    # Store colors in adata.uns
    adata.obs[sample_key] = adata.obs[sample_key].astype('category')
    adata.uns[sample_key+'_colors'] = np.array([sample_colors[d] for d in adata.obs[sample_key].cat.categories])
    #
    adata.obs[tissue_key] = adata.obs[tissue_key].astype('category')
    adata.uns[tissue_key+'_colors'] = np.array([tissue_colors[d] for d in adata.obs[tissue_key].cat.categories])
    #
    adata.obs[donor_key] = adata.obs[donor_key].astype('category')
    adata.uns[donor_key+'_colors'] = np.array([donor_colors[d] for d in adata.obs[donor_key].cat.categories])

    return {
        'sample_colors': sample_colors,
        'sample_tissue_colors': sample_tissue_colors,
        'tissue_colors': tissue_colors,
        'donor_colors': donor_colors
    }
    
def save_merged_anndata(adataset, h5ad_dir, suffix, runtag, sample_tissue_dict,
                        config, add_colors=True):
    """Convert AnnDataSet to AnnData and save as .h5ad."""

    adata_merged = adataset.to_adata()
    adataset.close()

    # Parse donor IDs for plotting
    new_names = [bc.split(":")[-1] for bc in adata_merged.obs_names]
    adata_merged.obs['donorID'] = [bc.split("_")[0] for bc in new_names]
    adata_merged.obs_names = new_names
    adata_merged.obs['sampleID'] = adata_merged.obs['sample']
    adata_merged.obs['tissue'] = adata_merged.obs['sampleID'].map(sample_tissue_dict)

    # Assign colors and save in adata.uns
    logging.info("Assigning colors...")
    if add_colors:
        all_colors = assign_obs_colors(adata_merged, sample_tissue_dict=sample_tissue_dict,
                                       config=config,
                                       tissue_key='tissue', donor_key='donorID', sample_key='sampleID')
    logging.info("Finished assigning colors...")

    merged_h5ad_path = os.path.join(h5ad_dir, f"{suffix}.ATAC.{runtag}.h5ad")
    adata_merged.write(merged_h5ad_path)
    logging.info(f"Wrote merged h5ad: {merged_h5ad_path}")
    return adata_merged, all_colors

def run_per_tissue(tissue, working_df, output_h5ad_dir, config, runtag, outdir, gencode_gtf, n_threads, n_features_merge):
    """
    Run joint AnnData analysis workflow for a single tissue.
    """
    tissue_std = standardize_tissue_name(tissue)

    sample_list = working_df['atacID'].unique()
    logging.info(f"tissue {tissue} with samples: {list(sample_list)}")
    sample_tissue_dict = dict(zip(working_df['atacID'], working_df['tissue']))

    # 1. Collect AnnData objects
    logging.info("Collecting anndata objects...")
    sample_ids, adatas_list = collect_input_adatas(sample_list, output_h5ad_dir, runtag)
    if not adatas_list:
        logging.error(f"No valid AnnData files found for tissue {tissue}. Exiting tissue run.")
        return

    # 2. Make joint AnnDataSet, update metrics, and embed
    logging.info(f"Joint all anndata objects for tissue {tissue_std}...")
    adataset = create_joint_anndataset(sample_ids, adatas_list, tissue_std, runtag, output_h5ad_dir)
    update_dataset_metrics(adataset, gencode_gtf, n_threads)
    logging.info(f'Joint AnnDataSet created for tissue {tissue_std}: cells={adataset.n_obs} samples={len(adatas_list)}')
    logging.info("Feature selection and dimentional reduction...")
    snap.pp.select_features(adataset, n_features=n_features_merge)
    snap.tl.spectral(adataset)
    snap.tl.umap(adataset, random_state=0)

    # 4. Save merged AnnData and plot UMAPs
    logging.info("Save AnnDataSet to regular AnnData obejct...")
    adata_merged, all_colors = save_merged_anndata(adataset, output_h5ad_dir, tissue_std, runtag,
                                                   sample_tissue_dict,config, add_colors=True)

    logging.info("Plotting UMAPs...")
    plot_umap_by_sample(adata_merged, tissue_std, runtag, outdir, sample_colors=all_colors['sample_colors'])
    plot_umap_by_donor(adata_merged, tissue_std, runtag, outdir, donor_colors=all_colors['donor_colors'])
    plot_umap_per_tissue_by_sample_all(adata_merged, tissue_std, runtag, outdir,
                                       sample_colors=all_colors['sample_tissue_colors'])
    plot_umap_with_QC(adata_merged, tissue_std, runtag, outdir, sample_colors=all_colors['sample_tissue_colors'])

    logging.info("Plotting stat figures...")
    plot_cells_per_tissue_by_donor(adata_merged, tissue_std, runtag, outdir, "tissue", "sampleID")
    plot_cells_per_donor_per_tissue(adata_merged, tissue_std, runtag, outdir, "tissue", "sampleID")

    logging.info(f"Finished mergeing for tissue: {tissue}")

def main(config_path, runtag):
    setup_logging()
    config = load_config(config_path)
    require_keys(config, [
        "paths.workdir", "paths.output_h5ad_dir", "paths.sample_metadata",
        "params.tissue", "references.gencode_gtf",
    ], context=config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    outdir = os.path.join(workdir, f'postprocessing_check.{runtag}')
    os.makedirs(outdir, exist_ok=True)
    os.chdir(workdir)
    
    tissue = config['params']['tissue']

    n_threads = config['params'].get('n_threads', 16)
    
    n_features_merge = config.get('merge', {}).get('n_features', 50000)
    gencode_gtf = config['references']['gencode_gtf']

    # Load sample info
    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None, index_col=False,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        
    if tissue == "---":
        sample_metadata = config['paths']['sample_metadata']
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
                         names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        
        tissues = sorted(df["tissue"].unique())
        logging.info(f"Concatenating & Clustering for MULTIPLE tissues: {tissues}")

        failed_tissues = []
        for idx, tissue_name in enumerate(tissues, 1):
            logging.info(f"========== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")

            try:
                working_df = df[df["tissue"] == tissue_name]
                run_per_tissue(tissue_name, working_df, output_h5ad_dir, config, runtag, outdir,
                               gencode_gtf, n_threads, n_features_merge)
            except Exception as e:
                logging.error(f"Concatenating & Clustering failed for {tissue_name}: {e}")
                failed_tissues.append(tissue_name)

        if failed_tissues:
            logging.error(f"concatenate_and_cluster failed for {len(failed_tissues)} tissue(s): {failed_tissues}")
            sys.exit(1)
    else:
        logging.info(f"========== Processing tissue: {tissue} ==========")

        try:
            working_df = df[df["tissue"] == tissue]
            run_per_tissue(tissue, working_df, output_h5ad_dir, config, runtag, outdir,
                           gencode_gtf, n_threads, n_features_merge)
        except Exception as e:
            logging.error(f"Concatenating & Clustering failed for {tissue}: {e}")
            sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate AnnData, embed, and plot joint ATAC UMAPs.")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("runtag", help="Run tag (e.g. round3)")
    args = parser.parse_args()
    try:
        main(args.config, args.runtag)
    except (FileNotFoundError, KeyError, ValueError) as e:
        logging.error(f"{type(e).__name__}: {e}")
        sys.exit(1)