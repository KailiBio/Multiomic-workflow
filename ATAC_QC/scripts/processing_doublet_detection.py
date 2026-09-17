#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Perform doublet detection for snATAC-seq samples, with summary plots.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message="Transforming to str index.")
warnings.filterwarnings("ignore", message="n_jobs value .* overridden to 1 by setting random_state. Use no seed for parallelism.")

import os

# set limits for resource usage
os.environ["OPENBLAS_NUM_THREADS"] = "64"
os.environ["OMP_NUM_THREADS"] = "64"
os.environ["MKL_NUM_THREADS"] = "64"

import sys
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import load_config, standardize_tissue_name, setup_logging, require_keys
from atac_qc.atac_plots import plot_doublet_score_probability

# Compatibility Hack for Python 3.12 + SciPy 1.11+
# This adds the 'nonzero' method back to Pandas Series so SciPy indexing doesn't crash.
if not hasattr(pd.Series, 'nonzero'):
    pd.Series.nonzero = lambda self: self.to_numpy().nonzero()

def main(config_path, runtag):
    setup_logging()
    config = load_config(config_path)
    require_keys(config, [
        "paths.workdir", "paths.output_h5ad_dir", "paths.sample_metadata", "params.tissue",
    ], context=config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    outdir = os.path.join(workdir, f'doublet_detection.{runtag}')
    os.makedirs(outdir, exist_ok=True)
    os.chdir(workdir)

    tissue = config['params']['tissue']
    
    n_threads = config['params'].get('n_threads', 16)

    # Load sample info
    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None, index_col=False,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    
    if tissue == "---":
        working_df = df
    else:
        working_df = df[df["tissue"] == tissue]

    tissues = sorted(working_df["tissue"].unique())
    logging.info(f"Working tissue: {', '.join(tissues)}")
    sample_list = working_df['atacID'].unique()
    logging.info(f"Samples: {list(sample_list)}")

    sample_tissue_dict = dict(zip(working_df['atacID'], working_df['tissue']))

    for i, fileID in enumerate(sample_list, 1):
        logging.info(f"========== Processing {fileID} ({i}/{len(sample_list)}) ==========")

        tissue_std = standardize_tissue_name(sample_tissue_dict[fileID])

        h5ad_path = os.path.join(output_h5ad_dir, f'{fileID}.filtered.{runtag}.h5ad')
        if not os.path.exists(h5ad_path):
            logging.error(f"No h5ad for {fileID} at {h5ad_path}.")
            continue

        logging.info("Loading anndata object...")
        adata = ad.read_h5ad(h5ad_path)

        # 1. tile matrix
        logging.info("Generating Tile matrix...")
        snap.pp.add_tile_matrix(adata, n_jobs=n_threads,
                                bin_size=config['params'].get('genomic_bin_size', 500))

        # 2. feature selection
        logging.info("Feature selection...")
        snap.pp.select_features(adata, n_features=config['process'].get('n_features', 250000))

        # 3. doublet detection
        logging.info("Running Scrublet...")
        snap.pp.scrublet(adata)

        # 4. get midpoint and plot
        logging.info("Calculate doublet midpoint and generate distribution plots...")
        doublet_probability = adata.obs['doublet_probability']
        probability_midpoint = ((doublet_probability.max() + doublet_probability.min()) / 2).round(1)
        plot_doublet_score_probability(
            adata.obs['doublet_score'], doublet_probability, probability_midpoint,
            tissue_std, fileID, outdir=outdir, runtag=runtag, show=False
        )
        logging.info(f"Default doublet probability cutoff is {probability_midpoint}")

        # 5. save
        h5ad_out_path = os.path.join(output_h5ad_dir, f"{fileID}.processed.{runtag}.h5ad")
        adata.write(h5ad_out_path)
        logging.info(f"wrote {os.path.relpath(h5ad_out_path)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ATAC doublet detection and visualization.")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("runtag", help="Tag for this run (e.g. round3 or v1)")
    args = parser.parse_args()
    try:
        main(args.config, args.runtag)
    except (FileNotFoundError, KeyError, ValueError) as e:
        logging.error(f"{type(e).__name__}: {e}")
        sys.exit(1)
