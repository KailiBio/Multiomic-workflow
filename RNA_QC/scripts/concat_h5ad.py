#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Concatenate single-cell h5ad files for a given tissue, integrating sample- and Scrinvex-based information, and return a unified h5ad file.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import os
import sys
import argparse
import logging
import pandas as pd
import scanpy as sc
import anndata as ad
import re

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from rna_qc.utils import load_config, standardize_tissue_name, setup_logging, require_keys
from rna_qc.rna_plots import assign_donor_colors

def load_cellranger_h5(input_dir, sampleID):
    """Load filtered CellRanger h5 (10X Genomics) for a sample."""
    adata_path = os.path.join(input_dir, sampleID, 'filtered_feature_bc_matrix.h5')
    if not os.path.isfile(adata_path):
        raise FileNotFoundError(f"[ERROR] CellRanger file not found: {adata_path}")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Variable names are not unique")
        adata = sc.read_10x_h5(adata_path)

    adata.var_names_make_unique()
    return adata

def add_scrinvex_info(scrinvex_dir, sampleID, adata):
    """Add 'pct_exon_reads' to adata.obs from Scrinvex output if available (prefer pct_exonic if present).
    """
    # 1. Try to read precomputed pct_exonic
    pct_exonic_file = os.path.join(scrinvex_dir, sampleID, f"{sampleID}.pct_exonic.tsv")
    if os.path.isfile(pct_exonic_file):
        try:
            pct_df = pd.read_csv(pct_exonic_file, sep=r'\s+')
            # Ensure barcode alignment
            adata.obs['pct_exon_reads'] = pct_df.set_index('barcode').reindex(adata.obs.index)['pct_exonic']
            logging.info(f"Loaded pct_exonic from {pct_exonic_file} for sample {sampleID}.")
            return adata
        except Exception as e:
            logging.warning(f"Could not parse {pct_exonic_file} for pct_exon_reads: {e}")
            # Fall through to legacy mode

    # 2. Fallback to legacy scrinvex if available
    scrinvex_file = os.path.join(scrinvex_dir, sampleID, f'{sampleID}.scrinvex.tsv')
    if os.path.isfile(scrinvex_file):
        try:
            scrinvex_all = pd.read_csv(scrinvex_file, sep='\t')
            scrinvex_count = scrinvex_all.groupby("barcode").sum(numeric_only=True)
            denom = scrinvex_count[['introns', 'junctions', 'exons']].sum(axis=1).replace(0, pd.NA)
            scrinvex_count['pct_exon_reads'] = (scrinvex_count['exons'] / denom * 100).round(2)
            adata.obs['pct_exon_reads'] = scrinvex_count['pct_exon_reads'].reindex(adata.obs.index)
            logging.info(f"Loaded legacy scrinvex for sample {sampleID}.")
            return adata
        except Exception as e:
            logging.warning(f"Could not parse legacy scrinvex for {sampleID}: {e}")

    # 3. No info found
    logging.warning(f"No pct_exon_reads found for sample {sampleID} - skipping.")
    return adata

def extract_batch_number(sampleID):
    """Extract batch number from sampleID with pattern 'EXP<digits>'; returns 'unknown' if not found."""
    match = re.search(r'EXP(\d+)', str(sampleID))
    return match.group(1) if match else "unknown"

def extract_chanel_number(sampleID):
    """Extract chancel number from sampleID with pattern 'EXP<digits>'; returns 'unknown' if not found."""
    match = re.search(r'-(\d+)$', str(sampleID))
    return match.group(1) if match else "unknown"
    
def reindex_obs_names(adata, donorID, batch_number, chanel_number):
    """Update cell barcodes for global uniqueness as donorID_batch_chancel_cellbarcode."""
    adata.obs['donorID'] = donorID
    adata.obs['cellbarcode'] = adata.obs_names
    adata.obs_names = [f"{donorID}_{batch_number}_{chanel_number}_{bc}" for bc in adata.obs_names]

def run_per_tissue(working_df, tissue, input_dir, scrinvex_dir, output_h5ad_dir, donor_colors, tissue_color):
    """Process all samples for a single tissue and concatenate h5ad files."""
    tissue_std = standardize_tissue_name(tissue)
    anndata_list = []

    # get colors
    all_donor_colors = assign_donor_colors(working_df, donor_colors, key='donorID')

    for _, row in working_df.iterrows():
        donorID = row["donorID"]
        sampleID = row["rnaID"]
        batch_number = extract_batch_number(sampleID)
        chanel_number = extract_chanel_number(sampleID)

        if sampleID == "---":
            logging.info("missing data, skip.")
        else:
            logging.info(f"Processing sample: {sampleID} (donor: {donorID}, batch: {batch_number}, chanel: {chanel_number})...")
            try:
                adata = load_cellranger_h5(input_dir, sampleID)
            except Exception as e:
                logging.error(f"Failed to load 10X data for {sampleID}: {e}")
                continue

            logging.info("Loading Scrinvex file for exon_reads% ...")
            if scrinvex_dir:
                adata = add_scrinvex_info(scrinvex_dir, sampleID, adata)
            else:
                logging.warning("No Scrinvex directory provided. Skipping exon_reads% info.")

            reindex_obs_names(adata, donorID, batch_number, chanel_number)

            adata.obs['sampleID'] = sampleID

            # add color
            logging.info(f"donor color for {donorID} is {all_donor_colors[donorID]}")
            adata.uns['donorID_colors'] = all_donor_colors[donorID]
            logging.info(f"tissue color for {tissue} is {tissue_color}")
            adata.uns['tissue_colors'] = tissue_color

            anndata_list.append(adata)

    if not anndata_list:
        logging.error(f"No AnnData objects for tissue '{tissue_std}'. Skipping concatenation.")
        return

    logging.info(f"Concatenating AnnData objects for tissue {tissue_std}...")
    concatenated_adata = ad.concat(anndata_list, join='inner', label=None, index_unique=None)
    concatenated_adata.obs['tissue'] = tissue_std
    concatenated_adata.var = anndata_list[0].var.copy()

    out_path = os.path.join(output_h5ad_dir, f'{tissue_std}_GEX.raw.h5ad')
    logging.info(f"Saving concatenated AnnData to: {out_path}")
    concatenated_adata.write_h5ad(out_path)

    logging.info("Tissue-level concatenation complete.")

def main(config_path):
    setup_logging()
    config = load_config(config_path)
    require_keys(config, [
        "paths.input_dir", "paths.output_h5ad_dir", "paths.sample_metadata",
        "params.tissue", "color",
    ], context=config_path)

    input_dir = config['paths']['input_dir']
    scrinvex_dir = config['paths'].get('scrinvex_dir', None)
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    os.makedirs(output_h5ad_dir, exist_ok=True)
    
    tissue = config['params']['tissue']

    donor_colors = config['color'].get("donor_colors")
    tissue_colors = config['color'].get("tissue_colors")

    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None,  index_col=False,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])

    failed_tissues = []
    if tissue == "---":  # Multiple tissues mode
        tissues = sorted(df["tissue"].unique())
        logging.info(f"Running analysis for MULTIPLE tissues: {tissues}")

        for idx, tissue_name in enumerate(tissues, 1):
            logging.info(f"========== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")

            try:
                tissue_col = tissue_colors.get(tissue_name, "#bdbdbd") if tissue_colors else "#bdbdbd"
                working_df = df[df["tissue"] == tissue_name]
                run_per_tissue(working_df, tissue_name, input_dir, scrinvex_dir, output_h5ad_dir, donor_colors, tissue_col)
            except Exception as e:
                logging.error(f"Encountered error for tissue {tissue_name}: {e}")
                failed_tissues.append(tissue_name)
    else:
        logging.info(f"========== Processing tissue: {tissue} ==========")

        try:
            tissue_col = tissue_colors.get(tissue, "#bdbdbd") if tissue_colors else "#bdbdbd"
            working_df = df[df["tissue"] == tissue]
            run_per_tissue(working_df, tissue, input_dir, scrinvex_dir, output_h5ad_dir, donor_colors, tissue_col)
        except Exception as e:
            logging.error(f"Encountered error for tissue {tissue}: {e}")
            failed_tissues.append(tissue)

    if failed_tissues:
        logging.error(f"concat_h5ad failed for {len(failed_tissues)} tissue(s): {failed_tissues}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Concatenate h5ad files across donors/samples for one tissue."
    )
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    try:
        main(args.config)
    except (FileNotFoundError, KeyError, ValueError) as e:
        logging.error(f"{type(e).__name__}: {e}")
        sys.exit(1)
