#!/usr/bin/env python3

"""
Author: Kaili Fan
Description:
    Process ATAC fragment files and save as h5ad (with QC and figures).
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import time
import argparse
import re
import traceback
import pandas as pd
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import load_config, standardize_tissue_name
from atac_qc.atac_plots import assign_donor_colors

def load_barcode_dicts(barcode_whitelist):
    """Load RNA/ATAC barcode mapping from file."""
    print("[INFO] Loading paired barcodes...")
    barcode_dic_rna = {}
    barcode_dic_atac = {}
    with open(barcode_whitelist) as f:
        for line in f:
            rna_barcode, atac_barcode = line.strip().split("\t")
            barcode_dic_rna[atac_barcode] = rna_barcode
            barcode_dic_atac[rna_barcode] = atac_barcode
    return barcode_dic_rna, barcode_dic_atac

def process_fragments(row, config, barcode_dic_rna, overwrite=False):
    """Convert barcodes, generate QC, save h5ad for one sample."""
    
    atacID = row['atacID']
    donorID = row['donorID']
    tissue = row['tissue']
    tissue_std = standardize_tissue_name(tissue)
    
    workdir = config['paths']['workdir']
    fragment_dir = config['paths']['fragment_dir']
    gencode_gtf = config['references']['gencode_gtf']
    n_threads = config['params'].get('n_threads', 16)
    min_fragments = config['params']['min_fragments']
    out_h5ad_dir = config['paths']['output_h5ad_dir']

    outdir = os.path.join(workdir, "fragment")
    os.makedirs(outdir, exist_ok=True)
    
    fragment_file = os.path.join(fragment_dir, atacID, "fragments.rmPCRchimeric.tsv.gz")
    output_h5ad = os.path.join(out_h5ad_dir, f"{atacID}.raw.h5ad")
    output_fig = os.path.join(outdir, f"{atacID}.fragment_size_distribution.pdf")

    if os.path.exists(output_h5ad):
        if overwrite:
            os.remove(output_h5ad)
            print(f"Overwriting: removed existing {output_h5ad}")
        else:
            print(f"Skipping {atacID}: output already exists.")
            return

    # Load chrom sizes: use file from config if provided, otherwise fall back to hg38
    chrom_sizes_file = config['references'].get('chrom_sizes', None)
    if chrom_sizes_file and os.path.isfile(chrom_sizes_file):
        print(f"[INFO] Loading chrom sizes from {chrom_sizes_file}")
        chrom_sizes = {}
        with open(chrom_sizes_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])
    else:
        if chrom_sizes_file:
            print(f"[WARNING] Configured chrom_sizes file not found: {chrom_sizes_file}")
        print("[INFO] Using default hg38 chrom sizes")
        chrom_sizes = snap.genome.hg38

    data = snap.pp.import_data(
        fragment_file,
        chrom_sizes = chrom_sizes,
        file = output_h5ad,
        min_num_fragments = min_fragments,
        sorted_by_barcode = False,
        n_jobs = n_threads
    )

    # Barcode conversion
    data.obs['ATAC_cellbarcode'] = data.obs_names
    print("\n[INFO] Converting barcodes...")
    
    new_bc = []
    for barcode in data.obs_names:
        base = barcode.split('-')[0]
        if base not in barcode_dic_rna:
            print(f"[ERROR] No matching RNA barcode for ATAC barcode {base}")
        new_bc.append(barcode_dic_rna.get(base, base)+'-1')
    data.obs['cellbarcode'] = new_bc

    channel_search = re.search(r'-(\d+)$', atacID)
    batch_search = re.search(r'EXP(\d+)', atacID)
    channel_number = channel_search.group(1) if channel_search else 'unknown'
    batch_number = batch_search.group(1) if batch_search else 'unknown'
    data.obs_names = [f"{donorID}_{batch_number}_{channel_number}_{bc}" for bc in new_bc ]

    # Figure and metrics
    print("[INFO] Plotting fragment size...")
    try:
        snap.pl.frag_size_distr(
            data, interactive=False,
            out_file=output_fig
        )
    except Exception as e:
        print(f"[WARNING] Fragment size plot failed (non-fatal): {e}")

    print("[INFO] Calculating TSS enrichment score...")
    snap.metrics.tsse(data, gene_anno=gencode_gtf, n_jobs = n_threads)

    # Verify tsse was computed
    if 'tsse' not in data.obs:
        data.close()
        os.remove(output_h5ad)
        raise RuntimeError(
            f"snap.metrics.tsse() did not produce 'tsse' column for {atacID}. "
            f"Available obs columns: {list(data.obs.keys())}"
        )

    print(f"[INFO] Adding sampleID: {atacID} to anndata object")
    data.obs['sampleID'] = [str(atacID) for bc in data.obs_names]

    print(f"[INFO] Adding tissue: {tissue_std} to anndata object")
    data.obs['tissue'] = [tissue_std for bc in data.obs_names]

    data.close()
    print(f"Saved raw .h5ad to {output_h5ad}")

def main(config_path, overwrite=False):
    config = load_config(config_path)
    
    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    os.makedirs(output_h5ad_dir, exist_ok=True)
    os.chdir(workdir)

    tissue = config['params']['tissue']
    suffix = config['params']['suffix']

    barcode_dic_rna, _ = load_barcode_dicts(config['references']['barcode_whitelist'])

    # Load sample info
    sample_metadata = config['paths']['sample_metadata']
    df = pd.read_csv(sample_metadata, sep='\t', header=None, index_col=False,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    
    if tissue == "---":
        working_df = df
    else:
        working_df = df[df["tissue"] == tissue]

    tissues = sorted(working_df["tissue"].unique())
    print(f"Working tissue: {', '.join(tissues)}")
    sample_list = working_df['atacID'].unique()
    print(f"Samples: {list(sample_list)}")

    for i, (_, row) in enumerate(working_df.iterrows(), 1):
        print(f"\n========== Processing {row['atacID']} ({i}/{len(working_df)}) ==========")

        try:
            process_fragments(row, config, barcode_dic_rna, overwrite=overwrite)
        except Exception as e:
            print(f"[ERROR] Encountered error for {row['atacID']}: {e}")
            traceback.print_exc()
            # Clean up partial h5ad so downstream steps don't find incomplete files
            partial_h5ad = os.path.join(output_h5ad_dir, f"{row['atacID']}.raw.h5ad")
            if os.path.exists(partial_h5ad):
                os.remove(partial_h5ad)
                print(f"[INFO] Removed partial h5ad: {partial_h5ad}")
            continue

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process ATAC fragment files and save as h5ad.")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing h5ad outputs")
    args = parser.parse_args()
    main(args.config, overwrite=args.overwrite)