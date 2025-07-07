#!/usr/bin/env python3

"""
Author: Kaili Fan
Description:
    Process ATAC fragment files and save as h5ad (with QC and figures).
"""

import os
import sys
import time
import argparse
import re
import pandas as pd
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from atac_qc.utils import load_config, standardize_tissue_name, print_elapsed_time

def load_barcode_dicts(barcode_whitelist):
    """Load RNA/ATAC barcode mapping from file."""
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
    fragment_dir = config['paths']['fragment_dir']
    gencode_gtf = config['references']['gencode_gtf']
    n_threads = config['params'].get('n_threads', 16)
    min_fragments = config['params']['min_fragments']
    out_h5ad_dir = config['paths']['output_h5ad_dir']
    out_fig_dir = config['paths']['output_figures_dir']
    fragment_file = os.path.join(fragment_dir, atacID, "fragments.rmPCRchimeric.tsv.gz")
    output_h5ad = os.path.join(out_h5ad_dir, f"{atacID}.raw.h5ad")
    output_fig = os.path.join(out_fig_dir, f"{atacID}.fragment_size_distribution.pdf")

    if os.path.exists(output_h5ad):
        if overwrite:
            os.remove(output_h5ad)
            print(f"Overwriting: removed existing {output_h5ad}")
        else:
            print(f"Skipping {atacID}: output already exists.")
            return

    if not os.path.exists(fragment_file):
        print(f"Warning: Fragment file missing for {atacID}. Skipping.")
        return

    print(f"Processing {atacID}")
    data = snap.pp.import_data(
        fragment_file,
        chrom_sizes = snap.genome.hg38,
        file = output_h5ad,
        min_num_fragments = min_fragments,
        sorted_by_barcode = False,
        n_jobs = n_threads
    )

    # Barcode conversion
    new_bc = []
    for barcode in data.obs_names:
        base = barcode.split('-')[0]
        if base not in barcode_dic_rna:
            print(f"Warning: No matching RNA barcode for ATAC barcode {base}")
        new_bc.append(barcode_dic_rna.get(base, base)+'-1')

    channel_search = re.search(r'-(\d+)$', atacID)
    batch_search = re.search(r'EXP(\d+)', atacID)
    channel_number = channel_search.group(1) if channel_search else 'NA'
    batch_number = batch_search.group(1) if batch_search else 'NA'
    data.obs_names = [
        f"{donorID}_{batch_number}_{channel_number}_{bc}"
        for bc in new_bc
    ]

    # Figure and metrics
    snap.pl.frag_size_distr(
        data, interactive=False,
        out_file=output_fig
    )
    snap.metrics.tsse(data, gene_anno=gencode_gtf, n_jobs = n_threads)
    
    data.close()
    print(f"Saved: {output_h5ad}, {output_fig}")

def main(config_path, overwrite=False):
    config = load_config(config_path)
    workdir = config['paths']['workdir']
    os.chdir(workdir)
    for path in [config['paths']['output_h5ad_dir'], config['paths']['output_figures_dir']]:
        os.makedirs(path, exist_ok=True)

    tissue = config['params']['tissue']
    suffix = config['params']['suffix']

    barcode_dic_rna, _ = load_barcode_dicts(config['references']['barcode_whitelist'])
    
    df = pd.read_csv(config['paths']['sample_metadata'], sep='\t', header=None,
                         names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    
    if tissue == "---":
        working_df = df
    else:
        working_df = df[df["tissue"] == tissue]
        
    sample_list = working_df['atacID'].unique()
    print(f"Samples: {list(sample_list)}")
    tissues = sorted(working_df["tissue"].unique())
    print(f"Working tissue: {', '.join(tissues)}")

    start_time = time.time()
    for i, (_, row) in enumerate(working_df.iterrows(), 1):
        print(f"[{i}/{len(working_df)}] Processing {row['atacID']}...")
        try:
            process_fragments(row, config, barcode_dic_rna, overwrite=overwrite)
        except Exception as e:
            print(f"[ERROR] {row['atacID']}: {e}", file=sys.stderr)
            continue
    end_time = time.time()
    print_elapsed_time(start_time, end_time)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process ATAC fragment files and save as h5ad.")
    parser.add_argument("config", help="YAML config file")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing h5ad outputs")
    args = parser.parse_args()
    main(args.config, overwrite=args.overwrite)