#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Parallel launcher for remove_pcr_chimeric_reads.py on all samples in a metadata table.
"""

import os
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import argparse

def process_atac_id(args, atac_id):
    inp = os.path.join(args.fragdir, atac_id, "fragments.tsv.gz")
    out = os.path.join(args.fragdir, atac_id, "fragments.rmPCRchimeric.tsv.gz")
    stat = os.path.join(args.fragdir, atac_id, "rmPCRchimeric_stat.txt")
    if os.path.exists(out):
        print(f"[SKIP] {atac_id} (output exists)")
        return
    cmd = ['python', args.script, inp, out, stat]
    try:
        subprocess.run(cmd, check=True)
        print(f"[OK]   {atac_id}")
    except subprocess.CalledProcessError as e:
        print(f"[FAIL] {atac_id}: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Parallel launcher for remove_pcr_chimeric_reads.py on all samples in a metadata table."
    )
    parser.add_argument('--fragdir', required=True, help='Fragment files directory')
    parser.add_argument('--script', required=True, help='Path to remove_pcr_chimeric_reads.py script')
    parser.add_argument('--metadata', required=True, help='Metadata file (with atacID in column 2)')
    parser.add_argument('--threads', type=int, default=8, help='Number of parallel workers (default: 8)')
    args = parser.parse_args()

    df = pd.read_csv(args.metadata, sep='\t', header=None)
    atac_ids = df[1].unique()
    print(f"Launching PCR chimera removal for {len(atac_ids)} samples using {args.threads} workers...")
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        for atac_id in atac_ids:
            pool.submit(process_atac_id, args, atac_id)
    print("Done.")

if __name__ == "__main__":
    main()