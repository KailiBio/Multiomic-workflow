#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Remove PCR chimeric reads from a 10x fragment file, keeping only the most supported fragment per locus/barcode group.

Usage:
    python scripts/remove_pcr_chimeric_reads.py input_file output_file stats_file [--force]
"""

import sys
import os
import pandas as pd
import gzip
import argparse

def remove_pcr_chimeric_reads(input_file, output_file):
    with gzip.open(input_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t', comment='#', header=None)
        df.columns = ['chrom', 'start', 'end', 'barcode', 'count']
        
        initial_fragments = len(df)
        initial_total_reads = df['count'].sum()

        # For each fragment (chrom, start, end, barcode), keep row with highest count
        idx = df.groupby(['chrom', 'start', 'end', 'barcode'])['count'].idxmax()
        filtered_df = df.loc[idx]

        filtered_fragments = len(filtered_df)
        filtered_total_reads = filtered_df['count'].sum()

        filtered_df.to_csv(output_file, sep='\t', index=False, header=False, compression='gzip')

    return initial_fragments, filtered_fragments, initial_total_reads, filtered_total_reads


def main():
    parser = argparse.ArgumentParser(
        description="Remove PCR chimeric reads from a 10x fragment file. Output is a filtered .tsv.gz. "
                    "For large datasets, parallelize externally for best performance."
    )
    parser.add_argument("input_file", help="Input fragments.tsv.gz")
    parser.add_argument("output_file", help="Output fragments.rmPCRchimeric.tsv.gz")
    parser.add_argument("stats_file", help="Output stats summary txt")
    parser.add_argument("--force", action="store_true", help="Overwrite output if exists")

    args = parser.parse_args()

    # Check input file
    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file {args.input_file} does not exist.", file=sys.stderr)
        sys.exit(1)

    # Optionally skip if output exists
    if os.path.exists(args.output_file) and not args.force:
        print(f"Output file {args.output_file} exists, skipping (use --force to overwrite).")
        sys.exit(0)
    try:
        initial_frags, filtered_frags, initial_reads, filtered_reads = remove_pcr_chimeric_reads(
            args.input_file, args.output_file
        )
        # Write stats
        with open(args.stats_file, 'w') as f:
            f.write(f"Initial unique fragments:\t{initial_frags}\n")
            f.write(f"Filtered unique fragments:\t{filtered_frags}\n")
            f.write(f"Initial total reads:\t{initial_reads}\n")
            f.write(f"Filtered total reads:\t{filtered_reads}\n")
        print(
            f"Processed sample:\n"
            f"  Input:    {args.input_file}\n"
            f"  Output:   {args.output_file}\n"
            f"  Stats:    {args.stats_file}\n"
            f"  Initial unique fragments: {initial_frags}\n"
            f"  Filtered unique fragments: {filtered_frags}\n"
            f"  Initial total reads: {initial_reads}\n"
            f"  Filtered total reads: {filtered_reads}\n"
        )
    except Exception as e:
        print(f"ERROR processing {args.input_file}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()