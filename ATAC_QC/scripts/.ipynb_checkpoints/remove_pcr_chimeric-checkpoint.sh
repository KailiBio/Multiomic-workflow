#!/usr/bin/env bash

# Author: Kaili Fan
# Description: Remove PCR chimeric reads from a 10x fragment file, keeping only the most supported fragment per locus/barcode group.

set -euo pipefail

usage() {
    echo "Usage: $0 input_file output_file stats_file [--force]" >&2
    exit 1
}

if [[ $# -lt 3 ]]; then
    usage
fi

input_file="$1"
output_file="$2"
stats_file="$3"
force="${4:-}"

if [[ ! -f "$input_file" ]]; then
    echo "ERROR: Input file $input_file does not exist." >&2
    exit 1
fi

if [[ -f "$output_file" && "$force" != "--force" ]]; then
    echo "Output file $output_file exists, skipping (use --force to overwrite)."
    exit 0
fi

# Temporary files
tmpfile1=$(mktemp)
tmpfile2=$(mktemp)
trap "rm -f '$tmpfile1' '$tmpfile2'" EXIT

# Read and preprocess the input — skip comments
zcat "$input_file" | grep -v '^#' > "$tmpfile1"

# Stats: initial unique reads and total reads
initial_unique_reads=$(wc -l < "$tmpfile1")
initial_total_reads=$(awk -F'\t' '{c += $5} END{print c}' "$tmpfile1")

# Sort by chrom, start, end, and (descending count) to prepare grouping
sort -k1,1 -k2,2n -k3,3n -k5,5nr "$tmpfile1" > "$tmpfile2"

# Now keep the FIRST occurrence per (chrom, start, end)
# This keeps the line with the highest count (because of sort above)
awk -F'\t' '
    {
        key = $1 "\t" $2 "\t" $3
        if (!(key in seen)) {
            seen[key] = 1
            print $0
        }
    }
' "$tmpfile2" > "$tmpfile1"

filtered_unique_reads=$(wc -l < "$tmpfile1")
filtered_total_reads=$(awk -F'\t' '{c += $5} END{print c}' "$tmpfile1")

# Write gzipped output
cat "$tmpfile1" | gzip > "$output_file"

# Write stats
{
    echo -e "Initial unique reads:\t$initial_unique_reads"
    echo -e "Filtered unique reads:\t$filtered_unique_reads"
    echo -e "Initial total reads:\t$initial_total_reads"
    echo -e "Filtered total reads:\t$filtered_total_reads"
} > "$stats_file"

# Display summary
echo "Processed sample:"
echo "  Input:    $input_file"
echo "  Output:   $output_file"
echo "  Stats:    $stats_file"
echo "  Initial unique reads: $initial_unique_reads"
echo "  Filtered unique reads: $filtered_unique_reads"
echo "  Initial total reads: $initial_total_reads"
echo "  Filtered total reads: $filtered_total_reads"