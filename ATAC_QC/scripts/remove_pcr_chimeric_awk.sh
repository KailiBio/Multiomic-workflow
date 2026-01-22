#!/usr/bin/env bash
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

# Streaming deduplication assuming input is sorted by chr, start, end
gunzip -c "$input_file" | awk -F'\t' '
BEGIN {
    OFS="\t"
    prev_key = ""
    max_count = 0
    best_line = ""
    initial_unique = 0
    initial_total = 0
    filtered_unique = 0
    filtered_total = 0
}
/^#/ { next }
{
    initial_unique++
    initial_total += $5
    
    key = $1 "\t" $2 "\t" $3
    
    if (key != prev_key) {
        # New group - output previous best
        if (prev_key != "") {
            print best_line
            filtered_unique++
            filtered_total += max_count
        }
        prev_key = key
        max_count = $5
        best_line = $0
    } else if ($5 > max_count) {
        # Better fragment in same group
        max_count = $5
        best_line = $0
    }
}
END {
    # Output last group
    if (prev_key != "") {
        print best_line
        filtered_unique++
        filtered_total += max_count
    }
    
    # Write stats to stderr for capture
    print initial_unique, filtered_unique, initial_total, filtered_total > "/dev/stderr"
}
' 2> >(read iu fu it ft; echo -e "Initial unique reads:\t$iu\nFiltered unique reads:\t$fu\nInitial total reads:\t$it\nFiltered total reads:\t$ft" > "$stats_file") | gzip > "$output_file"

# Read stats and display
echo "Processed sample:"
echo "  Input:    $input_file"
echo "  Output:   $output_file"
echo "  Stats:    $stats_file"
cat "$stats_file"
