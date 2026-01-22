#!/bin/bash

# Author: Kaili Fan
# Description: Run the full ATAC QC workflow, with parallel PCR chimera removal,
# followed by all sequential ATAC processing steps.
#
# Usage: bash run_atac_qc_workflow.sh [runtag] [nthread]
#    - default runtag: demo_run
#    - nthread: Number of samples to process in parallel for PCR chimeric removal; default = 4

set -euo pipefail

# ===========================
# PARAMETERS
# ===========================
metatable="/Users/neva/Documents/Kaili-Multiomic-workflow/data/human.txt"
runtag="${1:-v1}"
nthread="${2:-4}"

config="config/atac_qc_config.yaml"
fragdir="/Users/neva/Documents/Kaili-Multiomic-workflow/data/"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# ===========================
# PIPELINE STEPS
# ===========================

printf "\nStep 1: Remove PCR chimeric reads (in parallel across samples)\n"
awk '$2!="---"{print $2}' "$metatable" | xargs -n 1 -P "$nthread" -I{} \
  bash scripts/remove_pcr_chimeric_awk.sh \
    "$fragdir"/{}/fragments.tsv.gz \
    "$fragdir"/{}/fragments.rmPCRchimeric.tsv.gz \
    "$fragdir"/{}/rmPCRchimeric_stat.txt \
    --force

printf "\nStep 2: Convert fragments to h5ad\n"
python scripts/load_fragments.py "$config"

printf "\nStep 3: QC filtering and visualization\n"
python scripts/qc_kde_filter.py "$config" "$runtag"

printf "\nStep 4: Doublet detection\n"
python scripts/processing_doublet_detection.py "$config" "$runtag"

printf "\nStep 5: Doublet filtering, embedding, clustering\n"
python scripts/downstream_processing.py "$config" "$runtag"

printf "\nStep 6: Integrate samples and joint UMAP visualization\n"
python scripts/concatenate_and_cluster.py "$config" "$runtag"

printf "\nAll ATAC QC workflow steps completed. Cheers!\n"
