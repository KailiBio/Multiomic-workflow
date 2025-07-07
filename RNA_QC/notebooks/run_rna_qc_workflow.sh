#!/bin/bash

# Author: Kaili Fan
# Description: Run the full RNA QC workflow, with parallel PCR chimera removal,
# followed by all sequential RNA processing steps.
#
# Usage: bash run_rna_qc_workflow.sh [runtag]
#    - default runtag: demo_run

set -euo pipefail

# ===========================
# PARAMETERS
# ===========================
metatable="../ATAC_QC/test_data/test_metatable.txt"
runtag="${1:-demo_run}"

config="config/rna_qc_config.yaml"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

# ===========================
# PIPELINE STEPS
# ===========================

printf "\nStep 0: Check CellRanger QC Summary (Optional)\n"
python scripts/check_cellranger_qc_summary.py "$config"

printf "\nStep 1: Concatenate h5ad files by tissue\n"
python scripts/concat_h5ad.py "$config"

printf "\nStep 2: Pre-QC assessment\n"
python scripts/pre_qc_assessment.py "$config"

printf "\nStep 3: QC filtering\n"
python scripts/filter_qc.py "$config" "$runtag"

printf "\nStep 4: Batch correction (Optional)\n"
python scripts/batch_correction.py "$config" "$runtag"

printf "\nStep 5: Clustering and QC re-assessment\n"
python scripts/clustering_and_qc_reassessment.py "$config" "$runtag"

printf "\nAll RNA QC workflow steps completed. Cheers!\n"