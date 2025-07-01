#!/bin/bash

# Author: Kaili Fan
# Description: Run the full multiome_integration workflow.
#
# Usage: bash run_multiome_integration_workflow.sh 

set -euo pipefail

config="config/joint_multiomic_config.yaml"
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if [[ ! -f "$config" ]]; then
    echo "Error: Config file '$config' not found!"
    exit 1
fi

if ! command -v python &> /dev/null; then
    echo "Error: python not found in \$PATH. Activate your environment!"
    exit 1
fi

echo "====================="
echo " Joint Multiome Integration Workflow"
echo "====================="
date

printf "\nStep 1: Remove remaining doublets and check overlapping cells\n"
python scripts/check_overlapping_cells.py "$config"

printf "\nStep 2: Annotate ATAC cells\n"
python scripts/annotate_ATAC_cells.py "$config"

printf "\nStep 3: Identify cis-Regulatory Elements\n"
python scripts/identify_regulatory_elements.py "$config"

printf "\nStep 4: Data imputation\n"
python scripts/data_imputation.py "$config"

printf "\nStep 5: Final integration\n"
python scripts/data_integration.py "$config"

date
printf "\nAll multiome integration workflow steps completed. Cheers!\n"