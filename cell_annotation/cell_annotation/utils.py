#!/usr/bin/env python3

"""
Author: Kaili Fan
"""

import os
import re
import time
import yaml
import shutil


def print_elapsed_time(start_time, end_time):
    elapsed_time = end_time - start_time
    hours, rem = divmod(int(elapsed_time), 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"Elapsed time: {hours}h {minutes}m {seconds}s")

def load_config(config_path):
    """Load YAML config file to Python dict."""
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def standardize_tissue_name(tissue):
    """
    Standardize tissue names:
    - Remove parentheses (but keep words inside)
    - Replace hyphens and whitespace with underscores
    - Collapse multiple underscores into one
    - Strip leading/trailing underscores
    """
    tissue = re.sub(r'[()]', '', tissue)                # Remove ( and )
    tissue = tissue.replace('-', '_')                   # Hyphens to underscores
    tissue = re.sub(r'\s+', '_', tissue)                # Spaces to underscores
    tissue = re.sub(r'_+', '_', tissue)                 # Remove double/multi underscores
    return tissue.strip('_')

def move_figures_to_newdir(output_figures_dir, tissue_std, old, new):
    old_path = os.path.join(output_figures_dir, tissue_std, old)
    new_path = os.path.join(output_figures_dir, tissue_std, new)
    if os.path.exists(new_path):
        shutil.rmtree(new_path)
    if os.path.exists(old_path):
        os.rename(old_path, new_path)