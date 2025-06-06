#!/usr/bin/env python3

"""
Author: Kaili Fan
"""

import os
import re
import time
import yaml


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
    Remove parentheses but keep words inside, replace hyphens and whitespace with underscores
    """
    tissue = re.sub(r'\((.*?)\)', r'_\1', tissue)
    tissue = tissue.replace('-', '_')
    tissue = re.sub(r'\s+', '_', tissue)
    tissue = tissue.replace('__', '_')
    return tissue.strip('_')
