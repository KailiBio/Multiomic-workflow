#!/usr/bin/env python3

"""
Author: Kaili Fan
"""

import os
import re
import time
import logging
import yaml


def setup_logging(log_path=None, level=logging.INFO):
    """Configure the root logger: always logs to stderr, optionally also to a file."""
    handlers = [logging.StreamHandler()]
    if log_path:
        os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
        handlers.append(logging.FileHandler(log_path))
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )


def require_keys(config, dotted_keys, context="config"):
    """Raise a clear KeyError if any dotted-path key (e.g. 'paths.input_dir') is missing."""
    for dotted_key in dotted_keys:
        node = config
        parts = dotted_key.split(".")
        for i, part in enumerate(parts):
            if not isinstance(node, dict) or part not in node:
                raise KeyError(f"Missing required key '{'.'.join(parts[:i + 1])}' in {context}")
            node = node[part]


def print_elapsed_time(start_time, end_time):
    elapsed_time = end_time - start_time
    hours, rem = divmod(int(elapsed_time), 3600)
    minutes, seconds = divmod(rem, 60)
    logging.info(f"Elapsed time: {hours}h {minutes}m {seconds}s")

def load_config(config_path):
    """Load YAML config file to Python dict.

    Resolves every path under 'paths'/'references' to an absolute path
    (relative to the current working directory at load time). Scripts
    commonly os.chdir() into workdir later on, and since config paths are
    themselves relative to that same workdir, leaving them unresolved would
    make them silently re-interpreted against the new cwd after the chdir.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    for section in ("paths", "references"):
        for key, value in (config.get(section) or {}).items():
            if isinstance(value, str):
                config[section][key] = os.path.abspath(value)
    return config

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
