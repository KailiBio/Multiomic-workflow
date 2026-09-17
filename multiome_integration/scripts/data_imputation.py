#!/usr/bin/env python3
"""
Author: Kaili Fan
Description: Joint scRNA/scATAC imputation with MultiVI
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import snapatac2 as snap
import scvi
print("scvi-tools version:", scvi.__version__)
import time
scvi.settings.seed = 0
import torch
from scipy import sparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name, print_elapsed_time, setup_logging, require_keys

def process_overlap(rna, atac):
    rna_cells = set(rna.obs_names)
    atac_cells = set(atac.obs_names)
    shared_cells = rna_cells & atac_cells
    rna_only_cells = rna_cells - atac_cells
    atac_only_cells = atac_cells - rna_cells
    logging.info(f"n(RNA): {len(rna_cells)}, n(ATAC): {len(atac_cells)}, "
                 f"n(shared): {len(shared_cells)}, n(rna_only): {len(rna_only_cells)}, n(atac_only): {len(atac_only_cells)}")
    return rna_cells, atac_cells, shared_cells, rna_only_cells, atac_only_cells

def extract_atac_object(anndata_org, celllist):
    logging.info(f"Extracting ATAC for {len(celllist)} cells")
    anndata_extract = anndata_org[list(celllist), :].copy()
    obs_col = list(anndata_extract.obs.columns)
    if 'celltype_glue' in obs_col:
        obs_col.remove('celltype_glue')
    anndata_extract.var['modality'] = "Peaks"
    anndata_extract.obs.drop(columns=obs_col, inplace=True)
    return anndata_extract

def extract_rna_object(anndata_org, celllist, celltype_obs):
    logging.info(f"Extracting RNA for {len(celllist)} cells (with {celltype_obs})")
    anndata_extract = anndata_org[list(celllist), :].copy()
    anndata_extract.obs['celltype_glue'] = anndata_extract.obs[celltype_obs]
    obs_col = list(anndata_extract.obs.columns)
    if 'celltype_glue' in obs_col:
        obs_col.remove('celltype_glue')
    var_col = list(anndata_extract.var.columns)
    anndata_extract.var['modality'] = "Gene Expression"
    anndata_extract.obs.drop(columns=obs_col, inplace=True)
    anndata_extract.var.drop(columns=var_col, inplace=True, errors='ignore')
    anndata_extract.uns.clear()
    anndata_extract.obsm.clear()
    anndata_extract.varm.clear()
    anndata_extract.layers.clear()
    anndata_extract.obsp.clear()
    return anndata_extract

def generate_multivi_input(rna, atac, shared_cells, rna_only_cells, atac_only_cells, celltype_obs, output_h5ad_dir, tissue_std):
    logging.info("Assembling MultiVI input AnnData ...")
    atac_shared = extract_atac_object(atac, shared_cells)
    rna_shared = extract_rna_object(rna, shared_cells, celltype_obs)
    data_shared = ad.concat([rna_shared, atac_shared], axis=1, join='inner')
    data_shared.obs['modality'] = 'paired'
    data_shared.obs['celltype_glue'] = list(rna[data_shared.obs_names, :].obs[celltype_obs])
    
    atac_only = extract_atac_object(atac, atac_only_cells)
    atac_only.obs['modality'] = 'ATAC'
    
    rna_only = extract_rna_object(rna, rna_only_cells, celltype_obs)
    rna_only.obs['modality'] = 'GEX'

    adata_mvi = scvi.data.organize_multiome_anndatas(data_shared, rna_only, atac_only)
    adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
    if not sparse.issparse(adata_mvi.X):
        adata_mvi.X = sparse.csr_matrix(adata_mvi.X)
        logging.info("Converted adata_mvi.X to sparse CSR matrix.")
    logging.info(f"MultiVI input shape: {adata_mvi.shape}")
    adata_mvi.write(os.path.join(output_h5ad_dir,f'MultiVI_input.{tissue_std}.h5ad'))
    return adata_mvi

def train_multivi(adata_mvi, tissue_std):
    logging.info("Setting up MULTIVI anndata ...")
    logging.info(f"tissue: {tissue_std}")
    logging.info(f"shape before copy: {adata_mvi.shape} | is_view: {adata_mvi.is_view}")

    if adata_mvi.is_view:
        adata_mvi = adata_mvi.copy()
        logging.info(f"made a copy; is_view now: {adata_mvi.is_view}")

    logging.info(f"n_genes: {(adata_mvi.var['modality'] == 'Gene expression').sum()}")
    logging.info(f"n_regions: {(adata_mvi.var['modality'] == 'Peaks').sum()}")

    logging.info(f"cell with both modality: {(adata_mvi.obs['modality'] == 'paired').sum()}")
    logging.info(f"cell with RNA-only: {(adata_mvi.obs['modality'] == 'expression').sum()}")
    logging.info(f"cell with ATAC-only: {(adata_mvi.obs['modality'] == 'accessibility').sum()}")

    scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="modality")
    model = scvi.model.MULTIVI(
        adata_mvi,
        n_genes=(adata_mvi.var["modality"] == "Gene Expression").sum(),
        n_regions=(adata_mvi.var["modality"] == "Peaks").sum(),
    )
    
    logging.info("Training MULTIVI model ...")
    start_time = time.time()
    scvi.settings.seed = 0
    model.train()
    end_time = time.time()
    logging.info("MULTIVI model training complete.")
    print_elapsed_time(start_time, end_time)

    logging.info("Saving MULTIVI model ...")
    model.save(f"MultiVI.{tissue_std}.model", overwrite=True)
    return model

def get_imputation_values(model, adata_mvi, output_h5ad_dir, tissue_std):
    logging.info("Extracting latent space and embedding")
    adata_mvi.obsm['X_multivi'] = model.get_latent_representation()
    sc.pp.neighbors(adata_mvi, use_rep='X_multivi')
    sc.tl.leiden(adata_mvi, random_state=0)
    sc.tl.umap(adata_mvi, random_state=0)
    sc.pl.umap(adata_mvi, color=["modality", "celltype_glue"], frameon=False,
               save=f'_multiVI.{tissue_std}.png')
    adata_mvi.write(os.path.join(output_h5ad_dir, f'MultiVI_merged.{tissue_std}.h5ad'))

    logging.info("Imputing expression matrix ...")
    imputed_expression = model.get_normalized_expression(adata_mvi)
    imputed_expression.index = imputed_expression.index.str.rsplit('_', n=1).str[0]
    imputed_rna_adata = ad.AnnData(
        X=sparse.csr_matrix(imputed_expression.values.astype(np.float32)),
        obs=pd.DataFrame(index=imputed_expression.index),
        var=pd.DataFrame(index=imputed_expression.columns)
    )
    imputed_rna_adata.write(os.path.join(output_h5ad_dir,f'MultiVI_impute_RNA.{tissue_std}.h5ad'))
    logging.info(f"Wrote imputed RNA: MultiVI_impute_RNA.{tissue_std}.h5ad")

    logging.info("Imputing accessibility matrix ...")
    imputed_accessibility = model.get_accessibility_estimates(adata_mvi)
    imputed_accessibility.index = imputed_accessibility.index.str.rsplit('_', n=1).str[0]
    imputed_atac_adata = ad.AnnData(
        X=sparse.csr_matrix(imputed_accessibility.values.astype(np.float32)),
        obs=pd.DataFrame(index=imputed_accessibility.index),
        var=pd.DataFrame(index=imputed_accessibility.columns)
    )
    imputed_atac_adata.write(os.path.join(output_h5ad_dir,f'MultiVI_impute_ATAC.{tissue_std}.h5ad'))
    logging.info(f"Wrote imputed ATAC: MultiVI_impute_ATAC.{tissue_std}.h5ad")

def run_per_tissue(tissue, output_h5ad_dir, celltype_obs):
    tissue_std = standardize_tissue_name(tissue)
    
    logging.info(f"Loading RNA and ATAC AnnData for tissue: {tissue}")
    rna_path = os.path.join(output_h5ad_dir, f'RNA_removeDoublet.{tissue_std}.h5ad')
    atac_path = os.path.join(output_h5ad_dir, f'ATAC_cellByPeak.{tissue_std}.h5ad')
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)

    logging.info("Integrating all data...")
    rna_cells, atac_cells, shared_cells, rna_only_cells, atac_only_cells = process_overlap(rna, atac)
    adata_mvi = generate_multivi_input(rna, atac, shared_cells, rna_only_cells, atac_only_cells, 
                                       celltype_obs, output_h5ad_dir, tissue_std)

    logging.info(f"Tissue: {tissue_std}")
    logging.info(f"adata_mvi shape: {adata_mvi.shape}")
    logging.info(f"obs columns: {adata_mvi.obs.columns.tolist()}")
    logging.info(f"var columns: {adata_mvi.var.columns.tolist()}")
    logging.info(f"unique adata_mvi.var['modality']: {adata_mvi.var['modality'].unique()}")

    model = train_multivi(adata_mvi, tissue_std)
    get_imputation_values(model, adata_mvi, output_h5ad_dir, tissue_std)
    logging.info(f"Done integration with {tissue}\n")

def main(config_path):
    setup_logging()
    config = load_config(config_path)
    require_keys(config, [
        "paths.workdir", "paths.output_h5ad_dir", "params.tissue", "params.celltype_obs",
    ], context=config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    multivi_output_dir = os.path.join(workdir, 'multivi')
    os.makedirs(multivi_output_dir, exist_ok=True)
    os.makedirs(os.path.join(multivi_output_dir, 'figures'), exist_ok=True)
    os.chdir(multivi_output_dir)
    logging.info(f"Output path set to: {multivi_output_dir}")

    tissue = config['params']['tissue']
    
    celltype_obs = config['params']['celltype_obs']

    if tissue == "---":
        sample_metadata = config['paths']['sample_metadata']
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        
        tissues = sorted(df["tissue"].unique())
        logging.info(f"Data imputation for MULTIPLE tissues: {tissues}")
        
        failed_tissues = []
        for idx, tissue_name in enumerate(tissues, 1):
            logging.info(f"========== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")

            try:
                run_per_tissue(tissue_name, output_h5ad_dir, celltype_obs)
            except Exception as e:
                logging.error(f"Encountered error for tissue {tissue_name}: {str(e)}")
                failed_tissues.append(tissue_name)
    else:
        logging.info(f"========== Processing tissue: {tissue} ==========")

        failed_tissues = []
        try:
            run_per_tissue(tissue, output_h5ad_dir, celltype_obs)
        except Exception as e:
            logging.error(f"Encountered error for tissue {tissue}: {str(e)}")
            failed_tissues.append(tissue)

    if failed_tissues:
        logging.error(f"data_imputation failed for {len(failed_tissues)} tissue(s): {failed_tissues}")
        sys.exit(1)

if __name__ == "__main__":
    setup_logging()
    logging.info(f"torch.cuda.is_available(): {torch.cuda.is_available()}")
    logging.info(f"torch.cuda.device_count(): {torch.cuda.device_count()}")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f"Using device: {device}")

    parser = argparse.ArgumentParser(
        description="Joint imputation using MultiVI"
    )
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    try:
        main(args.config)
    except (FileNotFoundError, KeyError, ValueError) as e:
        logging.error(f"{type(e).__name__}: {e}")
        sys.exit(1)