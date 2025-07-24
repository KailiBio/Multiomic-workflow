#!/usr/bin/env python3
"""
Author: Kaili Fan
Description: Joint scRNA/scATAC imputation with MultiVI
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import snapatac2 as snap
import scvi
import time
scvi.settings.seed = 0
import torch
from scipy import sparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name, print_elapsed_time

def process_overlap(rna, atac):
    rna_cells = set(rna.obs_names)
    atac_cells = set(atac.obs_names)
    shared_cells = rna_cells & atac_cells
    rna_only_cells = rna_cells - atac_cells
    atac_only_cells = atac_cells - rna_cells
    print(f"[INFO] n(RNA): {len(rna_cells)}, n(ATAC): {len(atac_cells)}, \
n(shared): {len(shared_cells)}, n(rna_only): {len(rna_only_cells)}, n(atac_only): {len(atac_only_cells)}")
    return rna_cells, atac_cells, shared_cells, rna_only_cells, atac_only_cells

def extract_atac_object(anndata_org, celllist):
    print(f"[INFO] Extracting ATAC for {len(celllist)} cells")
    anndata_extract = anndata_org[list(celllist), :].copy()
    obs_col = list(anndata_extract.obs.columns)
    if 'celltype_glue' in obs_col:
        obs_col.remove('celltype_glue')
    anndata_extract.var['modality'] = "Peaks"
    anndata_extract.obs.drop(columns=obs_col, inplace=True)
    return anndata_extract

def extract_rna_object(anndata_org, celllist, celltype_col):
    print(f"[INFO] Extracting RNA for {len(celllist)} cells (celltype_col={celltype_col})")
    anndata_extract = anndata_org[list(celllist), :].copy()
    anndata_extract.obs['celltype_glue'] = anndata_extract.obs[celltype_col]
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

def generate_multivi_input(rna, atac, shared_cells, rna_only_cells, atac_only_cells, celltype_col):
    print("[INFO] Assembling MultiVI input AnnData ...")
    atac_shared = extract_atac_object(atac, shared_cells)
    rna_shared = extract_rna_object(rna, shared_cells, celltype_col)
    data_shared = ad.concat([rna_shared, atac_shared], axis=1, join='inner')
    data_shared.obs['modality'] = 'paired'
    data_shared.obs['celltype_glue'] = list(rna[data_shared.obs_names, :].obs[celltype_col])
    
    atac_only = extract_atac_object(atac, atac_only_cells)
    atac_only.obs['modality'] = 'ATAC'
    
    rna_only = extract_rna_object(rna, rna_only_cells, celltype_col)
    rna_only.obs['modality'] = 'GEX'

    adata_mvi = scvi.data.organize_multiome_anndatas(data_shared, rna_only, atac_only)
    adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
    if not sparse.issparse(adata_mvi.X):
        adata_mvi.X = sparse.csr_matrix(adata_mvi.X)
        print("[INFO] Converted adata_mvi.X to sparse CSR matrix.")
    print(f"[INFO] MultiVI input shape: {adata_mvi.shape}")
    return adata_mvi

def train_multivi(adata_mvi, tissue_std):
    print("[INFO] Setting up MULTIVI anndata ...")
    scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key="modality")
    model = scvi.model.MULTIVI(
        adata_mvi,
        n_genes=(adata_mvi.var["modality"] == "Gene Expression").sum(),
        n_regions=(adata_mvi.var["modality"] == "Peaks").sum(),
    )
    
    print("[INFO] Training MULTIVI model ...")
    start_time = time.time()
    scvi.settings.seed = 0
    model.train()
    end_time = time.time()
    print("[INFO] MULTIVI model training complete.")
    print_elapsed_time(start_time, end_time)

    print("[INFO] Saving MULTIVI model ...")
    model.save(f"MultiVI.{tissue_std}.model", overwrite=True)
    return model

def get_imputation_values(model, adata_mvi, output_h5ad_dir, tissue_std):
    print("[INFO] Extracting latent space and embedding")
    adata_mvi.obsm['X_multivi'] = model.get_latent_representation()
    sc.pp.neighbors(adata_mvi, use_rep='X_multivi')
    sc.tl.leiden(adata_mvi, random_state=0)
    sc.tl.umap(adata_mvi, random_state=0)
    sc.pl.umap(adata_mvi, color=["modality", "celltype_glue"], frameon=False,
               save=f'_multiVI.{tissue_std}.png')
    adata_mvi.write(os.path.join(output_h5ad_dir, f'MultiVI_merged.{tissue_std}.h5ad'), compression="gzip")

    print("[INFO] Imputing expression matrix ...")
    imputed_expression = model.get_normalized_expression(adata_mvi)
    imputed_expression.index = imputed_expression.index.str.rsplit('_', n=1).str[0]
    imputed_rna_adata = ad.AnnData(
        X=sparse.csr_matrix(imputed_expression.values.astype(np.float32)),
        obs=pd.DataFrame(index=imputed_expression.index),
        var=pd.DataFrame(index=imputed_expression.columns)
    )
    imputed_rna_adata.write(os.path.join(output_h5ad_dir,f'MultiVI_impute_RNA.{tissue_std}.h5ad'), compression="gzip")
    print(f"[INFO] Wrote imputed RNA: MultiVI_impute_RNA.{tissue_std}.h5ad")

    print("[INFO] Imputing accessibility matrix ...")
    imputed_accessibility = model.get_accessibility_estimates(adata_mvi)
    imputed_accessibility.index = imputed_accessibility.index.str.rsplit('_', n=1).str[0]
    imputed_atac_adata = ad.AnnData(
        X=sparse.csr_matrix(imputed_accessibility.values.astype(np.float32)),
        obs=pd.DataFrame(index=imputed_accessibility.index),
        var=pd.DataFrame(index=imputed_accessibility.columns)
    )
    imputed_atac_adata.write(os.path.join(output_h5ad_dir,f'MultiVI_impute_ATAC.{tissue_std}.h5ad'), compression="gzip")
    print(f"[INFO] Wrote imputed ATAC: MultiVI_impute_ATAC.{tissue_std}.h5ad")

def run_per_tissue(tissue, output_h5ad_dir, celltype_col):
    tissue_std = standardize_tissue_name(tissue)
    print(f"[INFO] Loading RNA and ATAC AnnData for tissue: {tissue}")
    rna_path = os.path.join(output_h5ad_dir, f'RNA_removeDoublet.{tissue_std}.h5ad')
    atac_path = os.path.join(output_h5ad_dir, f'ATAC_cellByPeak.{tissue_std}.h5ad')
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)
    rna_cells, atac_cells, shared_cells, rna_only_cells, atac_only_cells = process_overlap(rna, atac)
    adata_mvi = generate_multivi_input(
        rna, atac, shared_cells, rna_only_cells, atac_only_cells, celltype_col)
    model = train_multivi(adata_mvi, tissue_std)
    get_imputation_values(model, adata_mvi, output_h5ad_dir, tissue_std)
    print(f"[INFO] Done with {tissue}\n")

def main(config_path):
    config = load_config(config_path)
    
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    multivi_output_dir = os.path.join(config['paths']['workdir'], 'multivi')
    os.makedirs(multivi_output_dir, exist_ok=True)
    os.makedirs(os.path.join(multivi_output_dir, 'figures'), exist_ok=True)
    os.chdir(multivi_output_dir)
    tissue = config['params']['tissue']
    celltype_col = config['params']['celltype_col']

    if tissue == "---":
        working_df = pd.read_csv(config['paths']['sample_metadata'], sep='\t', header=None,
                                 names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        tissues = sorted(working_df["tissue"].unique())
        for idx, working_tissue in enumerate(tissues, 1):
            print(f"\n============== [Main] Processing tissue: {working_tissue} ({idx}/{len(tissues)}) ==============")
            try:
                run_per_tissue(working_tissue, output_h5ad_dir, celltype_col)
            except Exception as e:
                print(f"[ERROR] Encountered error for tissue {working_tissue}: {str(e)}")
    else:
        print(f"\n============== [Main] Processing tissue: {tissue} ==============")
        try:
            run_per_tissue(tissue, output_h5ad_dir, celltype_col)
        except Exception as e:
            print(f"[ERROR] Encountered error for tissue {tissue}: {str(e)}")

if __name__ == "__main__":
    print(f"[INFO] torch.cuda.is_available(): {torch.cuda.is_available()}")
    print(f"[INFO] torch.cuda.device_count(): {torch.cuda.device_count()}")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"[INFO] Using device: {device}")
    
    parser = argparse.ArgumentParser(
        description="Joint imputation using MultiVI"
    )
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    main(args.config)