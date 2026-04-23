#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Finalize MultiVI multiome integration; propagate metadata and layers, add imputed signal, split output.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name

def propagate_obs(adata, refdata, features, new_prefix):
    for feature in features:
        print(f"[INFO] Propagating {feature}")
        if new_prefix != "":
            adata.obs[new_prefix+"."+feature] = refdata.obs[feature].reindex(adata.obs_names).values
        else:
            adata.obs[feature] = refdata.obs[feature].reindex(adata.obs_names).values

def propagate_obs_with_fallback(adata, ref1, ref2, features, new_prefix=''):
    for feature in features:
        print(f"[INFO] Propagating {feature}")
        # 1. Map from first reference (ref1)
        tmp = ref1.obs[feature].reindex(adata.obs_names)
        # 2. For missing, fallback to second reference (ref2)
        if ref2 is not None:
            missing = tmp.isna()
            tmp.loc[missing] = ref2.obs[feature].reindex(adata.obs_names[missing]).values
        adata.obs[new_prefix+feature] = tmp.values
        
def propagate_ann_structures(adata, rna, atac, celllineage_obs, celltype_obs):
    """
    Copy AnnData var, uns, obsm, varm keys from both rna and atac into adata if they don't exist yet.
    Only features/cells in intersection are copied (other entries are set NaN if required).
    """
    # --- VAR columns (gene related) ---
    rna_var_cols = [
        'gene_ids', 'feature_types', 'genome', 'mt', 'ribo', 'hb', 'n_cells_by_counts',
        'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'n_cells', 'highly_variable',
        'means', 'dispersions', 'dispersions_norm'
    ]
    for col in rna_var_cols:
        print(f"[INFO] Adding var column (rna): {col}")
        values = pd.Series(index=adata.var_names, dtype=rna.var[col].dtype)
        values.loc[rna.var_names] = rna.var[col].values
        adata.var[col] = values.values

    # --- UNS keys (RNA related) ---
    uns_copy_keys = [
        'donorID_colors','hvg','leiden','leiden_before_harmony','leiden_colors',
        celltype_obs+'_colors', celllineage_obs+'_colors','neighbors','pca',
        'rank_genes_groups','sampleID_colors','scrublet','umap'
    ]
    for key in uns_copy_keys:
        if key in rna.uns and key not in adata.uns:
            print(f"[INFO] Adding uns key (rna): {key}")
            adata.uns[key] = rna.uns[key]

    # --- OBSM keys ---
    obsm_copy_keys = ['X_pca','X_pca_before_harmony','X_pca_harmony']
    for key in obsm_copy_keys:
        if key in rna.obsm:
            print(f"[INFO] Adding obsm key (rna): {key}")
            arr = rna.obsm[key]             # (n_rna, d)
            n_dim = arr.shape[1]
            arr_c = np.full((adata.n_obs, n_dim), np.nan, dtype=arr.dtype)
            common = adata.obs_names.isin(rna.obs_names)
            idx_glue = np.where(common)[0]
            idx_rna  = rna.obs_names.get_indexer(adata.obs_names[common])
            arr_c[idx_glue, :] = arr[idx_rna, :]
            adata.obsm[key] = arr_c

    for key in ['X_lsi']:
        if key in atac.obsm:
            print(f"[INFO] Adding obsm key (atac): {key}")
            arr = atac.obsm[key]             # (n_rna, d)
            n_dim = arr.shape[1]
            arr_c = np.full((adata.n_obs, n_dim), np.nan, dtype=arr.dtype)
            common = adata.obs_names.isin(atac.obs_names)
            idx_glue = np.where(common)[0]
            idx_atac  = atac.obs_names.get_indexer(adata.obs_names[common])
            arr_c[idx_glue, :] = arr[idx_atac, :]
            adata.obsm[key] = arr_c
            
        
    # --- OBSP keys ---
    for key in ['connectivities','distances']:
        if key in rna.obsp and key not in adata.obsp:
            print(f"[INFO] Adding obsp key: {key}")
            idx = adata.obs_names.intersection(rna.obs_names)
            array = rna.obsp[key]
            arr_c = np.full((adata.n_obs, adata.n_obs), np.nan)
            adata_idx = adata.obs_names.get_indexer(idx)
            arr_c[np.ix_(adata_idx, adata_idx)] = array[np.ix_(rna.obs_names.get_indexer(idx), rna.obs_names.get_indexer(idx))]
            adata.obsp[key] = sparse.csr_matrix(arr_c)
    print("[INFO] Finished copying all relevant AnnData structures to adata.")

def propagate_layers_from_rna(adata, rna, rna_layer_names):
    print("[INFO] Propagating RNA layers for gene expression features...")

    # genes that are RNA features in adata and present in rna
    mask_genes = (adata.var["is_rna_feature"] == 1) & adata.var_names.isin(rna.var_names)
    # cells present in both
    mask_cells = adata.obs_names.isin(rna.obs_names)

    adata_row = np.where(mask_cells)[0]
    adata_col = np.where(mask_genes)[0]

    cells_common = adata.obs_names[mask_cells]
    genes_common = adata.var_names[mask_genes]

    # these will now work if rna.obs_names / var_names are unique
    rna_row = rna.obs_names.get_indexer(cells_common)
    rna_col = rna.var_names.get_indexer(genes_common)

    for layer_name in rna_layer_names:
        print(f"[INFO] Propagating layer: {layer_name}")
        if layer_name not in rna.layers:
            print(f"  [WARN] Layer '{layer_name}' not found in RNA; skipping.")
            continue

        block = rna.layers[layer_name][np.ix_(rna_row, rna_col)]
        if not isinstance(block, np.ndarray):
            block = block.toarray()

        adata_layer = np.zeros(adata.shape, dtype=block.dtype)
        adata_layer[np.ix_(adata_row, adata_col)] = block
        adata.layers[layer_name] = sparse.csr_matrix(adata_layer)

def add_imputed_signal(adata, rna_impute, atac_impute):
    print("[INFO] Integrating imputed expression and accessibility signal into combined AnnData ...")
    genes = adata.var_names[adata.var['modality'] == "Gene Expression"]
    peaks = adata.var_names[adata.var['modality'] == "Peaks"]

    imp = np.zeros(adata.shape, dtype=float)

    # Imputed genes
    imp[:, adata.var_names.get_indexer(genes)] = (
        rna_impute[:, genes].X.toarray()
        if hasattr(rna_impute.X, "toarray")
        else rna_impute[:, genes].X
    )

    # Imputed peaks
    imp[:, adata.var_names.get_indexer(peaks)] = (
        atac_impute[:, peaks].X.toarray()
        if hasattr(atac_impute.X, "toarray")
        else atac_impute[:, peaks].X
    )

    adata.layers['imputation_signal'] = sparse.csr_matrix(imp)
    print("[INFO] Added layer: imputation_signal.")

    
def split_and_save_multiome(adata, output_h5ad_dir, tissue_std, celllineage_obs, celltype_obs):
    # RNA-only
    print("[INFO] Extracting and writing RNA AnnData ...")
    rna_final = adata[adata.obs['modality']!='accessibility-only', adata.var['is_rna_feature']==1].copy()
    rna_keep = [
        'rna.sampleID', 'rna.leiden', 'donorID', 'cellbarcode', 'total_counts', 'n_genes_by_counts', 
        'n_genes', 'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_exon_reads', 'MALAT1_CPM', 
        'doublet_score', 'doublet_probabilities', celllineage_obs, celltype_obs]
    drop_cols = [col for col in rna_final.obs.columns if col not in rna_keep]
    rna_final.obs.drop(columns=drop_cols, inplace=True)
    rna_file = os.path.join(output_h5ad_dir, f'multiome_final.RNA.{tissue_std}.h5ad')
    rna_final.write(rna_file, compression="gzip")
    print(f"[INFO] Wrote RNA AnnData: {rna_file}")
    sc.pl.umap(rna_final, color=[celllineage_obs], save=f"_multiome_final.RNA.{celllineage_obs}.{tissue_std}.png")
    sc.pl.umap(rna_final, color=[celltype_obs], save=f"_multiome_final.RNA.{celltype_obs}.{tissue_std}.png")
    
    # ATAC-only
    print("[INFO] Extracting and writing ATAC AnnData ...")
    atac_final = adata[adata.obs['modality']!='expression-only', adata.var['is_atac_feature']==1].copy()
    atac_keep = [
        'atac.sampleID', 'atac.leiden', 'modality', 'celltype_glue', 'n_fragment', 'frac_dup', 
        'frac_mito', 'tsse', 'donorID', 'cellbarcode', celllineage_obs, celltype_obs]
    drop_cols2 = [col for col in atac_final.obs.columns if col not in atac_keep]
    atac_final.obs.drop(columns=drop_cols2, inplace=True)
    if 'X_pca_before_harmony' in atac_final.obsm: del atac_final.obsm['X_pca_before_harmony']
    if 'X_pca_harmony' in atac_final.obsm: del atac_final.obsm['X_pca_harmony']
    if 'X_spectral' in atac_final.obsm: del atac_final.obsm['X_spectral']
    if 'X_umap_before_harmony' in atac_final.obsm: del atac_final.obsm['X_umap_before_harmony']
    if 'CPM' in atac_final.layers: del atac_final.layers['CPM']
    if 'rawcounts' in atac_final.layers: del atac_final.layers['rawcounts']
    atac_file = os.path.join(output_h5ad_dir, f'multiome_final.ATAC.{tissue_std}.h5ad')
    atac_final.write(atac_file, compression="gzip")
    print(f"[INFO] Wrote ATAC AnnData: {atac_file}")
    sc.pl.umap(atac_final, color=[celllineage_obs], save=f"_multiome_final.ATAC.{celllineage_obs}.{tissue_std}.png")
    sc.pl.umap(atac_final, color=[celltype_obs], save=f"_multiome_final.ATAC.{celltype_obs}.{tissue_std}.png")
    
def run_per_tissue(tissue, output_h5ad_dir, celllineage_obs, celltype_obs, imputation):
    
    tissue_std = standardize_tissue_name(tissue)
    
    print(f"[INFO] Loading AnnData for tissue: {tissue}")
    adata = ad.read_h5ad(os.path.join(output_h5ad_dir, f'Multiome_merged.GLUE.{tissue_std}.h5ad'))
    adata.obs['tissue'] = tissue
    rna = ad.read_h5ad(os.path.join(output_h5ad_dir, f'RNA.GLUE.{tissue_std}.h5ad'))
    atac = ad.read_h5ad(os.path.join(output_h5ad_dir, f'ATAC.GLUE.{tissue_std}.h5ad'))

    # propagate obs/metadata
    print("[INFO] Propagating RNA/ATAC obs fields ...")
    propagate_obs(adata, rna, ['sampleID','leiden'], 'rna')
    propagate_obs(adata, rna, ['cellbarcode', 'total_counts', 'n_genes_by_counts', 'n_genes',
                               'pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb', 'pct_exon_reads',
                               'MALAT1_CPM', 'doublet_score', 'doublet_probabilities'], '')
    propagate_obs_with_fallback(adata, rna, atac, ['donorID'], new_prefix='')
    propagate_obs(adata, atac, ['sampleID', 'leiden'], 'atac')
    propagate_obs(adata, atac, ['n_fragment', 'frac_dup', 'frac_mito', 'tsse'], '')

    # add cell lineage
    celltype_to_lineage = dict(zip(rna.obs[celltype_obs], rna.obs[celllineage_obs]))
    adata.obs['cell_lineage'] = adata.obs['celltype_glue'].map(celltype_to_lineage)
    adata.obs['celltype_broad'] = adata.obs['celltype_glue']

    propagate_ann_structures(adata, rna, atac, celllineage_obs, celltype_obs)

    propagate_layers_from_rna(adata, rna, ['CPM', 'rawcounts'])

    # Add imputation layers
    if imputation:
        rna_impute = ad.read_h5ad(os.path.join(output_h5ad_dir, f"MultiVI_impute_RNA.{tissue_std}.h5ad"))
        atac_impute = ad.read_h5ad(os.path.join(output_h5ad_dir, f"MultiVI_impute_ATAC.{tissue_std}.h5ad"))
        add_imputed_signal(adata, rna_impute, atac_impute)

    # Save main complete file
    final_file = os.path.join(output_h5ad_dir, f'multiome_final.{tissue_std}.h5ad')
    adata.write(final_file)
    print(f"[INFO] Wrote: {final_file}")

    # plot UMAP
    print(f"[INFO] plotting UMAPs")
    sc.pl.umap(adata, color=[celllineage_obs], save=f"_multiome_final.cell_lineage.{tissue_std}.png")
    sc.pl.umap(adata, color=[celltype_obs], save=f"_multiome_final.celltype_broad.{tissue_std}.png")
    sc.pl.umap(adata, color=["modality"], save=f"_multiome_final.modality.{tissue_std}.png")

    split_and_save_multiome(adata, output_h5ad_dir, tissue_std, celllineage_obs, celltype_obs)

    print(f"[DONE] Finish the whole workflow! Hooray!!!")

def main(config_path):
    config = load_config(config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    final_output_dir = os.path.join(workdir, 'final')
    os.makedirs(final_output_dir, exist_ok=True)
    os.chdir(final_output_dir)
    
    tissue = config['params']['tissue']
    
    celltype_obs = config['params']['celltype_obs']
    celllineage_obs = config['params']['celllineage_obs']
    imputation = config['params']['imputation']
    
    if tissue == "---":
        sample_metadata = config['paths']['sample_metadata']
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        
        tissues = sorted(df["tissue"].unique())
        print(f"[INFO] Data integration for MULTIPLE tissues: {tissues}")
        
        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n========== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")
            
            try:
                run_per_tissue(working_tissue, output_h5ad_dir, celllineage_obs, celltype_obs, imputation)
            except Exception as e:
                print(f"[ERROR] Encountered error for tissue {working_tissue}: {str(e)}")
    else:
        print(f"\n========== Processing tissue: {tissue} ==========")
        
        try:
            run_per_tissue(tissue, output_h5ad_dir, celllineage_obs, celltype_obs, imputation)
        except Exception as e:
            print(f"[ERROR] Encountered error for tissue {tissue}: {str(e)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Finalize MultiVI multiome outputs, propagate metadata/layers, add imputed signal, and split output."
    )
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    main(args.config)