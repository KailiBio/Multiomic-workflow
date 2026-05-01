#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Annotate ATAC cells via GLUE integrating scRNA and scATAC data.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

import os
import sys
import argparse
import time
from itertools import chain
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
import scglue
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name, print_elapsed_time
from multiome_integration.integration_plots import plot_multiome_celltype_count_bar

def prepare_rna(rna, gtf, output_h5ad_dir, suffix, celltype_obs):
    print(f"[INFO] Preparing RNA AnnData for '{suffix}' using cell type column '{celltype_obs}' ...")
    rna.obs['celltype_glue'] = rna.obs[celltype_obs]
    
    print(f"[INFO] Computing highly variable genes (HVGs) ...")
    sc.pp.highly_variable_genes(rna, flavor='seurat')
    print(f"Getting {rna.var['highly_variable'].value_counts()} highly variable genes")
    
    print(f"[INFO] Plotting RNA UMAP before GLUE integration ...")
    sc.pl.umap(rna, color="celltype_glue", save=f'_RNA.beforeGLUE.{suffix}.png')

    print(f"[INFO] Annotating genes using GTF: {gtf} ...")
    scglue.data.get_gene_annotation(rna, gtf=gtf, gtf_by="gene_name")
    good = (
        (~rna.var['chromStart'].isnull()) &
        (~rna.var['chromEnd'].isnull()) &
        (~np.isinf(rna.var['chromStart'])) &
        (~np.isinf(rna.var['chromEnd']))
    )
    num_good = good.sum()
    print(f"Retaining {num_good} genes with usable positions.")

    rna = rna[:, good].copy()
    rna.var['chromStart'] = rna.var['chromStart'].astype(int)
    rna.var['chromEnd'] = rna.var['chromEnd'].astype(int)
    if "artif_dupl" in rna.var:
        rna.var.drop(columns=['artif_dupl'], inplace=True)
        
    out_path = os.path.join(output_h5ad_dir, f"RNA.beforeGLUE.{suffix}.h5ad")
    rna.write(out_path)
    print(f"[INFO] Wrote prepared RNA AnnData: {out_path}")
    return rna

def prepare_atac(atac_all, gtf, output_h5ad_dir, suffix, n_features=50000):
    print(f"[INFO] Preparing ATAC AnnData for '{suffix}' ...")
    
    print(f"[INFO] Selecting {n_features} top features ...")
    snap.pp.select_features(atac_all, n_features=n_features)
    atac = atac_all[:, atac_all.var['selected']].copy()

    print(f"[INFO] Calculating LSI for ATAC ...")
    scglue.data.lsi(atac, n_components=100, n_iter=15)
    
    print(f"[INFO] Computing neighbors, and Leiden clustering for ATAC ...")
    sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
    sc.tl.umap(atac)
    sc.tl.leiden(atac, resolution=0.5, flavor="igraph")
    print(f"[INFO] Plotting ATAC UMAP before GLUE integration ...")
    sc.pl.umap(atac, color="leiden", save=f'_ATAC.beforeGLUE.{suffix}.png')

    print(f"[INFO] Assigning peak coordinates from var_names ...")
    split = atac.var_names.str.split(r"[:-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)

    out_path = os.path.join(output_h5ad_dir, f"ATAC.beforeGLUE.{suffix}.h5ad")
    atac.write(out_path)
    print(f"[INFO] Wrote prepared ATAC AnnData: {out_path}")
    return atac

def train_glue(rna, atac, glue_output_dir, tissue_std):
    print(f"[INFO] Building GLUE guidance graph for '{tissue_std}' ...")
    guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
    scglue.graph.check_graph(guidance, [rna, atac])
    print("[INFO] Guidance graph created and checked.")

    print(f"[INFO] Configuring GLUE datasets ...")
    scglue.models.configure_dataset(
        rna, "NB", use_highly_variable=True,
        use_layer="rawcounts", use_rep="X_pca",
        use_obs_names=True,
        use_cell_type="celltype_glue"
    )
    scglue.models.configure_dataset(
        atac, "NB", use_highly_variable=True,
        use_rep="X_lsi",
        use_obs_names=True
    )

    print(f"[INFO] Subgraph with highly variable features ...")
    guidance_hvf = guidance.subgraph(chain(
        rna.var.query("highly_variable").index,
        atac.var.query("highly_variable").index
    )).copy()

    print(f"[INFO] Training GLUE model for '{tissue_std}' ... This may take a while.")
    start_time = time.time()
    glue = scglue.models.fit_SCGLUE(
        {"rna": rna, "atac": atac}, guidance_hvf,
        model=scglue.models.PairedSCGLUEModel,
        fit_kws={"directory": glue_output_dir}, 
        skip_balance=True
    )
    end_time = time.time()
    print(f"[INFO] GLUE model training complete.")
    print_elapsed_time(start_time, end_time)

    glue_file = os.path.join(glue_output_dir, f"glue.{tissue_std}.dill")
    glue.save(glue_file)
    print(f"[INFO] Saved GLUE model: {glue_file}")

    print(f"[INFO] Checking integration consistency ...")
    dx = scglue.models.integration_consistency(glue, {"rna": rna, "atac": atac}, guidance_hvf)

    # Save the training curve
    fig_dir = os.path.join(glue_output_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    ax = sns.lineplot(x="n_meta", y="consistency", data=dx)
    ax.axhline(y=0.05, c="darkred", ls="--")
    fig = ax.get_figure()
    training_curve = os.path.join(fig_dir, f"GLUE_trainingCurve.{tissue_std}.png")
    fig.savefig(training_curve)
    plt.close(fig)
    print(f"[INFO] GLUE integration consistency curve saved: {training_curve}")

    print(f"[INFO] GLUE model training completed for tissue: {tissue_std}.")
    return glue

def rebuild_obs_from_modalities(
    adata_glue: AnnData,
    rna: AnnData,
    atac: AnnData,
    keys,
    priority: str = "rna",
    as_category: bool = True,
):
    """
    Rebuild selected .obs columns in adata_glue from rna and atac.

    Parameters
    ----------
    adata_glue : AnnData
        Concatenated object with union of cells (axis=1 concat).
    rna, atac : AnnData
        Original modality-specific objects.
    keys : list or str
        obs column names to rebuild (e.g. ["donorID", "tissue", "celltype_glue"]).
    priority : {"rna", "atac"}
        Which modality to prefer when both have a value for the same cell.
    as_category : bool
        If True, cast merged column back to categorical.
    """
    if isinstance(keys, str):
        keys = [keys]

    for key in keys:
        base = pd.Series(index=adata_glue.obs_names, dtype=object)

        # RNA column
        if key in rna.obs:
            s_rna = rna.obs[key].reindex(adata_glue.obs_names)
            if pd.api.types.is_categorical_dtype(s_rna):
                s_rna = s_rna.astype(object)
        else:
            s_rna = base.copy()

        # ATAC column
        if key in atac.obs:
            s_atac = atac.obs[key].reindex(adata_glue.obs_names)
            if pd.api.types.is_categorical_dtype(s_atac):
                s_atac = s_atac.astype(object)
        else:
            s_atac = base.copy()

        # Merge
        if priority == "rna":
            merged = s_rna.fillna(s_atac)
        elif priority == "atac":
            merged = s_atac.fillna(s_rna)
        else:
            raise ValueError("priority must be 'rna' or 'atac'")

        if as_category:
            merged = merged.astype("category")

        if key in adata_glue.obs:
            del adata_glue.obs[key]
        adata_glue.obs[key] = merged

    return adata_glue


def add_modality_label(adata_glue: AnnData, rna: AnnData, atac: AnnData):
    """
    Add .obs['modality'] with values:
    - 'expression-only'
    - 'accessibility-only'
    - 'multi-modality'
    """
    import pandas as pd

    in_rna  = adata_glue.obs_names.isin(rna.obs_names)
    in_atac = adata_glue.obs_names.isin(atac.obs_names)

    modality = pd.Series("unknown", index=adata_glue.obs_names)
    modality[in_rna & ~in_atac] = "expression-only"
    modality[~in_rna & in_atac] = "accessibility-only"
    modality[in_rna & in_atac]  = "multi-modality"

    adata_glue.obs["modality"] = modality.astype("category")
    return adata_glue


def transfer_umap_from_axis0(
    adata_glue_all: AnnData,  # axis=0 concat of [rna, atac]
    adata_glue: AnnData,      # axis=1 concat of {'rna': rna, 'atac': atac}
    rna: AnnData,
    atac: AnnData,
    key: str = "X_umap",
):
    """
    Transfer UMAP coordinates from adata_glue_all (axis=0 concat) to adata_glue (axis=1 concat),
    using RNA UMAP for overlapping cells and ATAC UMAP only for ATAC-only cells.

    Assumes adata_glue_all was built as:
        adata_glue_all = ad.concat([rna, atac])
    so rows 0..n_rna-1 == rna, rows n_rna..n_rna+n_atac-1 == atac.
    """
    X_all = adata_glue_all.obsm[key]

    n_rna  = rna.n_obs
    n_atac = atac.n_obs

    X_rna  = X_all[:n_rna, :]
    X_atac = X_all[n_rna:n_rna + n_atac, :]

    n_cells = adata_glue.n_obs
    n_dim   = X_rna.shape[1]
    X_new   = np.full((n_cells, n_dim), np.nan, dtype=X_rna.dtype)

    # RNA cells (including overlaps)
    mask_rna      = adata_glue.obs_names.isin(rna.obs_names)
    idx_glue_rna  = np.where(mask_rna)[0]
    idx_rna       = rna.obs_names.get_indexer(adata_glue.obs_names[mask_rna])
    X_new[idx_glue_rna, :] = X_rna[idx_rna, :]

    # ATAC-only cells
    mask_atac_only = adata_glue.obs_names.isin(atac.obs_names) & ~mask_rna
    idx_glue_atac  = np.where(mask_atac_only)[0]
    idx_atac       = atac.obs_names.get_indexer(adata_glue.obs_names[mask_atac_only])
    X_new[idx_glue_atac, :] = X_atac[idx_atac, :]

    adata_glue.obsm[key] = X_new
    return adata_glue

def annotate_and_merge(rna, atac, glue, output_h5ad_dir, tissue_std):
    print(f"[INFO] Annotating and merging RNA/ATAC for '{tissue_std}' ...")

    #------
    # RNA
    print("  [RNA] Encoding data ...")
    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    rna.obs['modality'] = 'expression'
    sc.pp.neighbors(rna, use_rep="X_glue", metric="cosine")
    sc.tl.umap(rna)
    print("  [RNA] Saving UMAP and AnnData ...")
    sc.pl.umap(rna, color=["celltype_glue"], save=f"_GLUEembed.RNA.{tissue_std}.png")
    rna.write(os.path.join(output_h5ad_dir, f'RNA.GLUE.{tissue_std}.h5ad'))

    #------
    # ATAC
    print("  [ATAC] Encoding data, running neighbors/UMAP/leiden, annotating cells ...")
    atac.obsm["X_glue"] = glue.encode_data("atac", atac)
    atac.obs['modality'] = 'accessibility'
    atac_cell_type = glue.classify_data("atac", atac)

    df = pd.DataFrame(atac_cell_type)
    df['Max_Column'] = df.idxmax(axis=1)
    df.to_csv(f'ATAC_celltype.{tissue_std}.tsv', sep='\t', index=False)

    atac.obs['celltype_glue'] = df.loc[atac.obs_names.tolist(),'Max_Column']
    sc.pp.neighbors(atac, use_rep="X_glue", metric="cosine")
    sc.tl.umap(atac)
    print("  [ATAC] Saving UMAP and AnnData ...")
    sc.pl.umap(atac, color=["celltype_glue"], save=f"_GLUEembed.ATAC.{tissue_std}.png")
    atac.write(os.path.join(output_h5ad_dir, f'ATAC.GLUE.{tissue_std}.h5ad'))

    #------
    # Merged AnnData
    print("  [MERGE] Concatenating and writing AnnData, embedding UMAP ...")
    
    #clean PCs from individual objects (to avoid varm shape issues)
    for adata in (rna, atac):
        if "PCs" in adata.varm:
            del adata.varm["PCs"]

    # 1) axis=0 concat for integrated neighbors/UMAP
    print("  [MERGE] Generate axis=0 UMAP")
    adata_glue_all = ad.concat([rna, atac]) 
    sc.pp.neighbors(adata_glue_all, use_rep="X_glue", metric="cosine")
    sc.tl.umap(adata_glue_all)
            
    # 2) axis=1 concat for merged features (union of cells)
    print("  [MERGE] concat cells")
    adata_glue = ad.concat({'rna': rna, 'atac': atac},
                           axis=1,          # stack features
                           join='outer',    # union of cells
                           merge='first',   # how to merge obs/var annotations
                           fill_value=0     # missing entries become 0 (keeps it sparse)
                          )

    # 3) transfer UMAP from adata_glue_all to adata_glue
    print("  [MERGE] add UMAP")
    adata_glue = transfer_umap_from_axis0(
        adata_glue_all=adata_glue_all,
        adata_glue=adata_glue,
        rna=rna,
        atac=atac,
        key="X_umap",
    )

    # 4) build a combined X_glue on adata_glue, preferring RNA on overlaps
    print("  [MERGE] add X_glue")
    X_rna  = rna.obsm["X_glue"]
    X_atac = atac.obsm["X_glue"]
    n_dim  = X_rna.shape[1]
    
    X_glue_new = np.zeros((adata_glue.n_obs, n_dim), dtype=X_rna.dtype)
    
    mask_rna       = adata_glue.obs_names.isin(rna.obs_names)
    idx_glue_rna   = np.where(mask_rna)[0]
    idx_rna        = rna.obs_names.get_indexer(adata_glue.obs_names[mask_rna])
    X_glue_new[idx_glue_rna, :] = X_rna[idx_rna, :]
    
    mask_atac_only = adata_glue.obs_names.isin(atac.obs_names) & ~mask_rna
    idx_glue_atac  = np.where(mask_atac_only)[0]
    idx_atac       = atac.obs_names.get_indexer(adata_glue.obs_names[mask_atac_only])
    X_glue_new[idx_glue_atac, :] = X_atac[idx_atac, :]
    adata_glue.obsm["X_glue"] = X_glue_new
    
    # 5) add modality label 
    print("  [MERGE] add modality")
    # (expression-only / accessibility-only / multi-modality)
    adata_glue = add_modality_label(adata_glue, rna, atac)
    # rna/atac feature
    in_rna_var  = adata_glue.var_names.isin(rna.var_names)
    in_atac_var = adata_glue.var_names.isin(atac.var_names)
    adata_glue.var["is_rna_feature"]  = in_rna_var.astype(int) 
    adata_glue.var["is_atac_feature"] = in_atac_var.astype(int)

    # 6) rebuild key obs fields from modalities
    print("  [MERGE] fix obs fields")
    adata_glue = rebuild_obs_from_modalities(adata_glue, rna, atac,
                                             keys=["sampleID", "donorID", "tissue", "celltype_glue"],
                                             priority="rna",)

    # 7) visualize final merged object
    print("  [MERGE] plot and save anndata")
    sc.pl.umap(adata_glue, color=["celltype_glue", "modality"], wspace=0.95, ncols=1, 
               save=f"_GLUEembed.merged.{tissue_std}.png")

    if "predicted_doublet" in adata_glue.obs: 
        del adata_glue.obs["predicted_doublet"]
    adata_glue.write(os.path.join(output_h5ad_dir, f'Multiome_merged.GLUE.{tissue_std}.h5ad'))

    print(f"[INFO] Annotation and merging completed for '{tissue_std}'.")

def run_per_tissue(tissue_std, rna, atac, output_h5ad_dir, glue_output_dir, gtf, celltype_obs, n_features=50000):
    
    # Prepare
    print("\nRNA...")
    rna_prepared = prepare_rna(rna, gtf, output_h5ad_dir, tissue_std, celltype_obs)
    print("\nATAC...")
    atac_prepared = prepare_atac(atac, gtf, output_h5ad_dir, tissue_std, n_features=n_features)

    # Train GLUE & annotate
    glue = train_glue(rna_prepared, atac_prepared, glue_output_dir, tissue_std)
    annotate_and_merge(rna_prepared, atac_prepared, glue, output_h5ad_dir, tissue_std)
    plot_multiome_celltype_count_bar(rna_prepared, atac_prepared, os.path.join(glue_output_dir, 'figures'), tissue_std)
    print(f"[INFO] Finished GLUE for tissue '{tissue_std}'.")

def main(config_path):
    config = load_config(config_path)

    workdir = config['paths']['workdir']
    output_h5ad_dir = config['paths']['output_h5ad_dir']
    glue_output_dir = os.path.join(workdir, 'glue')
    os.makedirs(glue_output_dir, exist_ok=True)
    os.makedirs(os.path.join(glue_output_dir, 'figures'), exist_ok=True)
    os.chdir(glue_output_dir)
    print(f"[INFO] Output path set to: {glue_output_dir}")

    tissue = config['params']['tissue']
    
    celltype_obs = config['params']['celltype_obs']
    gtf = config['references']['gencode_gtf']
    n_features = config['params'].get('n_features', 50000)

    if tissue == "---":
        sample_metadata = config['paths']['sample_metadata']
        df = pd.read_csv(sample_metadata, sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    
        tissues = sorted(df["tissue"].unique())
        print(f"[INFO] Annotating ATAC cells for MULTIPLE tissues: {tissues}")

        suffix = config['params']['suffix']
        rna_path = os.path.join(output_h5ad_dir, f'RNA_removeDoublet.{suffix}.h5ad')
        rna = ad.read_h5ad(rna_path)
        atac_path = os.path.join(output_h5ad_dir, f'ATAC_removeDoublet.{suffix}.h5ad')
        atac = ad.read_h5ad(atac_path)
    
        for idx, tissue_name in enumerate(tissues, 1):
            print(f"\n========== Processing tissue: {tissue_name} ({idx}/{len(tissues)}) ==========")
    
            try:
                tissue_std = standardize_tissue_name(tissue_name)
                
                rna_sel = rna[rna.obs['tissue']==tissue_std,:].copy()
                atac_sel = atac[atac.obs['tissue']==tissue_std,:].copy()

                run_per_tissue(tissue_std, rna_sel, atac_sel, output_h5ad_dir, glue_output_dir, gtf, 
                               celltype_obs, n_features=n_features)
            except Exception as e:
                print(f"[ERROR] Encountered error for tissue {tissue_name}: {str(e)}")
    else:
        print(f"\n========== Processing tissue: {tissue} ==========")
        
        try:
            tissue_std = standardize_tissue_name(tissue)

            rna_path = os.path.join(output_h5ad_dir, f'RNA_removeDoublet.{tissue_std}.h5ad')
            rna = ad.read_h5ad(rna_path)
            atac_path = os.path.join(output_h5ad_dir, f'ATAC_removeDoublet.{tissue_std}.h5ad')
            atac = ad.read_h5ad(atac_path)
            
            run_per_tissue(
                tissue_std, rna, atac, output_h5ad_dir, glue_output_dir, gtf, celltype_obs, n_features=n_features
            )
        except Exception as e:
            print(f"[ERROR] Encountered error for tissue {tissue}: {str(e)}")


if __name__ == "__main__":
    import torch
    print(f"[INFO] torch.cuda.is_available(): {torch.cuda.is_available()}")
    print(f"[INFO] torch.cuda.device_count(): {torch.cuda.device_count()}")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"[INFO] Using device: {device}\n")

    parser = argparse.ArgumentParser(
        description="Annotate ATAC cells via GLUE integrating scRNA and scATAC data."
    )
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    main(args.config)
