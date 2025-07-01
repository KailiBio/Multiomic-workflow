#!/usr/bin/env python3

"""
Author: Kaili Fan
Annotate ATAC cells via GLUE integrating scRNA and scATAC data.
"""

import os
import sys
import argparse
import time
from itertools import chain
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import scglue
import snapatac2 as snap

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name
from multiome_integration.integration_plots import plot_multiome_celltype_count_bar

def prepare_rna(rna, gtf, output_h5ad_dir, suffix, celltype_col):
    print(f"[INFO] Preparing RNA AnnData for '{suffix}' using cell type column '{celltype_col}' ...")
    rna.obs['celltype_glue'] = rna.obs[celltype_col]
    print(f"[INFO] Computing highly variable genes (HVGs) ...")
    sc.pp.highly_variable_genes(rna, flavor='seurat')
    print(rna.var['highly_variable'].value_counts())
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
    print(f"[INFO] Retaining {num_good} genes with usable positions.")

    rna = rna[:, good].copy()
    rna.var['chromStart'] = rna.var['chromStart'].astype(int)
    rna.var['chromEnd'] = rna.var['chromEnd'].astype(int)
    if "artif_dupl" in rna.var:
        print(f"[INFO] Dropping artif_dupl column from var ...")
        rna.var.drop(columns=['artif_dupl'], inplace=True)
    out_path = os.path.join(output_h5ad_dir, f"RNA.beforeGLUE.{suffix}.h5ad")
    rna.write(out_path, compression="gzip")
    print(f"[INFO] Wrote prepared RNA AnnData: {out_path}")
    return rna

def prepare_atac(atac_all, gtf, output_h5ad_dir, suffix, n_features=50000):
    print(f"[INFO] Preparing ATAC AnnData for '{suffix}' ...")
    print(f"[INFO] Selecting {n_features} top features using snapatac2 ...")
    snap.pp.select_features(atac_all, n_features=n_features)
    atac = atac_all[:, atac_all.var['selected']].copy()

    print(f"[INFO] Calculating LSI for ATAC ...")
    scglue.data.lsi(atac, n_components=100, n_iter=15)
    print(f"[INFO] Computing neighbors, UMAP, and Leiden clustering for ATAC ...")
    sc.pp.neighbors(atac, use_rep="X_lsi", metric="cosine")
    sc.tl.umap(atac)
    sc.tl.leiden(atac, resolution=0.5, flavor="igraph")
    sc.pl.umap(atac, color="leiden", save=f'_ATAC.beforeGLUE.{suffix}.png')

    print(f"[INFO] Assigning peak coordinates from var_names ...")
    split = atac.var_names.str.split(r"[:-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)

    out_path = os.path.join(output_h5ad_dir, f"ATAC.beforeGLUE.{suffix}.h5ad")
    atac.write(out_path, compression="gzip")
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
        fit_kws={"directory": glue_output_dir},  # Directory for checkpoints/logs
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

def annotate_and_merge(rna, atac, glue, output_h5ad_dir, tissue_std):
    print(f"[INFO] Annotating and merging RNA/ATAC for '{tissue_std}' ...")
    # RNA
    print("  [RNA] Encoding data ...")
    rna.obsm["X_glue"] = glue.encode_data("rna", rna)
    rna.obs['modality'] = 'expression'
    sc.pp.neighbors(rna, use_rep="X_glue", metric="cosine")
    sc.tl.umap(rna)
    print("  [RNA] Saving UMAP and AnnData ...")
    sc.pl.umap(rna, color=["celltype_glue"], save=f"_GLUEembed.RNA.{tissue_std}.png")
    rna.write(os.path.join(output_h5ad_dir, f'RNA.GLUE.{tissue_std}.h5ad'), compression="gzip")

    # ATAC
    print("  [ATAC] Encoding data, running neighbors/UMAP/leiden ...")
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
    atac.write(os.path.join(output_h5ad_dir, f'ATAC.GLUE.{tissue_std}.h5ad'), compression="gzip")

    # Merged AnnData
    print("  [MERGE] Concatenating and writing AnnData, embedding UMAP ...")
    adata_glue = ad.concat([rna, atac])
    sc.pp.neighbors(adata_glue, use_rep="X_glue", metric="cosine")
    sc.tl.umap(adata_glue)
    sc.pl.umap(adata_glue, color=["celltype_glue", "modality"], wspace=0.95, ncols=1, 
               save=f"_GLUEembed.merged.{tissue_std}.png")
    adata_glue.write(os.path.join(output_h5ad_dir, f'Multiome_merged.GLUE.{tissue_std}.h5ad'), compression="gzip")

    print(f"[INFO] Annotation and merging completed for '{tissue_std}'.")

def run_per_tissue(
    tissue,
    output_h5ad_dir,
    glue_output_dir,
    gtf,
    celltype_col,
    n_features=50000,
):
    tissue_std = standardize_tissue_name(tissue)
    print(f"\n============== [GLUE] Processing tissue: {tissue} ({tissue_std}) ==============")
    # Load filtered AnnData for the current tissue
    print(f"[INFO] Loading RNA and ATAC for tissue: {tissue_std}")
    rna_path = os.path.join(output_h5ad_dir, f'RNA_removeDoublet.{tissue_std}.h5ad')
    atac_path = os.path.join(output_h5ad_dir, f'ATAC_removeDoublet.{tissue_std}.h5ad')
    rna = ad.read_h5ad(rna_path)
    atac = ad.read_h5ad(atac_path)

    # Prepare
    rna_prepared = prepare_rna(rna, gtf, output_h5ad_dir, tissue_std, celltype_col)
    atac_prepared = prepare_atac(atac, gtf, output_h5ad_dir, tissue_std, n_features=n_features)

    # Train GLUE & annotate
    glue = train_glue(rna_prepared, atac_prepared, glue_output_dir, tissue_std)
    annotate_and_merge(rna_prepared, atac_prepared, glue, output_h5ad_dir, tissue_std)
    plot_multiome_celltype_count_bar(rna_prepared, atac_prepared, os.path.join(glue_output_dir, 'figures'), tissue_std)
    print(f"[INFO] Finished GLUE for tissue '{tissue}'.")

def main(config_path):
    config = load_config(config_path)

    output_h5ad_dir = config['paths']['output_h5ad_dir']
    glue_output_dir = os.path.join(config['paths']['workdir'], 'glue')
    os.makedirs(glue_output_dir, exist_ok=True)
    os.makedirs(os.path.join(glue_output_dir, 'figures'), exist_ok=True)
    os.chdir(glue_output_dir)
    print(f"[INFO] Output path set to: {glue_output_dir}")

    tissue = config['params']['tissue']
    celltype_col = config['params']['celltype_col']
    gtf = config['references']['gencode_gtf']
    n_features = config['params'].get('n_features', 50000)

    if tissue == "---":
        working_df = pd.read_csv(config['paths']['sample_metadata'], sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
        tissues = sorted(working_df["tissue"].unique())
        print(f"[INFO] Running analysis for MULTIPLE tissues: {tissues}")
        for idx, working_tissue in enumerate(tissues, 1):
            print(f"\n============== [Main] Processing tissue: {working_tissue} ({idx}/{len(tissues)}) ==============")
            try:
                run_per_tissue(
                    working_tissue, output_h5ad_dir, glue_output_dir, gtf, celltype_col, n_features=n_features
                )
            except Exception as e:
                print(f"[ERROR] Encountered error for tissue {working_tissue}: {str(e)}")
    else:
        print(f"\n============== [Main] Processing tissue: {tissue} ==============")
        try:
            run_per_tissue(
                tissue, output_h5ad_dir, glue_output_dir, gtf, celltype_col, n_features=n_features
            )
        except Exception as e:
            print(f"[ERROR] Encountered error for tissue {tissue}: {str(e)}")


if __name__ == "__main__":
    import torch
    print(f"[INFO] torch.cuda.is_available(): {torch.cuda.is_available()}")
    print(f"[INFO] torch.cuda.device_count(): {torch.cuda.device_count()}")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"[INFO] Using device: {device}")

    parser = argparse.ArgumentParser(
        description="Annotate ATAC cells via GLUE integrating scRNA and scATAC data."
    )
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    main(args.config)