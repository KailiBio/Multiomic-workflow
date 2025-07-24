#!/usr/bin/env python3

"""
Author: Kaili Fan
Description: Check overlapping cells between scRNA and scATAC for all samples in a tissue and generate Venn/bar plots.
"""

import os
import sys
import argparse
import pandas as pd
import anndata as ad
import csv

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from multiome_integration.utils import load_config, standardize_tissue_name
from multiome_integration.integration_plots import plot_venn_and_save, plot_overlap_bar

def remove_remaining_doublet(adata, doublet_dir, id_convert_dic, assay):
    doublet_barcodes_to_remove = set()
    sample_list = adata.obs['sampleID'].unique().tolist()
    
    for sample in sample_list:
        
        if sample in id_convert_dic:
            
            sample_edit = id_convert_dic[sample]
            tissue = adata[adata.obs['sampleID']==sample, :].obs['tissue'].unique()[0]
            print(f"[INFO] checking sample ${sample}: ${tissue}")
            
            filename = f"{assay}_doublet_results.{sample_edit}.tsv"
            doublet_file_path = os.path.join(doublet_dir, filename)
            if not os.path.exists(doublet_file_path):
                print(f"[WARNING] Doublet file not found for sample {sample}: {filename}")
                continue
            doublet_file = pd.read_csv(doublet_file_path, sep='\t', header=0)
            
            if 'cell_barcode' not in doublet_file.columns:
                raise ValueError(f"'cell_barcode' column not found in {doublet_file_path}")
            barcodes = doublet_file.loc[doublet_file['doublet_call'] == 'yes', 'cell_barcode']
            barcodes = list(barcodes)
    
            # edit barcodes for CZI
            if assay == "ATAC":
                barcodes_update = [f"{parts[0]}_{parts[1]}_{parts[3]}-{tissue}" 
                            if len(parts := n.split('_')) >= 4 else n for n in barcodes]
            if assay == "RNA":
                barcodes_update = [f'{b}-{tissue}' for b in barcodes]

            # special for skin samples
            if sample in ['ENC-1LGRB-310-SM-LZDVO-EXP01-GEX-01', 'ENC-1LGRB-310-SM-LZDVO-EXP01-ATAC-01', 'ENC-1K2DA-219-SM-LZDQI-EXP02-GEX-01', 'ENC-1K2DA-219-SM-LZDQI-EXP02-ATAC-01']:
                barcodes_update = ["Sun_Exposed_"+bc for bc in barcodes]
            elif sample in ['ENC-1LGRB-129-SM-LZDVH-EXP01-GEX-Merge', 'ENC-1LGRB-129-SM-LZDVH-EXP01-ATAC-Merge', 'ENC-1K2DA-120-SM-LZDV6-EXP01-GEX-Merge', 'ENC-1K2DA-120-SM-LZDV6-EXP01-ATAC-Merge']:
                barcodes_update = ["Not_Sun_Exposed_"+bc for bc in barcodes]
                                
            doublet_barcodes_to_remove.update(barcodes_update)

        else:
            print(f"[Warning] sample {sample} not on metatable")

    print(f"Removing {len(doublet_barcodes_to_remove)} remaining doublet cells for {assay}.")
    mask = ~adata.obs_names.isin(doublet_barcodes_to_remove)
    adata_filtered = adata[mask].copy()
    print(f"adata reduced from {adata.shape[0]} to {adata_filtered.shape[0]} cells.")
    return adata_filtered

def process_overlap(rna_sample, atac_sample, rna_df, atac_df, tissue_std, figures_dir):
    rna_cells = set(rna_df.loc[rna_df['sampleID']==rna_sample, 'cellbarcode'])
    atac_cells = set(atac_df.loc[atac_df['sampleID']==atac_sample, 'cellbarcode'])
    shared_cells = rna_cells & atac_cells
    rna_only = len(rna_cells - shared_cells)
    overlap = len(shared_cells)
    atac_only = len(atac_cells - shared_cells)
    plot_venn_and_save(list(rna_cells), list(atac_cells), rna_sample, atac_sample, figures_dir, tissue_std)
    return [rna_only, overlap, atac_only]

def main(config_path):
    config = load_config(config_path)

    output_h5ad_dir = config['paths']['output_h5ad_dir']
    figures_dir = os.path.join(config['paths']['workdir'], 'overlap_figures')
    os.makedirs(output_h5ad_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    df = pd.read_csv(config['paths']['sample_metadata'], sep='\t', header=None,
                     names=["rnaID", "atacID", "species", "donorID", "ageGroup", "gender", "tissue"])
    rna_to_atac = dict(zip(df["rnaID"], df["atacID"]))
    atac_to_rna = dict(zip(df["atacID"], df["rnaID"]))

    tissue = config['params']['tissue']
    suffix = config['params']['suffix']
    if tissue == "---":
        working_df = df
    else:
        working_df = df[df["tissue"] == tissue]

    rna_sample_list = working_df["rnaID"].unique().tolist()
    atac_sample_list = working_df["atacID"].unique().tolist()
    sample_tissue_dict = dict(zip(working_df['atacID'], working_df['tissue']))
    print(f"RNA Samples: {list(rna_sample_list)}")
    print(f"ATAC Samples: {list(atac_sample_list)}")
    tissues = sorted(working_df["tissue"].unique())
    print(f"Working tissue: {', '.join(tissues)}")

    # remove remaining doublets and save final h5ad
    rna_doublet_dir = config['paths']['rna_doublet_dir']
    atac_doublet_dir = config['paths']['atac_doublet_dir']

    #rna_all = ad.read_h5ad(config['paths']['rna_post_qc_h5ad'])
    #atac_all = ad.read_h5ad(config['paths']['atac_post_qc_h5ad'])

    # Remove potential doublets in RNA (based on ATAC doublet directory)
    #if os.path.exists(atac_doublet_dir):
     #   print(f"[INFO] ATAC doublet directory found: {atac_doublet_dir}")
      #  rna = remove_remaining_doublet(rna_all, atac_doublet_dir, rna_to_atac, assay="ATAC")
      #  rna.write(os.path.join(output_h5ad_dir, f'RNA_removeDoublet.{suffix}.h5ad'))
    #else:
    #    print(f"[WARNING] ATAC doublet directory not found: {atac_doublet_dir}. Skipping RNA doublet removal.")
    #    rna = rna_all
    #    rna.write(os.path.join(output_h5ad_dir, f'RNA_removeDoublet.{suffix}.h5ad'))

    # Remove potential doublets in ATAC (based on RNA doublet directory)
    #if os.path.exists(rna_doublet_dir):
    #    print(f"[INFO] RNA doublet directory found: {rna_doublet_dir}")
    #    atac = remove_remaining_doublet(atac_all, rna_doublet_dir, atac_to_rna, assay="RNA")
    #    atac.write(os.path.join(output_h5ad_dir, f'ATAC_removeDoublet.{suffix}.h5ad'))
    #else:
    #    print(f"[WARNING] RNA doublet directory not found: {rna_doublet_dir}. Skipping ATAC doublet removal.")
    #    atac = atac_all
    #    atac.write(os.path.join(output_h5ad_dir, f'ATAC_removeDoublet.{suffix}.h5ad'))

    rna = ad.read_h5ad("/data/fankaili/project/CZI/integration/h5ad/RNA_removeDoublet.CZI.h5ad")
    atac = ad.read_h5ad("/data/fankaili/project/CZI/integration/h5ad/ATAC_removeDoublet.CZI.h5ad")

    # For RNA
    rna_df = rna.obs[['sampleID']].copy()
    rna_df['cellbarcode'] = rna.obs_names.values
    # For ATAC
    atac_df = atac.obs[['sampleID']].copy()
    atac_df['cellbarcode'] = atac.obs_names.values
    #
    cell_counts = []
    for idx, row in enumerate(working_df.itertuples(index=False), 1):
        rna_sample = row.rnaID
        atac_sample = row.atacID
        print(f"[{idx}/{len(working_df)}] Processing sample: {rna_sample} / {atac_sample} ...")
        tissue_std = standardize_tissue_name(sample_tissue_dict[atac_sample])
        
        res = process_overlap(rna_sample, atac_sample, rna_df, atac_df, tissue_std, figures_dir)
        if res:
            cell_counts.append({
                "rna_sample": rna_sample,
                "atac_sample": atac_sample,
                "rna_only": res[0],
                "overlap": res[1],
                "atac_only": res[2]
            })
    
    if cell_counts:
        df_cell_counts = pd.DataFrame(cell_counts)
        df_cell_counts.to_csv(
            os.path.join(config['paths']['workdir'], f'multiome_cell_counts.{suffix}.tsv'),
            sep='\t', index=False
        )
        print("Cell counts saved to tsv")
        plot_overlap_bar(cell_counts, figures_dir, suffix)
        print(f"Overlap statistics complete. Figures saved in {figures_dir}")
    else:
        print("No overlaps found or no valid samples to plot.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check overlapping cells between scRNA and scATAC for all samples in a tissue and plot Venn/bar graphs."
    )
    parser.add_argument("config", help="YAML config file")
    args = parser.parse_args()
    main(args.config)