# multiome_integration

Joins the RNA and ATAC modalities per cell barcode, integrates them (GLUE,
MultiVI), identifies cis-regulatory elements, imputes cross-modality signal,
and produces the final annotated joint multiome AnnData objects. The last
step in the overall workflow — consumes [`RNA_QC`](../RNA_QC/),
[`ATAC_QC`](../ATAC_QC/), and [`cell_annotation`](../cell_annotation/)
output.

## Inputs

- **`cell_annotation` output** at `paths.rna_post_qc_h5ad` — the final
  annotated RNA `.h5ad` (`<tissue>.GEX.final.h5ad`).
- **`ATAC_QC` output** at `paths.atac_post_qc_h5ad` — the final clustered
  ATAC `.h5ad` for the tissue.
- **Doublet-detection intermediates** at `paths.rna_doublet_dir`
  (`RNA_QC`'s `filter_qc.<runtag>/`) and `paths.atac_doublet_dir`
  (`ATAC_QC`'s `doublet_filter_processing.<runtag>/`).
- **Sample metadata table** at `paths.sample_metadata` (shared across all
  modules) — only needed when `params.tissue` is `"---"`.
- **Reference files** under `references:` — the same GENCODE GTF used by
  `ATAC_QC` (`references.gencode_gtf`), and a `<chrom><TAB><length>`
  chromosome-sizes file for your reference genome (`references.chrom_size`,
  e.g. `hg38.chrom.sizes` from the UCSC goldenPath). The GTF is not shipped in
  this repo; supply your own. `references.chrom_size` is optional and defaults
  to the built-in hg38 sizes; set it for non-human genomes.

Config: `config/joint_multiomic_config.yaml`. Paths are relative to
`multiome_integration/`; paths starting with `../` reach into sibling module
directories. Key params: `params.celltype_obs` / `params.celllineage_obs`
(the `.obs` column names holding cell-type/cell-lineage annotation, as
written by `cell_annotation`), and `params.imputation` (whether to run the
imputation step).

## Pipeline Steps

Run in order (mirrors `notebooks/run_multiome_integration_workflow.sh`):

| # | Script | What it does |
|---|---|---|
| 1 | `scripts/check_overlapping_cells.py` | Remove remaining doublets, check RNA/ATAC barcode overlap per sample. |
| 2 | `scripts/annotate_ATAC_cells.py` | Transfer RNA-derived cell-type/lineage annotation onto matching ATAC cells. |
| 3 | `scripts/identify_regulatory_elements.py` | Build a consensus peak set per cell type and identify candidate cis-regulatory elements. |
| 4 | `scripts/data_imputation.py` | MultiVI-based cross-modality signal imputation (RNA→ATAC cells and vice versa). |
| 5 | `scripts/data_integration.py` | Finalize the joint AnnData: propagate metadata/layers from both modalities, add imputed signal, split into RNA-only/ATAC-only/joint outputs. |

## Running

```bash
cd multiome_integration
bash notebooks/run_multiome_integration_workflow.sh
```
or open `notebooks/multiome_integration_demo.ipynb` for the same steps
interactively.

Each script can also be run standalone, e.g.:
```bash
python scripts/data_integration.py config/joint_multiomic_config.yaml
```

## Outputs

Under `paths.workdir` (default `test_output/`):
- `h5ad/` (`paths.output_h5ad_dir`): intermediate `.h5ad` at each
  integration stage (overlap-checked, annotated, imputed, GLUE-merged).
- `glue/`, `multivi/`, `atac_peak/`, `overlap_figures/`: per-step figures,
  model checkpoints, and consensus-peak tables.
- `h5ad/` also holds the final joint `multiome_final.<tissue>.h5ad`, plus the
  split `multiome_final.RNA.<tissue>.h5ad` / `multiome_final.ATAC.<tissue>.h5ad`.
- `final/`: summary UMAPs of the final objects.

**Note on `imputation_signal`.** The final object's `imputation_signal` layer
is filled by matching feature names against the MultiVI imputed matrices.
Imputed expression covers all genes. Imputed accessibility is defined on the
consensus peaks (cREs) from step 3, so it only fills ATAC features that use
the same peak names; in the demo the combined object uses genomic tiles, so
those columns stay 0 (the script logs how many features matched, and the
counts are stored in `adata.uns["imputation_signal_features"]`). The complete
imputed matrices are in `MultiVI_impute_RNA/ATAC.<tissue>.h5ad`.
