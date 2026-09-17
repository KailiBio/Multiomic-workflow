# RNA_QC

QC, filtering, and clustering for snRNA-seq data (CellRanger output), per
tissue. First step in the overall workflow — its output feeds
[`cell_annotation`](../cell_annotation/) and, indirectly (via doublet
results), [`multiome_integration`](../multiome_integration/).

## Inputs

- **CellRanger output**, one directory per sample under `paths.input_dir`:
  ```
  <input_dir>/<rnaID>/
    ├── filtered_feature_bc_matrix.h5
    └── metrics_summary.csv
  ```
- **`scrinvex` intronic-read tables** (optional), one per sample under
  `paths.scrinvex_dir`, named `<rnaID>/<rnaID>.scrinvex.tsv`. Skip this
  directory entirely if you don't need exonic-read-ratio QC.
- **Sample metadata table** at `paths.sample_metadata`: tab-separated, no
  header, columns `rnaID, atacID, species, donorID, ageGroup, gender, tissue`
  (see `test_data/test_metatable.txt`).
- **QC cutoff table** at `qc.rna_qc_cutoff_table` (Excel or TSV, sheet given
  by `qc.sheet_name`): per-sample thresholds used for pre-QC plot cutoff
  lines and, in `filter_qc.py`, actual cell filtering (see
  `test_data/test_QC_cutoffs.xlsx` for the expected columns, e.g.
  `Min_genes_in_cells`, `Max_genes_in_cells`, `Max_counts_in_cells`,
  `Max_percent_mt_in_cells`, `Max_percent_ribo_in_cells`).

Config: `config/rna_qc_config.yaml`. Paths are relative to `RNA_QC/`; set
`params.tissue` to `"---"` to process every tissue in the metadata table, or
a single tissue name (matching the `tissue` column) to run just one.

## Pipeline Steps

Run in order (mirrors `notebooks/run_rna_qc_workflow.sh`):

| # | Script | What it does |
|---|---|---|
| 0 | `scripts/check_cellranger_qc_summary.py` | (Optional) Summarize CellRanger's own `metrics_summary.csv` across samples. |
| 1 | `scripts/concat_h5ad.py` | Load each sample's CellRanger matrix, concatenate into one raw `.h5ad` per tissue. |
| 2 | `scripts/pre_qc_assessment.py` | Compute QC metrics (mito/ribo/hb %, MALAT1 CPM, exon ratio), run doublet detection, generate pre-filtering diagnostic plots. |
| 3 | `scripts/filter_qc.py` | Apply cutoff-table-based cell filtering, normalize, HVG selection, PCA/clustering. |
| 4 | `scripts/batch_correction.py` | (Optional) Harmony batch correction across samples/donors. |
| 5 | `scripts/clustering_and_qc_reassessment.py` | Re-cluster post-filtering and re-check QC metrics by cluster. |

`scripts/rnaqc_presentation.py` is a standalone extra that regenerates a
polished set of presentation-ready QC figures from already-processed output;
it's not part of the sequential pipeline above.

## Running

```bash
cd RNA_QC
bash notebooks/run_rna_qc_workflow.sh [runtag]   # runtag defaults to demo_run
```
or open `notebooks/RNA_QC_demo.ipynb` for the same steps interactively.

Each script can also be run standalone, e.g.:
```bash
python scripts/pre_qc_assessment.py config/rna_qc_config.yaml
python scripts/filter_qc.py config/rna_qc_config.yaml demo_run
```

## Outputs

Under `paths.workdir` (default `test_output/`):
- `rna_h5ad/` (`paths.output_h5ad_dir`): raw, filtered, and processed `.h5ad`
  per tissue, at each pipeline stage.
- `pre_qc_assessment/`, `filter_qc.<runtag>/`,
  `clustering_and_qc_reassessment.<runtag>/`, `batch_correction.<runtag>/`:
  per-step figures and summary tables.
