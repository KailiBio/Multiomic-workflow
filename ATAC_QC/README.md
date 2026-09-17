# ATAC_QC

QC, filtering, and clustering for snATAC-seq data (fragment files), per
tissue. Runs independently of `RNA_QC`; its output feeds
[`multiome_integration`](../multiome_integration/) alongside `RNA_QC` and
`cell_annotation`.

## Inputs

- **Fragment files**, one directory per sample under `paths.fragment_dir`:
  ```
  <fragment_dir>/<atacID>/fragments.tsv.gz
  ```
- **Sample metadata table** at `paths.sample_metadata`: tab-separated, no
  header, columns `rnaID, atacID, species, donorID, ageGroup, gender, tissue`
  (see `test_data/test_metatable.txt`) — shared with `RNA_QC`.
- **QC cutoff table** at `qc.atac_qc_cutoff_table` (Excel or TSV, sheet given
  by `qc.sheet_name`): per-sample `num_fragment` / `TSS_enrichment_score`
  thresholds (see `test_data/test_QC_cutoffs.xlsx`).
- **Reference files** under `references:` — a 10x multiome barcode whitelist
  (`references.barcode_whitelist`, bundled at `ref/`) and a GENCODE GTF
  (`references.gencode_gtf`) for gene annotation; the GTF itself is too
  large to ship (see `.gitignore`) — download the matching GENCODE release
  and place it at that path.

Config: `config/atac_qc_config.yaml`. Paths are relative to `ATAC_QC/`; set
`params.tissue` to `"---"` to process every tissue in the metadata table, or
a single tissue name to run just one.

## Pipeline Steps

Run in order (mirrors `notebooks/run_atac_qc_workflow.sh`):

| # | Script | What it does |
|---|---|---|
| 1 | `scripts/remove_pcr_chimeric_awk.sh` (via the shell runner, in parallel across samples) | Strip PCR-chimeric read pairs from each sample's raw fragment file. |
| 2 | `scripts/load_fragments.py` | Load cleaned fragments with `snapatac2`, compute TSS enrichment, produce raw per-sample `.h5ad`. |
| 3 | `scripts/qc_kde_filter.py` | Filter cells by fragment count / TSS enrichment cutoffs; KDE before/after plots. |
| 4 | `scripts/processing_doublet_detection.py` | Doublet detection on filtered cells. |
| 5 | `scripts/downstream_processing.py` | Doublet filtering, embedding, clustering per sample. |
| 6 | `scripts/concatenate_and_cluster.py` | Merge samples for the tissue, joint UMAP/clustering. |

`scripts/atacqc_presentation.py` is a standalone extra that regenerates a
polished set of presentation-ready QC figures from already-processed output;
it's not part of the sequential pipeline above.

## Running

```bash
cd ATAC_QC
bash notebooks/run_atac_qc_workflow.sh [runtag] [nthread]   # runtag defaults to demo_run, nthread to 4
```
or open `notebooks/ATAC_QC_demo.ipynb` for the same steps interactively.

Each Python script can also be run standalone, e.g.:
```bash
python scripts/load_fragments.py config/atac_qc_config.yaml
python scripts/qc_kde_filter.py config/atac_qc_config.yaml demo_run
```

## Outputs

Under `paths.workdir` (default `test_output/`):
- `atac_h5ad/` (`paths.output_h5ad_dir`): raw, filtered, and processed
  `.h5ad` per sample/tissue, at each pipeline stage.
- `qc_filtering.<runtag>/`, `doublet_detection.<runtag>/`,
  `doublet_filter_processing.<runtag>/`, `postprocessing_check.<runtag>/`:
  per-step figures and summary tables.
- `logs/`: log output from individual steps (when redirected to a file).
