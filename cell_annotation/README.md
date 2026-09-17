# cell_annotation

Marker-gene-based DEG/dotplot analysis and automatic cell-type annotation
(via [ScType](https://github.com/kris-nader/sc-type-py)) on filtered,
clustered snRNA-seq data. Consumes [`RNA_QC`](../RNA_QC/) output; its output
feeds [`multiome_integration`](../multiome_integration/).

## Inputs

- **RNA_QC output `.h5ad`** at `paths.input_h5ad_dir` — the filtered,
  processed, clustered file produced by `RNA_QC/scripts/filter_qc.py` /
  `clustering_and_qc_reassessment.py` (e.g.
  `<tissue>_GEX.filtered.processed.<runtag>.h5ad`).
- **Marker gene table** at `annotation.marker_gene_file` (Excel or TSV, sheet
  given by `annotation.sheet_name`): two columns, `CellType` and `Gene` (see
  `test_data/test_marker_gene.xlsx`) — one row per marker gene, grouped by
  cell type.
- **Sample metadata table** at `paths.sample_metadata` (shared with
  `RNA_QC`/`ATAC_QC`) — only needed when `params.tissue` is `"---"`, to
  enumerate tissues.

Config: `config/cell_annotation_config.yaml`. Paths are relative to
`cell_annotation/`; paths starting with `../` reach into `RNA_QC/`. Set
`params.tissue` to `"---"` to process every tissue, or a single tissue name
to run just one.

## Pipeline Steps

There's no shell runner for this module (run via notebook or the two
scripts directly, in order):

| # | Script | What it does |
|---|---|---|
| 1 | `scripts/annotate_cells.py` | Rank DEGs per Leiden cluster, dotplots against the marker gene table, run ScType automatic annotation, save an annotated `.h5ad` + annotation table. |
| 2 | `scripts/generate_final_h5ad.py` | Merge the ScType annotation back onto the full RNA object and write the final per-tissue annotated `.h5ad`. |

## Running

```bash
cd cell_annotation
python scripts/annotate_cells.py config/cell_annotation_config.yaml
python scripts/generate_final_h5ad.py config/cell_annotation_config.yaml
```
or open `notebooks/cell_annotation_demo.ipynb` for the same steps
interactively.

## Outputs

Under `paths.workdir` (default `test_output/`):
- `cell_annotation_auto/`: DEG/marker dotplots, per-cluster top-DEG UMAPs,
  ScType annotation UMAP, `<tissue>_GEX.autoAnnotation.tsv`, and the
  intermediate annotated `.h5ad`.
- `final/`: the final annotated `.h5ad` (consumed by
  `multiome_integration`'s `paths.rna_post_qc_h5ad`) and summary UMAPs/cell
  counts.
