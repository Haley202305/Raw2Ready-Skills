---
name: geo-matrix-merge
description: Merge GEO (Gene Expression Omnibus) supplementary matrix files from multiple samples into a single AnnData object. Use when the user needs to combine GSM*.matrix.txt.gz files, merge single-cell count matrices downloaded from GEO, or build a unified AnnData (.h5ad) from multiple GEO samples.
---

# GEO Matrix Merge

Merge multiple GEO supplementary matrix files (e.g. `GSM*.filtered.matrix.txt.gz`) into a single AnnData object ready for downstream single-cell analysis.

## When to Use

- User has downloaded GEO supplementary files containing expression matrices
- User wants to combine multiple GSM samples into one `.h5ad` file
- User needs to inspect the structure of a GEO matrix file before merging

## Dependencies

Install before running:

```bash
pip install pandas numpy anndata scipy
```

## Workflow

### Step 1 — Inspect (optional)

Before merging, preview a matrix file to confirm its structure:

```bash
python scripts/geo_matrix_merge.py inspect /path/to/GSM1234567.filtered.matrix.txt.gz
```

Output includes shape, gene/cell name examples, data type, value range, and a head preview.

### Step 2 — Merge

Merge all matrix files in a directory:

```bash
python scripts/geo_matrix_merge.py merge /path/to/samples_dir/ -o merged_output.h5ad
```

If `-o` is omitted, the output defaults to `merged_geo_samples.h5ad` inside the input directory.

## What the Script Does

1. **File discovery** — scans the directory for files matching common matrix patterns (`*.matrix.txt.gz`, `*.mtx.gz`, `*.tsv.gz`, etc.), prioritising GSM-prefixed files.
2. **Format detection** — auto-detects plain text (tab/comma/space-delimited) vs Matrix Market (`.mtx`) format, with or without gzip compression.
3. **Gene alignment** — collects the union of all gene names across samples and fills missing genes with zeros so every sample shares the same gene axis.
4. **Cell barcode deduplication** — prepends a `{GSM}_{sampleName}_` prefix to each cell barcode to guarantee uniqueness across samples.
5. **AnnData assembly** — transposes the combined (gene x cell) matrix into (cell x gene) layout, stores sample identity in `adata.obs['sample']`, and writes a `.h5ad` file.

## Output

The resulting `.h5ad` file contains:

| Component | Content |
|-----------|---------|
| `adata.X` | Merged expression matrix (cells x genes) |
| `adata.obs` | Cell barcodes with a `sample` column |
| `adata.var` | Gene names (union of all samples, deduplicated) |

## Example

```
$ python scripts/geo_matrix_merge.py merge /data/geo_downloads/

Found 3 matrix files: ['GSM3589420_PP019swap.filtered.matrix.txt.gz', ...]
Loaded: GSM3589420_PP019swap.filtered.matrix.txt.gz, shape: (27654, 4521), sample: GSM3589420_PP019swap
...
Merged DataFrame shape: (27654, 12483)

GEO Matrix Merge Complete
==================================================
Input directory: /data/geo_downloads/
Matrix files: 3
Merged shape: 12483 cells x 27654 genes
Sample count: 3
Output file: /data/geo_downloads/merged_geo_samples.h5ad
```

## Pitfalls

- Files without `GSM` in the name are still merged if no GSM-prefixed files are found — verify this matches your intent.
- Matrix Market format requires companion `_genes.txt` / `_barcodes.txt` files; if missing, generic names (`GENE_0`, `CELL_0`) are used.
- Very large datasets (100k+ cells) may need significant RAM; consider sparse matrices for production pipelines.
