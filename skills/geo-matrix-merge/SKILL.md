---
name: geo-matrix-merge
description: Merge multiple GEO single-cell matrix files (txt.gz, tsv, mtx) from a directory into a single AnnData h5ad file. Use when the user wants to combine GEO expression matrices, convert GEO count tables to h5ad, or merge GSM sample files into one scRNA-seq object. Also supports inspecting individual matrix file structure.
argument-hint: <samples_directory> [-o output.h5ad]
---

# GEO Matrix Merge

Merge GEO matrix files into a single AnnData (h5ad) object, or inspect a single matrix file.

## Usage

### Merge mode

```bash
python3 ~/.qoder-cn/skills/geo-matrix-merge/scripts/geo_matrix_merge.py <samples_directory> [-o output.h5ad]
```

### Inspect mode

```bash
python3 ~/.qoder-cn/skills/geo-matrix-merge/scripts/geo_matrix_merge.py --inspect <file_path> [-n 10]
```

## Auto-detected Formats

The script auto-detects the data format from directory contents and dispatches to the appropriate reader:

| Format | Detection Rule | Handler |
|--------|---------------|---------|
| **10X** | Directory contains one `matrix.mtx.gz` (+ `features.tsv.gz` / `barcodes.tsv.gz`) | `merge_10x` |
| **Multi-10X** | Directory contains multiple `*_matrix.mtx.gz` + `*_features.tsv.gz` + `*_barcodes.tsv.gz` triplets | `merge_multi_10x` |
| **counts+cellname** | Directory contains `*.counts.*` + `*.cellname.*` paired files | `merge_counts_cellname` |
| **Single GEO** | Only one `.txt.gz` / `.tsv.gz` file | `merge_single_geo` |
| **Multi GEO** | Multiple `GSM*.txt.gz` / `*.filtered.matrix.txt.gz` files | `merge_multi_geo` |

### Format details

**10X format**: Standard CellRanger output with `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`. Read via `scipy.io.mmread` with sparse storage.

**Multi-10X format**: Multiple 10X samples in one directory. Files follow the naming pattern `GSM*_sampleName_matrix.mtx.gz`, `GSM*_sampleName_features.tsv.gz`, `GSM*_sampleName_barcodes.tsv.gz`. The script groups files by sample prefix and merges all samples.

**counts+cellname format**: GEO-deposited paired files where `*.counts.*` contains the expression matrix (genes x cell indices) and `*.cellname.*` maps cell indices to real barcode names. The script auto-matches pairs by sample prefix.

**Single GEO**: A single tab/comma-separated matrix file (genes x cells, with optional `Gene` symbol column).

**Multi GEO**: Multiple GSM-prefixed matrix files. Each file is one sample. Cell barcodes are prefixed with sample name for uniqueness. Missing genes across samples are filled with 0.

## Output

- Default output: `<samples_directory>/merged_geo_samples.h5ad`
- Sparse CSR matrix stored in `adata.X` (cells x genes)
- `adata.obs['sample']`: sample identifier per cell
- `adata.var['gene_symbol']`: gene symbol mapping (when available)

## Dependencies

Requires: `pandas`, `numpy`, `anndata`, `scipy`.
