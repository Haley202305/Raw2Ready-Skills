---
name: geo-matrix-merge
description: Merge GEO single-cell matrix files from a directory into a single AnnData h5ad. Auto-detects 7 formats — 10X MTX (single/multi), 10X HDF5, 10X CSV bz2, counts+cellname pairs, single GEO text, and multi GEO text. Use when the user wants to combine GEO expression matrices, convert GEO count tables to h5ad, or merge GSM sample files into one scRNA-seq object. Also supports inspecting individual matrix file structure.
argument-hint: <samples_directory> [-o output.h5ad]
---

# GEO Matrix Merge

Merge GEO / 10X single-cell matrix files into a single AnnData (`.h5ad`) object, or inspect a single matrix file. The script auto-detects the data format and dispatches to the correct reader — no manual format selection needed.

## Usage

### Merge mode

```bash
python3 ~/.qoderworkcn/skills/geo-matrix-merge/scripts/geo_matrix_merge.py <samples_directory> [-o output.h5ad]
```

### Inspect mode

```bash
python3 ~/.qoderworkcn/skills/geo-matrix-merge/scripts/geo_matrix_merge.py --inspect <file_path> [-n 10]
```

## Auto-detected Formats (7 total)

The script scans the directory contents and picks the first matching format:

| # | Format | Detection Rule | Handler |
|---|--------|---------------|---------|
| 1 | **10X HDF5** | `*_feature_bc_matrix.h5.gz` or `*_filtered_feature_bc_matrix.h5.gz` | `merge_10x_h5` |
| 2 | **10X CSV bz2** | `*_filtered_feature_count_matrix.csv.bz2` | `merge_10x_csv_bz2` |
| 3 | **10X MTX (single)** | One `matrix.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz` | `merge_10x` |
| 4 | **10X MTX (multi)** | Multiple `GSM*_matrix.mtx.gz` + `*_features.tsv.gz` + `*_barcodes.tsv.gz` triplets | `merge_multi_10x` |
| 5 | **counts+cellname** | `*.counts.*` + `*.cellname.*` paired files | `merge_counts_cellname` |
| 6 | **Single GEO** | One `.txt.gz` / `.tsv.gz` / `.csv.gz` text matrix | `merge_single_geo` |
| 7 | **Multi GEO** | Multiple `GSM*.txt.gz` / `*.filtered.matrix.txt.gz` text matrices | `merge_multi_geo` |

Detection priority: HDF5 → CSV bz2 → MTX → counts+cellname → GEO text.

### Format details

**10X HDF5** (`*.h5.gz`): CellRanger HDF5 output containing `matrix/features` and `matrix/barcodes` groups. Read via `h5py`, parsed as CSC sparse matrix and transposed. Supports multiple h5 files in one directory.

**10X CSV bz2** (`*.csv.bz2`): CSV-compressed count matrix (genes × cells). Read via `bz2.open`, transposed to cells × genes. Each file is one sample; sample name extracted from filename prefix.

**10X MTX (single)**: Standard CellRanger directory with `matrix.mtx.gz`, `features.tsv.gz`, `barcodes.tsv.gz`. Read via `scipy.io.mmread` into sparse CSR format.

**10X MTX (multi)**: Multiple 10X samples in one directory. Files follow the naming pattern `GSM*_sampleName_matrix.mtx.gz`, `GSM*_sampleName_features.tsv.gz`, `GSM*_sampleName_barcodes.tsv.gz`. The script groups files by sample prefix and merges all samples.

**counts+cellname**: GEO-deposited paired files where `*.counts.*` contains the expression matrix (genes × cell indices) and `*.cellname.*` maps cell indices to real barcode names. The script auto-matches pairs by sample prefix.

**Single GEO**: A single tab/comma-separated matrix file (genes × cells, with optional `Gene` symbol column that is extracted and dropped).

**Multi GEO**: Multiple GSM-prefixed matrix files. Each file is one sample. Cell barcodes are prefixed with sample name for uniqueness. Missing genes across samples are filled with 0.

## What the Script Does

1. **Format detection** — scans directory, dispatches to the matching handler.
2. **Gene ID normalization** — strips genome reference prefixes (`hg19_`, `hg38_`, `mm10_`, etc.) for consistent gene naming across samples.
3. **Sparse merge engine** — builds a global gene index, reindexes each sample's sparse matrix, transposes to cells × genes, and `vstack`s all blocks into one CSR matrix.
4. **AnnData assembly** — creates AnnData with sparse `adata.X`, stores sample identity in `adata.obs['sample']`, and gene symbol mapping in `adata.var['gene_symbol']` (when available).
5. **Output** — writes `.h5ad` with deduplicated cell and gene names.

## Output

| Component | Content |
|-----------|---------|
| `adata.X` | Sparse CSR matrix (cells × genes) |
| `adata.obs` | Cell barcodes with `sample` column |
| `adata.obs_names` | Unique cell barcodes (prefixed with sample name) |
| `adata.var` | Gene IDs with optional `gene_symbol` column |

Default output path: `<samples_directory>/merged_geo_samples.h5ad`

## Dependencies

```bash
pip install pandas numpy anndata scipy h5py
```

`h5py` is only required when processing 10X HDF5 files.
