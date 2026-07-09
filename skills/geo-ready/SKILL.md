---
name: geo-ready
description: End-to-end GEO single-cell data download and preprocessing pipeline. Use when the user wants to download a GEO dataset (by GSE ID) and preprocess it into an analysis-ready h5ad file. Handles downloading RAW matrix and SOFT metadata from GEO FTP, parsing metadata into CSV, filtering non-human samples, extracting and merging matrix files into h5ad, and registering the dataset in an inventory TSV. Calls geo-download, parse-geo-soft, and geo-matrix-merge skills internally.
argument-hint: <GSE_ID> -d <download_base_dir>
---

# GEO-ready Pipeline

Download a GEO dataset and preprocess it into an analysis-ready h5ad file with metadata.

## Quick Start (Full Pipeline)

Run the complete pipeline in one command:

```bash
python3 ~/.qoder-cn/skills/geo-ready/scripts/geo_ready.py <GSE_ID> -d <download_base_dir>
```

Example:
```bash
python3 ~/.qoder-cn/skills/geo-ready/scripts/geo_ready.py GSE213835 -d /public/home/jihongyi/download_data
```

Options:
- `--skip-download` — skip download, use existing files in the dataset directory
- `--skip-merge` — skip matrix merge step (useful if you only need metadata)

## What It Does (6 Steps)

### Step 1: Download
- Checks if `<base_dir>/<GSE_ID>/` already contains the data files
- If not, downloads from GEO FTP:
  - `<GSE_ID>_family.soft.gz` — metadata in SOFT format
  - `<GSE_ID>_RAW.tar` — raw expression matrix files
- Skips download if files already exist (idempotent)

### Step 2: Parse Metadata
- Decompresses `.soft.gz` to `.soft`
- Calls the `parse-geo-soft` skill to extract structured metadata
- Output: `<GSE_ID>_metadata.csv` with columns:
  - series_id, study_title, platform_id, platform_title
  - Sample_geo_accession, Sample_title
  - tissue, genotype, age (dynamic characteristics)
  - description, source_name, organism, molecule
  - extract_protocol, data_processing, instrument_model
  - library_selection, library_source, library_strategy

### Step 3: Extract Matrix Files
- Extracts `<GSE_ID>_RAW.tar` into `matrix_files/` subdirectory
- Uses a marker file to avoid re-extraction

### Step 3b: Filter by Organism
- Reads the metadata CSV to check the `organism` column for each sample (GSM)
- Matches GSM IDs in matrix filenames against metadata
- Non-human samples (not Homo sapiens) are moved to `other_samples/` directory
- Human samples remain in `matrix_files/` for downstream merging
- If all samples are human, no files are moved

### Step 4: Merge Matrices
- Calls the `geo-matrix-merge` skill to merge all matrix files
- Auto-detects data format: standard matrix or counts+cellname paired format
- For counts+cellname format: reads cellname files to replace generic column names (C1, C2...) with real cell barcodes
- Output: `output/<GSE_ID>_merged.h5ad`
- Cells x genes format, with `adata.obs['sample']` tracking sample origin

### Step 5: Update Inventory
- Reads h5ad to get cell count, gene count, sample count
- Reads metadata CSV for tissue/disease description
- Appends (or updates) a row in `<base_dir>/dataset_inventory.tsv`:
  ```
  dataset_id  data_type  disease_tissue_coverage  n_samples_or_cells  metadata_status  current_status
  GSE213835   single_cell  GBM organoid...  70.6K cells (8 samples)  collected  merged_h5ad
  ```

## Directory Structure

After running, the dataset directory looks like:
```
<base_dir>/
├── dataset_inventory.tsv
└── <GSE_ID>/
    ├── <GSE_ID>_family.soft.gz
    ├── <GSE_ID>_family.soft
    ├── <GSE_ID>_RAW.tar
    ├── <GSE_ID>_metadata.csv
    ├── matrix_files/          # human samples only
    │   ├── GSM*.counts.tsv.gz
    │   ├── GSM*.cellname.list.txt.gz
    │   └── .extracted
    ├── other_samples/         # non-human samples (if any)
    │   └── GSM*.counts.tsv.gz
    └── output/
        └── <GSE_ID>_merged.h5ad
```

## Step-by-Step (Manual / Using Individual Skills)

If you prefer to run steps individually:

### 1. Download with geo-download skill
```bash
python3 ~/.qoder-cn/skills/geo-download/scripts/geo_download.py soft <GSE_ID> -o <base_dir>/<GSE_ID>
curl -L -o <base_dir>/<GSE_ID>/<GSE_ID>_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/<GSEnnn>/<GSE_ID>/suppl/<GSE_ID>_RAW.tar"
```

### 2. Parse metadata with parse-geo-soft skill
```bash
gunzip -k <base_dir>/<GSE_ID>/<GSE_ID>_family.soft.gz
python3 ~/.qoder-cn/skills/parse-geo-soft/scripts/parse_soft.py \
  <base_dir>/<GSE_ID>/<GSE_ID>_family.soft \
  -o <base_dir>/<GSE_ID>/<GSE_ID>_metadata.csv
```

### 3. Extract RAW tar
```bash
mkdir -p <base_dir>/<GSE_ID>/matrix_files
tar xf <base_dir>/<GSE_ID>/<GSE_ID>_RAW.tar -C <base_dir>/<GSE_ID>/matrix_files
```

### 3b. Filter non-human samples
Check the metadata CSV organism column and move non-human files to other_samples/.

### 4. Merge matrices with geo-matrix-merge skill
```bash
python3 ~/.qoder-cn/skills/geo-matrix-merge/scripts/geo_matrix_merge.py \
  <base_dir>/<GSE_ID>/matrix_files \
  -o <base_dir>/<GSE_ID>/output/<GSE_ID>_merged.h5ad
```

### 5. Update inventory
Read the h5ad and metadata CSV, then append a row to `dataset_inventory.tsv`.

## Dependencies

- Python packages: `pandas`, `anndata`, `numpy`
- System tools: `curl`, `tar`, `gunzip`
- Other skills: `parse-geo-soft`, `geo-matrix-merge`

## Notes

- The pipeline is idempotent — re-running skips steps that are already done
- Large datasets may take significant time for download and merge steps
- The `--skip-download` flag is useful when files were downloaded manually
- The inventory TSV uses tab separation to match the existing format
- Non-human samples are preserved in `other_samples/` for reference
