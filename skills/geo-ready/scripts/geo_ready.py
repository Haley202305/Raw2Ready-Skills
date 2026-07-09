#!/usr/bin/env python3
"""GEO-ready: End-to-end GEO data download and preprocessing pipeline.

Usage:
    python geo_ready.py <GSE_ID> -d <download_base_dir> [--skip-download] [--skip-merge]

Steps:
    1. Check if data already exists; if not, download SOFT.gz and RAW.tar from GEO FTP
    2. Parse SOFT file -> metadata CSV
    3. Extract RAW.tar -> matrix_files/
    3b. Filter by organism: move non-human sample files to other_samples/
    4. Merge matrix files -> h5ad (via geo-matrix-merge skill)
    5. Update dataset_inventory.tsv with dataset summary
"""

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import tarfile
from collections import OrderedDict
from pathlib import Path

import anndata as ad
import pandas as pd

SKILLS_DIR = os.path.expanduser("~/.qoder-cn/skills")
PARSE_SOFT_SCRIPT = os.path.join(SKILLS_DIR, "parse-geo-soft/scripts/parse_soft.py")
MATRIX_MERGE_SCRIPT = os.path.join(SKILLS_DIR, "geo-matrix-merge/scripts/geo_matrix_merge.py")

GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def get_series_dir_name(geo_id):
    num = int(re.search(r'\d+', geo_id).group())
    s = str(num)
    if len(s) >= 4:
        return f"GSE{s[:-3]}nnn"
    else:
        return f"GSE{s[:-2]}nn"


def download_file(url, dest_path):
    print(f"  Downloading: {url}")
    print(f"  Saving to:   {dest_path}")
    result = subprocess.run(
        ["curl", "-L", "-o", str(dest_path), "--retry", "3", "--connect-timeout", "30", url],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"Download failed (curl exit {result.returncode}): {result.stderr}")
    size_mb = os.path.getsize(dest_path) / (1024 * 1024)
    print(f"  Downloaded: {size_mb:.2f} MB")
    return dest_path


def step1_download(geo_id, dataset_dir):
    series_dir = get_series_dir_name(geo_id)
    soft_url = f"{GEO_FTP_BASE}/{series_dir}/{geo_id}/soft/{geo_id}_family.soft.gz"
    raw_url = f"{GEO_FTP_BASE}/{series_dir}/{geo_id}/suppl/{geo_id}_RAW.tar"

    soft_gz = dataset_dir / f"{geo_id}_family.soft.gz"
    soft_file = dataset_dir / f"{geo_id}_family.soft"
    raw_tar = dataset_dir / f"{geo_id}_RAW.tar"

    if not soft_gz.exists() and not soft_file.exists():
        print("[Step 1a] Downloading SOFT metadata file...")
        download_file(soft_url, soft_gz)
    elif soft_gz.exists():
        print(f"[Step 1a] SOFT file already exists: {soft_gz}")
    else:
        print(f"[Step 1a] Decompressed SOFT already exists: {soft_file}")

    if not raw_tar.exists():
        print("[Step 1b] Downloading RAW data tar file...")
        download_file(raw_url, raw_tar)
    else:
        print(f"[Step 1b] RAW tar already exists: {raw_tar}")

    return soft_gz, raw_tar


def step2_parse_soft(dataset_dir, geo_id):
    soft_gz = dataset_dir / f"{geo_id}_family.soft.gz"
    soft_file = dataset_dir / f"{geo_id}_family.soft"

    if not soft_file.exists():
        if soft_gz.exists():
            print("[Step 2a] Decompressing SOFT file...")
            with gzip.open(str(soft_gz), 'rb') as f_in:
                with open(str(soft_file), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            raise FileNotFoundError(f"No SOFT file found in {dataset_dir}")
    else:
        print(f"[Step 2a] Decompressed SOFT already exists: {soft_file}")

    metadata_csv = dataset_dir / f"{geo_id}_metadata.csv"
    print(f"[Step 2b] Parsing SOFT -> {metadata_csv}")

    result = subprocess.run(
        [sys.executable, PARSE_SOFT_SCRIPT, str(soft_file), "-o", str(metadata_csv)],
        capture_output=True, text=True
    )
    print(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"parse-geo-soft failed: {result.stderr}")

    return metadata_csv


def step3_extract_raw(raw_tar, dataset_dir):
    matrix_dir = dataset_dir / "matrix_files"
    matrix_dir.mkdir(exist_ok=True)

    marker = matrix_dir / ".extracted"
    if marker.exists():
        n_files = len([f for f in matrix_dir.iterdir() if f.name != ".extracted"])
        print(f"[Step 3] RAW already extracted to {matrix_dir} ({n_files} files)")
        return matrix_dir

    print(f"[Step 3] Extracting RAW tar -> {matrix_dir}")
    with tarfile.open(str(raw_tar), 'r') as tar:
        members = tar.getmembers()
        print(f"  Archive contains {len(members)} files")
        tar.extractall(path=str(matrix_dir))

    extracted_files = [f for f in matrix_dir.iterdir() if f.name != ".extracted"]
    print(f"  Extracted {len(extracted_files)} files")

    marker.touch()
    return matrix_dir


HUMAN_KEYWORDS = ['homo sapiens', 'human', 'homo_sapiens']


def step3b_filter_organism(matrix_dir, metadata_csv):
    """Check organism for each sample and move non-human files to other_samples/."""
    other_dir = matrix_dir.parent / "other_samples"

    if not metadata_csv.exists():
        print(f"[Step 3b] 警告: metadata CSV 不存在，跳过物种过滤")
        return matrix_dir

    df_meta = pd.read_csv(metadata_csv)
    if 'organism' not in df_meta.columns or 'Sample_geo_accession' not in df_meta.columns:
        print(f"[Step 3b] 警告: metadata 缺少 organism 或 Sample_geo_accession 列，跳过物种过滤")
        return matrix_dir

    non_human_gsms = set()
    human_gsms = set()
    for _, row in df_meta.iterrows():
        gsm = str(row.get('Sample_geo_accession', '')).strip()
        organism = str(row.get('organism', '')).strip().lower()
        if not gsm:
            continue
        is_human = any(kw in organism for kw in HUMAN_KEYWORDS)
        if is_human:
            human_gsms.add(gsm)
        else:
            non_human_gsms.add(gsm)

    if not non_human_gsms:
        print(f"[Step 3b] 所有 {len(human_gsms)} 个样本均为人类 (Homo sapiens)，无需过滤")
        return matrix_dir

    print(f"[Step 3b] 物种过滤: {len(human_gsms)} 个人类样本, {len(non_human_gsms)} 个非人类样本")
    print(f"  非人类 GSM: {sorted(non_human_gsms)}")

    other_dir.mkdir(exist_ok=True)
    moved_count = 0
    for f in sorted(matrix_dir.iterdir()):
        if f.name.startswith('.'):
            continue
        for gsm in non_human_gsms:
            if f.name.startswith(gsm):
                dest = other_dir / f.name
                shutil.move(str(f), str(dest))
                moved_count += 1
                break

    print(f"  已移动 {moved_count} 个文件到 {other_dir}")

    remaining = len([f for f in matrix_dir.iterdir() if not f.name.startswith('.')])
    print(f"  matrix_files 剩余 {remaining} 个文件")

    return matrix_dir


def step4_merge_matrix(matrix_dir, geo_id, dataset_dir):
    output_dir = dataset_dir / "output"
    output_dir.mkdir(exist_ok=True)
    h5ad_file = output_dir / f"{geo_id}_merged.h5ad"

    if h5ad_file.exists():
        size_mb = h5ad_file.stat().st_size / (1024 * 1024)
        print(f"[Step 4] h5ad already exists: {h5ad_file} ({size_mb:.1f} MB)")
        return h5ad_file

    print(f"[Step 4] Merging matrix files -> {h5ad_file}")
    result = subprocess.run(
        [sys.executable, MATRIX_MERGE_SCRIPT, str(matrix_dir), "-o", str(h5ad_file)],
        capture_output=True, text=True
    )
    print(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"geo-matrix-merge failed: {result.stderr}")

    return h5ad_file


def step5_update_inventory(geo_id, h5ad_file, metadata_csv, base_dir):
    inventory_path = base_dir / "dataset_inventory.tsv"

    print(f"[Step 5] Updating inventory: {inventory_path}")

    adata = ad.read_h5ad(str(h5ad_file), backed='r')
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    n_samples = adata.obs['sample'].nunique() if 'sample' in adata.obs else 0
    adata.file.close()

    df_meta = pd.read_csv(metadata_csv)
    study_title = df_meta['study_title'].iloc[0] if 'study_title' in df_meta.columns else ""

    tissue_vals = []
    for col in ['tissue', 'source_name', 'genotype']:
        if col in df_meta.columns:
            unique_vals = df_meta[col].dropna().unique()
            tissue_vals.extend([str(v) for v in unique_vals if str(v).strip()])
    disease_tissue = "; ".join(dict.fromkeys(tissue_vals)) if tissue_vals else study_title

    if n_cells >= 1000:
        cells_str = f"{n_cells / 1000:.1f}K cells ({n_samples} samples)"
    else:
        cells_str = f"{n_cells} cells ({n_samples} samples)"

    new_row = {
        "dataset_id": geo_id,
        "data_type": "single_cell",
        "disease_tissue_coverage": disease_tissue,
        "n_samples_or_cells": cells_str,
        "metadata_status": "collected",
        "current_status": "merged_h5ad",
    }

    if inventory_path.exists():
        existing = pd.read_csv(inventory_path, sep='\t')
        if geo_id in existing['dataset_id'].values:
            print(f"  {geo_id} already in inventory, updating row...")
            idx = existing[existing['dataset_id'] == geo_id].index
            for col, val in new_row.items():
                existing.loc[idx, col] = val
            existing.to_csv(inventory_path, sep='\t', index=False)
        else:
            existing = pd.concat([existing, pd.DataFrame([new_row])], ignore_index=True)
            existing.to_csv(inventory_path, sep='\t', index=False)
    else:
        pd.DataFrame([new_row]).to_csv(inventory_path, sep='\t', index=False)

    print(f"  Inventory updated: {geo_id} | {cells_str} | {n_genes} genes")
    return inventory_path


def run_pipeline(geo_id, base_dir, skip_download=False, skip_merge=False):
    geo_id = geo_id.upper()
    if not geo_id.startswith('GSE'):
        print(f"Error: Invalid GEO ID '{geo_id}', must start with GSE", file=sys.stderr)
        sys.exit(1)

    base_dir = Path(base_dir).resolve()
    dataset_dir = base_dir / geo_id
    dataset_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"GEO-ready Pipeline: {geo_id}")
    print(f"Base directory: {base_dir}")
    print(f"Dataset directory: {dataset_dir}")
    print("=" * 60)

    if skip_download:
        soft_gz = dataset_dir / f"{geo_id}_family.soft.gz"
        soft_file = dataset_dir / f"{geo_id}_family.soft"
        raw_tar = dataset_dir / f"{geo_id}_RAW.tar"
        if not soft_gz.exists() and not soft_file.exists():
            print(f"Error: --skip-download but no SOFT file found in {dataset_dir}", file=sys.stderr)
            sys.exit(1)
        if not raw_tar.exists():
            print(f"Error: --skip-download but RAW tar not found: {raw_tar}", file=sys.stderr)
            sys.exit(1)
        print("[Step 1] Skipped (--skip-download)")
    else:
        step1_download(geo_id, dataset_dir)

    metadata_csv = step2_parse_soft(dataset_dir, geo_id)

    matrix_dir = step3_extract_raw(
        dataset_dir / f"{geo_id}_RAW.tar", dataset_dir
    )

    step3b_filter_organism(matrix_dir, metadata_csv)

    if skip_merge:
        h5ad_file = dataset_dir / "output" / f"{geo_id}_merged.h5ad"
        if not h5ad_file.exists():
            print(f"Warning: --skip-merge but h5ad not found: {h5ad_file}", file=sys.stderr)
        print("[Step 4] Skipped (--skip-merge)")
    else:
        h5ad_file = step4_merge_matrix(matrix_dir, geo_id, dataset_dir)

    if h5ad_file.exists():
        step5_update_inventory(geo_id, h5ad_file, metadata_csv, base_dir)
    else:
        print("[Step 5] Skipped (h5ad file not found)")

    print("\n" + "=" * 60)
    print(f"GEO-ready pipeline completed for {geo_id}")
    print(f"  Metadata CSV: {metadata_csv}")
    print(f"  Matrix files: {matrix_dir}")
    print(f"  Merged h5ad:  {h5ad_file}")
    print(f"  Inventory:    {base_dir / 'dataset_inventory.tsv'}")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description="GEO-ready: Download GEO data and preprocess into analysis-ready h5ad"
    )
    parser.add_argument("geo_id", help="GEO Series ID (e.g. GSE213835)")
    parser.add_argument("-d", "--base-dir", required=True,
                        help="Base download directory (e.g. ~/download_data)")
    parser.add_argument("--skip-download", action="store_true",
                        help="Skip download step (use existing files)")
    parser.add_argument("--skip-merge", action="store_true",
                        help="Skip matrix merge step")
    args = parser.parse_args()

    run_pipeline(args.geo_id, args.base_dir, args.skip_download, args.skip_merge)


if __name__ == "__main__":
    main()
