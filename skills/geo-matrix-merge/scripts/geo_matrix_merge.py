#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geo_matrix_merge.py - Merge multiple GEO matrix files into one AnnData object

Commands:
  merge   - Scan a directory for matrix files and merge into .h5ad
  inspect - Preview the structure of a single matrix file

Dependencies: pandas, numpy, anndata, scipy
Usage: python geo_matrix_merge.py {merge,inspect} [options]
"""

import os
import sys
import glob
import gzip
import re
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import anndata as ad


# --------------------------------------------------
# Internal helper functions (内部工具函数)
# --------------------------------------------------

def _read_matrix_file(file_path: str):
    """
    Read a matrix file; supports .gz, .txt, .tsv, .mtx and more.
    (读取矩阵文件，支持多种格式)
    """
    try:
        filename = os.path.basename(file_path)

        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                first_line = f.readline()
                f.seek(0)
                if first_line.startswith('%%MatrixMarket'):
                    return _read_matrix_market_format(file_path)
                else:
                    return pd.read_csv(f, sep=None, engine='python', index_col=0)
        else:
            if filename.endswith(('.mtx', 'matrix')):
                return _read_matrix_market_format(file_path)
            else:
                with open(file_path, 'r') as f:
                    first_line = f.readline()
                if '\t' in first_line:
                    sep = '\t'
                elif ',' in first_line:
                    sep = ','
                else:
                    sep = r'\s+'
                return pd.read_csv(file_path, sep=sep, engine='python', index_col=0)

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None


def _read_matrix_market_format(file_path: str):
    """
    Read Matrix Market format files (.mtx / .mtx.gz).
    (读取 Matrix Market 格式的矩阵文件)
    """
    try:
        import scipy.io
        import tempfile

        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rb') as f_in:
                with tempfile.NamedTemporaryFile(delete=False) as f_out:
                    f_out.write(f_in.read())
                    temp_path = f_out.name
            matrix = scipy.io.mmread(temp_path)
            os.unlink(temp_path)
        else:
            matrix = scipy.io.mmread(file_path)

        dense_matrix = matrix.toarray()

        base_path = file_path.replace('.mtx', '').replace('.gz', '')

        genes = None
        barcodes = None

        # Look for companion gene/feature files (查找基因文件)
        possible_gene_files = [
            base_path + '_genes.txt',
            base_path + '_rows.txt',
            base_path.rsplit('/', 1)[0] + '/genes.tsv',
            base_path.rsplit('/', 1)[0] + '/features.tsv',
        ]
        for gf in possible_gene_files:
            if os.path.exists(gf):
                genes = pd.read_csv(gf, header=None, squeeze=True).values
                break

        # Look for companion barcode files (查找细胞条码文件)
        possible_barcode_files = [
            base_path + '_barcodes.txt',
            base_path.rsplit('/', 1)[0] + '/barcodes.tsv',
        ]
        for bf in possible_barcode_files:
            if os.path.exists(bf):
                barcodes = pd.read_csv(bf, header=None, squeeze=True).values
                break

        if genes is None:
            genes = [f"GENE_{i}" for i in range(dense_matrix.shape[0])]
        if barcodes is None:
            barcodes = [f"CELL_{i}" for i in range(dense_matrix.shape[1])]

        return pd.DataFrame(dense_matrix, index=genes, columns=barcodes)

    except Exception as e:
        print(f"Error reading Matrix Market file {file_path}: {e}")
        return None


# --------------------------------------------------
# Public API (公开功能)
# --------------------------------------------------

def merge_geo_matrix_files(samples_directory: str, output_file: str = None) -> str:
    """
    Merge multiple GEO matrix files in a directory into one AnnData object.
    (将目录中的多个 GEO 矩阵文件合并成一个 AnnData 对象)

    Args:
        samples_directory: Path to directory containing GEO matrix files
        output_file: Output .h5ad path (default: merged_geo_samples.h5ad)

    Returns:
        Summary report string
    """
    # Validate input directory (验证输入目录)
    if not os.path.exists(samples_directory):
        return f"Error: directory does not exist - {samples_directory}"
    if not os.path.isdir(samples_directory):
        return f"Error: path is not a directory - {samples_directory}"

    # Discover matrix files (查找所有符合条件的文件)
    patterns = [
        "*.filtered.matrix.txt.gz",
        "*.filtered.matrix.txt",
        "*matrix.mtx.gz",
        "*matrix.mtx",
        "*.txt.gz",
        "*.txt",
        "*.tsv.gz",
        "*.tsv",
    ]

    matrix_files = []
    for pattern in patterns:
        matrix_files.extend(glob.glob(os.path.join(samples_directory, pattern)))

    # Prefer GSM-prefixed files (优先使用 GSM 格式文件)
    gsm_files = [f for f in matrix_files if 'GSM' in os.path.basename(f)]
    if not gsm_files:
        gsm_files = matrix_files

    if not gsm_files:
        return f"Error: no GEO matrix files found in {samples_directory}"

    print(f"Found {len(gsm_files)} matrix files: {[os.path.basename(f) for f in gsm_files]}")

    dataframes = {}
    sample_mapping = {}   # cell -> sample mapping (细胞 -> 样本 映射)
    all_genes = set()     # union of all gene names (所有基因名的并集)

    for file_path in gsm_files:
        filename = os.path.basename(file_path)
        # Extract sample name from filename, e.g. GSM3589420_PP019swap
        parts = re.split(r'[._]', filename)
        sample_name = parts[0] + '_' + parts[1] if len(parts) > 1 else parts[0]

        df = _read_matrix_file(file_path)
        if df is None:
            print(f"Warning: could not read {file_path}, skipping")
            continue

        # Prefix cell barcodes for uniqueness (为细胞名添加样本前缀以确保唯一性)
        df.columns = [f"{sample_name}_{col}" for col in df.columns]
        dataframes[sample_name] = df
        all_genes.update(df.index)

        for cell_barcode in df.columns:
            sample_mapping[cell_barcode] = sample_name

        print(f"Loaded: {filename}, shape: {df.shape}, sample: {sample_name}")

    if not dataframes:
        return "Error: no matrix files were successfully read"

    # Fill missing genes with zeros (为每个 DataFrame 补齐缺失基因)
    for sample_name, df in dataframes.items():
        missing_genes = all_genes - set(df.index)
        if missing_genes:
            missing_df = pd.DataFrame(0, index=list(missing_genes), columns=df.columns)
            dataframes[sample_name] = pd.concat([df, missing_df])

    # Concatenate by columns = more cells (按列合并 = 增加细胞数)
    combined_df = pd.concat(dataframes.values(), axis=1)
    combined_df = combined_df.reindex(sorted(all_genes))
    combined_df = combined_df[~combined_df.index.duplicated(keep='first')]

    print(f"Merged DataFrame shape: {combined_df.shape}")

    # Build AnnData: transpose (gene x cell) -> (cell x gene)
    # (创建 AnnData 对象 — 将 (基因, 细胞) 转置为 (细胞, 基因))
    adata = ad.AnnData(
        X=combined_df.values.T,
        obs=pd.DataFrame(index=combined_df.columns),
        var=pd.DataFrame(index=combined_df.index),
    )

    # Store sample identity in obs (添加样本信息到 obs)
    adata.obs['sample'] = [sample_mapping.get(cell, 'Unknown') for cell in combined_df.columns]
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # Save output (保存)
    if output_file is None:
        output_file = os.path.join(samples_directory, 'merged_geo_samples.h5ad')
    adata.write(output_file)

    # Generate report (生成报告)
    report = "GEO Matrix Merge Complete\n"
    report += "=" * 50 + "\n"
    report += f"Input directory: {samples_directory}\n"
    report += f"Matrix files: {len(gsm_files)}\n"
    report += f"Merged shape: {adata.shape[0]} cells x {adata.shape[1]} genes\n"
    report += f"Sample count: {len(set(adata.obs['sample']))}\n"
    report += f"Sample info stored in adata.obs['sample']\n"
    report += f"Output file: {output_file}\n"

    sample_counts = adata.obs['sample'].value_counts()
    report += "\nSample distribution:\n"
    for sample, count in sample_counts.items():
        report += f"  {sample}: {count} cells\n"

    return report


def inspect_geo_matrix_file(file_path: str, n_rows: int = 10) -> str:
    """
    Inspect the structure of a single GEO matrix file.
    (检查单个 GEO 矩阵文件的结构)

    Args:
        file_path: Path to the matrix file
        n_rows: Number of rows to preview

    Returns:
        File structure report string
    """
    df = _read_matrix_file(file_path)
    if df is None:
        return f"Error: could not read {file_path}"

    report = f"File inspection: {os.path.basename(file_path)}\n"
    report += "=" * 50 + "\n"
    report += f"Path: {file_path}\n"
    report += f"Shape: {df.shape[0]} rows x {df.shape[1]} cols\n"
    report += f"Row names (genes): {list(df.index[:5])}\n"
    report += f"Col names (cells): {list(df.columns[:5])}\n"
    report += f"Data type: {df.dtypes.iloc[0] if len(df.columns) > 0 else 'N/A'}\n"

    if not df.empty and df.select_dtypes(include=[np.number]).shape[1] > 0:
        report += f"Value range: {df.min().min():.2f} ~ {df.max().max():.2f}\n"
        report += f"Non-zero ratio: {(df != 0).sum().sum() / df.size:.2%}\n"

    report += f"\nFirst {n_rows} rows:\n"
    report += df.head(n_rows).to_string() + "\n"
    return report


# --------------------------------------------------
# CLI entry point (CLI 入口)
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="GEO Matrix Tool - merge samples or inspect file structure"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # merge sub-command
    merge_parser = subparsers.add_parser("merge", help="Merge matrix files into AnnData")
    merge_parser.add_argument("directory", help="Directory containing matrix files")
    merge_parser.add_argument("-o", "--output", default=None,
                              help="Output .h5ad path (default: merged_geo_samples.h5ad)")

    # inspect sub-command
    inspect_parser = subparsers.add_parser("inspect", help="Inspect a matrix file")
    inspect_parser.add_argument("file", help="Path to the matrix file")
    inspect_parser.add_argument("-n", "--rows", type=int, default=10,
                                help="Number of rows to preview (default: 10)")

    args = parser.parse_args()

    if args.command == "merge":
        result = merge_geo_matrix_files(args.directory, args.output)
    elif args.command == "inspect":
        result = inspect_geo_matrix_file(args.file, args.rows)
    else:
        parser.print_help()
        sys.exit(1)

    print(result)


if __name__ == "__main__":
    main()
