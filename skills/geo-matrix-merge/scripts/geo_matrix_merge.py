#!/usr/bin/env python3
"""Merge GEO matrix files into a single AnnData (h5ad) object.

Supports four auto-detected formats:
  1. 10X:           matrix.mtx.gz + features.tsv.gz + barcodes.tsv.gz
  2. counts+cellname: *.counts.* + *.cellname.*
  3. Single GEO:    one *.txt.gz / *.tsv.gz
  4. Multi GEO:     multiple GSM*.txt.gz / *.filtered.matrix.txt.gz

Usage:
    python geo_matrix_merge.py <samples_directory> [-o output.h5ad]
    python geo_matrix_merge.py --inspect <file_path> [-n 10]
"""

import argparse
import os
import sys
import gzip
import tempfile
import pandas as pd
import numpy as np
import anndata as ad
from scipy import sparse


# ──────────────────────────────────────────────
# Format detection
# ──────────────────────────────────────────────

def detect_format(directory):
    """Auto-detect data format from directory contents.

    Returns one of: '10x', 'counts_cellname', 'single_geo', 'multi_geo'
    """
    files = os.listdir(directory)
    basenames = [os.path.basename(f) for f in files]

    # 1. 10X format: matrix.mtx.gz (+ features/barcodes sidecar files)
    mtx_files = [f for f in basenames
                 if f.endswith('.mtx.gz') or f.endswith('.mtx')]
    if mtx_files:
        if len(mtx_files) == 1:
            print("[格式检测] 10X 格式 (matrix.mtx + features/barcodes)")
            return '10x'
        else:
            print(f"[格式检测] 多样本 10X 格式 ({len(mtx_files)} 个样本)")
            return 'multi_10x'

    # 2. counts + cellname paired format
    counts_files = [f for f in basenames if '.counts.' in f]
    cellname_files = [f for f in basenames if '.cellname.' in f]
    if counts_files and cellname_files:
        print(f"[格式检测] counts+cellname 配对格式")
        print(f"  counts 文件: {len(counts_files)} 个")
        print(f"  cellname 文件: {len(cellname_files)} 个")
        return 'counts_cellname'

    # 3 & 4. GEO text matrix files (*.txt.gz, *.tsv.gz, *.txt, *.tsv)
    geo_patterns = ('.txt.gz', '.tsv.gz', '.txt', '.tsv')
    geo_files = [f for f in basenames
                 if any(f.endswith(ext) for ext in geo_patterns)
                 and not f.startswith('.')]
    # Exclude cellname / counts files that might also match
    geo_files = [f for f in geo_files
                 if '.cellname.' not in f and '.counts.' not in f]

    if len(geo_files) == 1:
        print(f"[格式检测] 单个 GEO 矩阵文件: {geo_files[0]}")
        return 'single_geo'
    elif len(geo_files) > 1:
        print(f"[格式检测] 多个 GEO 矩阵文件: {len(geo_files)} 个")
        return 'multi_geo'

    print("错误：无法识别目录中的数据格式", file=sys.stderr)
    print(f"  目录: {directory}", file=sys.stderr)
    print(f"  文件列表: {basenames}", file=sys.stderr)
    sys.exit(1)


# ──────────────────────────────────────────────
# File readers
# ──────────────────────────────────────────────

def read_geo_text_matrix(file_path):
    """Read a GEO tab/comma-separated matrix file.

    Returns (DataFrame, gene_symbols_dict).
    DataFrame has gene IDs as index, cell barcodes as columns.
    """
    try:
        if file_path.endswith('.gz'):
            with gzip.open(file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t', index_col=0)
        else:
            with open(file_path, 'r') as f:
                first_line = f.readline()
            sep = '\t' if '\t' in first_line else (',' if ',' in first_line else r'\s+')
            df = pd.read_csv(file_path, sep=sep, engine='python', index_col=0)

        gene_symbols = {}
        if 'Gene' in df.columns:
            gene_symbols = dict(zip(df.index, df['Gene']))
            df = df.drop('Gene', axis=1)

        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) < df.shape[1]:
            df = df[numeric_cols]
        df = df.apply(pd.to_numeric, errors='coerce').fillna(0)

        return df, gene_symbols

    except Exception as e:
        print(f"读取文件 {file_path} 时出错: {e}", file=sys.stderr)
        return None, {}


def read_10x_directory(directory):
    """Read 10X-style format: matrix.mtx.gz + features.tsv.gz + barcodes.tsv.gz.

    Returns (sparse_matrix_genes_x_cells, gene_ids, gene_symbols, barcodes).
    """
    import scipy.io

    files = os.listdir(directory)

    # Find matrix file
    mtx_file = None
    for f in files:
        if f.endswith('.mtx.gz') or f.endswith('.mtx'):
            mtx_file = os.path.join(directory, f)
            break
    if mtx_file is None:
        print("错误：未找到 matrix.mtx.gz 文件", file=sys.stderr)
        sys.exit(1)

    # Find features/genes file
    features_file = None
    for f in files:
        if f in ('features.tsv.gz', 'features.tsv',
                 'genes.tsv.gz', 'genes.tsv'):
            features_file = os.path.join(directory, f)
            break
    if features_file is None:
        for f in files:
            if 'feature' in f.lower() or 'gene' in f.lower():
                if f.endswith('.tsv.gz') or f.endswith('.tsv') or f.endswith('.txt'):
                    features_file = os.path.join(directory, f)
                    break

    # Find barcodes file
    barcodes_file = None
    for f in files:
        if f in ('barcodes.tsv.gz', 'barcodes.tsv'):
            barcodes_file = os.path.join(directory, f)
            break
    if barcodes_file is None:
        for f in files:
            if 'barcode' in f.lower():
                if f.endswith('.tsv.gz') or f.endswith('.tsv') or f.endswith('.txt'):
                    barcodes_file = os.path.join(directory, f)
                    break

    # Read matrix
    print(f"  读取矩阵: {os.path.basename(mtx_file)}")
    if mtx_file.endswith('.gz'):
        with gzip.open(mtx_file, 'rb') as f_in:
            with tempfile.NamedTemporaryFile(delete=False, suffix='.mtx') as f_out:
                f_out.write(f_in.read())
                temp_path = f_out.name
        matrix = scipy.io.mmread(temp_path)
        os.unlink(temp_path)
    else:
        matrix = scipy.io.mmread(mtx_file)
    matrix = sparse.csr_matrix(matrix)

    # Read features/genes
    gene_ids = []
    gene_symbols = {}
    if features_file:
        print(f"  读取基因: {os.path.basename(features_file)}")
        feat_df = pd.read_csv(features_file, sep='\t', header=None)
        if feat_df.shape[1] >= 2:
            gene_ids = list(feat_df.iloc[:, 0].astype(str))
            gene_symbols = dict(zip(gene_ids, feat_df.iloc[:, 1].astype(str)))
        else:
            gene_ids = list(feat_df.iloc[:, 0].astype(str))
    else:
        gene_ids = [f"GENE_{i}" for i in range(matrix.shape[0])]
        print("  警告：未找到 features/genes 文件，使用占位基因名")

    # Read barcodes
    barcodes = []
    if barcodes_file:
        print(f"  读取barcode: {os.path.basename(barcodes_file)}")
        bc_df = pd.read_csv(barcodes_file, sep='\t', header=None)
        barcodes = list(bc_df.iloc[:, 0].astype(str))
    else:
        barcodes = [f"CELL_{i}" for i in range(matrix.shape[1])]
        print("  警告：未找到 barcodes 文件，使用占位barcode名")

    print(f"  矩阵形状: {matrix.shape} (基因 × 细胞), 非零: {matrix.nnz}")
    return matrix, gene_ids, gene_symbols, barcodes


def read_cellname_mapping(cellname_path):
    """Read a cellname.list.txt.gz file and return {CellIndex: CellName} mapping."""
    try:
        if cellname_path.endswith('.gz'):
            with gzip.open(cellname_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t')
        else:
            df = pd.read_csv(cellname_path, sep='\t')

        if 'CellName' in df.columns and 'CellIndex' in df.columns:
            return dict(zip(df['CellIndex'], df['CellName']))
        elif df.shape[1] >= 2:
            return dict(zip(df.iloc[:, 1], df.iloc[:, 0]))
        else:
            return {f"C{i+1}": name for i, name in enumerate(df.iloc[:, 0])}
    except Exception as e:
        print(f"读取 cellname 文件 {cellname_path} 时出错: {e}", file=sys.stderr)
        return None


def find_cellname_for_counts(counts_path, directory):
    """Find the matching cellname file for a counts file."""
    sample_prefix = os.path.basename(counts_path).split('.')[0]
    for f in os.listdir(directory):
        if '.cellname.' in f and f.startswith(sample_prefix):
            return os.path.join(directory, f)
    return None


# ──────────────────────────────────────────────
# Sparse merge engine
# ──────────────────────────────────────────────

def finalize_sparse_merge(sparse_blocks, all_gene_ids, cell_labels,
                          sample_labels, sample_names, gene_symbols,
                          directory, output_file, n_files):
    """Reindex genes, transpose, vstack, create AnnData, save."""
    if not sparse_blocks:
        print("错误：没有成功读取任何矩阵", file=sys.stderr)
        sys.exit(1)

    # Build global gene index
    print("\n构建全局基因集...")
    gene_to_idx = {}
    for gids in all_gene_ids:
        for g in gids:
            if g not in gene_to_idx:
                gene_to_idx[g] = len(gene_to_idx)
    all_genes = list(gene_to_idx.keys())
    n_genes = len(all_genes)
    print(f"唯一基因总数: {n_genes}")

    # Reindex + transpose each block
    print("重索引并转置矩阵...")
    final_blocks = []
    for idx, (sp, gids) in enumerate(zip(sparse_blocks, all_gene_ids)):
        n_cells = sp.shape[1]
        new_row_idx = np.array([gene_to_idx[g] for g in gids], dtype=np.int32)

        coo = sp.tocoo()
        new_coo = sparse.coo_matrix(
            (coo.data, (new_row_idx[coo.row], coo.col)),
            shape=(n_genes, n_cells),
            dtype=np.float32
        )
        final_blocks.append(new_coo.tocsr().T.tocsr())
        print(f"  {sample_names[idx]}: 完成")
        del coo, new_coo

    # Stack all blocks
    print("\n合并所有细胞...")
    combined = sparse.vstack(final_blocks, format='csr')
    print(f"合并后矩阵形状: {combined.shape} (细胞 × 基因)")
    total = combined.shape[0] * combined.shape[1]
    print(f"非零元素: {combined.nnz}, 密度: {combined.nnz / total:.4%}")

    # Create AnnData
    print("创建 AnnData 对象...")
    var_df = pd.DataFrame(index=all_genes)
    if gene_symbols:
        var_df['gene_symbol'] = [gene_symbols.get(g, '') for g in all_genes]

    adata = ad.AnnData(
        X=combined,
        obs=pd.DataFrame(index=cell_labels),
        var=var_df
    )
    adata.obs['sample'] = sample_labels
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if output_file is None:
        output_file = os.path.join(directory, 'merged_geo_samples.h5ad')

    adata.write(output_file)
    _print_summary(adata, directory, output_file, n_files)
    return adata


def _print_summary(adata, directory, output_file, n_files):
    """Print merge summary."""
    print(f"\nGEO矩阵文件合并完成")
    print("=" * 50)
    print(f"输入目录: {directory}")
    print(f"矩阵文件数: {n_files}")
    print(f"合并后数据形状: {adata.shape[0]} 细胞 × {adata.shape[1]} 基因")
    print(f"样本数量: {len(set(adata.obs['sample']))}")
    print(f"样本信息已添加到 adata.obs['sample']")
    print(f"输出文件: {output_file}")

    sample_counts = adata.obs['sample'].value_counts()
    print(f"\n样本分布:")
    for sample, count in sample_counts.items():
        print(f"  {sample}: {count} 个细胞")


# ──────────────────────────────────────────────
# Format-specific merge handlers
# ──────────────────────────────────────────────

def merge_10x(directory, output_file=None):
    """Handle 10X format: matrix.mtx.gz + features/barcodes."""
    matrix, gene_ids, gene_symbols, barcodes = read_10x_directory(directory)

    sample_name = os.path.basename(directory)
    prefixed = [f"{sample_name}_{bc}" for bc in barcodes]

    return finalize_sparse_merge(
        sparse_blocks=[matrix],
        all_gene_ids=[gene_ids],
        cell_labels=prefixed,
        sample_labels=[sample_name] * len(barcodes),
        sample_names=[sample_name],
        gene_symbols=gene_symbols,
        directory=directory,
        output_file=output_file,
        n_files=1
    )


def merge_multi_10x(directory, output_file=None):
    """Handle multiple 10X samples in one directory.

    Files are named like:
      GSM6063506_230093_primary_focus_matrix.mtx.gz
      GSM6063506_230093_primary_focus_features.tsv.gz
      GSM6063506_230093_primary_focus_barcodes.tsv.gz
    """
    import scipy.io

    all_files = sorted(os.listdir(directory))
    mtx_files = [f for f in all_files
                 if f.endswith('.mtx.gz') or f.endswith('.mtx')]

    # Group by sample prefix: strip _matrix.mtx.gz to get sample prefix
    samples = {}
    for mtx_f in mtx_files:
        prefix = mtx_f.replace('_matrix.mtx.gz', '').replace('_matrix.mtx', '')
        samples[prefix] = {'matrix': mtx_f}

    # Find matching features and barcodes for each sample
    for prefix in samples:
        for f in all_files:
            if f.startswith(prefix + '_'):
                if 'features' in f or 'genes' in f:
                    samples[prefix]['features'] = f
                elif 'barcode' in f:
                    samples[prefix]['barcodes'] = f

    print(f"\n找到 {len(samples)} 个 10X 样本")

    sparse_blocks = []
    all_gene_ids = []
    cell_labels = []
    sample_labels = []
    sample_names = []
    all_gene_symbols = {}

    for prefix in sorted(samples.keys()):
        info = samples[prefix]
        sample_name = prefix
        sample_names.append(sample_name)

        # Read matrix
        mtx_path = os.path.join(directory, info['matrix'])
        print(f"\n  样本: {sample_name}")
        print(f"    矩阵: {info['matrix']}")

        if mtx_path.endswith('.gz'):
            with gzip.open(mtx_path, 'rb') as f_in:
                with tempfile.NamedTemporaryFile(delete=False, suffix='.mtx') as f_out:
                    f_out.write(f_in.read())
                    temp_path = f_out.name
            matrix = scipy.io.mmread(temp_path)
            os.unlink(temp_path)
        else:
            matrix = scipy.io.mmread(mtx_path)
        matrix = sparse.csr_matrix(matrix)

        # Read features
        gene_ids = []
        gene_symbols = {}
        if 'features' in info:
            feat_path = os.path.join(directory, info['features'])
            print(f"    基因: {info['features']}")
            feat_df = pd.read_csv(feat_path, sep='\t', header=None)
            if feat_df.shape[1] >= 2:
                gene_ids = list(feat_df.iloc[:, 0].astype(str))
                gene_symbols = dict(zip(gene_ids, feat_df.iloc[:, 1].astype(str)))
            else:
                gene_ids = list(feat_df.iloc[:, 0].astype(str))
        else:
            gene_ids = [f"GENE_{i}" for i in range(matrix.shape[0])]

        all_gene_symbols.update(gene_symbols)

        # Read barcodes
        barcodes = []
        if 'barcodes' in info:
            bc_path = os.path.join(directory, info['barcodes'])
            print(f"    Barcode: {info['barcodes']}")
            bc_df = pd.read_csv(bc_path, sep='\t', header=None)
            barcodes = list(bc_df.iloc[:, 0].astype(str))
        else:
            barcodes = [f"CELL_{i}" for i in range(matrix.shape[1])]

        print(f"    形状: {matrix.shape} (基因×细胞), 非零: {matrix.nnz}")

        sparse_blocks.append(matrix)
        all_gene_ids.append(gene_ids)

        prefixed = [f"{sample_name}_{bc}" for bc in barcodes]
        cell_labels.extend(prefixed)
        sample_labels.extend([sample_name] * len(barcodes))

    return finalize_sparse_merge(
        sparse_blocks, all_gene_ids, cell_labels,
        sample_labels, sample_names, all_gene_symbols,
        directory, output_file, len(samples)
    )


def merge_counts_cellname(directory, output_file=None):
    """Handle counts+cellname paired format."""
    all_files = sorted(os.listdir(directory))
    counts_files = [f for f in all_files if '.counts.' in f]

    if not counts_files:
        print("错误：没有找到 counts 文件", file=sys.stderr)
        sys.exit(1)

    print(f"\n找到 {len(counts_files)} 个样本的 counts 文件")

    sparse_blocks = []
    all_gene_ids = []
    cell_labels = []
    sample_labels = []
    sample_names = []
    all_gene_symbols = {}

    for counts_name in counts_files:
        counts_path = os.path.join(directory, counts_name)
        sample_name = counts_name.split('.')[0]
        sample_names.append(sample_name)

        cellname_path = find_cellname_for_counts(counts_path, directory)
        barcode_map = None
        if cellname_path:
            barcode_map = read_cellname_mapping(cellname_path)
            print(f"  已匹配 cellname: {os.path.basename(cellname_path)}")

        df, gene_symbols = read_geo_text_matrix(counts_path)
        if df is None:
            print(f"警告：无法读取 {counts_path}，跳过", file=sys.stderr)
            continue

        if barcode_map:
            df.columns = [barcode_map.get(c, c) for c in df.columns]

        all_gene_symbols.update(gene_symbols)

        sp = sparse.csr_matrix(df.values.astype(np.float32))
        sparse_blocks.append(sp)
        all_gene_ids.append(list(df.index))

        prefixed = [f"{sample_name}_{c}" for c in df.columns]
        cell_labels.extend(prefixed)
        sample_labels.extend([sample_name] * df.shape[1])

        print(f"  已读取: {counts_name}, 形状: {df.shape}, 非零: {sp.nnz}")
        del df

    return finalize_sparse_merge(
        sparse_blocks, all_gene_ids, cell_labels,
        sample_labels, sample_names, all_gene_symbols,
        directory, output_file, len(counts_files)
    )


def merge_single_geo(directory, output_file=None):
    """Handle a single GEO text matrix file."""
    geo_files = [f for f in os.listdir(directory)
                 if (f.endswith('.txt.gz') or f.endswith('.tsv.gz')
                     or f.endswith('.txt') or f.endswith('.tsv'))
                 and not f.startswith('.')]
    file_path = os.path.join(directory, geo_files[0])
    # Extract sample name: GSM12345_SampleName from GSM12345_SampleName.filtered.matrix.txt.gz
    parts = geo_files[0].split('.')[0].split('_')
    sample_name = '_'.join(parts[:2]) if len(parts) >= 2 else parts[0]

    df, gene_symbols = read_geo_text_matrix(file_path)
    if df is None:
        print(f"错误：无法读取文件 {file_path}", file=sys.stderr)
        sys.exit(1)

    sp = sparse.csr_matrix(df.values.astype(np.float32))
    gene_ids = list(df.index)
    prefixed = [f"{sample_name}_{c}" for c in df.columns]

    print(f"已读取: {geo_files[0]}, 形状: {df.shape}, 非零: {sp.nnz}")
    del df

    return finalize_sparse_merge(
        [sp], [gene_ids],
        prefixed, [sample_name] * sp.shape[1],
        [sample_name], gene_symbols,
        directory, output_file, 1
    )


def merge_multi_geo(directory, output_file=None):
    """Handle multiple GEO text matrix files (GSM*.txt.gz etc.)."""
    geo_exts = ('.txt.gz', '.tsv.gz', '.txt', '.tsv')
    geo_files = sorted([
        os.path.join(directory, f) for f in os.listdir(directory)
        if any(f.endswith(ext) for ext in geo_exts)
        and not f.startswith('.')
        and '.cellname.' not in f
        and '.counts.' not in f
    ])

    if not geo_files:
        print("错误：没有找到 GEO 矩阵文件", file=sys.stderr)
        sys.exit(1)

    print(f"找到 {len(geo_files)} 个矩阵文件")

    sparse_blocks = []
    all_gene_ids = []
    cell_labels = []
    sample_labels = []
    sample_names = []
    all_gene_symbols = {}

    for file_path in geo_files:
        filename = os.path.basename(file_path)
        # Extract sample name: GSM12345_SampleName from GSM12345_SampleName.filtered.matrix.txt.gz
        parts = filename.split('.')[0].split('_')
        sample_name = '_'.join(parts[:2]) if len(parts) >= 2 else parts[0]
        sample_names.append(sample_name)

        df, gene_symbols = read_geo_text_matrix(file_path)
        if df is None:
            print(f"警告：无法读取 {file_path}，跳过", file=sys.stderr)
            continue

        all_gene_symbols.update(gene_symbols)

        sp = sparse.csr_matrix(df.values.astype(np.float32))
        sparse_blocks.append(sp)
        all_gene_ids.append(list(df.index))

        prefixed = [f"{sample_name}_{c}" for c in df.columns]
        cell_labels.extend(prefixed)
        sample_labels.extend([sample_name] * df.shape[1])

        print(f"  已读取: {filename}, 形状: {df.shape}, 非零: {sp.nnz}")
        del df

    return finalize_sparse_merge(
        sparse_blocks, all_gene_ids, cell_labels,
        sample_labels, sample_names, all_gene_symbols,
        directory, output_file, len(geo_files)
    )


# ──────────────────────────────────────────────
# Inspect mode
# ──────────────────────────────────────────────

def inspect_file(file_path, n_rows=10):
    """Preview a single matrix file's structure."""
    df, gene_symbols = read_geo_text_matrix(file_path)
    if df is None:
        print(f"错误：无法读取文件 {file_path}", file=sys.stderr)
        sys.exit(1)

    print(f"文件检查报告: {os.path.basename(file_path)}")
    print("=" * 50)
    print(f"文件路径: {file_path}")
    print(f"数据形状: {df.shape[0]} 行 × {df.shape[1]} 列")
    print(f"行名示例 (基因): {list(df.index[:5])}")
    print(f"列名示例 (细胞): {list(df.columns[:5])}")

    if not df.empty and df.select_dtypes(include=[np.number]).shape[1] > 0:
        print(f"数值范围: {df.min().min():.2f} ~ {df.max().max():.2f}")
        print(f"非零元素比例: {(df != 0).sum().sum() / df.size:.2%}")

    print(f"\n前{n_rows}行预览:")
    print(df.head(n_rows).to_string())


# ──────────────────────────────────────────────
# Main entry
# ──────────────────────────────────────────────

FORMAT_HANDLERS = {
    '10x': merge_10x,
    'multi_10x': merge_multi_10x,
    'counts_cellname': merge_counts_cellname,
    'single_geo': merge_single_geo,
    'multi_geo': merge_multi_geo,
}


def main():
    parser = argparse.ArgumentParser(
        description="Merge GEO matrix files into h5ad or inspect a single file")
    parser.add_argument("path",
                        help="Directory of matrix files (merge) or single file (inspect)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output h5ad file path")
    parser.add_argument("--inspect", action="store_true",
                        help="Inspect a single matrix file")
    parser.add_argument("-n", type=int, default=10,
                        help="Number of rows to show in inspect mode")
    args = parser.parse_args()

    if args.inspect:
        inspect_file(args.path, args.n)
        return

    if not os.path.isdir(args.path):
        print(f"错误：路径不存在或不是目录 - {args.path}", file=sys.stderr)
        sys.exit(1)

    fmt = detect_format(args.path)
    handler = FORMAT_HANDLERS[fmt]
    handler(args.path, args.output)


if __name__ == "__main__":
    main()
