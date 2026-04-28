"""
批量扫描 h5ad 文件，输出每个文件的候选 cell type 列信息到 CSV。
用户审核 CSV 后填入 selected_column，即可作为 extract_benchmark_labels.py 的配置。

用法:
  python batch_inspect.py <h5ad目录>
  python batch_inspect.py <h5ad目录> --output configs/dataset_config.csv
"""
from __future__ import annotations

import os
import sys
import argparse
import pandas as pd
import anndata as ad
from pathlib import Path

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ANNOTATION_KEYWORDS = [
    "cell_type", "celltype", "cell.type", "annotation", "cell_annotation",
    "cell_ontology", "cl_id", "cl_label", "cluster", "leiden", "louvain",
    "cell_identity", "author_cell_type", "free_annotation", "label",
    "most_likely", "immgen", "subtype",
]


def scan_one(file_path: str) -> list[dict]:
    """扫描单个 h5ad，返回所有候选 cell type 列的信息。"""
    fname = os.path.basename(file_path)
    try:
        adata = ad.read_h5ad(file_path, backed="r")
    except Exception as e:
        return [{"file": fname, "error": str(e)}]

    obs_cols = list(adata.obs.columns)
    rows = []

    for col in obs_cols:
        col_lower = col.lower()
        is_candidate = any(kw in col_lower for kw in ANNOTATION_KEYWORDS)
        is_numeric = pd.api.types.is_numeric_dtype(adata.obs[col])

        if not is_candidate:
            continue

        series = adata.obs[col].dropna()
        nunique = series.nunique()
        examples = sorted([str(v) for v in series.unique()[:15]])

        rows.append({
            "file": fname,
            "column": col,
            "dtype": str(adata.obs[col].dtype),
            "is_numeric": is_numeric,
            "n_unique": nunique,
            "n_cells": len(series),
            "examples": " | ".join(examples[:10]),
            "selected_column": "",
        })

    if not rows:
        rows.append({
            "file": fname,
            "column": "(无候选列)",
            "dtype": "",
            "is_numeric": "",
            "n_unique": "",
            "n_cells": adata.n_obs,
            "examples": f"所有obs列: {', '.join(obs_cols[:20])}",
            "selected_column": "",
        })

    if hasattr(adata, "file"):
        adata.file.close()

    return rows


def main():
    parser = argparse.ArgumentParser(description="批量扫描 h5ad 的候选 cell type 列")
    parser.add_argument("input_dir", help="h5ad 文件所在目录")
    default_output = os.path.join(_ROOT, "configs", "dataset_config.csv")
    parser.add_argument("--output", "-o", default=default_output,
                        help="输出 CSV 路径 (默认 configs/dataset_config.csv)")
    args = parser.parse_args()

    h5ad_dir = Path(args.input_dir)
    if not h5ad_dir.is_dir():
        print(f"错误: {args.input_dir} 不是一个目录")
        sys.exit(1)

    files = sorted(h5ad_dir.glob("**/*.h5ad"))
    if not files:
        print(f"在 {args.input_dir} 中未找到 .h5ad 文件")
        sys.exit(1)

    print(f"共找到 {len(files)} 个 h5ad 文件，开始扫描…\n")

    all_rows = []
    for i, fp in enumerate(files):
        print(f"  [{i+1}/{len(files)}] {fp.name}")
        rows = scan_one(str(fp))
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows)

    col_order = ["file", "column", "dtype", "is_numeric", "n_unique",
                 "n_cells", "examples", "selected_column"]
    for c in col_order:
        if c not in df.columns:
            df[c] = ""
    df = df[col_order]

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output, index=False, encoding="utf-8-sig")
    print(f"\n已写入: {args.output}")

    unique_files = df["file"].nunique() if "error" not in df.columns else df["file"].nunique()
    print(f"共扫描 {unique_files} 个文件，{len(df)} 条候选列记录")
    print(f"\n下一步: 打开 {args.output}，在 selected_column 列填入你选定的列名")
    print("  - 每个文件只填一个列名（主矩阵的 cell type 列）")
    print("  - 填好后运行 extract_benchmark_labels.py 即可精确提取")


if __name__ == "__main__":
    main()
