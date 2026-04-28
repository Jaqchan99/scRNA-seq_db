"""
统计 benchmark 在不同器官/组织上的覆盖情况。

优先级：
1) 从 h5ad 的 obs/uns 里找 tissue/organ 字段
2) 从文件名解析（规则字典）
3) 无法识别则标记为 UNKNOWN，输出待补充清单

用法：
  python summarize_benchmark_coverage.py --benchmark data/benchmark/benchmark_labels.csv --input_dir .
  python summarize_benchmark_coverage.py --benchmark data/benchmark/benchmark_labels.csv --file_list file_paths.txt
"""
from __future__ import annotations
import os
import re
import argparse
from collections import Counter, defaultdict

import pandas as pd


_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TISSUE_KEYS = [
    # obs 或 uns 可能出现的字段
    "tissue", "organ", "organism_part", "anatomical_region",
    "tissue_type", "source_tissue", "sample_tissue", "tissue_origin",
    "organism_part_ontology_term_id", "tissue_ontology_term_id",
]

# 文件名解析规则（可按需扩充）
FILENAME_TISSUE_MAP = {
    "bladder": "Bladder",
    "liver": "Liver",
    "kidney": "Kidney",
    "lung": "Lung",
    "heart": "Heart",
    "brain": "Brain",
    "spleen": "Spleen",
    "pancreas": "Pancreas",
    "intestine": "Intestine",
    "smallintestine": "Intestine",
    "colon": "Colon",
    "stomach": "Stomach",
    "muscle": "Muscle",
    "bonemarrow": "Bone marrow",
    "blood": "Blood",
    "skin": "Skin",
    "adipose": "Adipose",
    "thymus": "Thymus",
    "lymph": "Lymph",
    "retina": "Retina",
    "ovary": "Ovary",
    "testis": "Testis",
    "uterus": "Uterus",
    "embryo": "Embryo",
    "placenta": "Placenta",
}


def normalize_text(x: str) -> str:
    x = str(x).strip().lower()
    x = re.sub(r"[^a-z0-9]+", " ", x)
    return re.sub(r"\s+", " ", x).strip()


def infer_tissue_from_filename(fname: str) -> str:
    norm = normalize_text(fname)
    for key, tissue in FILENAME_TISSUE_MAP.items():
        if key in norm.replace(" ", "") or key in norm:
            return tissue
    return "UNKNOWN"


def infer_tissue_from_h5ad(path: str) -> tuple[str, str]:
    try:
        import anndata as ad
    except Exception:
        return "UNKNOWN", "no_anndata"

    if not os.path.exists(path):
        return "UNKNOWN", "missing_file"

    try:
        adata = ad.read_h5ad(path, backed="r")
    except Exception:
        return "UNKNOWN", "read_failed"

    # 1) obs 中找
    for key in TISSUE_KEYS:
        if key in adata.obs.columns:
            values = adata.obs[key].dropna().astype(str).tolist()
            if values:
                top = Counter(values).most_common(1)[0][0]
                return str(top), f"obs:{key}"

    # 2) uns 中找
    for key in TISSUE_KEYS:
        if key in adata.uns:
            val = adata.uns.get(key)
            if isinstance(val, (str, int, float)):
                return str(val), f"uns:{key}"

    return "UNKNOWN", "not_found"


def load_file_paths(input_dir: str, file_list: str | None) -> dict:
    """
    返回 {filename: full_path}，用于根据 benchmark 的 source_file 反查路径。
    """
    mapping = {}
    if file_list:
        with open(file_list, "r", encoding="utf-8") as f:
            for line in f:
                p = line.strip()
                if not p:
                    continue
                mapping[os.path.basename(p)] = p
        return mapping

    # 默认从目录中扫描
    if input_dir and os.path.isdir(input_dir):
        for fname in os.listdir(input_dir):
            if fname.endswith(".h5ad"):
                mapping[fname] = os.path.join(input_dir, fname)
    return mapping


def main():
    parser = argparse.ArgumentParser()
    default_benchmark = os.path.join(_ROOT, "data", "benchmark", "benchmark_labels.csv")
    parser.add_argument("--benchmark", default=default_benchmark)
    parser.add_argument("--input_dir", default=".")
    parser.add_argument("--file_list", default="")
    coverage_dir = os.path.join(_ROOT, "results", "tissue_coverage")
    parser.add_argument("--out_detail", default=os.path.join(coverage_dir, "tissue_coverage_detail.csv"))
    parser.add_argument("--out_summary", default=os.path.join(coverage_dir, "tissue_coverage_summary.csv"))
    parser.add_argument("--out_missing", default=os.path.join(coverage_dir, "tissue_coverage_missing.csv"))
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out_detail)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    bench = pd.read_csv(args.benchmark, encoding="utf-8-sig")
    files = bench["source_file"].dropna().unique().tolist()

    path_map = load_file_paths(args.input_dir, args.file_list or None)

    detail_rows = []
    missing = []

    # 统计每个文件在 benchmark 中的 label 数
    label_counts = bench.groupby("source_file")["label_text"].nunique().to_dict()
    cell_counts = bench.groupby("source_file")["cell_count"].sum().to_dict()

    for fname in files:
        full_path = path_map.get(fname, "")
        tissue, source = infer_tissue_from_h5ad(full_path) if full_path else ("UNKNOWN", "missing_path")
        if tissue == "UNKNOWN":
            tissue = infer_tissue_from_filename(fname)
            if tissue != "UNKNOWN":
                source = "filename"

        if tissue == "UNKNOWN":
            missing.append({"source_file": fname})

        detail_rows.append({
            "source_file": fname,
            "tissue": tissue,
            "tissue_source": source,
            "label_count": int(label_counts.get(fname, 0)),
            "cell_count": int(cell_counts.get(fname, 0)),
            "file_path": full_path,
        })

    detail = pd.DataFrame(detail_rows)
    detail.to_csv(args.out_detail, index=False, encoding="utf-8-sig")

    summary = detail.groupby("tissue").agg(
        datasets=("source_file", "nunique"),
        labels=("label_count", "sum"),
        cells=("cell_count", "sum"),
    ).reset_index().sort_values(by="datasets", ascending=False)
    summary.to_csv(args.out_summary, index=False, encoding="utf-8-sig")

    if missing:
        pd.DataFrame(missing).to_csv(args.out_missing, index=False, encoding="utf-8-sig")

    print(f"完成：{args.out_detail}, {args.out_summary}")
    if missing:
        print(f"未识别组织（需手动补充）：{args.out_missing}")


if __name__ == "__main__":
    main()
