"""
从多个 h5ad 文件中批量提取「原始细胞类型标签 → CL ID」用于构建 Benchmark。

工作流:
  1. 先跑 batch_inspect.py 生成 dataset_config.csv
  2. 人工在 CSV 的 selected_column 列填入每个文件要用的主矩阵 cell type 列名
     如需同时抽取 subtype，可额外新增 selected_subtype_column 列
  3. 再跑本脚本，精确提取

用法:
  python extract_benchmark_labels.py <h5ad目录> --config configs/dataset_config.csv
  python extract_benchmark_labels.py <h5ad目录> --config configs/dataset_config.csv --output data/benchmark/benchmark_labels.csv
"""
import os
import sys
import argparse
import pandas as pd
import anndata as ad
from pathlib import Path

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CL_ID_COLUMNS = [
    "cell_ontology_id", "cl_id", "ontology_id", "CL_id",
    "cell_type_ontology_term_id",          # CELLxGENE / Tabula Muris Senis 格式
]
CL_LABEL_COLUMNS = [
    "cell_ontology_class", "cl_label", "cell_ontology_class_label",
    "cell_ontology_label", "ontology_label", "CL_label",
    "cell_type",                            # CELLxGENE 格式（标准化 cell type 名）
]


def find_column(obs_columns, candidates):
    """在列名列表中找第一个匹配的候选列。"""
    for c in candidates:
        if c in obs_columns:
            return c
    return None


def load_config(config_path: str) -> dict:
    """
    读取 dataset_config.csv，返回 {filename: {"main": ..., "subtype": ...}} 映射。
    只保留 selected_column/selected_subtype_column 至少一列非空的行。
    """
    df = pd.read_csv(config_path, encoding="utf-8-sig")
    if "selected_column" not in df.columns:
        print("错误: config CSV 中缺少 selected_column 列")
        sys.exit(1)

    if "selected_subtype_column" not in df.columns:
        df["selected_subtype_column"] = ""

    df["selected_column"] = df["selected_column"].fillna("").astype(str).str.strip()
    df["selected_subtype_column"] = df["selected_subtype_column"].fillna("").astype(str).str.strip()

    # 同一 file 可能有多行（batch_inspect 每个候选列一行）。用户常把 selected_column
    # 填在「对应候选列」那一行，而不是第一行；不能只 keep="first"。
    def _first_nonempty(s: pd.Series) -> str:
        nonempty = s[s != ""]
        return nonempty.iloc[0] if not nonempty.empty else ""

    grouped = df.groupby("file", sort=False).agg(
        selected_column=("selected_column", _first_nonempty),
        selected_subtype_column=("selected_subtype_column", _first_nonempty),
    ).reset_index()

    mapping = {}
    for _, row in grouped.iterrows():
        fname = str(row["file"]).strip()
        if not fname:
            continue
        main_col = row["selected_column"]
        subtype_col = row["selected_subtype_column"]
        if not main_col and not subtype_col:
            continue
        mapping[fname] = {"main": main_col, "subtype": subtype_col}
    return mapping


def extract_from_one(file_path: str, cell_type_col: str, label_level: str) -> pd.DataFrame:
    """
    从单个 h5ad 精确提取：
      - cell_type_col 作为原始标签 (label_text)
      - 自动检测 CL ID 列作为 ground truth
      - 自动检测 CL label 列作为标准名
    """
    fname = os.path.basename(file_path)
    try:
        adata = ad.read_h5ad(file_path)
    except Exception as e:
        print(f"  [跳过] 无法读取: {fname} —— {e}")
        return pd.DataFrame()

    obs_cols = list(adata.obs.columns)

    if cell_type_col not in obs_cols:
        print(f"  [跳过] {fname} 中不存在列 '{cell_type_col}'")
        return pd.DataFrame()

    cl_id_col = find_column(obs_cols, CL_ID_COLUMNS)
    cl_label_col = find_column(obs_cols, CL_LABEL_COLUMNS)

    # 避免把 selected_column 自身也当作 cl_label_col
    if cl_label_col == cell_type_col:
        cl_label_col = None

    series_label = adata.obs[cell_type_col].dropna().astype(str)

    if cl_id_col and cl_id_col != cell_type_col:
        series_clid = adata.obs[cl_id_col].astype(str).reindex(series_label.index)
        tmp = pd.DataFrame({
            "label_text": series_label,
            "cl_id": series_clid,
        })
        if cl_label_col:
            tmp["cl_label"] = adata.obs[cl_label_col].astype(str).reindex(series_label.index)
        else:
            tmp["cl_label"] = ""

        tmp = tmp[tmp["cl_id"].notna() & (tmp["cl_id"] != "") & (tmp["cl_id"] != "nan")]

        agg = tmp.groupby(["label_text", "cl_id", "cl_label"], dropna=False).size().reset_index(name="cell_count")
        agg["has_ground_truth"] = True
    else:
        counts = series_label.value_counts()
        agg = pd.DataFrame({
            "label_text": counts.index,
            "cell_count": counts.values,
        })
        agg["cl_id"] = ""
        agg["cl_label"] = ""
        agg["has_ground_truth"] = False

    agg["source_file"] = fname
    agg["selected_column"] = cell_type_col
    agg["label_level"] = label_level
    return agg[["source_file", "selected_column", "label_level", "label_text", "cl_id", "cl_label", "cell_count", "has_ground_truth"]]


def collect_h5ad_paths(input_path: str):
    """输入可以是目录或一个每行一个路径的 .txt 文件。"""
    path = Path(input_path)
    if path.is_file() and path.suffix == ".txt":
        with open(path, "r", encoding="utf-8") as f:
            lines = [p.strip() for p in f if p.strip()]
        return [p for p in lines if p.endswith(".h5ad")]
    if path.is_dir():
        return sorted([str(p) for p in path.glob("**/*.h5ad")])
    return []


def main():
    parser = argparse.ArgumentParser(description="基于配置精确提取 cell type 标签用于 benchmark")
    parser.add_argument("input", help="h5ad 所在目录，或包含路径列表的 .txt 文件")
    parser.add_argument("--config", "-c", required=True,
                        help="dataset_config.csv 路径 (batch_inspect.py 生成，人工填好 selected_column)")
    default_output = os.path.join(_ROOT, "data", "benchmark", "benchmark_labels.csv")
    parser.add_argument("--output", "-o", default=default_output, help="输出 CSV 路径")
    args = parser.parse_args()

    config = load_config(args.config)
    if not config:
        print("错误: config 中没有任何文件有 selected_column/selected_subtype_column，请先填写")
        sys.exit(1)
    print(f"从 config 读取到 {len(config)} 个文件的列配置")

    paths = collect_h5ad_paths(args.input)
    if not paths:
        print(f"在 {args.input} 中未找到 .h5ad 文件")
        sys.exit(1)

    path_map = {os.path.basename(p): p for p in paths}

    all_dfs = []
    processed = 0
    skipped = 0
    for fname, cols in config.items():
        if fname not in path_map:
            print(f"  [未找到文件] {fname}")
            skipped += 1
            continue
        processed += 1
        main_col = cols.get("main", "").strip()
        subtype_col = cols.get("subtype", "").strip()
        print(f"  [{processed}/{len(config)}] {fname}")
        if main_col:
            print(f"     - main: {main_col}")
            df_main = extract_from_one(path_map[fname], main_col, "main")
            if df_main is not None and not df_main.empty:
                all_dfs.append(df_main)
        if subtype_col and subtype_col != main_col:
            print(f"     - subtype: {subtype_col}")
            df_sub = extract_from_one(path_map[fname], subtype_col, "subtype")
            if df_sub is not None and not df_sub.empty:
                all_dfs.append(df_sub)

    if not all_dfs:
        print("没有提取到任何记录")
        sys.exit(0)

    out = pd.concat(all_dfs, ignore_index=True)
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    out.to_csv(args.output, index=False, encoding="utf-8-sig")

    n_gt = out["has_ground_truth"].sum()
    n_unique_labels = out["label_text"].nunique()

    print(f"\n{'='*50}")
    print(f"已写入: {args.output}")
    print(f"总记录数: {len(out)}")
    print(f"唯一标签数 (label_text): {n_unique_labels}")
    print(f"已有 ground truth (cl_id 非空): {n_gt} 条")
    print(f"需人工标注: {len(out) - n_gt} 条")
    print(f"跳过文件: {skipped} 个（config 中有但目录中未找到）")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
