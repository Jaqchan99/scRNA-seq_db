"""
incremental_update_sop.py  —  Cell Type Mapping 增量更新一体化 SOP
=================================================================

将以下分散脚本整合为一条有序流水线（原脚本完全不动）：
  batch_inspect → extract_benchmark_labels → 基线评估 →
  merge_lookup_table → 更新后评估（对比） → coverage 统计

【快速开始】
  # 第一次运行（传入包含 h5ad 的目录）：
  python incremental_update_sop.py --h5ad_dir /path/to/h5ad_folder

  # 填完 dataset_config.csv 后续跑：
  python incremental_update_sop.py --resume sop_state_YYYYMMDD_HHMMSS.json

  # 跳过前面步骤，从指定阶段开始：
  python incremental_update_sop.py --h5ad_dir /path --start_phase merge

【流程阶段说明】
  scan         : 扫描新 h5ad 文件，生成 dataset_config_<run_id>.csv
  *** [手动] 打开 dataset_config csv，在 selected_column 列填入细胞类型列名 ***
  extract      : 提取 (label, CL ID) 对，合并进 benchmark_labels.csv
  eval_before  : 用当前 lookup table 跑评估，输出基线指标 + unmapped_labels.csv
  *** [可选] 在 unmapped_labels.csv 填写 cl_id / cl_label，可提升准确率 ***
  merge        : 把填好的 unmapped 合并进 cell_type_mapping.csv（自动备份）
  eval_after   : 重跑评估，打印前后对比表
  coverage     : 统计 benchmark 覆盖的器官/组织分布

【输出文件】
  results/logs/sop_state_<run_id>.json     : 断点状态文件（用于续跑）
  configs/dataset_config_<run_id>.csv      : 新文件的候选 cell type 列（需人工填写）
  data/benchmark/benchmark_labels.csv      : 持续追加的 benchmark（自动去重）
  data/labels/unmapped_labels.csv          : 未映射的标签（可选填写 cl_id）
  results/eval/eval_before_<run_id>.csv    : 更新前的逐条预测结果
  results/eval/eval_after_<run_id>.csv     : 更新后的逐条预测结果
  results/tissue_coverage/tissue_coverage_<run_id>.csv: 器官覆盖汇总
  results/logs/update_log_<run_id>.json    : 本轮更新的完整记录（供论文引用）
"""

import os
import sys
import json
import shutil
import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

# ── 将项目根目录加入 path，确保能导入同级脚本 ───────────────────────────────
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

# ── 从原脚本导入函数（只调用，绝不修改原脚本）──────────────────────────────
from batch_inspect import scan_one
from extract_benchmark_labels import load_config, extract_from_one, collect_h5ad_paths
from merge_lookup_table import merge as _merge_lookup
from evaluate_mapping import run_local_mapper, compute_metrics, normalize_cl_id, is_correct
from summarize_benchmark_coverage import (
    load_file_paths,
    infer_tissue_from_h5ad,
    infer_tissue_from_filename,
)

# ─── 路径常量（相对于项目根目录）──────────────────────────────────────────

BENCHMARK_PATH = os.path.join(_ROOT, "data", "benchmark", "benchmark_labels.csv")
LOOKUP_TABLE_PATH = os.path.join(_ROOT, "cell_type_mapping", "cell_type_mapping.csv")
UNMAPPED_OUTPUT = os.path.join(_ROOT, "data", "labels", "unmapped_labels.csv")
CONFIG_DIR = os.path.join(_ROOT, "configs")
EVAL_DIR = os.path.join(_ROOT, "results", "eval")
LOG_DIR = os.path.join(_ROOT, "results", "logs")
COVERAGE_DIR = os.path.join(_ROOT, "results", "tissue_coverage")

for _dir in (CONFIG_DIR, EVAL_DIR, LOG_DIR, COVERAGE_DIR):
    os.makedirs(_dir, exist_ok=True)

PHASES = ["scan", "extract", "resolve_cl_label", "eval_before", "merge", "eval_after", "coverage"]


# ═══════════════════════════════════════════════════════════════════════════
# 工具函数
# ═══════════════════════════════════════════════════════════════════════════

def _make_run_id() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _load_state(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _save_state(state: dict, path: str):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(state, f, ensure_ascii=False, indent=2)
    print(f"  [断点已保存 → {path}]")


def _backup(path: str) -> str:
    """备份文件，返回备份路径。文件不存在则静默跳过。"""
    if not os.path.exists(path):
        return ""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = f"{path}.bak_{ts}"
    shutil.copy2(path, backup)
    print(f"  [备份] {os.path.basename(path)} → {os.path.basename(backup)}")
    return backup


def _section(title: str):
    print()
    print("─" * 68)
    print(f"  {title}")
    print("─" * 68)


def _pause(heading: str, steps: list[str]):
    """打印暂停提示，等待用户按 Enter 继续。"""
    print()
    print("═" * 68)
    print(f"  ⏸  需要手动操作：{heading}")
    print("═" * 68)
    for s in steps:
        print(f"  {s}")
    print()
    input("  完成后按 [Enter] 继续...")
    print()


# ═══════════════════════════════════════════════════════════════════════════
# Phase 1 — 扫描新文件
# ═══════════════════════════════════════════════════════════════════════════

def phase_scan(state: dict) -> dict:
    _section("Phase 1 / 6 — 扫描 h5ad 文件（batch_inspect）")

    h5ad_dir = state["h5ad_dir"]
    run_id   = state["run_id"]

    # 收集所有 h5ad 路径
    all_paths = collect_h5ad_paths(h5ad_dir)
    if not all_paths:
        print(f"  ❌ 在 {h5ad_dir} 中未找到任何 .h5ad 文件，请检查路径。")
        sys.exit(1)
    print(f"  目录中共有 {len(all_paths)} 个 .h5ad 文件")

    # 找出 benchmark 中已有的文件
    existing_files: set[str] = set()
    if os.path.exists(BENCHMARK_PATH):
        bench = pd.read_csv(BENCHMARK_PATH, encoding="utf-8-sig")
        existing_files = set(bench["source_file"].dropna().unique())
        print(f"  benchmark_labels.csv 已覆盖 {len(existing_files)} 个数据集")

    new_paths = [p for p in all_paths if os.path.basename(p) not in existing_files]

    if not new_paths:
        print("  ✅ 所有文件已在 benchmark 中，无新文件需要处理。")
        print("     如需重新评估，可运行：python evaluate_mapping.py")
        state["new_files"] = []
        state["current_phase"] = "eval_before"
        return state

    print(f"  🆕 检测到 {len(new_paths)} 个新文件：")
    for p in new_paths:
        print(f"     • {os.path.basename(p)}")

    # 扫描候选 cell type 列
    config_path = os.path.join(CONFIG_DIR, f"dataset_config_{run_id}.csv")
    print(f"\n  正在扫描候选 cell type 列...")
    all_rows = []
    for i, fp in enumerate(new_paths):
        fname = os.path.basename(fp)
        print(f"    [{i+1}/{len(new_paths)}] {fname}")
        rows = scan_one(fp)
        all_rows.extend(rows)

    col_order = ["file", "column", "dtype", "is_numeric", "n_unique",
                 "n_cells", "examples", "selected_column"]
    df_config = pd.DataFrame(all_rows)
    for c in col_order:
        if c not in df_config.columns:
            df_config[c] = ""
    df_config[col_order].to_csv(config_path, index=False, encoding="utf-8-sig")
    print(f"\n  ✅ dataset_config 写入: {config_path}")

    state["new_files"]      = [os.path.basename(p) for p in new_paths]
    state["dataset_config"] = config_path
    state["current_phase"]  = "waiting_for_config"
    return state


# ═══════════════════════════════════════════════════════════════════════════
# Phase 2 — 提取标签并合并进 benchmark
# ═══════════════════════════════════════════════════════════════════════════

def phase_extract(state: dict) -> dict:
    _section("Phase 2 / 6 — 提取细胞类型标签（extract_benchmark_labels）")

    config_path = state.get("dataset_config", "")
    h5ad_dir    = state["h5ad_dir"]

    if not config_path or not os.path.exists(config_path):
        print(f"  ❌ 找不到 dataset_config: {config_path}")
        sys.exit(1)

    # 验证用户已填写 selected_column
    cfg_df = pd.read_csv(config_path, encoding="utf-8-sig")
    filled = (
        cfg_df["selected_column"]
        .dropna().astype(str).str.strip()
        .replace("", pd.NA).dropna()
    )
    if filled.empty:
        print("  ❌ selected_column 列全部为空！请先填写后重新运行。")
        sys.exit(1)

    # ── 修复 load_config 的 keep="first" 问题 ──────────────────────────
    # load_config 只读每个文件的第一行，若用户把 selected_column 填在
    # 非第一行（例如填在对应候选列的那一行），原函数会丢失填写内容。
    # 这里预处理：对每个文件，把任意行填写的 selected_column 值
    # 归并到该文件的第一行，再传给 load_config。
    if "selected_column" in cfg_df.columns:
        if "selected_subtype_column" not in cfg_df.columns:
            cfg_df["selected_subtype_column"] = ""
        cfg_df["selected_column"] = cfg_df["selected_column"].fillna("").astype(str).str.strip()
        cfg_df["selected_subtype_column"] = cfg_df["selected_subtype_column"].fillna("").astype(str).str.strip()

        # 对每个文件，取第一个非空的 selected_column / selected_subtype_column
        def _first_nonempty(series):
            nonempty = series[series != ""]
            return nonempty.iloc[0] if not nonempty.empty else ""

        file_main    = cfg_df.groupby("file")["selected_column"].apply(_first_nonempty)
        file_subtype = cfg_df.groupby("file")["selected_subtype_column"].apply(_first_nonempty)

        # 写回到每个文件的第一行，其余行清空，供 load_config 正确读取
        for fname in cfg_df["file"].unique():
            idx_all   = cfg_df.index[cfg_df["file"] == fname]
            first_idx = idx_all[0]
            cfg_df.at[first_idx, "selected_column"]         = file_main.get(fname, "")
            cfg_df.at[first_idx, "selected_subtype_column"] = file_subtype.get(fname, "")
            for i in idx_all[1:]:
                cfg_df.at[i, "selected_column"]         = ""
                cfg_df.at[i, "selected_subtype_column"] = ""

        fixed_config_path = config_path + ".fixed.tmp"
        cfg_df.to_csv(fixed_config_path, index=False, encoding="utf-8-sig")
        config_path = fixed_config_path

    # 提取标签
    config    = load_config(config_path)
    all_paths = collect_h5ad_paths(h5ad_dir)
    path_map  = {os.path.basename(p): p for p in all_paths}

    new_dfs = []
    for fname, cols in config.items():
        if fname not in path_map:
            print(f"  [未找到文件] {fname}")
            continue
        main_col    = cols.get("main", "").strip()
        subtype_col = cols.get("subtype", "").strip()
        if main_col:
            df = extract_from_one(path_map[fname], main_col, "main")
            if df is not None and not df.empty:
                new_dfs.append(df)
        if subtype_col and subtype_col != main_col:
            df = extract_from_one(path_map[fname], subtype_col, "subtype")
            if df is not None and not df.empty:
                new_dfs.append(df)

    # 清理临时文件
    if config_path.endswith(".fixed.tmp") and os.path.exists(config_path):
        os.remove(config_path)

    if not new_dfs:
        print("  ❌ 未提取到任何标签，请检查 selected_column 填写是否正确。")
        print("  提示：selected_column 可填在该文件任意一行（不必须是第一行）。")
        sys.exit(1)

    new_bench = pd.concat(new_dfs, ignore_index=True)
    n_new_labels = int(new_bench["label_text"].nunique())
    n_new_with_gt = int(new_bench["has_ground_truth"].sum())
    print(f"  新提取唯一标签数: {n_new_labels}  （已有 CL ID: {n_new_with_gt}）")

    # 与现有 benchmark 合并
    if os.path.exists(BENCHMARK_PATH):
        _backup(BENCHMARK_PATH)
        existing = pd.read_csv(BENCHMARK_PATH, encoding="utf-8-sig")
        combined = pd.concat([existing, new_bench], ignore_index=True)
    else:
        combined = new_bench

    # 去重
    before_n = len(combined)
    combined = combined.drop_duplicates(
        subset=["source_file", "label_text", "label_level"], keep="first"
    ).reset_index(drop=True)
    after_n = len(combined)
    print(f"  去重: {before_n} → {after_n} 条记录")

    combined.to_csv(BENCHMARK_PATH, index=False, encoding="utf-8-sig")
    print(f"  ✅ benchmark_labels.csv 更新完成（共 {after_n} 条）")

    state["n_new_labels"] = n_new_labels

    # 若新提取的标签中有无 ground truth 的（has_ground_truth=False），
    # 先走 resolve_cl_label 尝试通过 CL label 列补充 cl_id；
    # 否则直接进入 eval_before。
    n_no_gt = int((new_bench["has_ground_truth"] == False).sum())
    if n_no_gt > 0:
        print(f"  ℹ️  有 {n_no_gt} 条记录无 CL ID，将尝试通过 CL label 列自动补充映射")
        state["current_phase"] = "resolve_cl_label"
    else:
        state["current_phase"] = "eval_before"
    return state


# ═══════════════════════════════════════════════════════════════════════════
# Phase [可选] — 通过 CL label 反向补充 lookup table
# ═══════════════════════════════════════════════════════════════════════════

def phase_resolve_cl_label(state: dict) -> dict:
    """
    对没有 CL ID 但有 CL label 列的 h5ad 数据集，通过：
      label_text (free_annotation) —obs对齐→ cl_label (cell_ontology_class)
      cl_label —反查 cell_type_mapping.csv→ cl_id
    将解析成功的 (label_text, cl_id, cl_label) 追加进 unmapped_labels.csv，
    后续直接执行 merge 阶段即可写入 lookup table。
    """
    _section("Phase [可选] — 通过 CL label 反向解析 cl_id（resolve_cl_label）")

    config_path = state.get("dataset_config", "")
    h5ad_dir    = state["h5ad_dir"]
    run_id      = state["run_id"]

    if not config_path or not os.path.exists(config_path):
        print(f"  ❌ 找不到 dataset_config: {config_path}")
        sys.exit(1)

    try:
        import anndata as ad
    except ImportError:
        print("  ❌ 需要安装 anndata: pip install anndata")
        sys.exit(1)

    from extract_benchmark_labels import CL_LABEL_COLUMNS, find_column
    from difflib import SequenceMatcher

    # ── 构建 lookup table 的 cl_label → (cl_id, cl_label) 索引 ──────────
    if not os.path.exists(LOOKUP_TABLE_PATH):
        print(f"  ❌ 找不到 lookup table: {LOOKUP_TABLE_PATH}")
        sys.exit(1)

    lookup_df = pd.read_csv(LOOKUP_TABLE_PATH, encoding="utf-8-sig")

    def _norm(s):
        return str(s).strip().lower()

    cl_label_index = {}
    for _, row in lookup_df.iterrows():
        key = _norm(row.get("cl_label", ""))
        if key and key not in ("nan", ""):
            cl_label_index[key] = (str(row["cl_id"]).strip(), str(row["cl_label"]).strip())

    print(f"  lookup table 已索引 {len(cl_label_index)} 个 cl_label 条目")

    # ── 读取 dataset_config，找 label_col ────────────────────────────────
    cfg_df = pd.read_csv(config_path, encoding="utf-8-sig")
    cfg_df["selected_column"] = cfg_df["selected_column"].fillna("").astype(str).str.strip()

    def _first_nonempty(series):
        nonempty = series[series != ""]
        return nonempty.iloc[0] if not nonempty.empty else ""

    file_main = cfg_df.groupby("file")["selected_column"].apply(_first_nonempty)

    all_paths = collect_h5ad_paths(h5ad_dir)
    path_map  = {os.path.basename(p): p for p in all_paths}

    resolved_rows   = []
    unresolved_rows = []

    for fname, label_col in file_main.items():
        if not label_col:
            continue
        if fname not in path_map:
            print(f"  [未找到文件] {fname}")
            continue

        print(f"\n  处理: {fname}  (label列: {label_col})")

        try:
            adata = ad.read_h5ad(path_map[fname], backed="r")
        except Exception as e:
            print(f"    ❌ 无法读取: {e}")
            continue

        obs_cols     = list(adata.obs.columns)
        cl_label_col = find_column(obs_cols, CL_LABEL_COLUMNS)

        # 避免把 label_col 自身当作 cl_label_col
        if cl_label_col == label_col:
            cl_label_col = None

        if cl_label_col is None:
            print(f"    ℹ️  未找到 CL label 列（候选: {CL_LABEL_COLUMNS}），跳过")
            if hasattr(adata, "file"):
                adata.file.close()
            continue

        print(f"    CL label 列: {cl_label_col}")

        try:
            pair_df = adata.obs[[label_col, cl_label_col]].copy()
        except Exception as e:
            print(f"    ❌ 读取 obs 失败: {e}")
            if hasattr(adata, "file"):
                adata.file.close()
            continue

        if hasattr(adata, "file"):
            adata.file.close()

        pair_df.columns      = ["label_text", "cl_label_raw"]
        pair_df["label_text"]   = pair_df["label_text"].astype(str).str.strip()
        pair_df["cl_label_raw"] = pair_df["cl_label_raw"].astype(str).str.strip()
        pair_df = pair_df[pair_df["label_text"].ne("") & pair_df["cl_label_raw"].ne("")].copy()

        # 每个 label_text 取出现频率最高的 cl_label
        label_to_cllabel = (
            pair_df.groupby("label_text")["cl_label_raw"]
            .agg(lambda x: x.mode().iloc[0])
            .reset_index()
        )

        ok = 0
        fail = 0
        for _, row in label_to_cllabel.iterrows():
            label_text   = row["label_text"]
            cl_label_raw = row["cl_label_raw"]
            key          = _norm(cl_label_raw)

            if key in cl_label_index:
                cl_id, cl_label_canonical = cl_label_index[key]
                resolved_rows.append({
                    "label_text":     label_text,
                    "cl_id":          cl_id,
                    "cl_label":       cl_label_canonical,
                    "source_file":    fname,
                    "resolve_method": f"exact: {cl_label_raw}",
                })
                ok += 1
            else:
                # 模糊匹配 cl_label（阈值 0.85）
                best_score, best_key = 0.0, None
                for k in cl_label_index:
                    s = SequenceMatcher(None, key, k).ratio()
                    if s > best_score:
                        best_score, best_key = s, k

                if best_score >= 0.85 and best_key:
                    cl_id, cl_label_canonical = cl_label_index[best_key]
                    resolved_rows.append({
                        "label_text":     label_text,
                        "cl_id":          cl_id,
                        "cl_label":       cl_label_canonical,
                        "source_file":    fname,
                        "resolve_method": f"fuzzy({best_score:.2f}): {cl_label_raw}→{best_key}",
                    })
                    ok += 1
                else:
                    unresolved_rows.append({
                        "label_text":  label_text,
                        "cl_label_raw": cl_label_raw,
                        "source_file": fname,
                        "best_match":  best_key or "",
                        "best_score":  round(best_score, 3),
                    })
                    fail += 1

        print(f"    ✅ 解析成功: {ok}   ❓ 未解析: {fail}")

    # ── 保存无法解析的标签（供用户手动补充）────────────────────────────────
    if unresolved_rows:
        unresolved_path = os.path.join(_ROOT, "data", "labels", f"cl_label_unresolved_{run_id}.csv")
        pd.DataFrame(unresolved_rows).to_csv(unresolved_path, index=False, encoding="utf-8-sig")
        print(f"\n  ⚠️  {len(unresolved_rows)} 个标签 cl_label 未能匹配到 lookup table")
        print(f"  已写入供参考: {unresolved_path}")
        print("  → 可在该文件手动填写 cl_id 后，把内容复制追加到 unmapped_labels.csv 再执行 merge")

    if not resolved_rows:
        print("\n  ⚠️  没有成功解析任何映射，直接进入 eval_before。")
        state["current_phase"] = "eval_before"
        return state

    # ── 追加解析结果到 unmapped_labels.csv ──────────────────────────────
    resolved_df = pd.DataFrame(resolved_rows)[["label_text", "cl_id", "cl_label"]]

    if os.path.exists(UNMAPPED_OUTPUT):
        existing = pd.read_csv(UNMAPPED_OUTPUT, encoding="utf-8-sig")
        combined = pd.concat([existing, resolved_df], ignore_index=True)
    else:
        combined = resolved_df

    combined = combined.drop_duplicates(
        subset=["label_text", "cl_id"], keep="first"
    ).reset_index(drop=True)
    combined.to_csv(UNMAPPED_OUTPUT, index=False, encoding="utf-8-sig")

    print(f"\n  ✅ 通过 cl_label 解析出 {len(resolved_df)} 条映射")
    print(f"  已追加至 unmapped_labels.csv（共 {len(combined)} 条，含历史）")
    print("  下一步将自动执行 merge，把这些映射写入 lookup table")

    state["resolved_via_cl_label"] = len(resolved_df)
    state["current_phase"] = "merge"
    return state


# ═══════════════════════════════════════════════════════════════════════════
# Phase 3 — 更新前基线评估
# ═══════════════════════════════════════════════════════════════════════════

def phase_eval_before(state: dict) -> dict:
    _section("Phase 3 / 6 — 更新前基线评估")

    run_id = state["run_id"]

    bench = pd.read_csv(BENCHMARK_PATH, encoding="utf-8-sig")
    bench["cl_id"] = bench["cl_id"].apply(normalize_cl_id)

    # 只取有 ground truth 的去重标签
    bench_gt     = bench[bench["cl_id"] != ""].copy()
    bench_unique = bench_gt.drop_duplicates(
        subset=["label_text"], keep="first"
    ).reset_index(drop=True)
    labels = bench_unique["label_text"].tolist()

    print(f"  有 ground truth 的唯一标签数: {len(labels)}")
    print(f"  正在运行本地映射器（这可能需要几十秒）...")

    preds    = run_local_mapper(labels)
    pred_df  = pd.DataFrame(preds)
    merge_cols = ["label_text", "cl_id", "cl_label"]
    if "cell_count" in bench_gt.columns:
        merge_cols.append("cell_count")
    merged = bench_gt[merge_cols].merge(pred_df, on="label_text", how="left")
    merged["is_correct"] = merged.apply(
        lambda r: is_correct(r["pred_cl_id"], r["cl_id"]), axis=1
    )

    eval_before_path = os.path.join(EVAL_DIR, f"eval_before_{run_id}.csv")
    merged.to_csv(eval_before_path, index=False, encoding="utf-8-sig")

    metrics = compute_metrics(merged, "pred_cl_id", "cl_id")
    wnote = "细胞加权" if metrics.get("weighted_by_cells") else "按行"
    print(f"\n  ── 基线指标（更新 lookup table 前）──")
    print(f"  分母（{wnote}）: {metrics['total']}")
    print(f"  Coverage   : {metrics['coverage']}%   （预测出 CL ID 的比例）")
    print(f"  Accuracy   : {metrics['accuracy']}%   （正确映射数 / 分母）")
    print(f"  Precision  : {metrics['precision']}%  （正确 / 预测非空）")

    # ── 生成 unmapped_labels.csv ───────────────────────────────────────────
    # 1) 有 ground truth 但预测为空（按唯一标签列出）
    unmapped_gt = merged[
        merged["pred_cl_id"].apply(normalize_cl_id) == ""
    ][["label_text", "cl_id", "cl_label"]].drop_duplicates(subset=["label_text"]).copy()

    # 2) 无 ground truth 的标签（供参考填写）
    bench_no_gt = (
        bench[bench["cl_id"] == ""]
        .drop_duplicates(subset=["label_text"], keep="first")
        [["label_text"]]
        .copy()
    )
    bench_no_gt["cl_id"]    = ""
    bench_no_gt["cl_label"] = ""

    all_unmapped = pd.concat([unmapped_gt, bench_no_gt], ignore_index=True).drop_duplicates(
        subset=["label_text"], keep="first"
    )
    all_unmapped.to_csv(UNMAPPED_OUTPUT, index=False, encoding="utf-8-sig")

    print(f"\n  unmapped_labels.csv 写入: {UNMAPPED_OUTPUT}")
    print(f"    - 有 ground truth 但未映射: {len(unmapped_gt)} 条")
    print(f"    - 无 ground truth（可选填）: {len(bench_no_gt)} 条")

    state["eval_before"]      = metrics
    state["eval_before_path"] = eval_before_path
    state["current_phase"]    = "waiting_for_annotations"
    return state


# ═══════════════════════════════════════════════════════════════════════════
# Phase 4 — 合并标注数据到 lookup table
# ═══════════════════════════════════════════════════════════════════════════

def phase_merge(state: dict) -> dict:
    _section("Phase 4 / 6 — 合并标注到 Lookup Table（merge_lookup_table）")

    if not os.path.exists(UNMAPPED_OUTPUT):
        print(f"  ⚠️  未找到 {UNMAPPED_OUTPUT}，跳过合并。")
        state["current_phase"] = "eval_after"
        return state

    df = pd.read_csv(UNMAPPED_OUTPUT, encoding="utf-8-sig")
    annotated = df[
        df["cl_id"].notna()
        & (df["cl_id"].astype(str).str.strip() != "")
        & (df["cl_id"].astype(str).str.strip().str.lower() != "nan")
    ]

    if annotated.empty:
        print("  ℹ️  unmapped_labels.csv 中没有填写 cl_id 的行，跳过合并。")
        print("     （之后可手动填写后重新运行：--start_phase merge）")
        state["merged_count"]  = 0
        state["current_phase"] = "eval_after"
        return state

    print(f"  检测到 {len(annotated)} 条已标注记录，开始合并...")
    _backup(LOOKUP_TABLE_PATH)
    _merge_lookup(UNMAPPED_OUTPUT, LOOKUP_TABLE_PATH, LOOKUP_TABLE_PATH)

    state["merged_count"]  = int(len(annotated))
    state["current_phase"] = "eval_after"
    return state


# ═══════════════════════════════════════════════════════════════════════════
# Phase 5 — 更新后评估（对比）
# ═══════════════════════════════════════════════════════════════════════════

def phase_eval_after(state: dict) -> dict:
    _section("Phase 5 / 6 — 更新后评估（前后对比）")

    run_id = state["run_id"]

    bench = pd.read_csv(BENCHMARK_PATH, encoding="utf-8-sig")
    bench["cl_id"] = bench["cl_id"].apply(normalize_cl_id)
    bench_gt     = bench[bench["cl_id"] != ""].copy()
    bench_unique = bench_gt.drop_duplicates(
        subset=["label_text"], keep="first"
    ).reset_index(drop=True)
    labels = bench_unique["label_text"].tolist()

    print(f"  正在运行本地映射器（重新加载最新 lookup table）...")

    # 通过创建新实例保证读到最新的 cell_type_mapping.csv
    # run_local_mapper 每次调用都新建 MappingService 实例，无需手动重置缓存
    preds   = run_local_mapper(labels)
    pred_df = pd.DataFrame(preds)
    merge_cols = ["label_text", "cl_id", "cl_label"]
    if "cell_count" in bench_gt.columns:
        merge_cols.append("cell_count")
    merged = bench_gt[merge_cols].merge(pred_df, on="label_text", how="left")
    merged["is_correct"] = merged.apply(
        lambda r: is_correct(r["pred_cl_id"], r["cl_id"]), axis=1
    )

    eval_after_path = os.path.join(EVAL_DIR, f"eval_after_{run_id}.csv")
    merged.to_csv(eval_after_path, index=False, encoding="utf-8-sig")
    metrics_after  = compute_metrics(merged, "pred_cl_id", "cl_id")
    metrics_before = state.get("eval_before", {})

    # ── 打印对比表 ────────────────────────────────────────────────────────
    print()
    print("  ╔══════════════════════════════════════════════════╗")
    print("  ║         更新前  vs  更新后  —  指标对比           ║")
    print("  ╠══════════════════════════════════════════════════╣")
    print("  ║  指标       │   更新前   │   更新后   │   变化    ║")
    print("  ╠══════════════════════════════════════════════════╣")
    for key, label in [
        ("coverage",  "Coverage "),
        ("accuracy",  "Accuracy "),
        ("precision", "Precision"),
    ]:
        before = metrics_before.get(key, 0.0)
        after  = metrics_after.get(key, 0.0)
        delta  = after - before
        sign   = "+" if delta >= 0 else ""
        print(f"  ║  {label}  │  {before:6.2f}%   │  {after:6.2f}%   │  {sign}{delta:.2f}%  ║")
    print("  ╠══════════════════════════════════════════════════╣")
    print(f"  ║  分母(细胞等) │  {metrics_before.get('total', 0):>8}   │  {metrics_after['total']:>8}   │           ║")
    print("  ╚══════════════════════════════════════════════════╝")

    state["eval_after"]      = metrics_after
    state["eval_after_path"] = eval_after_path
    state["current_phase"]   = "coverage"
    return state


# ═══════════════════════════════════════════════════════════════════════════
# Phase 6 — 器官覆盖统计
# ═══════════════════════════════════════════════════════════════════════════

def phase_coverage(state: dict) -> dict:
    _section("Phase 6 / 6 — 器官/组织覆盖统计（summarize_benchmark_coverage）")

    run_id   = state["run_id"]
    h5ad_dir = state["h5ad_dir"]

    bench = pd.read_csv(BENCHMARK_PATH, encoding="utf-8-sig")
    files = bench["source_file"].dropna().unique().tolist()

    path_map    = load_file_paths(h5ad_dir, None)
    label_counts = bench.groupby("source_file")["label_text"].nunique().to_dict()
    cell_counts  = bench.groupby("source_file")["cell_count"].sum().to_dict()

    rows = []
    for fname in files:
        full_path = path_map.get(fname, "")
        if full_path:
            tissue, _ = infer_tissue_from_h5ad(full_path)
        else:
            tissue = "UNKNOWN"
        if tissue == "UNKNOWN":
            tissue = infer_tissue_from_filename(fname)
        rows.append({
            "source_file": fname,
            "tissue":      tissue,
            "label_count": int(label_counts.get(fname, 0)),
            "cell_count":  int(cell_counts.get(fname, 0)),
        })

    detail  = pd.DataFrame(rows)
    summary = (
        detail.groupby("tissue")
        .agg(datasets=("source_file", "nunique"),
             labels=("label_count", "sum"),
             cells=("cell_count", "sum"))
        .reset_index()
        .sort_values("datasets", ascending=False)
    )

    coverage_path = os.path.join(COVERAGE_DIR, f"tissue_coverage_{run_id}.csv")
    summary.to_csv(coverage_path, index=False, encoding="utf-8-sig")

    print(f"\n  器官覆盖摘要（共 {len(summary)} 个器官/组织，Top 10）：")
    for _, r in summary.head(10).iterrows():
        print(f"    {r['tissue']:<22}  {r['datasets']:>3} 数据集   {r['labels']:>5} 标签")
    if len(summary) > 10:
        print(f"    ... 其余 {len(summary)-10} 个组织详见 {os.path.basename(coverage_path)}")

    state["coverage_path"] = coverage_path
    state["current_phase"] = "done"
    return state


# ═══════════════════════════════════════════════════════════════════════════
# 保存 update log（供论文引用）
# ═══════════════════════════════════════════════════════════════════════════

def save_update_log(state: dict):
    run_id   = state["run_id"]
    log_path = os.path.join(LOG_DIR, f"update_log_{run_id}.json")

    eb = state.get("eval_before", {})
    ea = state.get("eval_after",  {})

    log = {
        "run_id":                   run_id,
        "timestamp":                state.get("timestamp", ""),
        "h5ad_dir":                 state.get("h5ad_dir", ""),
        "new_datasets_added":       state.get("new_files", []),
        "new_labels_extracted":     state.get("n_new_labels", 0),
        "labels_merged_to_lookup":  state.get("merged_count", 0),
        "eval_before": {
            "total":     eb.get("total", 0),
            "coverage":  eb.get("coverage", 0),
            "accuracy":  eb.get("accuracy", 0),
            "precision": eb.get("precision", 0),
        },
        "eval_after": {
            "total":     ea.get("total", 0),
            "coverage":  ea.get("coverage", 0),
            "accuracy":  ea.get("accuracy", 0),
            "precision": ea.get("precision", 0),
        },
        "delta_accuracy":   round(ea.get("accuracy", 0)  - eb.get("accuracy", 0),  2),
        "delta_coverage":   round(ea.get("coverage", 0)  - eb.get("coverage", 0),  2),
        "delta_precision":  round(ea.get("precision", 0) - eb.get("precision", 0), 2),
        "output_files": {
            "eval_before":   state.get("eval_before_path", ""),
            "eval_after":    state.get("eval_after_path", ""),
            "coverage":      state.get("coverage_path", ""),
        },
    }

    with open(log_path, "w", encoding="utf-8") as f:
        json.dump(log, f, ensure_ascii=False, indent=2)

    _section("✅  全流程完成")
    print(f"  run_id         : {run_id}")
    print(f"  新增数据集     : {len(log['new_datasets_added'])}")
    print(f"  新增标签数     : {log['new_labels_extracted']}")
    print(f"  合并至 lookup  : {log['labels_merged_to_lookup']} 条")
    if eb and ea:
        print(f"  Accuracy 变化  : {eb['accuracy']}%  →  {ea['accuracy']}%"
              f"  ({'+' if log['delta_accuracy'] >= 0 else ''}{log['delta_accuracy']}%)")
    print(f"  更新日志       : {log_path}")
    print()


# ═══════════════════════════════════════════════════════════════════════════
# 主流程
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Cell Type Mapping 增量更新 SOP",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python incremental_update_sop.py --h5ad_dir ./h5ad_data\n"
            "  python incremental_update_sop.py --resume sop_state_20260401_143022.json\n"
            "  python incremental_update_sop.py --h5ad_dir ./h5ad_data --start_phase merge\n"
        ),
    )
    parser.add_argument("--h5ad_dir", "-d",
                        help="包含 .h5ad 文件的目录（新建运行时必填）")
    parser.add_argument("--resume", "-r",
                        help="状态文件路径，用于断点续跑")
    parser.add_argument("--start_phase", "-p",
                        choices=PHASES,
                        help="直接跳到指定阶段开始")
    args = parser.parse_args()

    # ── 加载或初始化状态 ────────────────────────────────────────────────
    if args.resume:
        if not os.path.exists(args.resume):
            print(f"❌ 状态文件不存在: {args.resume}")
            sys.exit(1)
        state      = _load_state(args.resume)
        state_path = args.resume
        print(f"\n▶  续跑模式，读取状态: {args.resume}")
        print(f"   上次停在: {state.get('current_phase', '?')}")
    elif args.h5ad_dir:
        if not os.path.isdir(args.h5ad_dir):
            print(f"❌ 目录不存在: {args.h5ad_dir}")
            sys.exit(1)
        run_id     = _make_run_id()
        state_path = os.path.join(LOG_DIR, f"sop_state_{run_id}.json")
        state = {
            "run_id":        run_id,
            "timestamp":     datetime.now().isoformat(),
            "h5ad_dir":      str(Path(args.h5ad_dir).resolve()),
            "current_phase": "scan",
        }
        print(f"\n▶  新建运行  run_id = {run_id}")
    else:
        parser.print_help()
        sys.exit(1)

    # 手动指定起始阶段
    if args.start_phase:
        state["current_phase"] = args.start_phase
        print(f"   跳转至阶段: {args.start_phase}")

    # ── 阶段循环 ────────────────────────────────────────────────────────
    while state["current_phase"] != "done":
        phase = state["current_phase"]

        # ── Phase 1: scan ────────────────────────────────────────────
        if phase == "scan":
            state = phase_scan(state)
            _save_state(state, state_path)

            if state["current_phase"] == "waiting_for_config":
                cfg = state["dataset_config"]
                _pause(
                    "填写 dataset_config.csv 的 selected_column 列",
                    [
                        f"1. 用 Excel/VSCode 打开: {cfg}",
                        "2. 查看 'examples' 列，找到包含细胞类型名称的候选列",
                        "   （如 T cell、neuron、macrophage 等文本标签的列）",
                        "3. 在同一行的 'selected_column' 列填入该列名",
                        "4. 每个文件只填一个主列（subtype 列可选，新增 selected_subtype_column）",
                        "5. 保存文件",
                    ],
                )
                state["current_phase"] = "extract"
                _save_state(state, state_path)

        # ── Phase 2: extract ─────────────────────────────────────────
        elif phase == "extract":
            state = phase_extract(state)
            _save_state(state, state_path)

        # ── Phase [可选]: resolve_cl_label ───────────────────────────
        elif phase == "resolve_cl_label":
            state = phase_resolve_cl_label(state)
            _save_state(state, state_path)

        # ── Phase 3: eval_before ─────────────────────────────────────
        elif phase == "eval_before":
            state = phase_eval_before(state)
            _save_state(state, state_path)

            if state["current_phase"] == "waiting_for_annotations":
                # 检查 unmapped_labels.csv 是否已有内容（由 resolve_cl_label 自动填充，
                # 或用户手动填写）。有内容则自动进入 merge；否则暂停等待用户操作。
                has_unmapped_content = False
                if os.path.exists(UNMAPPED_OUTPUT):
                    try:
                        _um = pd.read_csv(UNMAPPED_OUTPUT, encoding="utf-8-sig")
                        has_unmapped_content = not _um.empty and _um["cl_id"].notna().any()
                    except Exception:
                        pass

                if has_unmapped_content:
                    print()
                    print("  ✅ 检测到 unmapped_labels.csv 已有内容，自动进入 merge 阶段")
                    state["current_phase"] = "merge"
                else:
                    print()
                    print("  ─────────────────────────────────────────────────────────")
                    print("  可选操作：为 unmapped 标签补充 CL ID（可提升映射准确率）")
                    print("  ─────────────────────────────────────────────────────────")
                    print(f"  打开: {UNMAPPED_OUTPUT}")
                    print("  在 cl_id 列填入 CL:XXXXXXX 格式的 ID，cl_label 填对应名称。")
                    print("  不填也没关系，直接回车跳过（本轮不更新 lookup table）。")
                    print()
                    ans = input("  是否已填写 unmapped_labels.csv？[y / 回车跳过] ").strip().lower()
                    if ans == "y":
                        state["current_phase"] = "merge"
                    else:
                        print("  ℹ️  跳过合并，eval_after 与 eval_before 指标相同。")
                        state["merged_count"]  = 0
                        state["current_phase"] = "eval_after"
                _save_state(state, state_path)

        # ── Phase 4: merge ───────────────────────────────────────────
        elif phase == "merge":
            state = phase_merge(state)
            _save_state(state, state_path)

        # ── Phase 5: eval_after ──────────────────────────────────────
        elif phase == "eval_after":
            state = phase_eval_after(state)
            _save_state(state, state_path)

        # ── Phase 6: coverage ────────────────────────────────────────
        elif phase == "coverage":
            state = phase_coverage(state)
            _save_state(state, state_path)

        else:
            print(f"❌ 未知阶段: {phase}，请检查状态文件。")
            sys.exit(1)

    save_update_log(state)


if __name__ == "__main__":
    main()
