"""
Cell Type Mapping 评估脚本

对比以下两种方法在 benchmark_labels.csv 上的表现：
  - Method A: 本平台的 MappingService（本地映射表 + SequenceMatcher 模糊匹配）
  - Method B: CellOntologyMapper（OmicVerse，NLP embedding + LLM）

用法:
  # 只跑 Method A（不需要安装 omicverse）
  python evaluate_mapping.py --benchmark data/benchmark/benchmark_labels.csv --method local

  # 只跑 Method B（需要 pip install omicverse）
  python evaluate_mapping.py --benchmark data/benchmark/benchmark_labels.csv --method cellontology

  # 两个都跑并对比
  python evaluate_mapping.py --benchmark data/benchmark/benchmark_labels.csv --method both

  # 指定只评估 main 层级标签
  python evaluate_mapping.py --benchmark data/benchmark/benchmark_labels.csv --method both --level main

输出:
  - 控制台打印汇总指标（Accuracy/Coverage/Precision 按 cell_count 对细胞加权，与 export 侧细胞级比例一致）
  - eval_results_<method>.csv：每条 benchmark 行的详细预测（含 cell_count 时一行对应一类细胞的汇总条数）
  - eval_summary.csv：两个方法的汇总对比表
"""
import os
import sys
import argparse
import time
import pandas as pd
from typing import Optional

# ─────────────────────────────────────────────
# 评估核心逻辑
# ─────────────────────────────────────────────

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SRC_DIR = os.path.join(_ROOT, "src")
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

INVALID_CL = {"", "nan", "not found", "none"}


def normalize_cl_id(cl_id) -> str:
    """统一 CL ID 格式，去除空白，转大写冒号前缀。"""
    if cl_id is None or (isinstance(cl_id, float)):
        return ""
    s = str(cl_id).strip()
    if s.lower() in INVALID_CL:
        return ""
    return s


def is_correct(pred_cl: str, true_cl: str) -> bool:
    """判断预测是否正确（精确匹配 CL ID）。"""
    p = normalize_cl_id(pred_cl)
    t = normalize_cl_id(true_cl)
    if not p or not t:
        return False
    return p == t


def compute_metrics(
    df: pd.DataFrame,
    pred_col: str,
    true_col: str = "cl_id",
    weight_col: Optional[str] = "cell_count",
) -> dict:
    """
    同时计算两套指标，供论文对比使用：

    Label-level（主指标，推荐写论文用）：
      每个唯一 label_text 权重相同（=1），不受高频类型主导。
      Accuracy_label  = 预测正确的唯一标签数 / 唯一标签总数
      Coverage_label  = 有预测的唯一标签数   / 唯一标签总数
      Precision_label = 预测正确的唯一标签数 / 有预测的唯一标签数

    Cell-weighted（辅助指标，反映实际数据覆盖质量）：
      按 cell_count 加权，若无 cell_count 列则退化为行级（与 label-level 相同）。
      Accuracy_cell   = 预测正确的细胞数 / 总细胞数
      Coverage_cell   = 有预测的细胞数   / 总细胞数
      Precision_cell  = 预测正确的细胞数 / 有预测的细胞数

    返回字段：
      label_total / label_predicted / label_correct
      coverage_label / accuracy_label / precision_label
      cell_total / cell_predicted / cell_correct
      coverage_cell / accuracy_cell / precision_cell
      weighted_by_cells   (True 表示 cell_count 列有效)
      # 向后兼容：total/predicted/correct/coverage/accuracy/precision 指向 label-level
    """
    pred_ok = df[pred_col].apply(normalize_cl_id).ne("")
    correct_ok = df.apply(lambda r: is_correct(r[pred_col], r[true_col]), axis=1)

    # ── Label-level：对 label_text 去重后计算（每种类型权重=1）
    if "label_text" in df.columns:
        dedup = df.drop_duplicates(subset=["label_text"], keep="first")
        pred_ok_lbl = dedup[pred_col].apply(normalize_cl_id).ne("")
        correct_ok_lbl = dedup.apply(lambda r: is_correct(r[pred_col], r[true_col]), axis=1)
        lbl_total = len(dedup)
        lbl_pred = int(pred_ok_lbl.sum())
        lbl_corr = int(correct_ok_lbl.sum())
    else:
        # 无 label_text 列时退化为行级
        lbl_total = len(df)
        lbl_pred = int(pred_ok.sum())
        lbl_corr = int(correct_ok.sum())

    cov_lbl = lbl_pred / lbl_total if lbl_total else 0.0
    acc_lbl = lbl_corr / lbl_total if lbl_total else 0.0
    prec_lbl = lbl_corr / lbl_pred if lbl_pred else 0.0

    # ── Cell-weighted：按 cell_count 加权
    has_cell_weight = False
    if weight_col and weight_col in df.columns:
        w = pd.to_numeric(df[weight_col], errors="coerce").fillna(0.0)
        if float(w.sum()) > 0:
            has_cell_weight = True
        else:
            w = pd.Series(1.0, index=df.index)
    else:
        w = pd.Series(1.0, index=df.index)

    total_w = float(w.sum())
    pred_w = float(w[pred_ok].sum())
    corr_w = float(w[correct_ok].sum())

    cov_cell = pred_w / total_w if total_w else 0.0
    acc_cell = corr_w / total_w if total_w else 0.0
    prec_cell = corr_w / pred_w if pred_w else 0.0

    return {
        # Label-level（主指标）
        "label_total": lbl_total,
        "label_predicted": lbl_pred,
        "label_correct": lbl_corr,
        "coverage_label": round(cov_lbl * 100, 2),
        "accuracy_label": round(acc_lbl * 100, 2),
        "precision_label": round(prec_lbl * 100, 2),
        # Cell-weighted（辅助指标）
        "cell_total": int(round(total_w)),
        "cell_predicted": int(round(pred_w)),
        "cell_correct": int(round(corr_w)),
        "coverage_cell": round(cov_cell * 100, 2),
        "accuracy_cell": round(acc_cell * 100, 2),
        "precision_cell": round(prec_cell * 100, 2),
        "weighted_by_cells": has_cell_weight,
        # 向后兼容旧调用方（指向 label-level）
        "total": lbl_total,
        "predicted": lbl_pred,
        "correct": lbl_corr,
        "coverage": round(cov_lbl * 100, 2),
        "accuracy": round(acc_lbl * 100, 2),
        "precision": round(prec_lbl * 100, 2),
    }


# ─────────────────────────────────────────────
# Method A: 本平台 MappingService
# ─────────────────────────────────────────────

def run_local_mapper(labels: list[str]) -> list[dict]:
    """调用 MappingService._match_cell_type 对每个标签做映射。"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mapping_service import MappingService

    svc = MappingService()
    svc._load_cell_type_mapping_data()

    results = []
    for label in labels:
        r = svc._match_cell_type(str(label))
        results.append({
            "label_text": label,
            "pred_cl_id": normalize_cl_id(r.get("cl_id")),
            "pred_cl_label": r.get("cl_label", ""),
            "pred_status": r.get("status", ""),
            "pred_confidence": r.get("confidence", ""),
            "pred_source": r.get("source", ""),
        })
    return results


# ─────────────────────────────────────────────
# Method B: CellOntologyMapper (OmicVerse)
# ─────────────────────────────────────────────

def run_cellontology_mapper(labels: list[str], cl_json: str = "new_ontology/cl.json",
                             model_name: str = "sentence-transformers/all-MiniLM-L6-v2",
                             model_dir: str = "./my_models") -> list[dict]:
    """
    调用 OmicVerse CellOntologyMapper 对每个标签做映射。
    需要先安装: pip install omicverse==1.7.2rc1
    """
    try:
        import omicverse as ov
    except ImportError:
        print("❌ omicverse 未安装，请先运行: pip install omicverse==1.7.2rc1")
        sys.exit(1)

    # 下载 cl.json（如果不存在）
    if not os.path.exists(cl_json):
        print(f"📥 cl.json 不存在，开始下载到 {os.path.dirname(cl_json)} ...")
        os.makedirs(os.path.dirname(cl_json), exist_ok=True)
        ov.single.download_cl(output_dir=os.path.dirname(cl_json),
                              filename=os.path.basename(cl_json))

    print(f"🔨 初始化 CellOntologyMapper (model={model_name}) ...")
    mapper = ov.single.CellOntologyMapper(
        cl_obo_file=cl_json,
        model_name=model_name,
        local_model_dir=model_dir,
    )

    results = []
    print(f"🔄 开始映射 {len(labels)} 个标签...")
    for i, label in enumerate(labels):
        if i % 50 == 0 and i > 0:
            print(f"  进度: {i}/{len(labels)}")
        try:
            # OmicVerse API: mapper.map(cell_type_name) 返回 dict
            # 参考: https://omicverse.readthedocs.io/en/latest/Tutorials-single/t_cellmatch/
            res = mapper.map(str(label))
            cl_id = normalize_cl_id(res.get("cl_id") or res.get("CL_id") or res.get("id"))
            cl_label = res.get("cl_label") or res.get("name") or res.get("label") or ""
        except Exception as e:
            cl_id = ""
            cl_label = ""
            print(f"  [警告] 映射失败: {label} — {e}")

        results.append({
            "label_text": label,
            "pred_cl_id": cl_id,
            "pred_cl_label": str(cl_label),
            "pred_status": "mapped" if cl_id else "unmapped",
            "pred_confidence": "",
            "pred_source": "cellontology_mapper",
        })
    return results


# ─────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Cell Type Mapping 评估脚本")
    default_benchmark = os.path.join(_ROOT, "data", "benchmark", "benchmark_labels.csv")
    parser.add_argument("--benchmark", "-b", default=default_benchmark,
                        help="benchmark CSV 路径 (默认 data/benchmark/benchmark_labels.csv)")
    parser.add_argument("--method", "-m", choices=["local", "cellontology", "both"],
                        default="local", help="评估方法")
    parser.add_argument("--level", "-l", choices=["main", "subtype", "all"],
                        default="all", help="只评估指定层级的标签 (默认 all)")
    default_output_dir = os.path.join(_ROOT, "results", "eval")
    parser.add_argument("--output_dir", "-o", default=default_output_dir,
                        help="结果输出目录 (默认 results/eval)")
    parser.add_argument("--cl_json", default="new_ontology/cl.json",
                        help="CellOntologyMapper 用的 cl.json 路径")
    parser.add_argument("--model", default="sentence-transformers/all-MiniLM-L6-v2",
                        help="CellOntologyMapper 用的 sentence-transformer 模型名")
    args = parser.parse_args()

    # ── 读取 benchmark ──
    bench = pd.read_csv(args.benchmark, encoding="utf-8-sig")
    bench["cl_id"] = bench["cl_id"].apply(normalize_cl_id)

    # 过滤掉没有 ground truth 的行
    bench = bench[bench["cl_id"] != ""].copy()

    if args.level != "all":
        bench = bench[bench["label_level"] == args.level].copy()

    # 映射器按唯一 label_text 调用；指标在完整 bench 行上按 cell_count 加权
    bench_unique = bench.drop_duplicates(subset=["label_text"], keep="first").reset_index(drop=True)
    labels = bench_unique["label_text"].tolist()
    if "cell_count" in bench.columns:
        cell_total = int(pd.to_numeric(bench["cell_count"], errors="coerce").fillna(0).sum())
    else:
        cell_total = len(bench)

    print(f"\n{'='*55}")
    print(f"Benchmark: {args.benchmark}")
    print(f"过滤层级: {args.level}  |  唯一标签数: {len(labels)}  |  有 GT 的细胞总数: {cell_total}")
    print(f"{'='*55}\n")

    summary_rows = []
    os.makedirs(args.output_dir, exist_ok=True)

    def run_and_save(method_name: str, pred_fn):
        t0 = time.time()
        preds = pred_fn(labels)
        elapsed = time.time() - t0

        pred_df = pd.DataFrame(preds)
        merge_cols = ["label_text", "cl_id", "cl_label"]
        if "cell_count" in bench.columns:
            merge_cols.append("cell_count")
        if "source_file" in bench.columns:
            merge_cols.append("source_file")
        if "label_level" in bench.columns:
            merge_cols.append("label_level")
        merged = bench[merge_cols].merge(pred_df, on="label_text", how="left")
        merged["is_correct"] = merged.apply(
            lambda r: is_correct(r["pred_cl_id"], r["cl_id"]), axis=1
        )

        out_path = os.path.join(args.output_dir, f"eval_results_{method_name}.csv")
        merged.to_csv(out_path, index=False, encoding="utf-8-sig")
        print(f"详细结果已写入: {out_path}")

        metrics = compute_metrics(merged, "pred_cl_id", "cl_id")
        metrics["method"] = method_name
        metrics["elapsed_sec"] = round(elapsed, 1)
        metrics["level"] = args.level

        print(f"\n── {method_name} 评估结果 ──")
        print(f"  ┌─────────────────────────────────────────────────────┐")
        print(f"  │  Label-level（主指标，推荐论文使用）                │")
        print(f"  │    唯一标签总数  : {metrics['label_total']:<6}                        │")
        print(f"  │    有预测的标签  : {metrics['label_predicted']:<6}                        │")
        print(f"  │    预测正确的标签: {metrics['label_correct']:<6}                        │")
        print(f"  │    Coverage      : {metrics['coverage_label']:>6.2f}%                      │")
        print(f"  │    Accuracy      : {metrics['accuracy_label']:>6.2f}%  (正确标签/唯一标签) │")
        print(f"  │    Precision     : {metrics['precision_label']:>6.2f}%  (正确/有预测)      │")
        print(f"  ├─────────────────────────────────────────────────────┤")
        if metrics["weighted_by_cells"]:
            print(f"  │  Cell-weighted（辅助指标，按细胞数加权）            │")
        else:
            print(f"  │  Cell-weighted（无 cell_count，退化为行级）          │")
        print(f"  │    细胞总数      : {metrics['cell_total']:<6}                        │")
        print(f"  │    有预测的细胞  : {metrics['cell_predicted']:<6}                        │")
        print(f"  │    预测正确的细胞: {metrics['cell_correct']:<6}                        │")
        print(f"  │    Coverage      : {metrics['coverage_cell']:>6.2f}%                      │")
        print(f"  │    Accuracy      : {metrics['accuracy_cell']:>6.2f}%  (正确细胞/总细胞)  │")
        print(f"  │    Precision     : {metrics['precision_cell']:>6.2f}%  (正确/有预测)      │")
        print(f"  └─────────────────────────────────────────────────────┘")
        print(f"  耗时: {metrics['elapsed_sec']}s")

        # 打印错误案例（前 20 条）
        wrong = merged[merged["is_correct"] == False].head(20)
        if not wrong.empty:
            print(f"\n  前 {len(wrong)} 条错误/未映射案例:")
            for _, row in wrong.iterrows():
                print(f"    [{row['label_text']}]  真实={row['cl_id']}  预测={row['pred_cl_id']}"
                      f"  ({row.get('pred_source', '')})")

        return metrics

    if args.method in ("local", "both"):
        m = run_and_save("local", run_local_mapper)
        summary_rows.append(m)

    if args.method in ("cellontology", "both"):
        m = run_and_save("cellontology",
                         lambda lbs: run_cellontology_mapper(lbs, args.cl_json, args.model))
        summary_rows.append(m)

    # ── 汇总对比表 ──
    if summary_rows:
        cols = [
            "method", "level",
            "label_total", "label_correct", "accuracy_label", "coverage_label", "precision_label",
            "cell_total", "cell_correct", "accuracy_cell", "coverage_cell", "precision_cell",
            "elapsed_sec",
        ]
        summary = pd.DataFrame(summary_rows)[[c for c in cols if c in pd.DataFrame(summary_rows).columns]]
        summary_path = os.path.join(args.output_dir, "eval_summary.csv")
        summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
        print(f"\n{'='*70}")
        print("汇总对比（Label-level 为主指标）:")
        print(summary.to_string(index=False))
        print(f"\n汇总已写入: {summary_path}")
        print(f"{'='*70}\n")


if __name__ == "__main__":
    main()
