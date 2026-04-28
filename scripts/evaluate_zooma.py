"""
ZOOMA benchmark 评估脚本

ZOOMA 是 EBI 提供的自动化本体注释服务，支持将 free-text 生物属性值
映射到多种生物医学本体（包括 Cell Ontology）。

API 文档: https://www.ebi.ac.uk/spot/zooma/v2/api
用法:
    python evaluate_zooma.py
    python evaluate_zooma.py --benchmark data/benchmark/benchmark_labels.csv --output results/eval/eval_results_zooma.csv
"""
import os
import sys
import time
import argparse
import requests
import pandas as pd

# 与 evaluate_mapping 共用按细胞加权的指标
from evaluate_mapping import compute_metrics, is_correct as em_is_correct

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ZOOMA_URL = "https://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate"
TIMEOUT = 15
# 每次请求之间的间隔（避免频繁请求被限速）
REQUEST_INTERVAL = 0.3

INVALID_CL = {"", "nan", "not found", "none"}


def normalize_cl_id(cl_id) -> str:
    if cl_id is None or isinstance(cl_id, float):
        return ""
    s = str(cl_id).strip()
    # ZOOMA 返回的是 URI 格式，如 http://purl.obolibrary.org/obo/CL_0000235
    # 转换为 CL:0000235
    if "CL_" in s:
        s = s.split("CL_")[-1]
        s = "CL:" + s
    if s.lower() in INVALID_CL:
        return ""
    return s


def query_zooma(cell_type: str) -> dict:
    """
    调用 ZOOMA API 将 cell type 名称映射到 CL ID。

    ZOOMA filter 说明:
      - required:[cl]    : 只接受 Cell Ontology
      - preferred:[cl]   : 优先使用 Cell Ontology
    confidence 参数过滤最低置信度（HIGH/GOOD/MEDIUM/LOW）
    """
    try:
        params = {
            "propertyValue": cell_type,
            # 只接受 CL，避免返回 NCIT 等非目标本体
            "filter": "required:[cl],preferred:[cl]",
        }
        resp = requests.get(ZOOMA_URL, params=params, timeout=TIMEOUT)

        if resp.status_code != 200:
            return {"cl_id": "", "cl_label": "", "status": "error",
                    "note": f"HTTP {resp.status_code}"}

        data = resp.json()
        if not data:
            return {"cl_id": "", "cl_label": "", "status": "unmapped", "note": "no results"}

        # 遍历结果，优先取 HIGH > GOOD > MEDIUM > LOW，且 semantic tag 含 CL_
        confidence_rank = {"HIGH": 4, "GOOD": 3, "MEDIUM": 2, "LOW": 1}
        best = None
        best_score = 0

        for item in data:
            conf = item.get("confidence", "LOW").upper()
            score = confidence_rank.get(conf, 0)
            # ZOOMA v2 返回顶层 semanticTags（字符串列表），部分旧结构在 _links.semanticTags
            tags = item.get("semanticTags", [])
            if not tags:
                tags = item.get("_links", {}).get("semanticTags", [])
            if not tags:
                # 兼容旧版 API 结构
                tags = []
                for ann in item.get("derivedFrom", {}).get("provenance", {}).get("evidence", []):
                    uri = ann.get("generatedFrom", "")
                    if uri:
                        tags.append({"href": uri})

            for tag in tags:
                href = tag.get("href", "") if isinstance(tag, dict) else str(tag)
                if "CL_" in href:
                    if score > best_score:
                        best_score = score
                        best = {
                            "cl_id": normalize_cl_id(href),
                            "cl_label": item.get("annotatedProperty", {}).get("propertyValue", ""),
                            "status": "mapped",
                            "confidence": conf,
                        }

        if best:
            return best

        # 没有 CL_tag，尝试从 annotationSummary 里找
        return {"cl_id": "", "cl_label": "", "status": "unmapped", "note": "no CL tag"}

    except requests.exceptions.Timeout:
        return {"cl_id": "", "cl_label": "", "status": "error", "note": "timeout"}
    except Exception as e:
        return {"cl_id": "", "cl_label": "", "status": "error", "note": str(e)}


def run_zooma_on_benchmark(benchmark_path: str, output_path: str):
    bench = pd.read_csv(benchmark_path, encoding="utf-8-sig")
    bench["cl_id"] = bench["cl_id"].apply(normalize_cl_id)
    bench = bench[bench["cl_id"] != ""].copy()
    bench_unique = bench.drop_duplicates(subset=["label_text"], keep="first").reset_index(drop=True)
    labels = bench_unique["label_text"].tolist()

    total = len(labels)
    print(f"\n{'='*55}")
    print(f"ZOOMA benchmark 评估")
    print(f"唯一标签数: {total}")
    print(f"{'='*55}\n")

    results = []
    t0 = time.time()

    for i, label in enumerate(labels):
        if i % 50 == 0:
            elapsed = time.time() - t0
            eta = (elapsed / (i + 1)) * (total - i - 1) if i > 0 else 0
            print(f"  [{i+1}/{total}]  已耗时 {elapsed:.0f}s  预计剩余 {eta:.0f}s")

        res = query_zooma(str(label))
        results.append({
            "label_text": label,
            "pred_cl_id": res.get("cl_id", ""),
            "pred_cl_label": res.get("cl_label", ""),
            "pred_status": res.get("status", ""),
            "pred_confidence": res.get("confidence", ""),
            "pred_source": "zooma",
            "note": res.get("note", ""),
        })
        time.sleep(REQUEST_INTERVAL)

    elapsed_total = time.time() - t0

    pred_df = pd.DataFrame(results)
    merge_cols = ["label_text", "cl_id", "cl_label"]
    if "cell_count" in bench.columns:
        merge_cols.append("cell_count")
    merged = bench[merge_cols].merge(pred_df, on="label_text", how="left")
    merged["is_correct"] = merged.apply(
        lambda r: em_is_correct(r["pred_cl_id"], r["cl_id"]),
        axis=1,
    )
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    merged.to_csv(output_path, index=False, encoding="utf-8-sig")

    metrics = compute_metrics(merged, "pred_cl_id", "cl_id")

    print(f"\n── ZOOMA 评估结果 ──")
    print(f"  唯一标签数（API 调用次数）: {total}")
    print(f"  ┌─────────────────────────────────────────────────────┐")
    print(f"  │  Label-level（主指标，推荐论文使用）                │")
    print(f"  │    唯一标签总数  : {metrics['label_total']:<6}                        │")
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
    print(f"  │    预测正确的细胞: {metrics['cell_correct']:<6}                        │")
    print(f"  │    Coverage      : {metrics['coverage_cell']:>6.2f}%                      │")
    print(f"  │    Accuracy      : {metrics['accuracy_cell']:>6.2f}%  (正确细胞/总细胞)  │")
    print(f"  │    Precision     : {metrics['precision_cell']:>6.2f}%  (正确/有预测)      │")
    print(f"  └─────────────────────────────────────────────────────┘")
    print(f"  总耗时: {elapsed_total:.1f}s")
    print(f"\n  详细结果 → {output_path}")

    # 错误案例前 20
    wrong = merged[merged["is_correct"] == False].head(20)
    if not wrong.empty:
        print(f"\n  前 {len(wrong)} 条错误/未映射:")
        for _, row in wrong.iterrows():
            print(f"    [{row['label_text']}]  真实={row['cl_id']}  "
                  f"预测={row['pred_cl_id']}  置信={row['pred_confidence']}")

    return {
        "method": "zooma",
        "label_total": metrics["label_total"],
        "label_correct": metrics["label_correct"],
        "coverage_label": metrics["coverage_label"],
        "accuracy_label": metrics["accuracy_label"],
        "precision_label": metrics["precision_label"],
        "cell_total": metrics["cell_total"],
        "cell_correct": metrics["cell_correct"],
        "coverage_cell": metrics["coverage_cell"],
        "accuracy_cell": metrics["accuracy_cell"],
        "precision_cell": metrics["precision_cell"],
        "elapsed": round(elapsed_total, 1),
        # 向后兼容
        "coverage": metrics["coverage_label"],
        "accuracy": metrics["accuracy_label"],
        "precision": metrics["precision_label"],
        "total": metrics["label_total"],
        "predicted": metrics["label_predicted"],
        "correct": metrics["label_correct"],
    }


def main():
    parser = argparse.ArgumentParser()
    default_benchmark = os.path.join(_ROOT, "data", "benchmark", "benchmark_labels.csv")
    default_output = os.path.join(_ROOT, "results", "eval", "eval_results_zooma.csv")
    parser.add_argument("--benchmark", default=default_benchmark)
    parser.add_argument("--output", default=default_output)
    args = parser.parse_args()
    run_zooma_on_benchmark(args.benchmark, args.output)


if __name__ == "__main__":
    main()
