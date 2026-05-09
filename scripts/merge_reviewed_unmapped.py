"""
将用户在平台下载的 unmapped CSV 中「已审核」行合并进 cell_type_mapping.csv。

闭环：平台 unmapped-labels 下载 → 用户填写 cl_id / cl_label → reviewed=yes →
维护者审核后运行本脚本 → 调用 merge_lookup_table.merge（自动备份建议由你手动执行）。

输入 CSV 列（与 GET /submissions/{id}/unmapped-labels 导出一致，可增列）：
  label_text  — 必填，原始细胞类型字符串（对应 merge 的 label_text）
  cl_id       — 审核通过必填，格式建议 CL:#####
  cl_label    — 建议填写标准 CL 名称
  reviewed    — yes / true / 1 / y 表示可合并（默认仅导出为 no）

用法:
  python scripts/merge_reviewed_unmapped.py --input data/labels/reviewed_batch.csv --dry-run
  python scripts/merge_reviewed_unmapped.py --input data/labels/reviewed_batch.csv

回滚稳定版本：在提交前 git tag / git commit；若要撤销合并，使用 merge_lookup_table 运行前备份的
cell_type_mapping.csv.bak_* 覆盖回 cell_type_mapping/cell_type_mapping.csv。
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
import tempfile

import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_merge():
    mpath = os.path.join(os.path.dirname(__file__), "merge_lookup_table.py")
    spec = importlib.util.spec_from_file_location("merge_lookup_table", mpath)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod.merge


def _truthy_reviewed(val) -> bool:
    if pd.isna(val):
        return False
    s = str(val).strip().lower()
    return s in ("1", "yes", "y", "true", "审核", "通过")


def main():
    parser = argparse.ArgumentParser(description="Merge reviewed unmapped labels into lookup table.")
    parser.add_argument("--input", "-i", required=True, help="填写后的 CSV 路径")
    parser.add_argument(
        "--table",
        default=os.path.join(_ROOT, "cell_type_mapping", "cell_type_mapping.csv"),
        help="lookup table 路径",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(_ROOT, "cell_type_mapping", "cell_type_mapping.csv"),
        help="输出路径（默认同 table 原地更新）",
    )
    parser.add_argument("--dry-run", action="store_true", help="只校验与打印摘要，不写文件")
    parser.add_argument(
        "--force-all-with-clid",
        action="store_true",
        help="CSV 无 reviewed 列时允许合并所有已填 cl_id 的行（慎用）",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input, encoding="utf-8-sig")
    col_label = None
    for c in ("label_text", "raw_cell_type", "label"):
        if c in df.columns:
            col_label = c
            break
    if not col_label:
        print("❌ 未找到 label_text / raw_cell_type 列")
        sys.exit(1)

    if "cl_id" not in df.columns:
        print("❌ 缺少 cl_id 列")
        sys.exit(1)

    has_review = "reviewed" in df.columns
    if not has_review and not args.force_all_with_clid:
        print("❌ CSV 缺少 reviewed 列。请为拟合并行填写 reviewed=yes，或显式使用 --force-all-with-clid")
        sys.exit(1)

    rows = []
    skipped = 0
    for _, r in df.iterrows():
        label = str(r[col_label]).strip()
        cl_id = str(r["cl_id"]).strip() if pd.notna(r.get("cl_id")) else ""
        cl_label = str(r["cl_label"]).strip() if "cl_label" in df.columns and pd.notna(r.get("cl_label")) else ""

        if not label or label.lower() == "nan":
            skipped += 1
            continue
        if not cl_id or cl_id.lower() == "nan":
            skipped += 1
            continue
        if not cl_id.upper().startswith("CL:"):
            print(f"⚠️ 跳过无效 cl_id（应以 CL: 开头）: {label!r} -> {cl_id!r}")
            skipped += 1
            continue
        if has_review and not _truthy_reviewed(r.get("reviewed")):
            skipped += 1
            continue

        if not cl_label or cl_label.lower() == "nan":
            cl_label = label

        rows.append({"label_text": label, "cl_id": cl_id, "cl_label": cl_label})

    if not rows:
        print("❌ 没有可合并的行（检查 reviewed / cl_id）")
        sys.exit(1)

    merge_df = pd.DataFrame(rows)
    print(f"准备合并 {len(merge_df)} 条标签（跳过 {skipped} 行）")
    print(merge_df.head(10).to_string(index=False))

    if args.dry_run:
        print("[dry-run] 未写入")
        return

    merge = _load_merge()
    ts = __import__("datetime").datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{args.table}.bak_reviewed_{ts}"
    if os.path.exists(args.table):
        shutil.copy2(args.table, bak)
        print(f"已备份: {bak}")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8-sig", newline="") as tf:
        merge_df.to_csv(tf.name, index=False)
        tmp = tf.name
    try:
        merge(tmp, args.table, args.output)
    finally:
        os.unlink(tmp)


if __name__ == "__main__":
    main()
