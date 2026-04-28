"""
把 unmapped_labels.csv 里的新标签合并进 cell_type_mapping.csv。

逻辑：
- 若 cl_id 已存在于 cell_type_mapping.csv → 把 label_text 追加到对应行的 synonyms 列
- 若 cl_id 不存在 → 新增一行，synonyms = label_text
- 同一 cl_id 可能有多个 label_text，全部 @-连接追加到同一行 synonyms
- 不会覆盖已有的 synonyms，只追加

用法:
  python merge_lookup_table.py
  python merge_lookup_table.py --new data/labels/unmapped_labels.csv --table cell_type_mapping/cell_type_mapping.csv
"""
import argparse
import os
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def merge(new_path: str, table_path: str, output_path: str):
    new_df = pd.read_csv(new_path, encoding="utf-8-sig")
    table_df = pd.read_csv(table_path, encoding="utf-8-sig")

    # 建立 cl_id -> 行索引 的快速查找
    clid_to_idx = {row["cl_id"]: i for i, row in table_df.iterrows()}

    added_to_existing = 0
    added_as_new_row = 0
    skipped_already_exists = 0

    for _, row in new_df.iterrows():
        label = str(row["label_text"]).strip()
        cl_id = str(row["cl_id"]).strip()
        cl_label = str(row["cl_label"]).strip()

        if not label or not cl_id:
            continue

        if cl_id in clid_to_idx:
            idx = clid_to_idx[cl_id]
            existing_syns = str(table_df.at[idx, "synonyms"])
            syn_list = [s.strip() for s in existing_syns.split("@") if s.strip()]

            if label.lower() in [s.lower() for s in syn_list]:
                skipped_already_exists += 1
            else:
                syn_list.append(label)
                table_df.at[idx, "synonyms"] = "@".join(syn_list)
                added_to_existing += 1
        else:
            new_row = {"synonyms": label, "cl_id": cl_id, "cl_label": cl_label}
            table_df = pd.concat([table_df, pd.DataFrame([new_row])], ignore_index=True)
            clid_to_idx[cl_id] = len(table_df) - 1
            added_as_new_row += 1

    table_df.to_csv(output_path, index=False, encoding="utf-8-sig")

    print(f"合并完成 → {output_path}")
    print(f"  追加到已有行的标签: {added_to_existing}")
    print(f"  新增行 (新 CL ID):  {added_as_new_row}")
    print(f"  跳过(已存在):       {skipped_already_exists}")
    print(f"  最终 lookup table 行数: {len(table_df)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--new", default=os.path.join(_ROOT, "data", "labels", "unmapped_labels.csv"))
    parser.add_argument("--table", default=os.path.join(_ROOT, "cell_type_mapping", "cell_type_mapping.csv"))
    parser.add_argument("--output", default=os.path.join(_ROOT, "cell_type_mapping", "cell_type_mapping.csv"),
                        help="输出路径（默认原地覆盖）")
    args = parser.parse_args()
    merge(args.new, args.table, args.output)


if __name__ == "__main__":
    main()
