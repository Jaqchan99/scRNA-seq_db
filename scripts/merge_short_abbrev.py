"""
把 short_abbrev_unmapped.csv 里的短缩写批量追加到 lookup table。
逻辑与 merge_lookup_table.py 相同：
  - cl_id 已存在 → 把 label_text 追加到对应行的 synonyms
  - cl_id 不存在 → 新增行
"""
import os
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEW_PATH = os.path.join(_ROOT, "data", "labels", "short_abbrev_unmapped.csv")
TABLE_PATH = os.path.join(_ROOT, "cell_type_mapping", "cell_type_mapping.csv")

new_df = pd.read_csv(NEW_PATH, encoding="utf-8-sig")
table_df = pd.read_csv(TABLE_PATH, encoding="utf-8-sig")

clid_to_idx = {str(row["cl_id"]).strip(): i for i, row in table_df.iterrows()}

added_to_existing = 0
added_as_new_row = 0
skipped = 0

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
            skipped += 1
        else:
            syn_list.append(label)
            table_df.at[idx, "synonyms"] = "@".join(syn_list)
            added_to_existing += 1
    else:
        new_row = {"synonyms": label, "cl_id": cl_id, "cl_label": cl_label}
        table_df = pd.concat([table_df, pd.DataFrame([new_row])], ignore_index=True)
        clid_to_idx[cl_id] = len(table_df) - 1
        added_as_new_row += 1

table_df.to_csv(TABLE_PATH, index=False, encoding="utf-8-sig")

print(f"完成 → {TABLE_PATH}")
print(f"  追加到已有行: {added_to_existing}")
print(f"  新增行:       {added_as_new_row}")
print(f"  跳过(已存在): {skipped}")
print(f"  最终行数:     {len(table_df)}")
