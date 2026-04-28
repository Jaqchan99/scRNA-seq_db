"""
修复 lookup table 中已知的错误映射。

分两类：
  redirect: synonym 目前挂在错误的 CL 行 → 把它从错误行删掉，追加到正确行
  add:      synonym 目前不在正确行 → 追加到正确行（fuzzy 打偏的情况）
"""
import os
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLE_PATH = os.path.join(_ROOT, "cell_type_mapping", "cell_type_mapping.csv")

# ── 修正清单 ──────────────────────────────────────────────────
# (label_text, correct_cl_id, correct_cl_label, action)
# action: "redirect" = 从错误行移走并写入正确行；"add" = 直接追加到正确行
CORRECTIONS = [
    # ── exact match 挂错行（redirect）──
    ("Lp",                  "CL:0002326", "luminal epithelial cell of mammary gland", "redirect"),
    ("delta",               "CL:0000173", "pancreatic D cell",                        "redirect"),
    ("K.cells",             "CL:0000164", "enteroendocrine cell",                     "redirect"),
    ("D.cells",             "CL:0002266", "type D cell of small intestine",           "redirect"),
    ("CPC",                 "CL:0000706", "choroid plexus epithelial cell",            "redirect"),
    ("EPC",                 "CL:0000065", "ependymal cell",                           "redirect"),
    ("PC",                  "CL:0000669", "pericyte cell",                            "redirect"),
    ("Basal",               "CL:1000348", "basal cell of epithelium of trachea",      "redirect"),
    # ── fuzzy 打偏（add synonym 到正确行）──
    ("Neutro",              "CL:0000775", "neutrophil",                               "add"),
    ("activated_stellate",  "CL:0002410", "pancreatic stellate cell",                 "add"),
    ("Enterocyte.dist",     "CL:1000334", "enterocyte of epithelium of small intestine", "add"),
    ("Enterocyte.prox",     "CL:1000334", "enterocyte of epithelium of small intestine", "add"),
    ("Paneth.progenitor",   "CL:0000510", "paneth cell",                              "add"),
    ("IL cell",             "CL:0001065", "innate lymphoid cell",                     "add"),
    ("Myoepithelial",       "CL:0000185", "myoepithelial cell",                       "add"),
    ("Mural Cell + Fibroblast",  "CL:0000057", "fibroblast",                          "add"),
    ("Mural Cell+ Fibroblast",   "CL:0000057", "fibroblast",                          "add"),
    ("Mural+Fibroblast",         "CL:0000057", "fibroblast",                          "add"),
    ("bipolar",             "CL:0000103", "retinal bipolar neuron",                   "add"),
    ("vascular_endothelium","CL:0000071", "blood vessel endothelial cell",             "add"),
    ("Bursal fibroblasts",  "CL:0000057", "fibroblast",                               "add"),
]


def normalize(name: str) -> str:
    import re
    name = str(name).strip().lower()
    name = name.replace('-', ' ').replace('_', ' ').replace('.', ' ')
    name = re.sub(r'\s+', ' ', name).strip()
    return name


def apply_corrections():
    df = pd.read_csv(TABLE_PATH, encoding="utf-8-sig")

    redirected = 0
    added = 0
    skipped = 0

    for label, correct_cl_id, correct_cl_label, action in CORRECTIONS:
        norm_label = normalize(label)

        # 找正确行的索引
        correct_idx = df.index[df['cl_id'] == correct_cl_id].tolist()
        if not correct_idx:
            # 正确 CL ID 不存在，新增行
            df = pd.concat([df, pd.DataFrame([{
                "synonyms": label,
                "cl_id": correct_cl_id,
                "cl_label": correct_cl_label
            }])], ignore_index=True)
            print(f"  [NEW ROW] {label} → {correct_cl_id} ({correct_cl_label})")
            added += 1
            continue

        correct_idx = correct_idx[0]

        if action == "redirect":
            # 1. 从所有行里删掉这个 synonym（防止 exact match 打到错误行）
            removed_from = []
            for i, row in df.iterrows():
                syns = [s.strip() for s in str(row['synonyms']).split('@') if s.strip()]
                norm_syns = [normalize(s) for s in syns]
                if norm_label in norm_syns and i != correct_idx:
                    # 删掉这个 synonym
                    new_syns = [s for s in syns if normalize(s) != norm_label]
                    df.at[i, 'synonyms'] = '@'.join(new_syns)
                    removed_from.append(str(row['cl_id']))
            if removed_from:
                print(f"  [REDIRECT] {label}: 从 {removed_from} 移走")
                redirected += 1
            else:
                print(f"  [REDIRECT-SKIP] {label}: 未找到错误行（可能已修复）")

        # 追加到正确行（无论 redirect 还是 add）
        existing_syns = [s.strip() for s in str(df.at[correct_idx, 'synonyms']).split('@') if s.strip()]
        norm_existing = [normalize(s) for s in existing_syns]
        if norm_label not in norm_existing:
            existing_syns.append(label)
            df.at[correct_idx, 'synonyms'] = '@'.join(existing_syns)
            if action == "add":
                print(f"  [ADD] {label} → {correct_cl_id} ({correct_cl_label})")
                added += 1
        else:
            print(f"  [SKIP] {label} 已在正确行 {correct_cl_id}")
            skipped += 1

    df.to_csv(TABLE_PATH, index=False, encoding="utf-8-sig")
    print()
    print(f"完成 → {TABLE_PATH}")
    print(f"  redirect: {redirected}  add: {added}  skip: {skipped}")
    print(f"  最终行数: {len(df)}")


if __name__ == "__main__":
    apply_corrections()
