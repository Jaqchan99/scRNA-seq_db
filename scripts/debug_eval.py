import os
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_EVAL_PATH = os.path.join(_ROOT, "results", "eval", "eval_results_local.csv")

res = pd.read_csv(_EVAL_PATH, encoding='utf-8-sig')
wrong = res[res['is_correct'] == False]
unmapped = wrong[wrong['pred_cl_id'].fillna('') == '']
mismatch = wrong[wrong['pred_cl_id'].fillna('') != '']

print(f"总标签: {len(res)}")
print(f"正确: {res['is_correct'].sum()}")
print(f"错误总数: {len(wrong)}  (未映射: {len(unmapped)}, 映射错误: {len(mismatch)})")

print("\n=== 映射错误（映射了但映射到错误CL ID）===")
for _, r in mismatch.iterrows():
    label = r['label_text']
    true_cl = r['cl_id']
    pred_cl = r['pred_cl_id']
    src = r['pred_source']
    print(f"  [{label}]  真实={true_cl}  预测={pred_cl}  来源={src}")

print("\n=== 未映射（前40条）===")
for _, r in unmapped.head(40).iterrows():
    label = r['label_text']
    true_cl = r['cl_id']
    print(f"  [{label}]  真实={true_cl}")

print("\n=== 正确映射的标签（前20条）===")
correct = res[res['is_correct'] == True]
for _, r in correct.head(20).iterrows():
    print(f"  [{r['label_text']}]  CL={r['cl_id']}  来源={r['pred_source']}")
