import os
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_EVAL_PATH = os.path.join(_ROOT, "results", "eval", "eval_results_local.csv")

res = pd.read_csv(_EVAL_PATH, encoding='utf-8-sig')

wrong = res[res['is_correct'] == False].copy()
unmapped = wrong[wrong['pred_cl_id'].fillna('') == '']
mismatched = wrong[wrong['pred_cl_id'].fillna('') != '']

print(f"总错误数: {len(wrong)}")
print(f"未映射 (unmapped): {len(unmapped)}")
print(f"映射错误 (mismatched): {len(mismatched)}")

print()
print("=== 未映射标签 (全部) ===")
for _, r in unmapped.iterrows():
    label = r['label_text']
    cl_id = r['cl_id']
    cl_label = r['cl_label']
    print(f"  [{label}]  真实={cl_id}  ({cl_label})")

print()
print("=== 映射错误标签 (全部) ===")
for _, r in mismatched.iterrows():
    label = r['label_text']
    true_id = r['cl_id']
    true_label = r['cl_label']
    pred_id = r['pred_cl_id']
    pred_label = r['pred_cl_label']
    source = r['pred_source']
    score = r.get('pred_confidence', '')
    print(f"  [{label}]  真实={true_id} ({true_label})  预测={pred_id} ({pred_label})  来源={source} 置信={score}")
