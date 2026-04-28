import os
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_EVAL_PATH = os.path.join(_ROOT, "results", "eval", "eval_results_local.csv")

res = pd.read_csv(_EVAL_PATH, encoding='utf-8-sig')
wrong = res[res['is_correct']==False].copy()
mismatched = wrong[wrong['pred_cl_id'].fillna('')!='']
print(f'当前 mismatched 总数: {len(mismatched)}')
print()
for _, r in mismatched.iterrows():
    label = r['label_text']
    true_id = r['cl_id']
    true_label = r['cl_label']
    pred_id = r['pred_cl_id']
    pred_label = r['pred_cl_label']
    source = r['pred_source']
    print(f'[{label}]')
    print(f'  真实: {true_id}  ({true_label})')
    print(f'  预测: {pred_id}  ({pred_label})  来源={source}')
