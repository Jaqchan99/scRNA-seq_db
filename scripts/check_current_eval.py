import os
import pandas as pd

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_EVAL_PATH = os.path.join(_ROOT, "results", "eval", "eval_results_local.csv")

res = pd.read_csv(_EVAL_PATH, encoding='utf-8-sig')
print('总行数:', len(res))
print('正确数:', int(res['is_correct'].sum()))
print('错误数:', int((res['is_correct']==False).sum()))
wrong = res[res['is_correct']==False]
unmapped = wrong[wrong['pred_cl_id'].fillna('')=='']
mismatched = wrong[wrong['pred_cl_id'].fillna('')!='']
print('  未映射 (unmapped):', len(unmapped))
print('  映射错误 (mismatched):', len(mismatched))
total = len(res)
covered = int((res['pred_cl_id'].fillna('')!='').sum())
correct = int(res['is_correct'].sum())
print()
print('Coverage:', str(round(covered/total*100, 1)) + '%')
print('Accuracy:', str(round(correct/total*100, 2)) + '%')
