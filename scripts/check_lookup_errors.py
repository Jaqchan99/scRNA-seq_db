import os
import pandas as pd
import re

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLE_PATH = os.path.join(_ROOT, "cell_type_mapping", "cell_type_mapping.csv")
table = pd.read_csv(TABLE_PATH, encoding='utf-8-sig')

def normalize(name):
    name = str(name).strip().lower()
    name = name.replace('-', ' ').replace('_', ' ').replace('.', ' ')
    name = re.sub(r'\s+', ' ', name).strip()
    return name

keywords = ['macrophage', 'EC', 'Basal', 'Enterocyte', 'AT2', 'delta', 'Lp',
            'K.cells', 'M.cell', 'D.cells', 'Red blood cells', 'COP', 'EPC', 'CPC']

for kw in keywords:
    nkw = normalize(kw)
    matches = []
    for _, row in table.iterrows():
        syns = [normalize(s) for s in str(row['synonyms']).split('@')]
        if nkw in syns:
            matches.append((row['cl_id'], row['cl_label']))
    if len(matches) == 1:
        print(f"[{kw}] -> {matches[0][0]}  {matches[0][1]}")
    elif len(matches) > 1:
        print(f"[{kw}] *** 多条匹配 ***")
        for m in matches:
            print(f"       {m[0]}  {m[1]}")
    else:
        print(f"[{kw}] NOT FOUND")
