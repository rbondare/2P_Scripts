"""Threshold analysis of the 30 stratified H.264 clips: user labels vs window-mean index, on_fraction, amplitude features."""
import os, csv, numpy as np, pandas as pd
CL = r'C:\Users\rbondarenko\DLC\whisking_check_h264_strat\clips'
key = pd.read_csv(os.path.join(CL, 'answer_key.csv'))
ans = {}
for r in csv.reader(open(os.path.join(CL, 'your_answers.csv'))):
    cells = [c.strip(' ;\t"').lower() for c in r[1:]]
    cells = [c for c in cells if c]
    if r and r[0].startswith('clip') and cells:
        ans[r[0]] = cells[0][0]
key['you'] = key['clip'].map(ans)
key['y'] = (key['you'] == 'y').astype(int)
print(key[['clip','session','start_s','mean_index','on_fraction','code_says','you']].sort_values('mean_index').to_string(index=False))
d = key[key.you.isin(['y','n'])]
from itertools import product
def auc(x, y):
    pos, neg = x[y == 1], x[y == 0]
    return np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg])
print('\nAUC mean_index:', round(auc(d.mean_index.values, d.y.values), 3), ' AUC on_fraction:', round(auc(d.on_fraction.values, d.y.values), 3))
print('\nthreshold on window-mean index -> accuracy / sens / spec')
for t in [1.2,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75]:
    p = (d.mean_index >= t).astype(int)
    tp = ((p==1)&(d.y==1)).sum(); tn = ((p==0)&(d.y==0)).sum(); fp = ((p==1)&(d.y==0)).sum(); fn = ((p==0)&(d.y==1)).sum()
    print('  >= %.2f: acc %.2f  sens %d/%d  spec %d/%d' % (t, (tp+tn)/len(d), tp, tp+fn, tn, tn+fp))
print('\nanswers by bin:')
print(pd.crosstab(key.bin, key.you))
print(pd.crosstab(key.code_says, key.you))
