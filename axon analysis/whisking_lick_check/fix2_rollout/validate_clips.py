"""Score the new index (front_top ME / REST) against the user's blind judgement of 24 clips (5 s windows, 23 decided,
all MPEG-4 sessions) and compare with the old index (mean of fan_front and fan_down).
usage: python validate_clips.py [rest_mpeg4]"""
import os, sys, csv
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter1d

here = os.path.dirname(os.path.abspath(__file__))
CL = r'C:\Users\rbondarenko\DLC\whisking_check\clips'
key = {r['clip']: r for r in csv.DictReader(open(os.path.join(CL, 'answer_key.csv')))}
ans = {}
for r in csv.reader(open(os.path.join(CL, 'your_answers.csv'))):
    cells = [c.strip(' ;\t"').lower() for c in r[1:]]
    cells = [c for c in cells if c]
    if r and r[0] in key and cells:
        ans[r[0]] = cells[0][0]
rest_new = float(sys.argv[1]) if len(sys.argv) > 1 else None
rows = []
for c, a in sorted(ans.items()):
    if a not in 'yn':
        continue
    s, st = key[c]['session'], float(key[c]['start_s'])
    new = pd.read_csv(os.path.join(here, 'pass1', s, s + '_whisker.csv'))
    old = pd.read_csv(os.path.join(r'C:\Users\rbondarenko\DLC\work', s, 'whisker_results_v1' if os.path.isdir(os.path.join(r'C:\Users\rbondarenko\DLC\work', s, 'whisker_results_v1')) else 'whisker_results', s + '_whisker.csv'))
    f0, f1 = int(round(st * 50)), int(round(st * 50)) + 250
    sm = gaussian_filter1d(new['me_front_top'].values, 3)[f0:f1]
    smo = gaussian_filter1d(old['me_fan_front'].values, 3)[f0:f1] / 0.447
    smd = gaussian_filter1d(old['me_fan_down'].values, 3)[f0:f1] / 0.450
    rows.append(dict(clip=c, session=s, user=a, code_old=key[c]['code_says'], top=sm.mean(), old_index=old['whisk_index'].values[f0:f1].mean(),
                     front_all=smo.mean(), jaw=smd.mean(), pad=old['me_pad'].values[f0:f1].mean()))
D = pd.DataFrame(rows)
D['y'] = (D.user == 'y').astype(int)


def auc(x, y):
    pos, neg = x[y == 1], x[y == 0]
    return np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg])


print('%d decided clips: %d whisking, %d still (user)' % (len(D), D.y.sum(), (1 - D.y).sum()))
print('AUC (window mean vs user y/n):')
for k in ('old_index', 'top', 'front_all', 'jaw', 'pad'):
    print('  %-10s %.3f' % (k, auc(D[k].values, D.y.values)))
print('old code agreement with user: %d of %d' % (((D.code_old == 'whisking') == (D.y == 1)).sum(), len(D)))
if rest_new:
    thr = 1.6 * rest_new
    pred = D.top > thr
    print('new index at threshold 1.6 x REST (%.3f): agreement %d of %d; whisking clips called still: %d; still clips called whisking: %d' % (
        rest_new, (pred == (D.y == 1)).sum(), len(D), ((D.y == 1) & ~pred).sum(), ((D.y == 0) & pred).sum()))
D['top_over_rest'] = D.top / rest_new if rest_new else np.nan
print(D.sort_values('top')[['clip', 'session', 'user', 'code_old', 'old_index', 'top', 'top_over_rest', 'jaw', 'pad']].round(2).to_string(index=False))
D.to_csv(os.path.join(here, 'validate_clips.csv'), index=False)
