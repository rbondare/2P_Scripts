"""Candidate on/off rules for H.264 scored against the user's labels on both blind clip sets (set1 5 s, strat 3 s)."""
import os, csv, numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter1d
W = r'C:\Users\rbondarenko\DLC\work'
FPS = 50.0
def labels(cl):
    key = pd.read_csv(os.path.join(cl, 'answer_key.csv'))
    ans = {}
    for r in csv.reader(open(os.path.join(cl, 'your_answers.csv'))):
        cells = [c.strip(' ;\t"').lower() for c in r[1:]]
        cells = [c for c in cells if c]
        if r and r[0].startswith('clip') and cells:
            ans[r[0]] = cells[0][0]
    key['you'] = key['clip'].map(ans)
    return key
rows = []
cache = {}
for tag, cl, clip_s in (('set1', r'C:\Users\rbondarenko\DLC\whisking_check_h264\clips', 5.0), ('strat', r'C:\Users\rbondarenko\DLC\whisking_check_h264_strat\clips', 3.0)):
    k = labels(cl)
    for _, r in k.iterrows():
        s = r.session
        if s not in cache:
            cache[s] = pd.read_csv(os.path.join(W, s, 'whisker_results', s + '_whisker.csv'))['whisk_index'].values
        a = int(r.start_s * FPS)
        seg = cache[s][a:a + int(clip_s * FPS)]
        rows.append(dict(set=tag, clip=r.clip, session=s, you=r.you, seg=seg))
D = pd.DataFrame(rows)
J = D[D.you.isin(['y', 'n'])].copy()
J['y'] = (J.you == 'y').astype(int)
print('judged clips: %d (%d whisking, %d still)' % (len(J), J.y.sum(), (1 - J.y).sum()))
J['mean'] = J.seg.map(np.mean); J['med'] = J.seg.map(np.median); J['p75'] = J.seg.map(lambda x: np.percentile(x, 75)); J['p90'] = J.seg.map(lambda x: np.percentile(x, 90))
def auc(x, y):
    pos, neg = x[y == 1], x[y == 0]
    return np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg])
for f in ('mean', 'med', 'p75', 'p90'):
    print('AUC %-4s %.3f' % (f, auc(J[f].values, J.y.values)))
def hyst(x, hi, lo):
    on = np.zeros(len(x), bool); st = False
    for i, v in enumerate(x):
        if not st and v > hi: st = True
        elif st and v < lo: st = False
        on[i] = st
    return on
print('\nwindow-mean cut: best thresholds')
for t in np.arange(1.30, 1.56, 0.025):
    p = (J['mean'] >= t); y = J.y.astype(bool)
    print('  mean >= %.3f: wrong %2d  (missed whisking %d, false whisking %d)' % (t, (p != y).sum(), (y & ~p).sum(), (~y & p).sum()))
print('\nper-frame hysteresis rule (on if above hi until below lo), clip called whisking if on >50%% of frames; current = hi 1.76 lo 1.44')
res = []
for hi in np.arange(1.40, 1.81, 0.05):
    for lo in (1.2, 1.3, 1.35, 1.4, 1.44):
        if lo > hi: continue
        p = J.seg.map(lambda x: hyst(x, hi, lo).mean() > 0.5); y = J.y.astype(bool)
        res.append((hi, lo, (p != y).sum(), (y & ~p).sum(), (~y & p).sum()))
res.sort(key=lambda r: (r[2], r[0]))
for r in res[:10]:
    print('  hi %.2f lo %.2f: wrong %d (missed %d, false %d)' % r)
cur = [r for r in res if abs(r[0] - 1.75) < 1e-6 and r[1] == 1.44] or [r for r in res if abs(r[0]-1.75) < 0.03 and r[1] == 1.44]
print('current-like:', cur)
wrong = J[(J['mean'] >= 1.4) != (J.y == 1)]
print('\nclips wrong at mean>=1.40:'); print(wrong[['set', 'clip', 'session', 'you', 'mean']].to_string(index=False))
J.drop(columns='seg').to_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'strat_rules_table.csv'), index=False)
