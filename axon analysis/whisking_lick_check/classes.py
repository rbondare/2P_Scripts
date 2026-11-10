"""Lick x locomotion classes of 1 s windows and the motion-energy boxes in each (whisker csv), plus the window lists for the pixel maps.
usage: python classes.py <session>"""
import sys, os, numpy as np, pandas as pd, scipy.io as sio

sess = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
W = r'C:\Users\rbondarenko\DLC\work\%s\whisker_results\%s_whisker.csv' % (sess, sess)
X = sio.loadmat(r'C:\Users\rbondarenko\2P_cache\xs\%s_xs.mat' % sess, squeeze_me=True)
T = sio.loadmat(r'C:\Users\rbondarenko\2P_cache\xs\%s_task.mat' % sess, squeeze_me=True)
lick = np.sort(np.atleast_1d(T['lick_t']).astype(float))
d = pd.read_csv(W)
fps = 50.0
REST = {'me_fan_front': 0.447, 'me_fan_down': 0.450}
tc = np.atleast_1d(X['tcam']).astype(float)
lt, lf, ls, lb = [np.atleast_1d(X[k]).astype(float) for k in ('loco_t', 'loco_fwd', 'loco_side', 'loco_blk')]
print(sess, 'loco_sync_ok', X['loco_sync_ok'], 'n_values/n_trig', X['loco_n_values'], X['loco_n_trig'])
n = min(len(d), len(tc))
tc = tc[:n]
speed = np.full(n, np.nan)
for k in np.unique(lb[np.isfinite(lb)]):
    m = lb == k
    tt = lt[m]
    if m.sum() < 100:
        continue
    dtl = np.median(np.diff(tt))
    inc = np.hypot(lf[m] - np.nanmedian(lf[m]), ls[m] - np.nanmedian(ls[m]))
    w = max(1, int(round(1.0 / dtl)))
    sp = np.convolve(np.nan_to_num(inc), np.ones(w) / w, 'same') / dtl
    mm = (tc >= tt[0]) & (tc <= tt[-1])
    speed[mm] = np.interp(tc[mm], tt, sp)
print('speed finite %.2f' % np.isfinite(speed).mean(), 'median speed %.1f' % np.nanmedian(speed))

wl = int(1.0 * fps)
nwin = n // wl
rows = []
for w in range(nwin):
    a, b = w * wl, (w + 1) * wl
    if tc[b - 1] - tc[a] > 1.1:
        continue
    t0, t1 = tc[a], tc[b - 1]
    nl = np.searchsorted(lick, t1) - np.searchsorted(lick, t0)
    nl_m = np.searchsorted(lick, t1 + 1.0) - np.searchsorted(lick, t0 - 1.0)
    sp = np.nanmedian(speed[a:b]) if np.isfinite(speed[a:b]).any() else np.nan
    art = d['artifact'].values[a:b].mean()
    rows.append((a, b, t0, nl, nl_m, sp, art, *[np.nanmean(d[c].values[a:b]) for c in ('me_fan_front', 'me_fan_down', 'me_pad', 'me_ctrl')]))
R = pd.DataFrame(rows, columns=['a', 'b', 't0', 'nl', 'nl_margin', 'speed', 'art', 'front', 'down', 'pad', 'ctrl'])
R = R[(R.art < 0.05) & np.isfinite(R.speed)]
R['lickcls'] = np.where(R.nl >= 3, 'lick', np.where(R.nl_margin == 0, 'nolick', 'other'))
R['runcls'] = np.where(R.speed > 25, 'run', np.where(R.speed < 5, 'still', 'mid'))
R['index'] = 0.5 * (R.front / REST['me_fan_front'] + R.down / REST['me_fan_down'])
print('windows', len(R))
g = R[(R.lickcls != 'other') & (R.runcls != 'mid')].groupby(['runcls', 'lickcls'])
print(g[['front', 'down', 'pad', 'ctrl', 'index']].median().round(2).join(g.size().rename('n')))
R.to_pickle(os.path.join(here, 'win_%s.pkl' % sess))
