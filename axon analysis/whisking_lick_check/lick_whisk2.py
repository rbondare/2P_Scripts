"""Second pass on RB24_261109_1506: (a) how much variance of each signal is lick-driven (held-out, unique R2 beyond speed + pupil), (b) is there whisking-band rhythm in front_top,
(c) effect sizes run vs still in lick-free frames, (d) stability of the matched lick effect (session halves, other lick definitions), (e) which region carries the old index' lick effect."""
import os, json, glob, numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter1d, uniform_filter1d, median_filter
from scipy.signal import welch

here = os.path.dirname(os.path.abspath(__file__))
sess = 'AnimalRB24_261109_1506'
FPS = 50.0
A = np.load(os.path.join(here, 'align_%s.npz' % sess))
tcam, tball = A['aux1'], A['aux2']
lk = np.load(os.path.join(here, 'lick_clean_%s.npy' % sess))
n = len(tcam)
cells = np.load(os.path.join(here, 'cells_%s.npy' % sess), mmap_mode='r')
meta = json.load(open(os.path.join(here, 'cells_%s.json' % sess)))
CS, GX0, GY0 = meta['CS'], meta['GX'][0], meta['GY'][0]
W = pd.read_csv(r'C:\Users\rbondarenko\DLC\work\%s\whisker_results\%s_whisker.csv' % (sess, sess))
art = W['artifact'].values.astype(bool)


def reg(x0, y0, x1, y1):
    ix0, ix1, iy0, iy1 = int((x0 - GX0) // CS), int((x1 - GX0) // CS), int((y0 - GY0) // CS), int((y1 - GY0) // CS)
    out = np.empty(n, np.float32)
    for k in range(0, n, 20000):
        out[k:k + 20000] = np.nanmean(cells[k:k + 20000, iy0:iy1, ix0:ix1], axis=(1, 2))
    return out


REG = {'front_all': (-260, -120, 60, 200), 'front_top': (-260, -120, 60, 0), 'front_low': (-260, 0, 60, 200), 'jaw': (60, 80, 300, 300), 'pad': (60, -100, 260, 60), 'noise': (-320, -160, -200, -60)}
ME = {k: reg(*v) for k, v in REG.items()}
REST = {'front_all': 0.447, 'jaw': 0.450}
ME['old_index'] = 0.5 * (ME['front_all'] / REST['front_all'] + ME['jaw'] / REST['jaw'])

ball = np.load(os.path.join(here, 'ball_%s.npy' % sess))
tb = tball[:len(ball)]
mag = np.hypot(ball[:, 0], ball[:, 1])
dtb = np.median(np.diff(tb))
speed = np.interp(tcam, tb, uniform_filter1d(mag, int(round(0.5 / dtb))) / dtb).astype(np.float32)
E = pd.read_csv(glob.glob(r'C:\Users\rbondarenko\DLC\work\%s\eye_results\*_eye.csv' % sess)[0], skiprows=2)
pup = pd.Series(np.pi / 4 * E['pupil_width'].values * E['pupil_height'].values).interpolate(limit_direction='both').values[:n]
pup_s = median_filter(pup, size=int(2 * FPS) | 1)
lick_f = np.searchsorted(tcam, lk)
lick_f = lick_f[(lick_f > 0) & (lick_f < n)]
imp = np.zeros(n, np.float32)
np.add.at(imp, lick_f, 1)
ili = np.diff(lk)
b0 = np.flatnonzero(np.r_[True, ili >= 1.0])
b1 = np.r_[b0[1:], len(lk)] - 1
bout = np.zeros(n, bool)
for a, b in zip(b0, b1):
    fa, fb = np.searchsorted(tcam, [lk[a] - 0.3, lk[b] + 0.5])
    bout[fa:fb] = True
near1 = np.convolve(imp, np.ones(int(2 * FPS) + 1), 'same') > 0
good = ~art
sz = (speed - speed.mean()) / speed.std()
pz = (pup_s - pup_s.mean()) / pup_s.std()
Zbase = np.stack([np.ones(n), sz, sz ** 2, pz, sz * pz], 1)
KL = np.arange(-10, 26)
cols = []
for k in KL:
    c = np.roll(imp, k).copy()
    if k > 0:
        c[:k] = 0
    elif k < 0:
        c[k:] = 0
    cols.append(c)
Xk = np.stack(cols, 1)
Zl = np.hstack([Xk, bout[:, None].astype(np.float32)])
blk = (np.arange(n) // int(300 * FPS)) % 2


def r2cv(y, X):
    """held-out R2 with 2-fold alternating 5-min blocks, on 0.2 s smoothed y, artifact frames excluded"""
    y = gaussian_filter1d(np.asarray(y, np.float64), 5)
    pred = np.zeros(n)
    for h in (0, 1):
        tr = (blk != h) & good
        D = X[tr]
        lam = 1e-3 * np.trace(D.T @ D) / D.shape[1]
        beta = np.linalg.solve(D.T @ D + lam * np.eye(D.shape[1]), D.T @ y[tr])
        pred[blk == h] = X[blk == h] @ beta
    m = good
    return 1 - np.sum((y[m] - pred[m]) ** 2) / np.sum((y[m] - y[m].mean()) ** 2)


print('== (a) held-out R2 of 0.2 s-smoothed ME: movement state (speed, pupil) vs + lick terms (kernel + bout)')
print('%-12s %10s %12s %12s %10s' % ('signal', 'state only', 'state+lick', 'lick only', 'unique lick'))
for k in ('old_index', 'jaw', 'pad', 'front_all', 'front_low', 'front_top', 'noise'):
    a = r2cv(ME[k], Zbase)
    b = r2cv(ME[k], np.hstack([Zbase, Zl]))
    c = r2cv(ME[k], np.hstack([Zbase[:, :1], Zl]))
    print('%-12s %10.3f %12.3f %12.3f %10.3f' % (k, a, b, c, b - a))

print('\n== (b) spectrum of front_top / front_all / noise ME in lick-free frames, running vs still (Welch, 2 s segments)')
cl = ~near1 & good
run = cl & (speed > np.percentile(speed, 80))
stl = cl & (speed < np.percentile(speed, 20))


def runs_of(mask, minlen=100):
    d = np.diff(np.r_[False, mask, False].astype(int))
    return [(a, b) for a, b in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)) if b - a >= minlen]


def spec(x, mask):
    P = []
    for a, b in runs_of(mask):
        f, p = welch(x[a:b] - x[a:b].mean(), fs=FPS, nperseg=100, noverlap=50)
        P.append(p)
    return f, np.mean(P, 0), len(P)


for k in ('front_top', 'front_all', 'pad', 'noise'):
    f, pr, nr = spec(ME[k], run)
    f, ps, ns = spec(ME[k], stl)
    bands = [(0.5, 3), (3, 8), (8, 15), (15, 25)]
    print('%-10s (%d run segs, %d still segs) power run/still by band: ' % (k, nr, ns) + '  '.join('%g-%g Hz: %.2f' % (lo, hi, pr[(f >= lo) & (f < hi)].sum() / ps[(f >= lo) & (f < hi)].sum()) for lo, hi in bands))

print('\n== (c) effect size (Cohen d) running - still, lick-free frames, 0.06 s smoothed')
for k in ('old_index', 'jaw', 'pad', 'front_all', 'front_top', 'noise'):
    s = gaussian_filter1d(np.asarray(ME[k], np.float64), 3)
    a, b = s[run], s[stl]
    d = (a.mean() - b.mean()) / np.sqrt(0.5 * (a.var() + b.var()))
    print('%-10s run %.3f  still %.3f  d = %.2f' % (k, a.mean(), b.mean(), d))

print('\n== (d) stability of the matched lick effect (windows of 1 s, strata = speed decile x pupil tercile)')
WL = int(FPS)
nw = n // WL


def wmean(x):
    return np.nanmean(np.asarray(x[:nw * WL], np.float64).reshape(nw, WL), axis=1)


def build(lick_min=2, margin=1.0):
    mg = int(margin * FPS)
    nearm = np.convolve(imp, np.ones(2 * mg + 1), 'same') > 0
    w = pd.DataFrame(dict(speed=wmean(speed), pupil=wmean(pup_s), nlick=imp[:nw * WL].reshape(nw, WL).sum(1),
                          clear=~(nearm[:nw * WL].reshape(nw, WL).any(1)), art=art[:nw * WL].reshape(nw, WL).any(1), half=(np.arange(nw) >= nw // 2).astype(int)))
    w['lick'] = w.nlick >= lick_min
    w['sbin'] = pd.qcut(w.speed.rank(method='first'), 10, labels=False)
    w['pbin'] = pd.qcut(w.pupil.rank(method='first'), 3, labels=False)
    w['stratum'] = w.sbin * 3 + w.pbin
    return w[~w.art]


def eff(w, sig, sel=None):
    d = w.assign(v=wmean(sig)[w.index])
    if sel is not None:
        d = d[sel(d)]
    num = den = 0.0
    for s, g in d.groupby('stratum'):
        a, b = g[g.lick].v, g[g.clear].v
        if len(a) >= 3 and len(b) >= 3:
            num += len(a) * (a.mean() - b.mean())
            den += len(a)
    c = d[d.clear]
    loc = c[c.sbin == 9].v.mean() - c[c.sbin == 0].v.mean()
    return 100 * (num / den) / loc if den else np.nan


w2 = build(2, 1.0)
print('%-10s %8s %8s %8s %10s %10s' % ('signal', 'all', '1st half', '2nd half', '>=5 licks', 'margin .5s') + '   (lick effect as % of the locomotion contrast)')
w5 = build(5, 1.0)
w05 = build(2, 0.5)
for k in ('old_index', 'jaw', 'pad', 'front_all', 'front_low', 'front_top', 'noise'):
    print('%-10s %8.1f %8.1f %8.1f %10.1f %10.1f' % (k, eff(w2, ME[k]), eff(w2, ME[k], lambda d: d.half == 0), eff(w2, ME[k], lambda d: d.half == 1), eff(w5, ME[k]), eff(w05, ME[k])))

print('\n== (e) where the old index lick effect comes from (matched DC, index units)')
e_all = eff(w2, ME['old_index'])
tot = {}
for k, wgt in (('front_all', 0.5 / REST['front_all']), ('jaw', 0.5 / REST['jaw'])):
    d = w2.assign(v=wmean(ME[k])[w2.index])
    num = den = 0.0
    for s, g in d.groupby('stratum'):
        a, b = g[g.lick].v, g[g.clear].v
        if len(a) >= 3 and len(b) >= 3:
            num += len(a) * (a.mean() - b.mean())
            den += len(a)
    tot[k] = wgt * num / den
print({k: round(v, 3) for k, v in tot.items()}, ' share of jaw %.0f%%' % (100 * tot['jaw'] / sum(tot.values())))
