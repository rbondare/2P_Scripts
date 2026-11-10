"""How does licking contaminate the face-motion-energy (whisking) index, and which correction works best?  Session RB24_261109_1506 (lick task, lick line aux3 = independent lick times).
usage: python lick_whisk.py          (needs align_*.npz, cells_*.npy, ball_*.npy, lick_clean_*.npy, whisker csv, eye csv)
Prints a results table and writes lw_cache.npz + figures LW1_*.png, LW2_*.png next to the script."""
import os, json, glob, numpy as np, pandas as pd, cv2
from scipy.ndimage import gaussian_filter1d, uniform_filter1d, median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
sess = 'AnimalRB24_261109_1506'
FPS = 50.0
rng = np.random.default_rng(1)

# ---------------------------------------------------------------- data
A = np.load(os.path.join(here, 'align_%s.npz' % sess))
tcam, tball = A['aux1'], A['aux2']
lk = np.load(os.path.join(here, 'lick_clean_%s.npy' % sess))
n = len(tcam)
cells = np.load(os.path.join(here, 'cells_%s.npy' % sess), mmap_mode='r')
meta = json.load(open(os.path.join(here, 'cells_%s.json' % sess)))
CS, GX0, GY0, NX, NY = meta['CS'], meta['GX'][0], meta['GY'][0], meta['NX'], meta['NY']
W = pd.read_csv(r'C:\Users\rbondarenko\DLC\work\%s\whisker_results\%s_whisker.csv' % (sess, sess))
art = W['artifact'].values.astype(bool)
assert len(W) == n == cells.shape[0], (len(W), n, cells.shape)


def reg(x0, y0, x1, y1):
    """mean ME over the cells of a region given in px relative to the nose anchor (multiples of the 20 px grid)"""
    ix0, ix1, iy0, iy1 = int((x0 - GX0) // CS), int((x1 - GX0) // CS), int((y0 - GY0) // CS), int((y1 - GY0) // CS)
    out = np.empty(n, np.float32)
    for k in range(0, n, 20000):
        out[k:k + 20000] = np.nanmean(cells[k:k + 20000, iy0:iy1, ix0:ix1], axis=(1, 2))
    return out


REG = {
    'front_all': (-260, -120, 60, 200),     # old fan_front
    'front_top': (-260, -120, 60, 0),       # above the spout
    'front_low': (-260, 0, 60, 200),        # spout shaft + tip
    'jaw': (60, 80, 300, 300),              # old fan_down: chin, lower jaw, tongue
    'pad': (60, -100, 260, 60),
    'noise': (-320, -160, -200, -60),       # empty background corner of the grid
}
ME = {k: reg(*v) for k, v in REG.items()}
print('region mean ME:', {k: round(float(np.nanmean(v)), 2) for k, v in ME.items()})
REST = {'front_all': 0.447, 'jaw': 0.450}
old_index = 0.5 * (gaussian_filter1d(ME['front_all'], 3) / REST['front_all'] + gaussian_filter1d(ME['jaw'], 3) / REST['jaw'])
print('own old-style index vs csv whisk_index: corr %.3f, mean %.2f vs %.2f' % (np.corrcoef(old_index, W.whisk_index.values)[0, 1], old_index.mean(), W.whisk_index.mean()))

# speed proxy (ball sensor 1 increments) on the camera clock
ball = np.load(os.path.join(here, 'ball_%s.npy' % sess))
nb = len(ball)
tb = tball[:nb]
mag = np.hypot(ball[:, 0], ball[:, 1])
dtb = np.median(np.diff(tb))
sp_b = uniform_filter1d(mag, int(round(0.5 / dtb))) / dtb
speed = np.interp(tcam, tb, sp_b).astype(np.float32)
# pupil (DLC eye ellipse, area ~ w * h) smoothed
E = pd.read_csv(glob.glob(r'C:\Users\rbondarenko\DLC\work\%s\eye_results\*_eye.csv' % sess)[0], skiprows=2)
pup = (np.pi / 4 * E['pupil_width'].values * E['pupil_height'].values).astype(float)
print('eye csv rows', len(E), 'camera frames', n)
pup = pup[:n] if len(pup) >= n else np.r_[pup, np.full(n - len(pup), np.nan)]
pup = pd.Series(pup).interpolate(limit_direction='both').values
pup_s = median_filter(pup, size=int(2 * FPS) | 1)
lick_f = np.searchsorted(tcam, lk)
lick_f = lick_f[(lick_f > 0) & (lick_f < n)]
imp = np.zeros(n, np.float32)
np.add.at(imp, lick_f, 1)
# bout indicator: first lick - 0.3 s ... last lick + 0.5 s of every bout (gap >= 1 s)
ili = np.diff(lk)
b0 = np.flatnonzero(np.r_[True, ili >= 1.0])
b1 = np.r_[b0[1:], len(lk)] - 1
bout = np.zeros(n, bool)
for a, b in zip(b0, b1):
    fa, fb = np.searchsorted(tcam, [lk[a] - 0.3, lk[b] + 0.5])
    bout[fa:fb] = True
near1 = np.convolve(imp, np.ones(int(2 * FPS) + 1), 'same') > 0     # lick within +-1 s
print('frames: bout %.3f, lick within +-1 s %.3f, artifact %.3f' % (bout.mean(), near1.mean(), art.mean()))

# ---------------------------------------------------------------- 1 s windows for matched comparisons
WL = int(FPS)
nw = n // WL


def wmean(x):
    return np.nanmean(x[:nw * WL].reshape(nw, WL), axis=1)


win = pd.DataFrame(dict(speed=wmean(speed), pupil=wmean(pup_s), nlick=imp[:nw * WL].reshape(nw, WL).sum(1), clear=~(near1[:nw * WL].reshape(nw, WL).any(1)), art=art[:nw * WL].reshape(nw, WL).any(1),
                        chunk=(np.arange(nw) // 60)))
win['lick'] = win.nlick >= 2
ok = ~win.art
win['sbin'] = pd.qcut(win.speed.rank(method='first'), 10, labels=False)
win['pbin'] = pd.qcut(win.pupil.rank(method='first'), 3, labels=False)
win['stratum'] = win.sbin * 3 + win.pbin
print('windows %d | lick (>= 2 licks) %d | lick-free (no lick within +-1 s) %d' % (nw, (win.lick & ok).sum(), (win.clear & ok).sum()))


def matched_effect(sig, boot=0):
    """mean(sig | lick window) - mean(sig | lick-free window), inside strata of ball speed decile x pupil tercile (weights = n lick windows); also the
    locomotion contrast = lick-free mean in the fastest speed decile minus the slowest"""
    w = wmean(sig)
    d = win.assign(v=w)
    d = d[ok]

    def calc(dd):
        num = den = 0.0
        for s, g in dd.groupby('stratum'):
            a, b = g[g.lick].v, g[g.clear].v
            if len(a) >= 3 and len(b) >= 3:
                num += len(a) * (a.mean() - b.mean())
                den += len(a)
        c = dd[dd.clear]
        loc = c[c.sbin == 9].v.mean() - c[c.sbin == 0].v.mean()
        return (num / den if den else np.nan), loc

    e, loc = calc(d)
    if not boot:
        return e, loc
    ch = d.chunk.unique()
    groups = {c: g for c, g in d.groupby('chunk')}
    es = []
    for _ in range(boot):
        pick = rng.choice(ch, len(ch))
        es.append(calc(pd.concat([groups[c] for c in pick]))[0])
    return e, loc, np.nanpercentile(es, [2.5, 97.5])


# ---------------------------------------------------------------- lick-triggered average + null
LAGS = np.arange(-int(0.4 * FPS), int(0.6 * FPS) + 1)


def lta(sig, events, frames_ok=None):
    ev = events[(events + LAGS[0] >= 0) & (events + LAGS[-1] < n)]
    if frames_ok is not None:
        ev = ev[frames_ok[ev]]
    M = sig[ev[:, None] + LAGS[None, :]]
    return np.nanmean(M, axis=0), len(ev)


def amp(w, sig):
    sel = (LAGS >= -int(0.2 * FPS)) & (LAGS <= int(0.4 * FPS))
    return (w[sel].max() - w[sel].min()) / np.nanmean(sig[~art])


def null_amp(sig, nsh=40, test_mask=None):
    a = []
    for _ in range(nsh):
        sh = int(rng.uniform(30, 600) * FPS) * rng.choice([-1, 1])
        ev = (lick_f + sh) % n
        w, _ = lta(sig, ev, test_mask)
        a.append(amp(w, sig))
    return np.percentile(a, 95)


# ---------------------------------------------------------------- kernel subtraction (phase-locked lick component), 2-fold block cross-validation
KL = np.arange(-10, 26)


def design(use_speed=True):
    cols = [np.roll(imp, k) for k in KL]
    for c, k in zip(cols, KL):          # np.roll wraps; blank the wrapped part
        if k > 0:
            c[:k] = 0
        elif k < 0:
            c[k:] = 0
    X = np.stack(cols, 1)
    Z = [bout.astype(np.float32), np.ones(n, np.float32)]
    if use_speed:
        sz = (speed - speed.mean()) / speed.std()
        pz = (pup_s - pup_s.mean()) / pup_s.std()
        Z += [sz, sz ** 2, pz]
    return X, np.stack(Z, 1)


Xk, Zn = design(True)
blk = (np.arange(n) // int(300 * FPS)) % 2        # alternate 5-min blocks
good = ~art


def fit_kernel(y, train):
    D = np.hstack([Xk, Zn])[train & good]
    yy = y[train & good]
    lam = 1e-3 * np.trace(D.T @ D) / D.shape[1]
    beta = np.linalg.solve(D.T @ D + lam * np.eye(D.shape[1]), D.T @ yy)
    return beta[:len(KL)]


def clean_cv(y):
    """lick-kernel part removed; each half is cleaned with the kernel fitted on the OTHER half"""
    out = y.copy()
    for h in (0, 1):
        k = fit_kernel(y, blk != h)
        out[blk == h] = y[blk == h] - (Xk[blk == h] @ k)
    return out


def clean_full(y):
    k = fit_kernel(y, np.ones(n, bool))
    return y - Xk @ k, k


# ---------------------------------------------------------------- signals to compare
SIG = {
    'old index (front_all + jaw)': 0.5 * (ME['front_all'] / REST['front_all'] + ME['jaw'] / REST['jaw']),
    'jaw': ME['jaw'], 'pad': ME['pad'], 'front_all': ME['front_all'], 'front_top': ME['front_top'], 'front_low': ME['front_low'], 'noise (control)': ME['noise'],
}
SIG['front_all, kernel-cleaned'] = clean_cv(ME['front_all'])
SIG['front_top, kernel-cleaned'] = clean_cv(ME['front_top'])
SIG['old index, kernel-cleaned'] = clean_cv(SIG['old index (front_all + jaw)'])
SIG['jaw, kernel-cleaned'] = clean_cv(ME['jaw'])
SIG['pad, kernel-cleaned'] = clean_cv(ME['pad'])

lick_ev = lick_f
test_all = np.ones(n, bool)
rows = []
LT = {}
for name, s in SIG.items():
    s = np.asarray(s, np.float64)
    w, ne = lta(s, lick_ev)
    LT[name] = w
    a = amp(w, s)
    na = null_amp(s)
    e, loc, ci = matched_effect(s, boot=150)
    sm = gaussian_filter1d(s, 3)
    cl = ~near1 & good
    run = cl & (speed > np.percentile(speed, 80))
    stl = cl & (speed < np.percentile(speed, 20))
    dr = np.nanmedian(sm[run]) - np.nanmedian(sm[stl])
    c_pad = np.corrcoef(gaussian_filter1d(ME['pad'], 3)[cl], sm[cl])[0, 1] if 'pad' not in name else np.nan
    rows.append(dict(signal=name, mean=np.nanmean(s[good]), lta_amp=a, null95=na, lick_dc=e, lick_dc_lo=ci[0], lick_dc_hi=ci[1], loco_contrast=loc,
                     lick_pct_of_loco=100 * e / loc if loc else np.nan, run_minus_still=dr, corr_pad_lickfree=c_pad))
R = pd.DataFrame(rows).set_index('signal')
pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 30)
print(R.round(3).to_string())
R.to_csv(os.path.join(here, 'lw_table.csv'))
np.savez(os.path.join(here, 'lw_cache.npz'), LAGS=LAGS, **{'lta_' + k.replace(' ', '_'): v for k, v in LT.items()})

# ---------------------------------------------------------------- lick-locked kernel of every region (full-data fit) and lick-sensitivity map
Kfull = {k: clean_full(np.asarray(ME[k], np.float64))[1] for k in ('front_all', 'front_top', 'front_low', 'jaw', 'pad', 'noise')}
# cell map: mean ME of cells in lick windows / lick-free windows within strata (log ratio) + lick-locked amplitude (LTA peak-to-peak / mean)
ratio = np.full((NY, NX), np.nan, np.float32)
ptp = np.full((NY, NX), np.nan, np.float32)
okw = ok.values
strat = win.stratum.values
lickw = win.lick.values
clearw = win.clear.values
w_l = np.zeros(nw, bool)
for r in range(NY):
    cr = np.asarray(cells[:, r, :], np.float32)      # [n, NX]
    for c in range(NX):
        v = cr[:, c]
        wv = np.nanmean(v[:nw * WL].reshape(nw, WL), axis=1)
        num = den = 0.0
        base = 0.0
        for s in np.unique(strat):
            m = (strat == s) & okw
            a, b = wv[m & lickw], wv[m & clearw]
            if len(a) >= 3 and len(b) >= 3 and b.mean() > 0:
                num += len(a) * np.log(a.mean() / b.mean())
                den += len(a)
        ratio[r, c] = num / den if den else np.nan
        w, _ = lta(v, lick_f)
        sel = (LAGS >= -int(0.2 * FPS)) & (LAGS <= int(0.4 * FPS))
        ptp[r, c] = (w[sel].max() - w[sel].min()) / max(np.nanmean(v[good]), 1e-3)
np.savez(os.path.join(here, 'lw_maps.npz'), ratio=ratio, ptp=ptp)

# ---------------------------------------------------------------- figures
vid = glob.glob(r'C:\Users\rbondarenko\DLC\work\%s\camera\*.mp4' % sess)[0]
cap = cv2.VideoCapture(vid)
fi = int(lick_f[len(lick_f) // 3] + 4)
cap.set(cv2.CAP_PROP_POS_FRAMES, fi)
okk, fr = cap.read()
axn, ayn = meta['anchor_median']
ex0, ey0 = axn + GX0, ayn + GY0
fig, ax = plt.subplots(2, 3, figsize=(19, 10.5))
ext = (GX0, GX0 + NX * CS, GY0 + NY * CS, GY0)
if okk:
    img = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY)[ey0:ey0 + NY * CS, ex0:ex0 + NX * CS]
    ax[0, 0].imshow(img, cmap='gray', extent=ext)
ax[0, 0].set_title('example frame (lick frame), regions')
cols = {'front_top': 'lime', 'front_low': 'yellow', 'jaw': 'magenta', 'pad': 'orange', 'noise': 'cyan'}
for k, c in cols.items():
    x0, y0, x1, y1 = REG[k]
    for a_ in (ax[0, 0], ax[0, 1], ax[0, 2]):
        a_.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=c, lw=1.8))
m1 = ax[0, 1].imshow(ratio, extent=ext, cmap='RdBu_r', vmin=-1, vmax=1)
ax[0, 1].set_title('licking vs lick-free, matched ball speed x pupil:\nlog ratio of cell ME (red = more motion with licks)')
plt.colorbar(m1, ax=ax[0, 1], fraction=0.046)
m2 = ax[0, 2].imshow(ptp, extent=ext, cmap='magma', vmin=0, vmax=np.nanpercentile(ptp, 99))
ax[0, 2].set_title('lick-LOCKED modulation: peak-to-peak of the lick-triggered\naverage / mean ME of the cell')
plt.colorbar(m2, ax=ax[0, 2], fraction=0.046)
for a_ in ax[0]:
    a_.plot(0, 0, 'r.')
t = LAGS / FPS * 1000
for k in ('jaw', 'pad', 'front_low', 'front_all', 'front_top', 'noise (control)'):
    s = np.asarray(SIG[k], np.float64)
    ax[1, 0].plot(t, LT[k] / np.nanmean(s[good]), label=k)
ax[1, 0].axvline(0, c='k', lw=.5)
ax[1, 0].set_xlabel('time from lick contact (ms)')
ax[1, 0].set_ylabel('ME / mean ME')
ax[1, 0].set_title('lick-triggered average (all %d licks)' % len(lick_f))
ax[1, 0].legend(fontsize=8)
for k in ('front_all', 'front_top', 'old index (front_all + jaw)'):
    s = np.asarray(SIG[k], np.float64)
    c = np.asarray(SIG[k + ', kernel-cleaned'] if k + ', kernel-cleaned' in SIG else SIG['old index, kernel-cleaned'], np.float64)
    l = ax[1, 1].plot(t, LT[k] / np.nanmean(s[good]), label=k)[0]
    kk = k + ', kernel-cleaned' if k + ', kernel-cleaned' in SIG else 'old index, kernel-cleaned'
    ax[1, 1].plot(t, LT[kk] / np.nanmean(s[good]), '--', c=l.get_color(), label=kk)
ax[1, 1].axvline(0, c='k', lw=.5)
ax[1, 1].set_xlabel('time from lick contact (ms)')
ax[1, 1].set_title('lick-triggered average before / after removing the lick kernel (held-out halves)')
ax[1, 1].legend(fontsize=7)
tk = np.array(KL) / FPS * 1000
for k in ('jaw', 'pad', 'front_all', 'front_top', 'noise'):
    ax[1, 2].plot(tk, Kfull[k], label=k)
ax[1, 2].axvline(0, c='k', lw=.5)
ax[1, 2].set_title('lick impulse kernels (gray levels per lick; fitted with bout, speed, pupil terms)')
ax[1, 2].legend(fontsize=8)
fig.suptitle('%s: how licking enters the motion-energy boxes' % sess)
fig.tight_layout()
fig.savefig(os.path.join(here, 'LW1_%s.png' % sess), dpi=65)
print('figures written')
