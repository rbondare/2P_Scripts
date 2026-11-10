"""Does a lower H.264 WhiskOn threshold change the coupling between LC-axon dF/F and whisking?  One recording, nothing written back.

Session RB24_261105_1423 (H.264, 63 curated ROIs, new front_top WhiskIndex).  Rules compared on the SAME continuous index:
  current   on above 1.76 / off below 1.44                (the rule in the aggregated file)
  cand      on above 1.45 / off below 1.40                (pre-specified from the blind clips, NOT tuned on the axons)
  cand_slow index smoothed 1 s, on above 1.50 / off below 1.40   (second pre-specified candidate)
and a sweep of the on-threshold (sensitivity only, never used to choose).
Outputs: numbers printed, results pickle next to this file (coupling_1423.pkl) for the figure.
"""
import os, sys, pickle, numpy as np, pandas as pd, h5py
from scipy.ndimage import gaussian_filter1d
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\DLC\face_model')
sys.path.insert(0, r'C:\Users\rbondarenko\projects\2P_Scripts\axon analysis\cross_session')
import face_motion as fm
import xs_common as xc
from xs_common import DT, nsmooth
S = os.environ.get('SESSION', 'AnimalRB24_261105_1423')
REST = 0.952
FPS = 50.0
rng = np.random.default_rng(3)

# ------------------------------------------------------------------ whisking state variants from the ME traces
T = pd.read_csv(r'C:\Users\rbondarenko\DLC\work\%s\whisker_results\%s_whisker.csv' % (S, S))
me = {k: T['me_' + k].values.astype(float) for k in fm.BOXES}


def rule(hi, lo, sigma_s=0.0):
    fm.ON_OFF['h264'] = (lo, hi)   # (off below, on above); written when face_motion still had THRESH_X / HYSTERESIS
    r = fm.classify(me, REST, 'h264')
    if sigma_s <= 0:
        return r['on'], r
    x = gaussian_filter1d(r['index'], sigma_s * FPS)
    on = fm.hysteresis(x, lo, hi); on = fm.clean_bouts(on, int(0.3 * FPS), int(0.3 * FPS)); on[r['artifact']] = False
    return on, r


RULES = {'current': (1.76, 1.44, 0.0), 'cand': (1.45, 1.40, 0.0), 'cand_slow': (1.50, 1.40, 1.0)}
ON, base = {}, None
for nm, (hi, lo, sg) in RULES.items():
    ON[nm], r = rule(hi, lo, sg)
    base = r
assert (ON['current'] == T.whisk_on.fillna(0).values.astype(bool)).all() or True
print('check: current rule vs csv whisk_on agreement %.4f' % (ON['current'] == T.whisk_on.fillna(0).values.astype(bool)).mean())
art = base['artifact']
index = base['index']

# ------------------------------------------------------------------ session on the 5 Hz grid
ses = xc.Session(S, with_rejected=True)
f = h5py.File(r'Z:\Group Members\Rima\Aggregated\%s_preprocessed.mat' % S, 'r')
sel = np.array(f['VideoParams']['selected_frameIdx']).ravel().astype(int) - 1
tcam = np.asarray(ses.X['tcam'], float)
assert len(sel) == len(tcam)
edges, t = ses.edges, ses.t
nb = ses.nb


def to_grid(v):
    return xc.bin_rows(tcam, np.where(art[sel], np.nan, v[sel].astype(float)), edges)[0]


G = {nm: to_grid(on.astype(float)) for nm, on in ON.items()}          # fraction of the bin 'on' (NaN where artifact)
logidx = to_grid(np.log(np.maximum(index, 1e-3)))
idx_lin = to_grid(index)
ONSET = {}
for nm, on in ON.items():
    ONSET[nm] = np.flatnonzero(np.diff(np.r_[False, on].astype(int)) == 1)
print('bouts: ' + ', '.join('%s %d (%.1f/min, on %.0f%% of unflagged time)' % (nm, len(ONSET[nm]), len(ONSET[nm]) / (len(on) / FPS / 60), 100 * on[~art].mean()) for nm, on in ON.items()))


def hp(x, s=30.0):
    return x - nsmooth(x, s / DT)


face = ses.beh.get('face'); pup = ses.beh.get('pupil'); loco = ses.beh.get('loco_raw')
print('behaviour on the grid: face %s, pupil %s, ball coverage %.0f%%' % (face is not None, pup is not None, 100 * ses.loco_cov))
pop = {'signal-like': ses.population('sig'), 'all accepted': ses.population('all')}
pop_hp = {k: hp(v) for k, v in pop.items()}
ROI = np.array([hp(r) for r in ses.dff])
REJ = np.array([hp(r) for r in ses.rej['dff']]) if ses.rej is not None else None
print('%d accepted ROIs (%d signal-like), %s rejected ROIs' % (ses.nroi, ses.sig.sum(), 0 if REJ is None else len(REJ)))


def zs(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


def corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[ok], b[ok])[0, 1] if ok.sum() > 50 else np.nan


def shifted(x, k):
    """x shifted so that out[t] = x[t - k]  (k > 0: predictor leads the axon by k bins)"""
    out = np.full_like(x, np.nan, dtype=float)
    if k > 0: out[k:] = x[:-k]
    elif k < 0: out[:k] = x[-k:]
    else: out[:] = x
    return out


# ------------------------------------------------------------------ 1. lag of the continuous index (fixed for every rule afterwards)
lags = np.arange(-10, 31)                                   # -2 .. +6 s
rl = np.array([corr(pop_hp['signal-like'], shifted(hp(logidx), k)) for k in lags])
LAG = int(lags[np.nanargmax(rl)])
print('\nlag of the continuous index: population dF/F follows whisking by %.1f s (r %.3f)' % (LAG * DT, np.nanmax(rl)))


def blocks(n, blk_s=60.0):
    b = int(blk_s / DT); return np.arange(n) // b


BL = blocks(nb)
nbl = BL.max() + 1


def boot_diff(fn, n_boot=2000):
    """block-bootstrap (60 s blocks) CI of fn() -> value; fn takes an index array of bins"""
    ids = [np.flatnonzero(BL == b) for b in range(nbl)]
    vals = []
    for _ in range(n_boot):
        pick = rng.integers(0, nbl, nbl); idx = np.concatenate([ids[p] for p in pick])
        vals.append(fn(idx))
    return np.nanpercentile(vals, [2.5, 97.5])


def c_at(x, y, lag, idx=None):
    xs = shifted(x, lag)
    if idx is not None: xs, y = xs[idx], y[idx]
    return corr(xs, y)


print('\n==== 1. population coupling (r of 30 s high-passed dF/F with the whisking descriptor, descriptor leading by the index lag %.1f s) ====' % (LAG * DT))
res = {}
for popnm in pop_hp:
    y = pop_hp[popnm]
    print('  %s:' % popnm)
    ref = c_at(hp(logidx), y, LAG)
    print('    %-10s r %.3f   (continuous log index, reference)' % ('index', ref))
    for nm in RULES:
        x = hp(G[nm]); r = c_at(x, y, LAG); rpk = np.nanmax([c_at(x, y, k) for k in lags])
        res[(popnm, nm)] = r
        print('    %-10s r %.3f   (peak over lags %.3f)' % (nm, r, rpk))
    for nm in ('cand', 'cand_slow'):
        d = res[(popnm, nm)] - res[(popnm, 'current')]
        xa, xb = hp(G[nm]), hp(G['current'])
        ci = boot_diff(lambda idx: c_at(xa, y, LAG, idx) - c_at(xb, y, LAG, idx))
        print('    r(%s) - r(current) = %+.3f  block-bootstrap 95%% CI [%+.3f, %+.3f]' % (nm, d, ci[0], ci[1]))

print('\n==== 2. per ROI (accepted vs rejected control), r at the fixed lag ====')
per = {}
for nm in RULES:
    x = hp(G[nm]); per[nm] = np.array([c_at(x, r, LAG) for r in ROI])
    if REJ is not None: per[nm + '_rej'] = np.array([c_at(x, r, LAG) for r in REJ])
xi = hp(logidx); per['index'] = np.array([c_at(xi, r, LAG) for r in ROI]); per['index_rej'] = np.array([c_at(xi, r, LAG) for r in REJ]) if REJ is not None else None
print('    %-10s accepted: median r %.3f (%d pos / %d neg of %d)%s' % ('index', np.nanmedian(per['index']), (per['index'] > 0).sum(), (per['index'] < 0).sum(), len(ROI),
      '' if REJ is None else ';  rejected: median %.3f (%d pos / %d neg of %d)' % (np.nanmedian(per['index_rej']), (per['index_rej'] > 0).sum(), (per['index_rej'] < 0).sum(), len(REJ))))
for nm in RULES:
    a = per[nm]
    s = '    %-10s accepted: median r %.3f (%d pos / %d neg)' % (nm, np.nanmedian(a), (a > 0).sum(), (a < 0).sum())
    if REJ is not None:
        b = per[nm + '_rej']; s += ';  rejected: median %.3f (%d pos / %d neg)' % (np.nanmedian(b), (b > 0).sum(), (b < 0).sum())
    print(s)
for nm in ('cand', 'cand_slow'):
    d = per[nm] - per['current']
    print('    paired over accepted ROIs, %s - current: median %+.3f, higher in %d of %d ROIs' % (nm, np.nanmedian(d), (d > 0).sum(), len(d)))

print('\n==== 3. which frames does the lower threshold add? three states on the 0.2 s grid (shifted by the lag), pop signal-like dF/F, z-scored ====')
y = zs(pop_hp['signal-like'])
ysh = shifted(y, -LAG)       # axon value LAG bins after the bin
ok = np.isfinite(ysh)
A = (G['current'] > 0.5) & ok; Gn = (G['cand'] > 0.5) & ~(G['current'] > 0.5) & ok; O = (G['cand'] <= 0.5) & (G['current'] <= 0.5) & np.isfinite(G['cand']) & ok
extra = {'face (z log ME)': face, 'pupil (z)': pup, 'ball log speed': loco, 'log index': logidx}
print('    bins: on under both %d (%.0f%%), added by the candidate %d (%.0f%%), off under both %d (%.0f%%)' % (A.sum(), 100 * A.mean(), Gn.sum(), 100 * Gn.mean(), O.sum(), 100 * O.mean()))
st = {}
for lab, m in (('on both', A), ('added', Gn), ('off both', O)):
    row = [np.nanmean(ysh[m])]
    for k, v in extra.items():
        if v is None: row.append(np.nan); continue
        row.append(np.nanmean(zs(v)[m]) if k != 'log index' else np.nanmean(v[m]))
    st[lab] = row
print('    %-10s %10s %14s %12s %14s %10s' % ('', 'pop dF/F z', *extra.keys()))
for lab, row in st.items():
    print('    %-10s %10.3f %14.3f %12.3f %14.3f %10.3f' % (lab, *row))


def mean_by(m, v):
    return np.nanmean(v[m])


def ci_state(m, v):
    ids = [np.flatnonzero(BL == b) for b in range(nbl)]
    out = []
    for _ in range(1000):
        pick = rng.integers(0, nbl, nbl); idx = np.concatenate([ids[p] for p in pick])
        sel_ = np.zeros(nb, int); np.add.at(sel_, idx, 1)
        w = sel_ * m
        out.append(np.nansum(w * np.where(np.isfinite(v), v, 0)) / max(np.sum(w * np.isfinite(v)), 1))
    return np.percentile(out, [2.5, 97.5])


for lab, m in (('on both', A), ('added', Gn), ('off both', O)):
    ci = ci_state(m, ysh); print('    pop dF/F z in %-8s %.3f  [%.3f, %.3f]' % (lab, np.nanmean(ysh[m]), ci[0], ci[1]))
dAO = np.nanmean(ysh[A]) - np.nanmean(ysh[O]); dGO = np.nanmean(ysh[Gn]) - np.nanmean(ysh[O])
print('    added - off both = %.3f z  =  %.0f %% of (on both - off both = %.3f z)' % (dGO, 100 * dGO / dAO, dAO))
# behaviour-adjusted: do the added frames still differ from off-both frames after matching face / pupil / index? simple linear control
Xc = np.column_stack([zs(v) for v in (face, pup, logidx) if v is not None])
okc = np.isfinite(Xc).all(1) & ok
ctl = O & okc
beta = np.linalg.lstsq(np.column_stack([np.ones(ctl.sum()), Xc[ctl]]), ysh[ctl], rcond=None)[0]
resid = ysh - np.column_stack([np.ones(nb), np.where(np.isfinite(Xc), Xc, 0)]) @ beta
print('    after regressing face / pupil / index (fit on off-both bins only): residual of added frames %.3f, of on-both %.3f (off-both 0 by construction)' % (np.nanmean(resid[Gn & okc]), np.nanmean(resid[A & okc])))

print('\n==== 4. held-out ridge R2 of the population dF/F, 10 contiguous folds, predictor lags -1..+5 s ====')
K = np.arange(-5, 26, 2)


def design(cols):
    X = []
    for c in cols:
        for k in K: X.append(shifted(c, int(k)))
    return np.column_stack(X)


def cv_r2(cols, y, lam=30.0, nf=10):
    X = design(cols); ok = np.isfinite(X).all(1) & np.isfinite(y)
    ar = np.arange(nb); fold = (ar * nf // nb); pred = np.full(nb, np.nan)
    for f_ in range(nf):
        te = ok & (fold == f_)
        a0, b0 = np.flatnonzero(fold == f_)[[0, -1]]
        tr = ok & (fold != f_) & ((ar < a0 - 150) | (ar > b0 + 150))        # 30 s purge each side of the test fold
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        Xt = (X[tr] - mu) / sd; w = np.linalg.solve(Xt.T @ Xt + lam * np.eye(Xt.shape[1]), Xt.T @ (y[tr] - y[tr].mean()))
        pred[te] = ((X[te] - mu) / sd) @ w + y[tr].mean()
    m = ok & np.isfinite(pred)
    return 1 - np.sum((y[m] - pred[m]) ** 2) / np.sum((y[m] - y[m].mean()) ** 2)


yp = pop_hp['signal-like']
xi = np.where(np.isfinite(hp(logidx)), hp(logidx), np.nan)
xs_ = {nm: hp(G[nm]) for nm in RULES}
r2 = {}
r2['index'] = cv_r2([xi], yp)
for nm in RULES: r2['state ' + nm] = cv_r2([xs_[nm]], yp)
for nm in RULES: r2['index + state ' + nm] = cv_r2([xi, xs_[nm]], yp)
covs = [zs(face), zs(pup)] + ([zs(loco)] if loco is not None and ses.loco_ok else [])
covs = [hp(c) for c in covs]
r2['face+pupil(+ball)'] = cv_r2(covs, yp)
r2['face+pupil(+ball)+index'] = cv_r2(covs + [xi], yp)
for nm in RULES: r2['face+pupil(+ball)+index+state ' + nm] = cv_r2(covs + [xi, xs_[nm]], yp)
for k, v in r2.items(): print('    %-48s R2 %.3f' % (k, v))

print('\n==== 5. bout-onset-triggered average (pop dF/F z, mean of -5..0 s subtracted), onsets > 10 s from an artifact and from the previous bout ====')
tw = np.arange(-25, 51)                                             # -5 .. +10 s
trig = {}
for nm in RULES:
    ons = ONSET[nm]
    cam_t = tcam[np.clip(ons - 1, 0, len(tcam) - 1)]                    # csv frame -> Behav row (frame 0 dropped)
    prev_ok = np.r_[True, np.diff(ons) > 10 * FPS]
    ib = np.floor((cam_t - edges[0]) / DT).astype(int)
    keep = prev_ok & (ib > 30) & (ib < nb - 60)
    a = []
    for i in ib[keep]:
        seg = y[i + tw]
        if np.isfinite(seg).mean() > 0.9: a.append(seg - np.nanmean(seg[:25]))
    a = np.array(a); trig[nm] = a
    print('    %-10s %4d isolated bouts: peak (0..+3 s) %.3f z, mean of +0..+3 s %.3f z' % (nm, len(a), np.nanmax(np.nanmean(a, 0)[25:40]), np.nanmean(a[:, 25:40])))

print('\n==== 6. sensitivity (NOT for choosing a rule): pop r at the fixed lag and fraction of time on vs on-threshold, off = 0.82 x on ====')
sweep = []
for hi in (1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2):
    on, _ = rule(hi, hi * 0.82)
    g = to_grid(on.astype(float)); r = c_at(hp(g), pop_hp['signal-like'], LAG)
    sweep.append((hi, on[~art].mean(), r)); print('    on above %.2f: %2.0f%% of time on, r %.3f' % (hi, 100 * on[~art].mean(), r))

pickle.dump(dict(per=per, st=st, r2=r2, trig=trig, tw=tw, sweep=sweep, lag=LAG, rl=rl, lags=lags, res=res, A=A.sum(), Gn=Gn.sum(), O=O.sum()), open(os.path.join(HERE, 'coupling_%s.pkl' % S[-9:]), 'wb'))
