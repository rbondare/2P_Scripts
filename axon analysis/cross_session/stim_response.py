"""Cross-session analysis 2: stimulus-evoked LC-axon responses and what is left of them after the concurrent behaviour is accounted for.

usage: python stim_response.py [session ...]      (default: every session of STIM_SESSIONS; output stim_response_results.pkl next to the figures)

Population signal P = mean best dF/F of the signal-like ROIs of a session (0.2 s grid, high-passed > 30 s), behaviour = pupil / face motion / whisking / locomotion (ball only where plausible).
Per session and stimulus type (flash, fff, moving bar, expanding / receding loom, task square):
  raw response       mean P in the response window minus the baseline window, per trial; null = the onsets circularly shifted inside their own block (same spacing, wrong phase)
  behaviour          the same for pupil, face, whisking, locomotion (an equally large behaviour change = arousal / startle)
  adjusted response  a behaviour-only model ('mv' = face + whisking + locomotion, primary; 'all' adds pupil, which carries the light reflex of visual stimuli and is only a
                     secondary check) (ridge, behaviour at lags -1..+3 s) is trained on the data OUTSIDE every stimulus window [-3, +12 s]; its prediction for the stimulus
                     windows is subtracted from P; the residual response is tested against the same shifted-onset null. A residual response = more than what the concurrent
                     pupil / face / whisker / locomotion change predicts (conservative: part of an arousal response driven by the stimulus is removed as well)
  trials             Spearman(dF/F amplitude, evoked behaviour) and the partial correlation controlling for the other behaviours
  contrast           expanding minus receding disc (paired within sessions that have both), raw and adjusted
"""
import sys, os, pickle, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xs_common as xc
from xs_common import DT
from scipy.stats import spearmanr

STIM_SESSIONS = ['AnimalRB24_261021_1630', 'AnimalRB24_261023_1527', 'AnimalRB24_261023_1538', 'AnimalRB24_261027_1456', 'AnimalRB24_261027_1534', 'AnimalRB24_261029_1752',
                 'AnimalRB24_261030_1744', 'AnimalRB24_261104_1640', 'AnimalRB24_261104_1651', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261106_1317', 'AnimalRB24_261106_1529',
                 'AnimalRB23_261023_1642', 'AnimalRB23_261027_1306', 'AnimalRB23_261027_1423', 'AnimalRB23_261029_1704', 'AnimalRB22_261023_1427', 'AnimalRB22_261027_1138']
OUT = os.environ.get('STIM_OUT', os.path.join(r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session', 'stim_response_results.pkl'))
KEYS = ['flash', 'fff', 'bar', 'loom_expand', 'loom_recede', 'square']
LAB = {'flash': 'local/global flashes', 'fff': 'full-field flash', 'bar': 'moving bar', 'loom_expand': 'expanding loom', 'loom_recede': 'receding loom', 'square': 'square (task / passive)'}
WIN = {'flash': (-2, 6), 'fff': (-5, 35), 'bar': (-2, 8), 'loom_expand': (-3, 10), 'loom_recede': (-3, 10), 'square': (-5, 12)}
BWIN = {'flash': (-1, 0), 'fff': (0, 2), 'bar': (-1, 0), 'loom_expand': (-2, 0), 'loom_recede': (-2, 0), 'square': (-3, -0.5)}
AWIN = {'flash': (0, 2), 'fff': (2, 32), 'bar': (0.5, 4), 'loom_expand': (0, 4), 'loom_recede': (0, 4), 'square': (0, 3)}
N_NULL = 1000; HP_S = 30.0
rng = np.random.default_rng(3)
BEH_ORDER = ['pupil', 'face', 'whisk', 'loco']


def hp(x): return x - xc.nsmooth(x, HP_S / DT)
def zs(v): return (v - np.nanmean(v)) / np.nanstd(v)


def cum_of(x):
    ok = np.isfinite(x); y = np.where(ok, x, 0.0); return np.concatenate([[0], np.cumsum(y)]), np.concatenate([[0], np.cumsum(ok.astype(float))])


def win_mean(c, edges0, ons, w):
    cs, cn = c; n = len(cs) - 1
    i1 = np.floor((ons + w[0] - edges0) / DT).astype(int); i2 = np.floor((ons + w[1] - edges0) / DT).astype(int) + 1
    ok = (i1 >= 0) & (i2 <= n) & (i2 > i1); i1c = np.clip(i1, 0, n); i2c = np.clip(i2, 0, n)
    with np.errstate(invalid='ignore', divide='ignore'):
        m = (cs[i2c] - cs[i1c]) / (cn[i2c] - cn[i1c])
    m[~ok] = np.nan; return m


def amp(c, edges0, ons, key):
    return win_mean(c, edges0, ons, AWIN[key]) - win_mean(c, edges0, ons, BWIN[key])


def shift_in_blocks(ons, spans):
    o2 = ons.copy()
    for sp in np.unique(spans, axis=0):
        ii = np.all(spans == sp, axis=1); dur = sp[1] - sp[0]
        if not np.isfinite(dur) or dur <= 0: continue
        d = (0.2 + 0.6 * rng.random()) * dur; o2[ii] = sp[0] + np.mod(ons[ii] - sp[0] + d, dur)
    return o2


def align(x, edges0, ons, w):
    rel = np.arange(w[0], w[1] + DT / 2, DT); idx = np.floor((ons[:, None] + rel[None, :] - edges0) / DT).astype(int)
    out = np.full(idx.shape, np.nan, np.float32); ok = (idx >= 0) & (idx < len(x)); out[ok] = x[idx[ok]]; return rel, out


def fit_ridge(X, y, folds=5, purge=int(5 / DT), lams=(0.1, 1, 10, 100, 1000)):
    ok = np.isfinite(X).all(1) & np.isfinite(y); idx = np.where(ok)[0]; n = len(idx)
    if n < 400: return None
    Xt, yt = X[idx], y[idx]; mu, sd = Xt.mean(0), Xt.std(0) + 1e-9; Z = (Xt - mu) / sd; ym = yt.mean(); e = np.linspace(0, n, folds + 1).astype(int)
    best, bl = np.inf, lams[0]
    for lam in lams:
        sse = 0
        for f in range(folds):
            te = np.arange(e[f], e[f + 1]); tr = np.setdiff1d(np.arange(n), np.arange(max(0, e[f] - purge), min(n, e[f + 1] + purge)))
            w = np.linalg.solve(Z[tr].T @ Z[tr] + lam * len(tr) * np.eye(Z.shape[1]) / 100.0, Z[tr].T @ (yt[tr] - yt[tr].mean())); sse += ((Z[te] @ w + yt[tr].mean() - yt[te]) ** 2).sum()
        if sse < best: best, bl = sse, lam
    w = np.linalg.solve(Z.T @ Z + bl * n * np.eye(Z.shape[1]) / 100.0, Z.T @ (yt - ym))
    return mu, sd, ym, w, bl


def predict(model, X):
    mu, sd, ym, w, _ = model; out = np.full(len(X), np.nan); ok = np.isfinite(X).all(1); out[ok] = ((X[ok] - mu) / sd) @ w + ym; return out


def analyse(S):
    X = S.X; out = dict(name=S.name, n_sig=int(S.sig.sum()), n_roi=S.nroi, dur_min=S.dur / 60, loco_ok=S.loco_ok, keys={})
    if S.sig.sum() < 3: out['why'] = 'fewer than 3 signal-like ROIs'; return out
    on_t = np.atleast_1d(np.asarray(X['on_t'], float)); on_key = np.atleast_1d(np.asarray(X['on_key'])).astype(str); on_blk = np.atleast_1d(np.asarray(X['on_blk'], int))
    bspan = np.atleast_2d(np.asarray(X['blk_span'], float)); binter = np.atleast_1d(np.asarray(X['blk_interrupted'], float)).astype(bool)
    if len(on_t) == 0: out['why'] = 'no stimulus onsets'; return out
    keep = np.isfinite(on_t) & ~binter[on_blk - 1]; on_t, on_key, on_blk = on_t[keep], on_key[keep], on_blk[keep]
    P = hp(S.population('sig')); e0 = S.edges[0]
    behs = {b: hp(zs(S.beh[b])) for b in BEH_ORDER if b in S.beh and (b != 'loco' or S.loco_ok) and np.isfinite(S.beh[b]).mean() > 0.3}
    out['behaviours'] = list(behs)
    # behaviour-only model trained outside every stimulus window
    mask_train = np.isfinite(P).copy()
    for o in on_t:
        i1 = int(np.floor((o - 3 - e0) / DT)); i2 = int(np.floor((o + 12 - e0) / DT)) + 1; mask_train[max(i1, 0):max(i2, 0)] = False
    lags = (np.arange(-1.0, 3.01, 0.4) / DT).round().astype(int); out['models'] = {}; extra = {}
    for tag, fam in (('mv', [b for b in ('face', 'whisk', 'loco') if b in behs]), ('all', [b for b in ('pupil', 'face', 'whisk', 'loco') if b in behs])):
        if not fam: continue
        Xd = np.concatenate([np.stack([np.roll(behs[b], -l) for l in lags], 1) for b in fam], 1)
        model = fit_ridge(np.where(mask_train[:, None], Xd, np.nan), np.where(mask_train, P, np.nan))
        if model is None: continue
        Phat = predict(model, Xd); ok = mask_train & np.isfinite(Phat)
        out['models'][tag] = dict(families=fam, lam=model[4], train_r=float(np.corrcoef(Phat[ok], P[ok])[0, 1]))
        extra['Phat_' + tag] = Phat; extra['res_' + tag] = P - Phat
    series = {'P': P}; series.update(extra)
    for b, v in behs.items(): series['b_' + b] = v
    cums = {k: cum_of(v) for k, v in series.items()}
    for key in KEYS:
        m = (on_key == key)
        if m.sum() < 5: continue
        ons = on_t[m]; spans = bspan[on_blk[m] - 1]
        d = dict(n=int(m.sum()), n_blocks=int(len(np.unique(on_blk[m]))))
        rel, C = align(P, e0, ons, WIN[key]); d['rel'] = rel; d['P'] = C - np.nanmean(C[:, (rel >= BWIN[key][0]) & (rel < BWIN[key][1])], 1, keepdims=True)
        for nm in [k2 for k2 in series if k2.startswith(('Phat_', 'res_'))]:
            _, Cx = align(series[nm], e0, ons, WIN[key]); d[nm] = Cx - np.nanmean(Cx[:, (rel >= BWIN[key][0]) & (rel < BWIN[key][1])], 1, keepdims=True)
        d['beh'] = {}
        for b, v in behs.items(): _, Cb = align(v, e0, ons, WIN[key]); d['beh'][b] = Cb - np.nanmean(Cb[:, (rel >= BWIN[key][0]) & (rel < BWIN[key][1])], 1, keepdims=True)
        d['amp'] = {}; d['null'] = {}
        for nm, cc in cums.items():
            d['amp'][nm] = amp(cc, e0, ons, key)
            d['null'][nm] = np.array([np.nanmean(amp(cc, e0, shift_in_blocks(ons, spans), key)) for _ in range(N_NULL)])
        d['z'] = {nm: float((np.nanmean(d['amp'][nm]) - d['null'][nm].mean()) / d['null'][nm].std()) for nm in cums}
        # trial-wise: dF/F amplitude vs evoked behaviour (Spearman) and partial correlation
        tr = {}
        a = d['amp']['P']
        for b in behs:
            ab = d['amp']['b_' + b]; ok = np.isfinite(a) & np.isfinite(ab)
            if ok.sum() >= 8: tr[b] = dict(rho=float(spearmanr(a[ok], ab[ok])[0]), p=float(spearmanr(a[ok], ab[ok])[1]), n=int(ok.sum()))
        d['trial_corr'] = tr
        out['keys'][key] = d
    return out


def pooled(results, key, series='P'):
    """mean over sessions of the session-mean amplitude, with the null of that mean (independent shifts per session)"""
    ss = [r for r in results.values() if key in r.get('keys', {}) and series in r['keys'][key]['amp']]
    if not ss: return None
    A = np.array([np.nanmean(r['keys'][key]['amp'][series]) for r in ss]); N = np.stack([r['keys'][key]['null'][series] for r in ss], 1)
    nul = N.mean(1); m = A.mean(); p = (1 + (np.abs(nul) >= abs(m)).sum()) / (1 + len(nul))
    return dict(n_sessions=len(ss), mean=float(m), null95=(float(np.percentile(nul, 2.5)), float(np.percentile(nul, 97.5))), p=float(p), sessions=[r['name'] for r in ss], A=A)


if __name__ == '__main__':
    names = sys.argv[1:] or STIM_SESSIONS
    os.makedirs(os.path.dirname(OUT), exist_ok=True); res = {}
    for nm in names:
        t0 = time.time()
        try:
            S = xc.Session(nm, with_rejected=False); res[nm] = analyse(S)
            r = res[nm]; print('  -> %s | keys %s | %.0f s' % (r.get('why', 'ok'), {k: v['n'] for k, v in r['keys'].items()}, time.time() - t0), flush=True)
        except Exception as e:
            import traceback; traceback.print_exc(); res[nm] = dict(name=nm, keys={}, why='error: %s' % e)
        pickle.dump(res, open(OUT, 'wb'))
    print('STIM_DONE')
