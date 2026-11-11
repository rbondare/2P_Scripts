"""Lick-class encoding model (FIR) of the LC-axon signal in a lick-task session.

usage: python dd_fir.py [session ...]      (one pickle per session in OUTDIR/cache)

Question: does the axon signal respond differently to licks that are ANTICIPATORY (wait licks shortly before the end of a wait: ant_to = the wait then ends in a time-out, ant_sq = the
square follows), plain wait licks (wait_early), licks that COLLECT the reward (0.1-3 s after the square onset) and POST-reward licks (3-6 s, > 6 s after the onset)?
Every lick is assigned to one class by its own time in the trial structure (dd_common.Task.classify); the model has one finite-impulse-response kernel (delays -1 .. +4 s) per lick
class, kernels for the square onset (all trials, hit trials only = hit-minus-miss kernel), the pump (free = no lick before it / after the first lick) and the silent time-out, and as
covariates the movement the licks come with (jaw, face, whisking, pupil, rigid shift, ball where plausible; delays -1 .. +3 s).  Ridge (lambda by 5-fold blocked CV).
Outcome: population mean of the signal-like ROIs (30 s high-passed best dF/F, 5 Hz) and of ROI groups; control = noise-like ROIs.
Statistics: block bootstrap (60 s blocks, refit) for kernel summaries; null = outcome circularly shifted against the regressors (>= 60 s) for the class DIFFERENCES.
Summary of a kernel = mean weight over delays 0 .. 2 s = dF/F the axon signal moves per lick (mean over the 2 s after one lick).
"""
import os, sys, pickle, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dd_common as dd
import xs_common as xc
from xs_common import DT

LICK_CLASSES = ['wait_early', 'ant_to', 'ant_sq', 'win', 'post_a', 'post_b']
LAG_Y = np.arange(int(round(-1.0 / DT)), int(round(4.0 / DT)) + 1)             # delays of the event kernels (bins): -1 .. +4 s
LAG_B = np.arange(int(round(-1.0 / DT)), int(round(3.0 / DT)) + 1, 2)          # delays of the covariates: -1 .. +3 s, every 0.4 s
SUM_WIN = (0.0, 2.0)
N_BOOT, N_SHIFT, BLOCK_S = 400, 500, 60.0


def zs(v): return (v - np.nanmean(v)) / np.nanstd(v)


def lagged(x, lags):
    """columns x(t - lag): the column of delay l holds the series delayed by l bins (zero filled)"""
    n = len(x); out = np.zeros((n, len(lags)))
    for j, l in enumerate(lags):
        if l >= 0: out[l:, j] = x[:n - l] if l > 0 else x
        else: out[:n + l, j] = x[-l:]
    return out


def event_series(times, edges, weight=None):
    c, _ = np.histogram(times, bins=edges)
    return c.astype(float)


def groups_from_log(name, S):
    """ROI groups of the task script: 1506 as printed by the user's run (V1 section 8); other sessions: Ward tree on 1 - r of the signal-like ROIs, cut at 4 groups"""
    if name == 'AnimalRB24_261109_1506':
        g = [[1, 2, 3, 4, 5, 42, 44, 47, 50, 53, 54, 59, 64, 76, 82, 83, 85, 87, 88, 95, 100], [6, 8, 12, 15, 16, 18, 27, 45, 46, 57, 63, 68, 89, 92, 96, 102, 105],
             [9, 10, 20, 35, 37, 56, 70, 75, 81], [13, 17, 19, 51, 52]]
        return [np.array(x) - 1 for x in g]
    from scipy.cluster.hierarchy import linkage, fcluster
    idx = np.where(S.sig)[0]
    if len(idx) < 8: return [idx]
    D = np.nan_to_num(S.dff[idx].astype(float)); D = D - D.mean(1, keepdims=True)
    C = np.corrcoef(D); Z = linkage(1 - C[np.triu_indices(len(idx), 1)], 'ward')
    lab = fcluster(Z, 4, 'maxclust'); gs = [idx[lab == k] for k in np.unique(lab)]
    return sorted(gs, key=lambda a: -len(a))


def build_design(S, T, covs=('jaw', 'face', 'whisk', 'pupil', 'shift', 'loco')):
    """design matrix X (bins x columns) and the column map {name: slice}; every series on the 5 Hz grid"""
    edges, nb = S.edges, S.nb
    cols, cmap, k = [], {}, 0
    def add(name, M):
        nonlocal k; cols.append(M); cmap[name] = slice(k, k + M.shape[1]); k += M.shape[1]
    cls, _ = T.classify(T.lick)
    cls = np.array(['win' if c in ('win_other', 'win_early') else c for c in cls], dtype=object)
    for c in LICK_CLASSES:
        add('lick_' + c, lagged(event_series(T.lick[cls == c], edges), LAG_Y))
    add('sq', lagged(event_series(T.on, edges), LAG_Y))
    add('sq_hit', lagged(event_series(T.on[T.hit], edges), LAG_Y))
    pk = np.isfinite(T.rew_t)
    add('pump_free', lagged(event_series(T.rew_t[pk & (T.rew_kind == 2)], edges), LAG_Y))
    add('pump_lick', lagged(event_series(T.rew_t[pk & (T.rew_kind == 4)], edges), LAG_Y))
    if T.has_timeouts:
        add('timeout', lagged(event_series(T.seg['t1'][T.seg['to']], edges), LAG_Y))
    used = []
    for b in covs:
        if b == 'loco' and not S.loco_ok: continue
        if b not in S.beh: continue
        v = S.beh[b].astype(float); v = zs(v); v = v - xc.nsmooth(v, 30.0 / DT); v = np.where(np.isfinite(v), v, 0.0)
        add('cov_' + b, lagged(v, LAG_B)); used.append(b)
    X = np.concatenate(cols, 1)
    return X, cmap, used, cls


def ridge_solve(XtX, Xty, lam, n):
    return np.linalg.solve(XtX + lam * n * np.eye(XtX.shape[0]), Xty)


def cv_lambda(X, Y, folds=5, purge=25, lams=(1e-4, 1e-3, 1e-2, 1e-1, 1.0)):
    """blocked CV R2 of the population outcome (Y: bins x outcomes, NaN rows dropped per outcome) -> best lambda and its held-out R2 per lambda"""
    n = X.shape[0]; edges = np.linspace(0, n, folds + 1).astype(int); res = np.zeros((len(lams), Y.shape[1]))
    for fo in range(folds):
        te = np.arange(edges[fo], edges[fo + 1]); tr = np.ones(n, bool); tr[max(0, edges[fo] - purge):min(n, edges[fo + 1] + purge)] = False
        for oi in range(Y.shape[1]):
            y = Y[:, oi]; ok = np.isfinite(y); trm = tr & ok; tem = np.zeros(n, bool); tem[te] = True; tem &= ok
            Xtr = X[trm]; ytr = y[trm]; mu = ytr.mean()
            XtX = Xtr.T @ Xtr; Xty = Xtr.T @ (ytr - mu)
            for li, lam in enumerate(lams):
                w = ridge_solve(XtX, Xty, lam, trm.sum()); pred = X[tem] @ w + mu
                res[li, oi] += ((y[tem] - pred) ** 2).sum() / ((y[tem] - y[trm].mean()) ** 2).sum() / folds
    r2 = 1 - res
    return lams[int(np.argmax(r2.mean(1)))], r2


def kernel_summary(w, cmap, name, lags=LAG_Y, win=SUM_WIN):
    k = w[cmap[name]]; m = (lags * DT >= win[0] - 1e-9) & (lags * DT <= win[1] + 1e-9); return float(k[m].mean())


def run_session(name, verbose=True):
    t0 = time.time(); S, T = dd.load_session(name)
    X, cmap, covs, cls = build_design(S, T)
    nb = S.nb; n_roi = S.nroi
    hp = lambda x: x - xc.nsmooth(x, 30.0 / DT)
    groups = groups_from_log(name, S)
    outcomes = {'P': hp(S.population('sig'))}
    for gi, g in enumerate(groups): outcomes['g%d' % (gi + 1)] = hp(np.nanmean(S.dff[g], 0))
    if S.noise.sum() >= 3: outcomes['control'] = hp(S.population('noise'))
    Y = np.stack(list(outcomes.values()), 1).astype(float)
    Yroi = np.stack([hp(r.astype(float)) for r in S.dff[S.sig]], 1)     # bins x signal-like ROIs
    # drop bins without behaviour / axon (NaN in P) from the fit
    okrow = np.isfinite(Y[:, 0]) & np.isfinite(X).all(1)
    lam, r2cv = cv_lambda(X[okrow], Y[okrow], lams=(1e-4, 1e-3, 1e-2, 1e-1, 1.0))
    if verbose: print('%s: design %d x %d (%s), lambda %.0e, held-out R2 P %.3f' % (name, X.shape[0], X.shape[1], ','.join(covs), lam, r2cv[list((1e-4, 1e-3, 1e-2, 1e-1, 1.0)).index(lam), 0]), flush=True)
    rng = np.random.default_rng(7)
    # per-block sufficient statistics for the bootstrap (blocks of BLOCK_S s)
    bl = int(BLOCK_S / DT); nblk = int(np.ceil(nb / bl)); bid = np.minimum(np.arange(nb) // bl, nblk - 1)
    res = dict(name=name, lam=lam, covs=covs, cmap={k: (v.start, v.stop) for k, v in cmap.items()}, lags=LAG_Y * DT, r2cv=r2cv, outcome_names=list(outcomes.keys()), groups=[g.tolist() for g in groups],
               n_licks={c: int((cls == c).sum()) for c in LICK_CLASSES}, nb=nb)
    names = list(cmap.keys())
    def summaries(w):
        return {nm: kernel_summary(w, cmap, nm) for nm in names if not nm.startswith('cov_')}
    def kern(w):
        return {nm: w[cmap[nm]].copy() for nm in names if not nm.startswith('cov_')}
    res['out'] = {}
    for oi, on in enumerate(outcomes):
        y = Y[:, oi]; ok = okrow & np.isfinite(y); Xo = X[ok]; yo = y[ok]; mu = yo.mean()
        w = ridge_solve(Xo.T @ Xo, Xo.T @ (yo - mu), lam, ok.sum())
        base = dict(w=w, sum=summaries(w), kern=kern(w))
        # block bootstrap
        XtX_b = np.zeros((nblk, X.shape[1], X.shape[1])); Xty_b = np.zeros((nblk, X.shape[1])); n_b = np.zeros(nblk)
        for b in range(nblk):
            m = ok & (bid == b)
            if m.sum() == 0: continue
            Xb = X[m]; yb = y[m] - mu; XtX_b[b] = Xb.T @ Xb; Xty_b[b] = Xb.T @ yb; n_b[b] = m.sum()
        bs = []; bk = []
        for _ in range(N_BOOT):
            pick = rng.integers(0, nblk, nblk); wb = ridge_solve(XtX_b[pick].sum(0), Xty_b[pick].sum(0), lam, n_b[pick].sum()); bs.append(summaries(wb)); bk.append(kern(wb))
        base['boot'] = {nm: np.array([b[nm] for b in bs]) for nm in bs[0]}
        base['boot_kern'] = {nm: np.stack([b[nm] for b in bk]) for nm in bk[0]}
        # null: shift the outcome against the design (circular, >= 60 s)
        nulls = []
        XtX = Xo.T @ Xo
        for _ in range(N_SHIFT):
            sh = int(rng.integers(int(60 / DT), nb - int(60 / DT))); ysh = np.roll(y, sh); oks = okrow & np.isfinite(ysh); Xs = X[oks]; ys = ysh[oks]; mus = ys.mean()
            ws = ridge_solve(Xs.T @ Xs, Xs.T @ (ys - mus), lam, oks.sum()); nulls.append(summaries(ws))
        base['null'] = {nm: np.array([b[nm] for b in nulls]) for nm in nulls[0]}
        # held-out R2 with and without the class kernels (all lick classes removed / one class removed)
        res['out'][on] = base
        if verbose: print('   %-8s' % on + ' '.join('%s %+.4f' % (c[5:] if c.startswith('lick_') else c, base['sum'][c]) for c in names if not c.startswith('cov_')), flush=True)
    # per ROI (signal-like) kernel summaries, no bootstrap: ridge with the same lambda, all ROIs at once
    okr = okrow & np.isfinite(Yroi).all(1)
    Xo = X[okr]; Yo = Yroi[okr]; mu = Yo.mean(0); W = np.linalg.solve(Xo.T @ Xo + lam * okr.sum() * np.eye(X.shape[1]), Xo.T @ (Yo - mu))
    res['roi_sum'] = {nm: np.array([float(W[cmap[nm], r][(LAG_Y * DT >= 0) & (LAG_Y * DT <= 2)].mean()) for r in range(W.shape[1])]) for nm in names if not nm.startswith('cov_')}
    res['roi_idx'] = np.where(S.sig)[0]
    # same for the noise-like control ROIs (heterogeneity reference) and the covariate summary of the population
    idx_n = np.where(S.noise)[0]
    Yn = np.stack([hp(r.astype(float)) for r in S.dff[idx_n]], 1); okn = okrow & np.isfinite(Yn).all(1); Xn = X[okn]; Wn = np.linalg.solve(Xn.T @ Xn + lam * okn.sum() * np.eye(X.shape[1]), Xn.T @ (Yn[okn] - Yn[okn].mean(0)))
    res['noise_sum'] = {nm: np.array([float(Wn[cmap[nm], r][(LAG_Y * DT >= 0) & (LAG_Y * DT <= 2)].mean()) for r in range(Wn.shape[1])]) for nm in names if not nm.startswith('cov_')}
    res['sec'] = time.time() - t0
    return res


if __name__ == '__main__':
    os.makedirs(os.path.join(dd.OUTDIR, 'cache'), exist_ok=True)
    for nm in (sys.argv[1:] or dd.TASK_SESSIONS):
        f = os.path.join(dd.OUTDIR, 'cache', 'fir_' + nm + '.pkl')
        if os.path.exists(f) and not os.environ.get('REFRESH'): print('cached', nm); continue
        try:
            r = run_session(nm); pickle.dump(r, open(f, 'wb')); print('  %.0f s' % r['sec'])
        except Exception:
            import traceback; traceback.print_exc()
    print('FIR_DONE')
