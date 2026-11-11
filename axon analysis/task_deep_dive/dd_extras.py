"""Extra tests of the lick-task deep dive (one pickle per session: OUTDIR/cache/extras_<session>.pkl).

usage: python dd_extras.py [session ...]

 E1 heterogeneity   lick-class kernels (dd_fir design, same ridge) of EVERY signal-like ROI and every noise-like ROI; z of each ROI from 300 circular shifts of ITS OWN trace (>= 60 s) against the
                    design -> is the spread of lick responses across the signal-like ROIs larger than what the noise-like ROIs (common-mode, no axon) show?  Reported: fraction with |z| > 3
                    (expected 0.3 % for a normal null) and SD of z across ROIs for signal-like vs noise-like, per lick class.
 E2 subsessions     one-class lick kernel (all licks, with square / pump regressors and the movement covariates) fitted separately inside every subsession (30 trials, ~25 min): does the lick
                    coupling of the population / groups change with time in the session (satiety, engagement)?  block bootstrap (60 s) per subsession.
 E3 first-lick      response aligned to the first lick of the square window vs aligned to the square onset in hit trials: which alignment gives the larger / tighter average (is the response to
                    the stimulus or to the lick?): variance of the trial-averaged trace, and the trial-by-trial correlation of the response size with the first-lick latency.
"""
import os, sys, pickle, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dd_common as dd, dd_fir as df
import xs_common as xc
from xs_common import DT
from scipy.stats import spearmanr


def run_session(name, verbose=True):
    t0 = time.time(); S, T = dd.load_session(name); rng = np.random.default_rng(3)
    X, cmap, covs, cls = df.build_design(S, T); nb = S.nb
    hp = lambda x: x - xc.nsmooth(x, 30.0 / DT)
    groups = df.groups_from_log(name, S)
    outcomes = {'P': hp(S.population('sig'))}
    for gi, g in enumerate(groups): outcomes['g%d' % (gi + 1)] = hp(np.nanmean(S.dff[g], 0))
    if S.noise.sum() >= 3: outcomes['control'] = hp(S.population('noise'))
    res = dict(name=name, covs=covs, groups=[g.tolist() for g in groups])
    okrow = np.isfinite(outcomes['P']) & np.isfinite(X).all(1)
    lam = 1e-2
    names = [k for k in cmap if not k.startswith('cov_')]
    win = (df.LAG_Y * DT >= 0) & (df.LAG_Y * DT <= 2)
    # ---------------------------------------------------------------- E1 per-ROI z of the lick-class kernels
    Xo = X[okrow]; A = np.linalg.inv(Xo.T @ Xo + lam * okrow.sum() * np.eye(X.shape[1])) @ Xo.T            # (p x n) ridge operator
    e1 = {}
    for kind, mask in (('signal', S.sig), ('noise', S.noise)):
        rows = np.where(mask)[0]; Y = np.stack([hp(S.dff[r].astype(float)) for r in rows], 1)[okrow]
        Yc = Y - Y.mean(0); W = A @ Yc
        obs = {nm: W[cmap[nm]][win].mean(0) for nm in names}
        nul = {nm: np.zeros((300, len(rows))) for nm in names}
        for i in range(300):
            sh = int(rng.integers(int(60 / DT), len(Yc) - int(60 / DT))); Ws = A @ np.roll(Yc, sh, axis=0)
            for nm in names: nul[nm][i] = Ws[cmap[nm]][win].mean(0)
        z = {nm: (obs[nm] - nul[nm].mean(0)) / (nul[nm].std(0) + 1e-12) for nm in names}
        e1[kind] = dict(roi=rows, obs=obs, z=z)
    res['e1'] = e1
    # ---------------------------------------------------------------- E2 subsession stability (all licks one class)
    edges = S.edges; c_all, _ = np.histogram(T.lick, bins=edges)
    cols = [df.lagged(c_all.astype(float), df.LAG_Y), df.lagged(df.event_series(T.on, edges), df.LAG_Y)]
    pk = np.isfinite(T.rew_t); cols.append(df.lagged(df.event_series(T.rew_t[pk], edges), df.LAG_Y))
    for b in covs:
        v = S.beh[b].astype(float); v = df.zs(v); v = v - xc.nsmooth(v, 30.0 / DT); cols.append(df.lagged(np.where(np.isfinite(v), v, 0.0), df.LAG_B))
    X2 = np.concatenate(cols, 1); n_l = len(df.LAG_Y)
    sub_bins = [(int(np.floor((T.sub_t0[i] - S.edges[0]) / DT)), int(np.floor(((T.sub_t1[i] if i < len(T.sub_t1) else T.tr1[-1]) - S.edges[0]) / DT))) for i in range(len(T.sub_t0))] if hasattr(T, 'sub_t0') else []
    e2 = {}
    for oi, (on, y) in enumerate(outcomes.items()):
        e2[on] = []
        for (a, b) in sub_bins:
            m = np.zeros(nb, bool); m[max(a, 0):min(b, nb)] = True; m &= np.isfinite(y) & np.isfinite(X2).all(1)
            if m.sum() < 600: e2[on].append((np.nan, np.nan, np.nan)); continue
            Xs = X2[m]; ys = y[m]; mu = ys.mean(); w = np.linalg.solve(Xs.T @ Xs + lam * m.sum() * np.eye(X2.shape[1]), Xs.T @ (ys - mu)); est = float(w[:n_l][win].mean())
            # block bootstrap (60 s blocks) of the subsession
            idx = np.where(m)[0]; bl = int(60 / DT); bid = idx // bl; ub = np.unique(bid); bs = []
            XtX_b = {u: Xs[bid == u].T @ Xs[bid == u] for u in ub}; Xty_b = {u: Xs[bid == u].T @ (ys[bid == u] - mu) for u in ub}; n_b = {u: (bid == u).sum() for u in ub}
            for _ in range(200):
                pick = rng.choice(ub, len(ub)); wb = np.linalg.solve(sum(XtX_b[u] for u in pick) + lam * sum(n_b[u] for u in pick) * np.eye(X2.shape[1]), sum(Xty_b[u] for u in pick)); bs.append(float(wb[:n_l][win].mean()))
            e2[on].append((est, np.percentile(bs, 2.5), np.percentile(bs, 97.5)))
    res['e2'] = e2
    # ---------------------------------------------------------------- E3 stimulus- vs first-lick-aligned response of hit trials
    e0 = S.edges[0]; tau = np.arange(int(-3 / DT), int(8 / DT) + 1) * DT
    hit = T.hit & np.isfinite(T.lat1); fl = T.on + T.lat1
    e3 = {}
    for on, y in outcomes.items():
        def al(ev):
            ii = np.floor((ev - e0) / DT).astype(int)[:, None] + np.arange(len(tau))[None, :]; v = (ii >= 0) & (ii < nb); M = np.full(ii.shape, np.nan); M[v] = y[ii[v]]; return M
        Mo = al(T.on[hit]); Mf = al(fl[hit]); bo = np.nanmean(Mo[:, (tau >= -3) & (tau <= -0.5)], 1); bf = np.nanmean(Mf[:, (tau >= -3) & (tau <= -0.5)], 1)
        resp = np.nanmean(Mf[:, (tau >= 0) & (tau <= 2)], 1) - bf
        lat = T.lat1[hit]
        e3[on] = dict(tau=tau, onset=np.nanmean(Mo - bo[:, None], 0), lick=np.nanmean(Mf - bf[:, None], 0), n=int(hit.sum()),
                      peak_onset=float(np.nanmax(np.abs(np.nanmean(Mo - bo[:, None], 0)[tau >= 0]))), peak_lick=float(np.nanmax(np.abs(np.nanmean(Mf - bf[:, None], 0)[tau >= 0]))),
                      rho_resp_lat=float(spearmanr(resp, lat, nan_policy='omit')[0]))
    res['e3'] = e3
    res['sec'] = time.time() - t0
    if verbose:
        for kind in ('signal', 'noise'):
            z = e1[kind]['z']
            print('  %-6s |z|>3 fraction per class: ' % kind + ' '.join('%s %.2f' % (nm.replace('lick_', ''), (np.abs(z[nm]) > 3).mean()) for nm in names if nm.startswith('lick_')) + '   SD of z: ' + ' '.join('%.1f' % z[nm].std() for nm in names if nm.startswith('lick_')))
        for on in ('P', 'g1', 'g2', 'g3', 'g4', 'control'):
            if on in e2: print('  E2 %-7s per subsession: ' % on + ' | '.join('%+.4f [%+.4f,%+.4f]' % v for v in e2[on]))
    return res


if __name__ == '__main__':
    os.makedirs(os.path.join(dd.OUTDIR, 'cache'), exist_ok=True)
    for nm in (sys.argv[1:] or dd.TASK_SESSIONS):
        f = os.path.join(dd.OUTDIR, 'cache', 'extras_' + nm + '.pkl')
        if os.path.exists(f) and not os.environ.get('REFRESH'): print('cached', nm); continue
        try:
            print(nm, flush=True); r = run_session(nm); pickle.dump(r, open(f, 'wb')); print('  %.0f s' % r['sec'], flush=True)
        except Exception:
            import traceback; traceback.print_exc()
    print('EXTRAS_DONE')
