"""Hit vs miss trials of the lick task: what is in the axon activity beyond the licking itself?

usage: python dd_hitmiss.py [session ...]        (one pickle per session in OUTDIR/cache/hitmiss_<session>.pkl)

hit = a lick 0.1-3 s after the square onset; miss = none.  Hits contain licking, misses do not, so the plain hit - miss difference is mostly movement.  Tests (all window means of the
30 s high-passed best dF/F, population P of the signal-like ROIs, the ROI groups, the noise-like control):
  1 time course        trial-aligned P(t), hit and miss separately, hit - miss with a trial bootstrap; windows pre (-3..0), early (0..0.8, hits with first lick > 0.9 s only: no licking in the
                       window in either class), mid (0.8..3), late (3..8 s)
  2 pre-trial state    can the axon level / movement in the 5 s before the square (-6..-1) predict the outcome?  cross-validated logistic regression (contiguous folds), movement + lick
                       history + level (M0) vs M0 + axon (M1); delta AUC against a null that permutes the axon rows inside subsessions
  3 per-ROI scan       hit - miss of every signal-like ROI and noise-like ROI in the four windows, raw and after regressing the window movement (face, whisking, jaw, lick rate, ball, pupil
                       state); null = outcome labels permuted inside subsessions; the noise-like ROIs calibrate how many ROIs 'pass' by chance / common-mode
"""
import os, sys, pickle, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dd_common as dd, dd_matched as dm
import xs_common as xc
from xs_common import DT
from dd_fir import groups_from_log
from matched_behaviour import _cum, _wstat
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

WINS = {'pre': (-3.0, 0.0), 'early': (0.0, 0.8), 'mid': (0.8, 3.0), 'late': (3.0, 8.0)}
TC = (-6.0, 10.0)
N_PERM = 2000; N_BOOT = 1000


def wm(series_cum, idx, win, mincov=0.8):
    a, b = int(round(win[0] / DT)), int(round(win[1] / DT)); return _wstat(series_cum, idx + a, idx + b, mincov)


def trace(x, idx, lo, hi):
    tau = np.arange(int(round(lo / DT)), int(round(hi / DT)) + 1); ii = idx[:, None] + tau[None, :]; o = np.full(ii.shape, np.nan); v = (ii >= 0) & (ii < len(x)); o[v] = x[ii[v]]; return tau * DT, o


def run_session(name, verbose=True):
    t0 = time.time(); R, T = dm.make_rec(name); S = R.S; rng = np.random.default_rng(5)
    hp = lambda x: x - xc.nsmooth(x, R.par['hp_s'] / DT)
    groups = groups_from_log(name, S)
    outs = {'P': R.P}
    for gi, g in enumerate(groups): outs['g%d' % (gi + 1)] = hp(np.nanmean(S.dff[g], 0))
    if S.noise.sum() >= 3: outs['control'] = hp(S.population('noise'))
    idx = np.floor((T.on - R.e0) / DT).astype(int); ok_idx = (idx >= 80) & (idx < R.nb - 80)
    hit = T.hit.copy(); n = T.n; sub = T.sub
    slow = hit & (T.lat1 > 0.9)
    res = dict(name=name, n=n, n_hit=int(hit.sum()), groups=[g.tolist() for g in groups], outcome_names=list(outs.keys()))
    # ---------------------------------------------------------------- 1 time course and window means
    res['tc'] = {}; res['win'] = {}
    for on, x in outs.items():
        tau, M = trace(x, idx, *TC); base = np.nanmean(M[:, (tau >= -3) & (tau <= -0.5)], 1)
        cum = _cum(x); W = {k: wm(cum, idx, w) for k, w in WINS.items()}
        def boot_diff(v, A, B):
            a = v[A & np.isfinite(v)]; b = v[B & np.isfinite(v)]
            if len(a) < 3 or len(b) < 3: return np.nan, (np.nan, np.nan)
            bs = [rng.choice(a, len(a)).mean() - rng.choice(b, len(b)).mean() for _ in range(N_BOOT)]; return a.mean() - b.mean(), tuple(np.percentile(bs, [2.5, 97.5]))
        res['tc'][on] = dict(tau=tau, hit=np.nanmean(M[hit], 0), miss=np.nanmean(M[~hit], 0), hit_bs=np.nanmean((M - base[:, None])[hit], 0), miss_bs=np.nanmean((M - base[:, None])[~hit], 0),
                             diff_bs=np.nanmean((M - base[:, None])[hit], 0) - np.nanmean((M - base[:, None])[~hit], 0))
        res['win'][on] = {k: dict(zip(('diff', 'ci'), boot_diff(W[k], hit if k != 'early' else slow, ~hit))) for k in WINS}
        res['win'][on]['W'] = W
    # ---------------------------------------------------------------- 2 pre-trial state -> outcome
    pre = (-6.0, -1.0); B = {}
    for b in ('face', 'whisk', 'jaw', 'lick', 'pupil') + (('loco',) if 'loco' in R.beh else ()):
        if b in R.beh: B[b] = wm(R.cum[b], idx, pre)
    prev = np.r_[np.nan, hit[:-1].astype(float)]; prev5 = np.array([np.nanmean(hit[max(0, t - 5):t]) if t > 0 else np.nan for t in range(n)])
    hist = np.stack([prev, prev5, (T.lvl == 3).astype(float), T.n_to.astype(float), np.arange(n) / n], 1)
    Xb = np.stack(list(B.values()), 1); X0 = np.concatenate([Xb, hist], 1)
    axon = np.stack([wm(_cum(outs[k]), idx, pre) for k in outs if k != 'control'], 1)
    ok = ok_idx & np.isfinite(X0).all(1) & np.isfinite(axon).all(1)
    def cv_auc(X, y, sub_ok, folds=10, C=0.3):
        Xs = (X - X.mean(0)) / (X.std(0) + 1e-9); nn = len(y); e = np.linspace(0, nn, folds + 1).astype(int); p = np.full(nn, np.nan)
        for f in range(folds):
            te = np.arange(e[f], e[f + 1]); tr = np.setdiff1d(np.arange(nn), te)
            if len(np.unique(y[tr])) < 2: continue
            m = LogisticRegression(C=C, max_iter=2000).fit(Xs[tr], y[tr]); p[te] = m.predict_proba(Xs[te])[:, 1]
        k = np.isfinite(p); return roc_auc_score(y[k], p[k]) if len(np.unique(y[k])) > 1 else np.nan
    y = hit[ok].astype(int); X0o = X0[ok]; axo = axon[ok]; sub_o = sub[ok]
    auc0 = cv_auc(X0o, y, None); auc1 = cv_auc(np.concatenate([X0o, axo], 1), y, None); aucA = cv_auc(axo, y, None); aucB = cv_auc(Xb[ok], y, None)
    null = []
    for _ in range(300):
        ax_p = axo.copy()
        for s_ in np.unique(sub_o):
            ii = np.where(sub_o == s_)[0]; ax_p[ii] = axo[rng.permutation(ii)]
        null.append(cv_auc(np.concatenate([X0o, ax_p], 1), y, None) - auc0)
    res['predict'] = dict(n=int(ok.sum()), auc_behaviour=float(aucB), auc_m0=float(auc0), auc_axon_only=float(aucA), auc_m1=float(auc1), d_auc=float(auc1 - auc0), null_d=np.array(null),
                          p=float((1 + (np.array(null) >= auc1 - auc0).sum()) / (1 + len(null))), features=list(B.keys()) + ['prev_hit', 'hit_rate_prev5', 'level3', 'n_timeouts', 'trial_pos'])
    # first-lick latency of hits vs axon pre-state (partial on movement)
    lat = np.log(T.lat1); h = hit & ok_idx & np.isfinite(lat) & np.isfinite(X0).all(1) & np.isfinite(axon).all(1)
    if h.sum() >= 15:
        from numpy.linalg import lstsq
        Z = np.concatenate([X0[h], np.ones((h.sum(), 1))], 1); rl = lat[h] - Z @ lstsq(Z, lat[h], rcond=None)[0]
        res['latency'] = dict(n=int(h.sum()), rho_P_raw=float(spearmanr(axon[h, 0], lat[h])[0]), rho_P_resid=float(spearmanr(axon[h, 0] - Z @ lstsq(Z, axon[h, 0], rcond=None)[0], rl)[0]))
    # ---------------------------------------------------------------- 3 per-ROI scan
    rows = []
    for kind, mask in (('signal', S.sig), ('noise', S.noise)):
        for r in np.where(mask)[0]:
            x = hp(S.dff[r].astype(float)); cum = _cum(x); rows.append((kind, r, {k: wm(cum, idx, w) for k, w in WINS.items()}))
    Bwin = {}   # movement features per window for the adjustment
    for k, w in WINS.items():
        cols = [wm(R.cum[b], idx, w) for b in ('face', 'whisk', 'jaw', 'lick') + (('loco',) if 'loco' in R.beh else ()) if b in R.beh]
        cols.append(wm(R.cum['pupil'], idx, (w[0] - 5.0, w[0] - 1.0)) if 'pupil' in R.cum else np.zeros(n)); Bwin[k] = np.stack(cols, 1)
    scan = {}
    for k in WINS:
        lab = (slow if k == 'early' else hit); labm = (~hit) | lab                     # early: slow hits vs misses only
        valid = ok_idx & labm & np.isfinite(Bwin[k]).all(1)
        W = np.stack([rw[2][k] for rw in rows]); valid &= np.isfinite(W).all(0)
        Wv = W[:, valid]; L = lab[valid]; sv = sub[valid]; Bv = Bwin[k][valid]
        Zb = np.concatenate([(Bv - Bv.mean(0)) / (Bv.std(0) + 1e-9), np.ones((valid.sum(), 1))], 1)
        Wadj = Wv - (Zb @ np.linalg.lstsq(Zb, Wv.T, rcond=None)[0]).T
        out = {}
        for tag, M in (('raw', Wv), ('adj', Wadj)):
            diff = M[:, L].mean(1) - M[:, ~L].mean(1)
            perm = np.zeros((N_PERM, M.shape[0]))
            for pi in range(N_PERM):
                lp = L.copy()
                for s_ in np.unique(sv): ii = np.where(sv == s_)[0]; lp[ii] = L[rng.permutation(ii)]
                perm[pi] = M[:, lp].mean(1) - M[:, ~lp].mean(1)
            z = (diff - perm.mean(0)) / (perm.std(0) + 1e-12); p = (1 + (np.abs(perm - perm.mean(0)) >= np.abs(diff - perm.mean(0))).sum(0)) / (1 + N_PERM)
            out[tag] = dict(diff=diff, z=z, p=p)
        scan[k] = dict(n_hit=int(L.sum()), n_miss=int((~L).sum()), **out)
    res['scan'] = dict(kind=np.array([r[0] for r in rows]), roi=np.array([r[1] for r in rows]), win=scan)
    res['sec'] = time.time() - t0
    if verbose:
        print('%s: %d trials (%d hit), P windows (hit - miss, 95%% CI):' % (name, n, hit.sum()), ' '.join('%s %+.3f [%+.3f,%+.3f]' % (k, res['win']['P'][k]['diff'], *res['win']['P'][k]['ci']) for k in WINS))
        pr = res['predict']; print('   predict hit from pre state: AUC behaviour %.2f, M0 %.2f, axon only %.2f, M1 %.2f, dAUC %+.3f (p %.3f, n %d)' % (pr['auc_behaviour'], pr['auc_m0'], pr['auc_axon_only'], pr['auc_m1'], pr['d_auc'], pr['p'], pr['n']))
        for k in WINS:
            for tag in ('raw', 'adj'):
                kk = res['scan']['kind']; sg = kk == 'signal'; nz = kk == 'noise'; z = scan[k][tag]['z']; p = scan[k][tag]['p']
                print('   ROI scan %-5s %-3s signal-like p<0.05 %2d/%2d (%.0f%%)  noise-like %2d/%2d (%.0f%%)' % (k, tag, (p[sg] < 0.05).sum(), sg.sum(), 100 * (p[sg] < 0.05).mean(), (p[nz] < 0.05).sum(), nz.sum(), 100 * (p[nz] < 0.05).mean()))
    return res


if __name__ == '__main__':
    os.makedirs(os.path.join(dd.OUTDIR, 'cache'), exist_ok=True)
    for nm in (sys.argv[1:] or dd.TASK_SESSIONS):
        f = os.path.join(dd.OUTDIR, 'cache', 'hitmiss_' + nm + '.pkl')
        if os.path.exists(f) and not os.environ.get('REFRESH'): print('cached', nm); continue
        try:
            r = run_session(nm); pickle.dump(r, open(f, 'wb')); print('  %.0f s' % r['sec'], flush=True)
        except Exception:
            import traceback; traceback.print_exc()
    print('HITMISS_DONE')
