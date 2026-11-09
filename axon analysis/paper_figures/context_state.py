"""Analysis for Figure 11: axon activity in the square sessions with a spout / reward (lick task) and without (passive), with the movement / arousal confounds regressed out,
and engaged vs not engaged trials inside the task sessions.   usage: python context_state.py [session ...]  -> cache/ctx_<session>.pkl

Signal: P = mean best dF/F of the signal-like ROIs on the 5 Hz grid, high-passed > 30 s (as in the stimulus analysis). Behaviour: pupil, face motion, whisking, ball (only where plausible), z-scored and
high-passed, at lags -1..+3 s (0.4 s steps) = nuisance design B (ridge, penalty fixed at 0.01 * n on the standardised columns).
Epoch regressors (square sessions): pre = -3..0 s, cue = 0-3 s, post = 3-10 s after square onset; the rest of the trial cycle (long wait + late grey) is the reference; reported contrasts: cue - pre (= the onset response of the stimulus analysis) and post - pre.
  raw   : P ~ epochs                       (difference of the epoch mean from the rest of the cycle)
  mv    : P ~ B(face, whisking, ball) + epochs    ('movement removed')
  all   : P ~ B(pupil, face, whisking, ball) + epochs   (pupil also removed: arousal; the cue is a light stimulus, so this is the conservative version)
Null for every epoch coefficient / contrast: the epoch pattern (pre + cue + post together) circularly shifted by >= 60 s, behaviours and P unchanged (500 draws).
Engagement (task sessions, rows inside trial cycles only): H = 1 during the cycle of a hit trial, 0 during a miss trial; E = hit rate of the surrounding 10 trials (z). Model: P ~ B + cue + post + (H or E);
null = hit vector circularly shifted by >= 5 trials (run structure kept).
Residual traces: res = P - B w (behaviour part of the 'all' model fitted with the epoch terms), so that responses of the epochs stay in res.
"""
import os, sys, pickle, numpy as np
from fig_common import *
import xs_common as xc
from xs_common import DT
from stim_response import hp, align
from task_engagement import load_task

SESS = ['AnimalRB24_261029_1752', 'AnimalRB24_261030_1744', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261106_1317', 'AnimalRB23_261029_1704',      # task (spout + reward)
        'AnimalRB24_261027_1456', 'AnimalRB23_261027_1423', 'AnimalRB24_261104_1651', 'AnimalRB24_261106_1529']                                                      # passive (no spout)
PRE, CUE, POST = (-3.0, 0.0), (0.0, 3.0), (3.0, 10.0); K = np.array([[-1.0, 1, 0], [-1.0, 0, 1]]); WIN = (-5.0, 20.0); BASE = (-3.0, -0.5)
LAGS = (np.arange(-1.0, 3.01, 0.4) / DT).round().astype(int); N_NULL = 500; LAM = 0.01
rng = np.random.default_rng(21)


def zs(v): return (v - np.nanmean(v)) / np.nanstd(v)


def design(behs, names):
    cols = []
    for b in names:
        for l in LAGS: cols.append(np.roll(behs[b], -l))
    return np.stack(cols, 1)


def epoch_vec(e0, nb, ons, w):
    v = np.zeros(nb)
    for o in ons:
        i1 = int(np.floor((o + w[0] - e0) / DT)); i2 = int(np.floor((o + w[1] - e0) / DT)); v[max(i1, 0):max(min(i2, nb), 0)] = 1.0
    return v


class Fit:
    """ridge on the nuisance design N (rows already selected, centred), unpenalised tested regressors T"""

    def __init__(self, N, y, pen):
        self.N, self.y = N, y
        self.G = np.linalg.inv(N.T @ N + np.diag(pen)) if N.shape[1] else np.zeros((0, 0)); self.Ny = N.T @ y

    def effect(self, T):
        T = T - T.mean(0)
        if self.N.shape[1] == 0: return np.linalg.solve(T.T @ T, T.T @ self.y)
        C = self.N.T @ T; S = T.T @ T - C.T @ self.G @ C; r = T.T @ self.y - C.T @ self.G @ self.Ny
        return np.linalg.solve(S + 1e-9 * np.eye(S.shape[0]), r)

    def weights(self, T):
        """nuisance weights given the tested regressors"""
        T = T - T.mean(0)
        if self.N.shape[1] == 0: return np.zeros(0)
        b = self.effect(T); return self.G @ (self.Ny - self.N.T @ T @ b)


def analyse(name):
    S = xc.Session(name, with_rejected=False, verbose=False); X = S.X; e0 = S.edges[0]; nb = S.nb; out = dict(name=name, n_sig=int(S.sig.sum()), n_roi=S.nroi, dur_min=S.dur / 60, loco_ok=S.loco_ok, task=name in TASK_SESSIONS)
    P = hp(S.population('sig')); behs = {b: hp(zs(S.beh[b])) for b in ('pupil', 'face', 'whisk', 'loco') if b in S.beh and (b != 'loco' or S.loco_ok) and np.isfinite(S.beh[b]).mean() > 0.3}
    out['behaviours'] = list(behs)
    on_t = np.atleast_1d(np.asarray(X['on_t'], float)); on_blk = np.atleast_1d(np.asarray(X['on_blk'], int)); binter = np.atleast_1d(np.asarray(X['blk_interrupted'], float)).astype(bool)
    keep = np.isfinite(on_t) & ~binter[on_blk - 1]; ons = on_t[keep]; T = None
    if out['task']:
        T = load_task(name); ons_t = np.asarray(T['on_t'], float); out['clock_diff_s'] = float(np.median(np.abs(ons_t - ons[:len(ons_t)]))) if len(ons) >= len(ons_t) else np.nan; ons = ons_t
    out['n_trials'] = int(len(ons))
    pre = epoch_vec(e0, nb, ons, PRE); cue = epoch_vec(e0, nb, ons, CUE); post = epoch_vec(e0, nb, ons, POST)
    fam = {'mv': [b for b in ('face', 'whisk', 'loco') if b in behs], 'all': [b for b in ('pupil', 'face', 'whisk', 'loco') if b in behs]}
    ok_b = np.all([np.isfinite(behs[b]) for b in fam['all']], 0); okrow = np.isfinite(P) & ok_b
    # rows: everything with finite P and behaviours; zero out the lag wrap-around edge
    edge = np.zeros(nb, bool); edge[:LAGS.max() + 1] = True; edge[-(abs(LAGS.min()) + 1):] = True; okrow &= ~edge
    okrow &= np.isfinite(design(behs, fam['all'])).all(1)                    # lagged copies must be finite as well
    out['frac_rows'] = float(okrow.mean())
    out['epochs'] = {}; E = np.stack([pre, cue, post], 1)
    models = {'raw': [], 'mv': fam['mv'], 'all': fam['all']}; weights = {}
    for tag, names in models.items():
        if tag != 'raw' and not names: continue
        Bm = design(behs, names)[okrow] if names else np.zeros((okrow.sum(), 0)); mu = Bm.mean(0) if Bm.size else 0; sd = Bm.std(0) + 1e-9 if Bm.size else 1
        Zb = (Bm - mu) / sd if Bm.size else Bm; y = P[okrow] - P[okrow].mean(); n = len(y)
        F = Fit(Zb, y, np.full(Zb.shape[1], LAM * n)); Eok = E[okrow]; beta = F.effect(Eok)
        nul = np.empty((N_NULL, 3)); cs = rng.integers(int(60 / DT), nb - int(60 / DT), N_NULL)
        for i, d in enumerate(cs): nul[i] = F.effect(np.roll(E, d, axis=0)[okrow])
        con = K @ beta; cn = nul @ K.T
        out['epochs'][tag] = dict(beta=beta, contrast=con, c_null_sd=cn.std(0), c_null_mean=cn.mean(0), c_z=(con - cn.mean(0)) / cn.std(0), c_p=(1 + (np.abs(cn) >= np.abs(con)).sum(0)) / (1 + N_NULL))
        if tag != 'raw':
            w = F.weights(Eok); hatB = np.full(nb, np.nan); hatB[okrow] = Zb @ w; weights[tag] = (hatB, mu, sd, w)
            Ec = Eok - Eok.mean(0); fit_all = Zb @ w + Ec @ beta
            r2_total = 1 - ((y - fit_all) ** 2).sum() / (y ** 2).sum(); r2_B_only = 1 - ((y - Zb @ w) ** 2).sum() / (y ** 2).sum(); r2_E_only = 1 - ((y - Ec @ beta) ** 2).sum() / (y ** 2).sum()
            out['epochs'][tag].update(r2_total=float(r2_total), r2_B=float(r2_B_only), r2_E=float(r2_E_only))
    hatB = weights['all'][0]; res = P - hatB; hatB_mv = weights['mv'][0] if 'mv' in weights else None
    out['sd_P'] = float(np.nanstd(P)); out['sd_res'] = float(np.nanstd(res))
    # aligned averages
    def al(x, sel=None):
        o = ons if sel is None else ons[sel]
        rel, C = align(x, e0, o, WIN); return rel, C - np.nanmean(C[:, (rel >= BASE[0]) & (rel < BASE[1])], 1, keepdims=True)
    rel, Cp = al(P); _, Cr = al(res); out['rel'] = rel; out['P_tr'] = Cp.astype(np.float32); out['res_tr'] = Cr.astype(np.float32)
    out['res_mv_tr'] = al(P - hatB_mv)[1].astype(np.float32) if hatB_mv is not None else None
    out['beh_tr'] = {b: al(v)[1].astype(np.float32) for b, v in behs.items()}
    # engagement inside the task sessions
    if out['task']:
        hit = np.asarray(T['hit'], float); t0 = np.asarray(T['tr_t0'], float); t1 = np.asarray(T['tr_t1'], float); out['hit'] = hit; out['hit_rate'] = float(np.nanmean(hit))
        ntr = len(hit); idx = np.full(nb, -1);
        for k in range(ntr):
            i1 = int(np.floor((t0[k] - e0) / DT)); i2 = int(np.floor((t1[k] - e0) / DT)); idx[max(i1, 0):max(min(i2, nb), 0)] = k
        inside = (idx >= 0) & okrow; Ein = E[inside]; ti = idx[inside]; y = P[inside] - P[inside].mean()
        def trial_series(h):  # rolling hit rate of 10 trials, z-scored
            k = np.ones(10) / 10; hh = np.where(np.isfinite(h), h, np.nanmean(h)); r = np.convolve(np.pad(hh, 5, mode='edge'), k, mode='same')[5:-5]; return (r - r.mean()) / (r.std() + 1e-9)
        eng = {}
        for tag, names in (('raw', []), ('mv', fam['mv']), ('all', fam['all'])):
            Bm = design(behs, names)[inside] if names else np.zeros((inside.sum(), 0)); Zb = (Bm - Bm.mean(0)) / (Bm.std(0) + 1e-9) if Bm.size else Bm; n = len(y)
            N = np.concatenate([Zb, Ein - Ein.mean(0)], 1) if Zb.size else Ein - Ein.mean(0); pen = np.concatenate([np.full(Zb.shape[1], LAM * n), np.zeros(3)]) if Zb.size else np.zeros(3)
            F = Fit(N, y, pen); r = {}
            for lab, vec in (('H', hit), ('E', trial_series(hit))):
                good = np.isfinite(vec); v = vec.copy(); v[~good] = np.nanmean(vec[good])
                if lab == 'H' and np.nanstd(v) == 0: continue
                beta = float(F.effect((v[ti])[:, None])[0]); nul = []
                for s in rng.integers(5, max(ntr - 5, 6), 300):
                    nul.append(F.effect((np.roll(v, s)[ti])[:, None])[0])
                nul = np.array(nul); r[lab] = dict(beta=beta, null_sd=float(nul.std()), null_mean=float(nul.mean()), z=float((beta - nul.mean()) / nul.std()), p=float((1 + (np.abs(nul) >= abs(beta)).sum()) / (1 + len(nul))))
            eng[tag] = r
        out['engage'] = eng; out['n_hit'] = int(np.nansum(hit)); out['n_miss'] = int(np.sum(hit == 0))
        # hit vs miss: cycle-averaged residual and raw (descriptive)
        cyc = np.array([np.nanmean(res[(idx == k)]) if (idx == k).any() else np.nan for k in range(ntr)]); cycP = np.array([np.nanmean(P[(idx == k)]) if (idx == k).any() else np.nan for k in range(ntr)])
        out['cyc_res'] = cyc; out['cyc_P'] = cycP; out['rate_E'] = trial_series(hit)
        out['licks'] = np.asarray(T['lick_t'], float); out['st_end'] = np.asarray(T['st_end'], float)
    # excerpt for the raw-trace panel
    t = S.t; out['t'] = t.astype(np.float32); out['ex'] = dict(P=P.astype(np.float32), fitB=hatB.astype(np.float32), res=res.astype(np.float32), ons=ons, behs={b: v.astype(np.float32) for b, v in behs.items()})
    return out


if __name__ == '__main__':
    names = sys.argv[1:] or SESS
    for nm in names:
        f = os.path.join(CACHE, 'ctx_%s.pkl' % nm)
        if os.path.exists(f) and not os.environ.get('REFRESH'): print('cached', nm); continue
        import time; t0 = time.time(); r = analyse(nm); pickle.dump(r, open(f, 'wb'))
        e = r['epochs']; print('%-26s %s ntr %3d | cue-pre raw %+.3f mv %+.3f all %+.3f (z %.1f / %.1f / %.1f) | post-pre raw %+.3f all %+.3f (z %.1f / %.1f) | R2(B) %.2f | %.0f s' % (nm, 'task' if r['task'] else 'pass', r['n_trials'], e['raw']['contrast'][0], e['mv']['contrast'][0], e['all']['contrast'][0], e['raw']['c_z'][0], e['mv']['c_z'][0], e['all']['c_z'][0], e['raw']['contrast'][1], e['all']['contrast'][1], e['raw']['c_z'][1], e['all']['c_z'][1], e['all']['r2_B'], time.time() - t0), flush=True)
        if r['task']: print('     engagement:', {k: {kk: (round(vv['beta'], 4), round(vv['z'], 1)) for kk, vv in v.items()} for k, v in r['engage'].items()}, 'hit rate %.2f' % r['hit_rate'], 'clock diff %.2f' % r.get('clock_diff_s', np.nan), flush=True)
