"""Cross-session analysis 1: the shared 'behavioural state' signal of the LC axons, standardised over sessions.

usage: python state_signal.py [session ...]      (default: all sessions of SESSIONS; results -> state_signal_results.pkl next to the figures)

Per session (everything on the 0.2 s grid of xs_common, dF/F = BestDFF recipe, signals band-passed 1 s .. 30 s):
  coherence   PC1 variance fraction, participation ratio, median pairwise r, same-sign fraction of the signal-like ROIs, against random sets of the same
              number of REJECTED suite2p ROIs (the control population)
  lead / lag  cross-correlation of the population signal (mean of the signal-like ROIs) with pupil, pupil dilation rate, face motion, whisking, eyelid, locomotion
              (only where the ball is plausible), brain shift; peak r / lag, null = the same statistic at lags >= 120 s away (circular), lag CI by block bootstrap
  ROI signs   r at the population lag for every ROI: signal-like / noise-like / rejected ROIs; significant positive vs negative (circular null p < 0.01)
  ridge       held-out R2 of the population signal from the behaviour families (lags -1..+3 s), unique contribution per family, same model on random rejected sets
  states      mean population signal in quiet / active (face motion terciles), running vs still (ball) with block-bootstrap CI
"""
import sys, os, pickle, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xs_common as xc
from xs_common import DT

SESSIONS = ['AnimalRB24_261021_1630', 'AnimalRB24_261023_1527', 'AnimalRB24_261023_1538', 'AnimalRB24_261027_1456', 'AnimalRB24_261027_1534', 'AnimalRB24_261029_1752',
            'AnimalRB24_261030_1744', 'AnimalRB24_261104_1640', 'AnimalRB24_261104_1651', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261106_1317', 'AnimalRB24_261106_1529',
            'AnimalRB23_261023_1642', 'AnimalRB23_261027_1306', 'AnimalRB23_261027_1423', 'AnimalRB23_261029_1704', 'AnimalRB22_261023_1427', 'AnimalRB22_261027_1138']
OUT = os.environ.get('STATE_OUT', os.path.join(r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session', 'state_signal_results.pkl'))
BEHS = ['pupil', 'dpupil', 'face', 'whisk', 'eyelid', 'loco', 'shift']
FAMS = ['pupil', 'face', 'whisk', 'eyelid', 'loco']            # ridge families
MAXLAG = int(10 / DT); SEP = int(120 / DT); BLOCK = int(60 / DT)
rng = np.random.default_rng(0)


def zs(v): return (v - np.nanmean(v)) / np.nanstd(v)


def coherence(D):
    """D: ROIs x bins (band-passed). PC1 fraction, participation ratio, median pairwise r, same-sign fraction"""
    ok = np.isfinite(D).all(0); D = D[:, ok]
    if D.shape[0] < 3 or D.shape[1] < 200: return None
    Z = (D - D.mean(1, keepdims=True)) / D.std(1, keepdims=True); C = np.corrcoef(Z); w, V = np.linalg.eigh(C); w = w[::-1]; v1 = V[:, ::-1][:, 0]
    iu = np.triu_indices(len(C), 1)
    return dict(pc1=float(w[0] / w.sum()), pr=float(w.sum() ** 2 / (w ** 2).sum()), r_med=float(np.median(C[iu])), same_sign=float(max((v1 > 0).mean(), (v1 < 0).mean())), n=int(len(C)))


def peak_and_ci(fx, fy, n_boot=200):
    lags, r, full = xc.xcorr_fft(fx, fy, MAXLAG)
    (lg, rr), nul = xc.peak_stats(full, MAXLAG, SEP)
    nb = len(fx); nblk = nb // BLOCK; lg_b = []
    if nblk >= 6:
        for _ in range(n_boot):
            ix = np.concatenate([np.arange(b * BLOCK, (b + 1) * BLOCK) for b in rng.integers(0, nblk, nblk)])
            _, rb, _ = xc.xcorr_fft(fx[ix], fy[ix], MAXLAG); lg_b.append(lags[np.argmax(np.abs(rb))])
    ci = (np.percentile(lg_b, 2.5) * DT, np.percentile(lg_b, 97.5) * DT) if lg_b else (np.nan, np.nan)
    p = float((1 + (nul >= abs(rr)).sum()) / (1 + len(nul)))
    return dict(lags_s=lags * DT, r=r, peak_r=float(rr), peak_lag_s=float(lg * DT), null95=float(np.percentile(nul, 95)), p=p, lag_ci=ci)


def roi_signs(F_roi, fy, lag_bins):
    """r at a fixed lag for every ROI + circular-null p; returns r, p arrays"""
    r = np.zeros(len(F_roi)); p = np.ones(len(F_roi))
    for i, row in enumerate(F_roi):
        _, _, full = xc.xcorr_fft(row, fy, 1); n = len(full)
        ro = full[lag_bins % n]; k = np.arange(SEP, n - SEP); nul = full[k]
        r[i] = ro; p[i] = (1 + (np.abs(nul) >= abs(ro)).sum()) / (1 + len(nul))
    return r, p


def design(beh_list, lags_bins):
    cols = []
    for b in beh_list:
        for l in lags_bins:
            cols.append(np.roll(b, -l))        # column = b[t + l]
    return np.stack(cols, 1)


def ridge_heldout(Xm, y, n_fold=10, purge=int(5 / DT), lams=(1, 10, 100, 1e3, 1e4, 1e5)):
    """contiguous-fold held-out R2 with inner-CV lambda; rows with NaN are dropped"""
    ok = np.isfinite(Xm).all(1) & np.isfinite(y); idx = np.where(ok)[0]
    if len(idx) < 400: return np.nan
    X = Xm[idx]; Y = y[idx]; n = len(idx); edges = np.linspace(0, n, n_fold + 1).astype(int); pred = np.full(n, np.nan)
    def fit(Xt, Yt, lam):
        mu, sd = Xt.mean(0), Xt.std(0) + 1e-9; ym = Yt.mean(); Z = (Xt - mu) / sd
        A = Z.T @ Z + lam * len(Z) * np.eye(Z.shape[1]) / 100.0; w = np.linalg.solve(A, Z.T @ (Yt - ym)); return mu, sd, ym, w
    def predict(m, Xs): mu, sd, ym, w = m; return ((Xs - mu) / sd) @ w + ym
    for f in range(n_fold):
        te = np.arange(edges[f], edges[f + 1]); tr = np.setdiff1d(np.arange(n), np.arange(max(0, edges[f] - purge), min(n, edges[f + 1] + purge)))
        e2 = np.linspace(0, len(tr), 5).astype(int); best, bl = np.inf, lams[0]
        for lam in lams:
            sse = 0
            for g in range(4):
                v = np.arange(e2[g], e2[g + 1]); t = np.setdiff1d(np.arange(len(tr)), np.arange(max(0, e2[g] - purge), min(len(tr), e2[g + 1] + purge)))
                m = fit(X[tr[t]], Y[tr[t]], lam); sse += ((predict(m, X[tr[v]]) - Y[tr[v]]) ** 2).sum()
            if sse < best: best, bl = sse, lam
        m = fit(X[tr], Y[tr], bl); pred[te] = predict(m, X[te])
    return float(1 - ((Y - pred) ** 2).sum() / ((Y - Y.mean()) ** 2).sum())


def ridge_partition(y, behs, lags_bins):
    """full / alone / unique R2 of the families in behs (dict name -> series), target y"""
    names = list(behs); cols = {n: design([behs[n]], lags_bins) for n in names}
    full = ridge_heldout(np.concatenate([cols[n] for n in names], 1), y)
    alone = {n: ridge_heldout(cols[n], y) for n in names}
    uniq = {n: full - ridge_heldout(np.concatenate([cols[m] for m in names if m != n], 1), y) if len(names) > 1 else np.nan for n in names}
    return dict(full=full, alone=alone, unique=uniq)


def block_boot_diff(x, mask_a, mask_b, n_boot=300):
    ok = np.isfinite(x); nblk = len(x) // BLOCK; est = np.nanmean(x[mask_a & ok]) - np.nanmean(x[mask_b & ok]); bs = []
    for _ in range(n_boot):
        ix = np.concatenate([np.arange(b * BLOCK, (b + 1) * BLOCK) for b in rng.integers(0, nblk, nblk)])
        a, b = mask_a[ix] & ok[ix], mask_b[ix] & ok[ix]
        if a.sum() > 20 and b.sum() > 20: bs.append(x[ix][a].mean() - x[ix][b].mean())
    return float(est), (float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))) if bs else (np.nan, np.nan)


def analyse(S):
    out = dict(name=S.name, n_roi=S.nroi, n_sig=int(S.sig.sum()), n_noise=int(S.noise.sum()), n_neuropil=int(S.neuropil.sum()), dur_min=S.dur / 60, fs=S.fs, n_rej=0 if S.rej is None else int(S.rej['dff'].shape[0]),
               loco_note=S.loco_note, loco_ok=S.loco_ok, behaviours=sorted(S.beh), usable=False, why='')
    if S.sig.sum() < 3: out['why'] = 'fewer than 3 signal-like ROIs'; return out
    if 'pupil' not in S.beh and 'face' not in S.beh: out['why'] = 'no camera behaviour'; return out
    if S.dur < 240: out['why'] = 'recording shorter than 4 min'; return out
    out['usable'] = True
    bp = xc.bandpass
    sigD = np.stack([bp(r) for r in S.dff[S.sig]]); P = np.nanmean(S.dff[S.sig], 0); fP = bp(P)
    noiseD = np.stack([bp(r) for r in S.dff[S.noise]]) if S.noise.sum() >= 3 else None
    out['t'] = S.t; out['P'] = P.astype(np.float32); out['P_pipe'] = S.population('sig', 'pipe').astype(np.float32)
    out['beh_grid'] = {k: v.astype(np.float32) for k, v in S.beh.items()}
    # random control sets of rejected ROIs (same size as the signal-like set)
    rej_sets = []
    if S.rej is not None and S.rej['dff'].shape[0] >= max(S.sig.sum(), 10):
        R = S.rej['dff']; k = int(S.sig.sum())
        for _ in range(20): rej_sets.append(np.sort(rng.choice(R.shape[0], min(k, R.shape[0]), replace=False)))
    out['coh'] = dict(sig=coherence(sigD), noise=coherence(noiseD) if noiseD is not None else None,
                      rej=[coherence(np.stack([bp(r) for r in S.rej['dff'][ix]])) for ix in rej_sets[:10]] if rej_sets else [])
    fbeh = {b: bp(S.beh[b]) for b in BEHS if b in S.beh and (b != 'loco' or S.loco_ok)}
    out['xc'] = {}; out['roi'] = {}
    rejP = [np.nanmean(S.rej['dff'][ix], 0) for ix in rej_sets[:20]] if rej_sets else []
    for b, fy in fbeh.items():
        if np.isfinite(fy).mean() < 0.3: continue
        d = peak_and_ci(fP, fy)
        d['r_pipe'] = xc.xcorr_fft(bp(out['P_pipe']), fy, MAXLAG)[1]
        d['r_noise'] = xc.xcorr_fft(np.nanmean(noiseD, 0), fy, MAXLAG)[1] if noiseD is not None else None
        d['rej_peaks'] = [xc.peak_stats(xc.xcorr_fft(bp(rp), fy, MAXLAG)[2], MAXLAG, SEP)[0] for rp in rejP]
        d['r_rej_curves'] = np.stack([xc.xcorr_fft(bp(rp), fy, MAXLAG)[1] for rp in rejP]) if rejP else None
        out['xc'][b] = d
        lag_bins = int(round(d['peak_lag_s'] / DT)); rr = {}
        for grp, rows in (('sig', sigD), ('noise', noiseD), ('rej', np.stack([bp(r) for r in S.rej['dff'][:200]]) if S.rej is not None else None)):
            if rows is None or len(rows) == 0: continue
            r_i, p_i = roi_signs(rows, fy, lag_bins)
            rr[grp] = dict(n=len(r_i), pos=int(((r_i > 0) & (p_i < 0.01)).sum()), neg=int(((r_i < 0) & (p_i < 0.01)).sum()), frac_pos=float((r_i > 0).mean()), r_med=float(np.median(r_i)), r=r_i, p=p_i)
        out['roi'][b] = rr
    # ridge (held-out R2) on the population and on random rejected sets
    fams = {f: zs(bp(S.beh[f])) for f in FAMS if f in S.beh and (f != 'loco' or S.loco_ok) and np.isfinite(S.beh[f]).mean() > 0.3}
    lags_bins = (np.arange(-1.0, 3.01, 0.4) / DT).round().astype(int)
    if len(fams) >= 1:
        fam_use = {f: np.where(np.isfinite(v), v, np.nan) for f, v in fams.items()}
        out['ridge'] = dict(sig=ridge_partition(zs(fP), fam_use, lags_bins), families=list(fam_use))
        out['ridge']['rej'] = [ridge_partition(zs(bp(rp)), fam_use, lags_bins)['full'] for rp in rejP[:6]] if rejP else []
        if noiseD is not None: out['ridge']['noise'] = ridge_partition(zs(bp(np.nanmean(noiseD, 0))), fam_use, lags_bins)['full']
        null = []
        for _ in range(8):
            sh = rng.integers(SEP, len(fP) - SEP); null.append(ridge_partition(zs(np.roll(fP, sh)), {'face': fam_use['face']} if 'face' in fam_use else fam_use, lags_bins)['full'])
        out['ridge']['null_full'] = null
    # states (face-motion terciles; ball where plausible)
    st = {}
    if 'face' in S.beh:
        f = nsmooth_state(S.beh['face']); q1, q2 = np.nanpercentile(f, [33.3, 66.7]); quiet, active = f <= q1, f >= q2
        st['face_active_minus_quiet'] = block_boot_diff(P, active, quiet); st['face_active_minus_quiet_rej'] = [block_boot_diff(rp, active, quiet)[0] for rp in rejP[:10]]
    if 'loco' in S.beh and S.loco_ok:
        sp = np.expm1(S.beh['loco']); run, still = sp > 50, sp < 10
        if run.sum() > 100 and still.sum() > 100: st['run_minus_still'] = block_boot_diff(P, run, still); st['run_frac'] = float(np.nanmean(run[np.isfinite(sp)]))
    out['states'] = st
    return out


def nsmooth_state(x): return xc.nsmooth(x, 2.0 / DT)


if __name__ == '__main__':
    names = sys.argv[1:] or SESSIONS
    res = pickle.load(open(OUT, 'rb')) if os.path.exists(OUT) and sys.argv[1:] else {}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    for nm in names:
        t0 = time.time()
        try:
            S = xc.Session(nm); res[nm] = analyse(S)
            r = res[nm]; print('  -> %s | %.0f s' % ('usable' if r['usable'] else 'EXCLUDED: ' + r['why'], time.time() - t0), flush=True)
        except Exception as e:
            import traceback; traceback.print_exc(); res[nm] = dict(name=nm, usable=False, why='error: %s' % e)
        pickle.dump(res, open(OUT, 'wb'))
    print('STATE_DONE')
