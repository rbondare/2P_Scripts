"""Heterogeneity across ROIs, PCA beyond PC1 and spatial structure (Figure 4).
usage: python roi_heterogeneity.py [session ...]   -> cache/roi_<session>.pkl
Per session (band-passed 1-30 s best dF/F on the 5 Hz grid):
  signal-like and noise-like ROIs: zero-lag r with face motion / whisking and r at +1.2 s with pupil, PCA of the signal-like ROIs (variance fractions, PC1 / PC2 scores, r of PC1 / PC2 with each behaviour),
  pairwise r between signal-like ROIs vs their distance in the FOV (Spearman, ROI-position permutation null), r with face motion vs position.
"""
import os, sys, pickle, numpy as np
from fig_common import *
import xs_common as xc
from scipy.stats import spearmanr
from event_triggered import BEH_LAG_S

DT = xc.DT


def corr_rows(R, y):
    ok = np.isfinite(y); out = np.full(len(R), np.nan)
    for i, r in enumerate(R):
        m = ok & np.isfinite(r)
        if m.sum() > 300: out[i] = np.corrcoef(r[m], y[m])[0, 1]
    return out


def shifted(y, lag_s):
    k = int(round(lag_s / DT)); b = np.full_like(y, np.nan)
    if k >= 0: b[:len(y) - k] = y[k:]
    else: b[-k:] = y[:len(y) + k]
    return b


def analyse(name, rng):
    S = xc.Session(name, with_rejected=False, verbose=False); out = dict(name=name, n_sig=int(S.sig.sum()), n_noise=int(S.noise.sum()))
    if S.sig.sum() < 3: out['why'] = 'fewer than 3 signal-like ROIs'; return out
    bp = xc.bandpass; cent = np.asarray(S.X['centroid'], float)[:, :2]; beh = {k: bp(S.beh[k]) for k in ('face', 'whisk', 'pupil') if k in S.beh}
    for grp, mask in (('sig', S.sig), ('noise', S.noise)):
        if mask.sum() == 0: continue
        R = np.array([bp(r) for r in S.dff[mask]]); out['r_' + grp] = {k: corr_rows(R, shifted(v, BEH_LAG_S.get(k, 0.0))) for k, v in beh.items()}
        out['xy_' + grp] = cent[mask]; out['idx_' + grp] = np.where(mask)[0]
        if grp == 'sig': Rsig = R
    ok = np.isfinite(Rsig).all(0); Z = Rsig[:, ok]; Z = (Z - Z.mean(1, keepdims=True)) / Z.std(1, keepdims=True)
    C = np.corrcoef(Z); w, V = np.linalg.eigh(C); order = np.argsort(w)[::-1]; w = w[order]; V = V[:, order]
    out['var_frac'] = (w / w.sum())[:8].astype(float); out['C'] = C.astype(np.float32); out['loadings'] = V[:, :2].astype(float)
    sc = (V[:, :2].T @ Z); full = np.full((2, Rsig.shape[1]), np.nan); full[:, ok] = sc; out['pc_scores'] = full.astype(np.float32); out['t'] = S.t
    out['pc_r'] = {k: [float(np.corrcoef(full[i][np.isfinite(full[i]) & np.isfinite(shifted(v, BEH_LAG_S.get(k, 0.0)))], shifted(v, BEH_LAG_S.get(k, 0.0))[np.isfinite(full[i]) & np.isfinite(shifted(v, BEH_LAG_S.get(k, 0.0)))])[0, 1]) for i in range(2)] for k, v in beh.items()}
    # after removing the behaviour-predicted part (ridge R2 on the first PC is in Fig 2): residual structure = PC2 / PC1 variance ratio
    out['pc2_over_pc1'] = float(w[1] / w[0])
    # space: pairwise r vs distance (signal-like ROIs)
    xy = out['xy_sig']; n = len(xy)
    if n >= 8:
        iu = np.triu_indices(n, 1); d = np.hypot(*(xy[iu[0]] - xy[iu[1]]).T); r = C[iu]
        rho = spearmanr(d, r)[0]; nul = []
        for _ in range(1000):
            p = rng.permutation(n); xyp = xy[p]; dd = np.hypot(*(xyp[iu[0]] - xyp[iu[1]]).T); nul.append(spearmanr(dd, r)[0])
        nul = np.array(nul); out['space'] = dict(rho=float(rho), p=float((1 + (np.abs(nul) >= abs(rho)).sum()) / (1 + len(nul))), null95=float(np.percentile(np.abs(nul), 95)), n=int(n), dist=d.astype(np.float32), r=r.astype(np.float32))
        # coupling vs position: ROI r_face vs distance to the FOV centre of mass of the best-coupled quartile
        rf = out['r_sig'].get('face')
        if rf is not None and np.isfinite(rf).sum() >= 8:
            q = rf >= np.nanpercentile(rf, 75); cm = xy[q].mean(0); dist_cm = np.hypot(*(xy - cm).T); rr = spearmanr(dist_cm[~q], rf[~q])[0]
            out['space']['rf_vs_dist_cm_rho'] = float(rr)
    return out


if __name__ == '__main__':
    names = sys.argv[1:] or ALL
    for nm in names:
        rng = np.random.default_rng(7)
        try:
            r = analyse(nm, rng); pickle.dump(r, open(os.path.join(CACHE, 'roi_%s.pkl' % nm), 'wb'))
            print(short(nm), 'sig %d noise %d' % (r['n_sig'], r['n_noise']), 'PC var', np.round(r.get('var_frac', [])[:3], 2), 'space', {k: round(v, 3) for k, v in r.get('space', {}).items() if k in ('rho', 'p')}, flush=True)
        except Exception as e:
            import traceback; traceback.print_exc()
