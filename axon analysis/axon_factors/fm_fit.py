"""Factor model of the accepted ROIs of one session:  z_i(t) = sum_k lambda_ik f_k(t) + e_i(t),  e_i ~ N(0, psi_i)   (factor analysis, ROI-specific noise variance psi_i).

Why a factor model instead of merging ROIs: every ROI keeps its OWN loading (amplitude), the combination of ROIs that estimates a factor is weighted by loading / noise
(optimal for 'same signal, different amplitude'), nothing is decided by thresholding a correlation, and the number of shared signals is a measured quantity (cross-validated rank).

Steps (analyse()):
  1 data          accepted ROIs, band-passed (1-30 s) best dF/F ('band', default) or the unfiltered 5 Hz dF/F ('dff'), valid frames only, every row z-scored (its SD is kept)
  2 rank          leave-block-out (8 contiguous blocks, 60 s buffer) AND leave-ROI-out: ROI i (and every ROI whose footprint is within 2 um of it) is predicted from the factor scores of the OTHER ROIs;
                  held-out R2 per ROI for k = 0..K_MAX; k* = first k whose gain in the median held-out R2 of the signal-like ROIs is < MIN_GAIN
  3 factors       full fit at k*, scores, share of the standardised variance per factor, varimax loadings (maps), sign convention (mean loading > 0)
  4 what are they correlation / ridge R2 with behaviour, and with nuisance signals (rigid shift, mean neuropil of the accepted ROIs, PC1 of the rejected ROIs not touching any accepted ROI)
  5 stability     split-half principal angles between the loading subspaces vs random subspaces
  6 residual      leave-one-ROI-out residual of every accepted ROI: does it still carry upward transients (V3 asymmetry test)? the same for rejected noise blobs projected on the factors
  7 rejected ROIs held-out R2 of every rejected ROI regressed on the factor scores (amplitude-aware candidate list)
"""
import os, sys, numpy as np
from scipy.optimize import linear_sum_assignment
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'axon_groups'))
import ag_data as ad, ag_pairs as ap, ag_noise as an, ag_group as ag          # noqa: E402
from ag_data import xc, v3_features, fill_rows, DT                              # noqa: E402
from sklearn.decomposition import FactorAnalysis
from scipy.stats import skew                                # noqa: E402

K_MAX, N_BLOCKS, BUF, MIN_GAIN, EX_UM = 8, 8, int(60 / DT), 0.01, 2.0
CACHE = r'C:\Users\rbondarenko\2P_cache\axon_factors'


def get_matrix(U, domain):
    if domain == 'band': return ap.bandpassed(U)
    Y = U.dff.astype(np.float64); v = (~U.bad) & (np.isfinite(Y).mean(0) >= 0.95); return Y, v


def prep(Y, v, rows):
    X = np.where(np.isfinite(Y[rows]), Y[rows], 0.0)[:, v]; X = X - X.mean(1, keepdims=True); sd = X.std(1, keepdims=True); sd[sd < 1e-9] = 1.0
    return X / sd, sd[:, 0]


def fit_fa(Z, k):
    fa = FactorAnalysis(n_components=k, max_iter=300, tol=5e-3, svd_method='lapack', random_state=0).fit(Z.T)
    return fa.components_.T.copy(), np.maximum(fa.noise_variance_, 1e-3)


def loo_predict(L, psi, Z, ex_sets):
    """leave-ROI-out prediction of every row of Z from the factor scores estimated with the OTHER rows (rows in ex_sets[i], i included, are left out)"""
    n, k = L.shape; W = L / psi[:, None]; A = W.T @ Z; M = np.eye(k) + L.T @ W; pred = np.empty_like(Z)
    for i in range(n):
        ex = ex_sets[i]; S = np.linalg.solve(M - L[ex].T @ W[ex], A - W[ex].T @ Z[ex]); pred[i] = L[i] @ S
    return pred


def scores(L, psi, Z):
    W = L / psi[:, None]; return np.linalg.solve(np.eye(L.shape[1]) + L.T @ W, W.T @ Z)


def varimax(L, q=100, tol=1e-6):
    p, k = L.shape; R = np.eye(k); d = 0.0
    for _ in range(q):
        Lr = L @ R; u, s, vh = np.linalg.svd(L.T @ (Lr ** 3 - Lr @ np.diag((Lr ** 2).sum(0)) / p)); R = u @ vh; dn = s.sum()
        if d != 0 and dn < d * (1 + tol): break
        d = dn
    return R


def cv_rank(Z, ex_sets, kmax=K_MAX):
    n, Tv = Z.shape; edges = np.linspace(0, Tv, N_BLOCKS + 1).astype(int); num = np.zeros((kmax + 1, n)); den = np.zeros(n)
    for b in range(N_BLOCKS):
        a, e = edges[b], edges[b + 1]; tr = np.r_[np.arange(0, max(a - BUF, 0)), np.arange(min(e + BUF, Tv), Tv)]; Zte = Z[:, a:e]
        den += (Zte ** 2).sum(1); num[0] += (Zte ** 2).sum(1)
        for k in range(1, kmax + 1):
            L, psi = fit_fa(Z[:, tr], k); num[k] += ((Zte - loo_predict(L, psi, Zte, ex_sets)) ** 2).sum(1)
    return 1 - num / den


def ridge_r2(X, y, halves, lam_rel=0.1):
    out = []
    for h in (0, 1):
        tr, te = halves[h], halves[1 - h]; G = X[:, tr] @ X[:, tr].T; lam = lam_rel * np.trace(G) / len(G); b = np.linalg.solve(G + lam * np.eye(len(G)), X[:, tr] @ y[tr]); p = b @ X[:, te]
        out.append(1 - ((y[te] - p) ** 2).sum() / max(((y[te] - y[te].mean()) ** 2).sum(), 1e-12))
    return float(np.mean(out))


def principal_cos(L1, L2):
    Q1, _ = np.linalg.qr(L1); Q2, _ = np.linalg.qr(L2); return np.sort(np.linalg.svd(Q1.T @ Q2, compute_uv=False))[::-1]


def analyse(name, domain='band', k_fixed=None, verbose=True):
    U = ad.load_universe(name, verbose=False); Y, v = get_matrix(U, domain); acc = np.where(U.src == 0)[0]; sigm = (U.cls[acc] == 'sig'); n = len(acc)
    Z, sd = prep(Y, v, acc); Tv = Z.shape[1]; gap = U.gap_um(ad.GAP_SCAN_UM); ex_sets = [np.where(gap[acc[i], acc] <= EX_UM)[0] for i in range(n)]
    R2 = cv_rank(Z, ex_sets); med = np.median(R2[:, sigm], 1); gain = np.diff(med)
    kstar = k_fixed or max(next((k for k in range(1, K_MAX + 1) if (med[k] - med[k - 1]) < MIN_GAIN), K_MAX + 1) - 1, 1)     # last k that still gained >= MIN_GAIN
    out = dict(name=name, domain=domain, n_roi=n, n_sig=int(sigm.sum()), Tv=Tv, R2=R2, med_R2_sig=med, med_R2_all=np.median(R2, 1), kstar=int(kstar), acc_idx=U.idx[acc], cls=U.cls[acc], sd=sd)
    # ---- factors at k*
    L, psi = fit_fa(Z, kstar); F = scores(L, psi, Z); order = np.argsort(-(L ** 2).sum(0)); L, F = L[:, order], F[order]
    sgn = np.sign(L.sum(0)); sgn[sgn == 0] = 1; L, F = L * sgn, F * sgn[:, None]; Rv = varimax(L) if kstar > 1 else np.eye(kstar); Lv, Fv = L @ Rv, Rv.T @ F
    share = (L ** 2).sum(0) / n; share_sig = (L[sigm] ** 2).sum(0) / max(sigm.sum(), 1); out.update(L=L, psi=psi, F=F, Lv=Lv, Fv=Fv, share=share, share_sig=share_sig, communality=(L ** 2).sum(1))
    # ---- behaviour and nuisance (frames of v)
    pos = np.where(v)[0]; Dm, vb, bnames = ap.behaviour_design(U, v); inb = vb[pos]; beh = {}
    for nm in ('face', 'whisk', 'pupil', 'loco'):
        if nm in U.S.beh: b = xc.bandpass(U.S.beh[nm], ap.BP_LO, ap.BP_HI)[pos] if domain == 'band' else U.S.beh[nm][pos]; beh[nm] = b
    nuis = {}
    if 'shift' in U.S.beh: nuis['rigid shift'] = (xc.bandpass(U.S.beh['shift'], ap.BP_LO, ap.BP_HI) if domain == 'band' else U.S.beh['shift'])[pos]
    npil = np.nanmean(U.Fng[acc], 0); nuis['mean neuropil'] = (xc.bandpass(npil[None], ap.BP_LO, ap.BP_HI)[0] if domain == 'band' else npil)[pos]
    geom = an.footprint_geometry(U); rej = np.where((U.src == 1) & (geom['adjacent'] == 0))[0]
    Zr, _ = prep(Y, v, rej) if len(rej) >= 10 else (None, None)
    if Zr is not None: u, s, vt = np.linalg.svd(Zr, full_matrices=False); nuis['rejected-ROI PC1'] = vt[0]
    def corr(a, b):
        ok = np.isfinite(a) & np.isfinite(b); return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 200 else np.nan
    out['corr_beh'] = {nm: [corr(F[k], b) for k in range(kstar)] for nm, b in beh.items()}; out['corr_nuis'] = {nm: [corr(F[k], b) for k in range(kstar)] for nm, b in nuis.items()}
    # ridge R2 of every factor on the lagged behaviours (frames with behaviour, halves) and PC1-of-accepted for orientation
    inbpos = np.where(inb)[0]; halves = [inbpos[:len(inbpos) // 2], inbpos[len(inbpos) // 2:]]; Dp = Dm[:, pos]
    out['ridge_r2_beh'] = [ridge_r2(Dp, F[k], halves) if len(inbpos) > 400 else np.nan for k in range(kstar)]; out['behaviours'] = bnames
    u, s, vt = np.linalg.svd(Z[sigm] if sigm.sum() >= 3 else Z, full_matrices=False); out['corr_pc1'] = [corr(F[k], vt[0]) for k in range(kstar)]; out['pc1_share_sig'] = float(s[0] ** 2 / (s ** 2).sum())
    # ---- split-half stability of the loading subspace
    h = Tv // 2; La, _ = fit_fa(Z[:, :h], kstar); Lb, _ = fit_fa(Z[:, h:], kstar); cos = principal_cos(La, Lb); rng = np.random.default_rng(0)
    nullc = np.array([principal_cos(rng.standard_normal((n, kstar)), rng.standard_normal((n, kstar))) for _ in range(200)])
    out['split_half_cos'] = cos; out['split_half_null_mean'] = nullc.mean(0); out['split_half_null_p95'] = np.quantile(nullc, 0.95, axis=0)
    # ---- leave-one-ROI-out residual of accepted ROIs, V3 asymmetry on it
    pred = loo_predict(L, psi, Z, ex_sets); E = Z - pred; v3 = v3_features(E); out['resid_skew'] = skew(E, axis=1); out['resid_pass'] = v3['passes']; out['resid_slow'] = v3['slow_ratio']; out['resid_ev'] = np.c_[v3['ev_pos'], v3['ev_neg']]
    out['loo_r2_full'] = 1 - (E ** 2).sum(1) / (Z ** 2).sum(1)
    # ---- rejected ROIs on the factor scores (ridge, halves) and the same residual test for the strict noise pool
    rr = np.where(U.src == 1)[0]; Zrj, _ = prep(Y, v, rr); half_t = [np.arange(Tv // 2), np.arange(Tv // 2, Tv)]; Fh = F
    r2rej = np.zeros(len(rr)); r2rej1 = np.zeros(len(rr)); slope = np.zeros((len(rr), kstar)); res_pass = np.zeros(len(rr), bool)
    for j in range(len(rr)):
        r2rej[j] = ridge_r2(Fh, Zrj[j], half_t, lam_rel=1e-3); r2rej1[j] = ridge_r2(Fh[:1], Zrj[j], half_t, lam_rel=1e-3)
        G = Fh @ Fh.T; b = np.linalg.solve(G + 1e-3 * np.trace(G) / kstar * np.eye(kstar), Fh @ Zrj[j]); slope[j] = b
    resid = Zrj - slope @ Fh; vr = v3_features(resid); res_pass = vr['passes']
    pool = an.noise_pool(U)[rr]; out.update(rej_idx=U.idx[rr], rej_r2=r2rej, rej_r2_f1=r2rej1, rej_r2_partial=np.where(kstar > 1, (r2rej - r2rej1) / np.maximum(1 - r2rej1, 1e-6), 0.0), rej_resid_skew=skew(resid, axis=1), rej_slope=slope, rej_resid_pass=res_pass, rej_in_noise_pool=pool, rej_inmerge=U.src[rr] == 2)
    if verbose:
        print('%s [%s]: %d accepted (%d signal-like), %d frames | held-out R2 (median, signal-like) by k: %s -> k* = %d' % (name[-12:], domain, n, sigm.sum(), Tv, ' '.join('%.2f' % x for x in med), kstar))
        print('   variance share of the factors (all ROIs / signal-like): %s / %s | PC1 share of signal-like ROIs %.2f' % (np.round(share, 2), np.round(share_sig, 2), out['pc1_share_sig']))
        for k in range(kstar):
            print('   factor %d: r with PC1 %+.2f | behaviour ridge R2 %.2f | max |r| beh %s | nuisance r %s' % (k + 1, out['corr_pc1'][k], out['ridge_r2_beh'][k], ' '.join('%s %+.2f' % (nm, c[k]) for nm, c in out['corr_beh'].items()), ' '.join('%s %+.2f' % (nm, c[k]) for nm, c in out['corr_nuis'].items())))
        print('   split-half subspace cosines %s (random %s, p95 %s)' % (np.round(cos, 2), np.round(out['split_half_null_mean'], 2), np.round(out['split_half_null_p95'], 2)))
        print('   leave-one-ROI-out residual still shows upward transients (V3 asymmetry test): signal-like %d / %d, noise-like %d / %d | rejected noise pool %d / %d, other rejected %d / %d' % (
            out['resid_pass'][sigm].sum(), sigm.sum(), out['resid_pass'][~sigm].sum(), (~sigm).sum(), res_pass[pool].sum(), pool.sum(), res_pass[~pool].sum(), (~pool).sum()))
        print('   residual skewness (positive = upward transients): signal-like median %.2f, noise-like %.2f | rejected noise pool %.2f, other rejected %.2f' % (np.median(out['resid_skew'][sigm]), np.median(out['resid_skew'][~sigm]) if (~sigm).any() else np.nan, np.median(out['rej_resid_skew'][pool]), np.median(out['rej_resid_skew'][~pool])))
        if kstar > 1: print('   rejected ROIs: R2 on factor 1 alone > 0.3: %d (noise pool %d) | extra R2 from factors >= 2 (partial) > 0.1: %d (noise pool %d)' % ((out['rej_r2_f1'] > 0.3).sum(), (out['rej_r2_f1'][pool] > 0.3).sum(), (out['rej_r2_partial'] > 0.1).sum(), (out['rej_r2_partial'][pool] > 0.1).sum()))
        print('   rejected ROIs held-out R2 on the factors: noise pool median %.2f, p95 %.2f | others median %.2f, p95 %.2f; > 0.3: %d of %d (noise pool %d)' % (
            np.median(r2rej[pool]), np.quantile(r2rej[pool], .95), np.median(r2rej[~pool]), np.quantile(r2rej[~pool], .95), (r2rej > 0.3).sum(), len(rr), (r2rej[pool] > 0.3).sum()))
    return out
