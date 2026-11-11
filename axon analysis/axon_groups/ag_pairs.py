"""Axon grouping layer, step 2: the pair statistic.

For every pair of ROIs: correlation of the band-passed dF/F AFTER the shared signal has been removed ("residual r"), a significance from a circular-shift null,
and the footprint gap in um (ag_data). Three removal levels, because what counts as 'the same axon' depends on what is removed first:
  raw     band-passed best dF/F only (shared arousal signal included: every axon correlates with every other, PC1 is 45-80 % of the variance)
  pc1     the top population component (PC1 of the signal-like accepted ROIs) projected out of every ROI (what merge_suite2p_rois.py does)
  pc1beh  PC1 AND the behaviour-predicted part (ridge on face / whisking / pupil / ball at lags -1..+3 s) removed: grouping is then blind to the readouts

Significance: z = r / sigma_ij with sigma_ij from Bartlett's formula (autocorrelations of both traces, i.e. the effective sample size), compared with a POOLED null of
standardised z values from circular shifts (>= 120 s) of one trace of every pair. Pooling gives a tail resolution of ~1e-7 (a per-pair shift null has ~6000 draws, too
coarse for Benjamini-Hochberg over 1e5 pairs), and keeps the heavy tails of transient-dominated traces.
"""
import numpy as np
from ag_data import xc, DT

BP_LO, BP_HI = 1.0, 30.0                 # s, as roi_heterogeneity.py
SEP = int(120 / DT)                      # minimum circular shift (bins), the convention of xs_common.peak_stats / state_signal
MAXLAG_ACF = int(90 / DT)                # autocorrelation lags used for Bartlett's variance (the 1-30 s band has no memory beyond that)
BEH_NAMES = ('face', 'whisk', 'pupil', 'loco')
BEH_LAGS = np.arange(-int(1 / DT), int(3 / DT) + 1)      # -1 .. +3 s, as axon_analysis_behaviour_V1
LEVELS = ('raw', 'pc1', 'pc1beh')


def zscore_valid(X, v):
    """rows z-scored over the valid frames; invalid frames set to 0 (keeps the length for FFT / circular shifts)"""
    Z = np.zeros_like(X, dtype=np.float32); Xv = X[:, v].astype(np.float64)
    m = np.nanmean(Xv, 1, keepdims=True); s = np.nanstd(Xv, 1, keepdims=True); s[s < 1e-9] = np.inf
    Zv = (Xv - m) / s; Zv[~np.isfinite(Zv)] = 0.0; Z[:, v] = Zv; return Z


def bandpassed(U):
    """N x T band-passed best dF/F (NaN where undefined) and the frame mask used for everything: not bad registration, finite in >= 95 % of the ROIs"""
    Y = xc.bandpass(U.dff, BP_LO, BP_HI)
    v = (~U.bad) & (np.isfinite(Y).mean(0) >= 0.95)
    return Y, v


def common_mode(Y, v, rows, k=1):
    """k leading temporal components (unit norm, zero at invalid frames) of the z-scored rows `rows` + the share of variance of each"""
    Z = zscore_valid(Y[rows], v)[:, v]; u, s, vt = np.linalg.svd(Z, full_matrices=False)
    Q = np.zeros((k, Y.shape[1]), np.float64); Q[:, v] = vt[:k]
    return Q, (s[:k] ** 2) / (s ** 2).sum()


def project_out(Y, v, Q):
    """remove the temporal components Q (k x T, orthonormal rows) from every row of Y (NaN -> 0 on valid frames)"""
    X = np.where(np.isfinite(Y), Y, 0.0).astype(np.float64); X[:, ~v] = 0.0
    return X - (X @ Q.T) @ Q


def behaviour_design(U, v):
    """lagged, band-passed, z-scored behaviour regressors on the grid; frames where any behaviour is missing (video gaps) are dropped from the mask"""
    S = U.S; cols = []; vb = v.copy(); names = []
    for nm in BEH_NAMES:
        if nm not in S.beh: continue
        b = xc.bandpass(S.beh[nm], BP_LO, BP_HI); vb &= np.isfinite(b); cols.append(b); names.append(nm)
    if not cols: return None, v, []
    B = np.array(cols); B = (B - np.nanmean(B[:, vb], 1, keepdims=True)) / np.nanstd(B[:, vb], 1, keepdims=True); B = np.where(np.isfinite(B), B, 0.0)
    X = np.concatenate([np.roll(B, -L, axis=1) for L in BEH_LAGS], axis=0)         # row: behaviour(t + L); circular wrap affects only the first / last 3 s
    return X, vb, names


def ridge_remove(Y, vb, X, lam_rel=0.1):
    """residual of every row of Y after a joint ridge fit on the regressors X (rows = regressors), on the valid frames vb; invalid frames 0"""
    Yv = np.where(np.isfinite(Y[:, vb]), Y[:, vb], 0.0).astype(np.float64); Xv = X[:, vb].astype(np.float64)
    G = Xv @ Xv.T; lam = lam_rel * np.trace(G) / len(G); B = np.linalg.solve(G + lam * np.eye(len(G)), Xv @ Yv.T)
    R = np.zeros(Y.shape, np.float64); R[:, vb] = Yv - (B.T @ Xv); return R


class PairStats:
    """Residual correlations of all ROIs of a Universe at one removal level"""

    def __init__(self, U, level='pc1', npc=1, n_null=150, seed=0, verbose=True):
        assert level in LEVELS
        self.U, self.level, self.npc = U, level, npc
        Y, v = bandpassed(U); self.Y = Y
        rows = np.where(U.src == 0)[0]; sigrows = np.where((U.src == 0) & (U.cls == 'sig'))[0]
        if len(sigrows) < 5: sigrows = rows
        self.Q, self.pc_share = common_mode(Y, v, sigrows, k=max(npc, 1))
        info = {}; self.v0 = v; self._D = None
        if level == 'raw': X = np.where(np.isfinite(Y), Y, 0.0); X[:, ~v] = 0.0; self.v = v
        elif level == 'pc1':
            X = project_out(Y, v, self.Q[:npc]); self.v = v
        else:
            X = project_out(Y, v, self.Q[:npc]); D, vb, names = behaviour_design(U, v)
            if D is None: raise ValueError('no behaviour on the grid for %s' % U.name)
            X = ridge_remove(X, vb, D); self.v = vb; info['beh'] = names; self._D = D
        v = self.v; self.n = int(v.sum()); self.info = info
        self.Zr = zscore_valid(X, v)                                         # residual, z-scored: input of all pair statistics
        self.res_frac = (np.var(X[:, v], 1) / np.maximum(np.var(np.where(np.isfinite(Y[:, v]), Y[:, v], 0.0), 1), 1e-12))     # variance left after removal
        Z = self.Zr; n = self.n; self.r, self.sigma, self.z = self._stats(Z)
        # pooled null of the standardised z from circular shifts of one trace
        rng = np.random.default_rng(seed); T = Z.shape[1]; shifts = rng.integers(SEP, T - SEP, n_null); nz = []
        iu = np.triu_indices(len(Z), 1)
        for s in shifts:
            R = (Z @ np.roll(Z, int(s), axis=1).T) / n; nz.append((R / self.sigma)[iu])
        self.null_abs = np.sort(np.abs(np.concatenate(nz)).astype(np.float32))
        if verbose:
            print('  [%s] n=%d frames, PC1 share %.2f, residual variance left median %.2f, null |z| 99.9%% %.2f (Gaussian 3.29), 1e-6 tail %.2f' % (
                level, n, self.pc_share[0], np.median(self.res_frac), np.quantile(self.null_abs, 0.999), self.null_abs[int(len(self.null_abs) * (1 - 1e-6))]))

    def _stats(self, Z):
        """r (N x N), Bartlett sigma_ij and z for the rows Z (zero-filled outside the valid frames)"""
        n = self.n; r = (Z @ Z.T / n).astype(np.float32); np.fill_diagonal(r, 1.0)
        # sigma_ij: Bartlett's variance of r for independent stationary series, sigma^2 = (1/n) sum_k rho_i(k) rho_j(k)
        F = np.fft.rfft(Z, n=2 * Z.shape[1], axis=1); acf = np.fft.irfft(F * np.conj(F), axis=1)[:, :MAXLAG_ACF + 1] / n
        w = np.r_[1.0, 2.0 * np.ones(MAXLAG_ACF)]; S2 = (acf * w) @ acf.T
        sigma = np.sqrt(np.maximum(S2, 1.0) / n).astype(np.float32)      # >= 1/sqrt(n): lag 0 term alone
        return r, sigma, r / sigma

    def with_rows(self, dff_rows):
        """copy of this PairStats whose ROI rows are replaced by `dff_rows` (N' x T unfiltered grid dF/F; same removal, same null): null runs and simulations"""
        import copy
        Q = copy.copy(self); Y, Z = self.process(dff_rows); Q.Y, Q.Zr = Y, Z; Q.r, Q.sigma, Q.z = Q._stats(Z); return Q

    def process(self, dff_rows):
        """run EXTRA unfiltered grid dF/F rows (k x T, e.g. simulated fragments) through exactly the same band-pass / removal / z-scoring: returns (band-passed rows, residual z-scored rows)"""
        Y = xc.bandpass(np.atleast_2d(dff_rows), BP_LO, BP_HI)
        if self.level == 'raw': X = np.where(np.isfinite(Y), Y, 0.0); X[:, ~self.v] = 0.0
        else:
            X = project_out(Y, self.v0, self.Q[:self.npc])
            if self.level == 'pc1beh': X = ridge_remove(X, self.v, self._D)
        return Y, zscore_valid(X, self.v)

    def p_of_z(self, z):
        """two-sided empirical p of |z| against the pooled shift null (resolution 1 / null size)"""
        a = np.abs(np.asarray(z)); k = len(self.null_abs) - np.searchsorted(self.null_abs, a, side='left'); return (1.0 + k) / (1.0 + len(self.null_abs))

    def vs_trace(self, x, rows=None):
        """r, z, p of residual rows against one external trace x (length T, same residualisation, z-scored over valid frames, 0 elsewhere)"""
        Z = self.Zr if rows is None else self.Zr[rows]; n = self.n
        r = Z @ x / n; F = np.fft.rfft(x, n=2 * len(x)); ax = (np.fft.irfft(F * np.conj(F))[:MAXLAG_ACF + 1]) / n
        Fz = np.fft.rfft(Z, n=2 * Z.shape[1], axis=1); az = np.fft.irfft(Fz * np.conj(Fz), axis=1)[:, :MAXLAG_ACF + 1] / n
        w = np.r_[1.0, 2.0 * np.ones(MAXLAG_ACF)]; sg = np.sqrt(np.maximum((az * w) @ ax, 1.0) / n); z = r / sg
        return r, z, self.p_of_z(z)


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values"""
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); q = np.empty(m); r = p[o] * m / (np.arange(m) + 1.0)
    q[o] = np.minimum.accumulate(r[::-1])[::-1]; return np.minimum(q, 1.0)


def candidate_pairs(U, dmax_um, rows=None):
    """(i, j) with i < j and footprint gap <= dmax_um; restricted to `rows` if given"""
    G = U.gap_um(max(dmax_um, 0)); ok = np.isfinite(G) & (G <= dmax_um); ok = np.triu(ok, 1)
    if rows is not None:
        m = np.zeros(U.N, bool); m[rows] = True; ok &= m[:, None] & m[None, :]
    return np.where(ok)


def gap_diagnostic(P, bins=(0, 0.01, 2, 5, 10, 20, 40, 60)):
    """residual r versus footprint gap for the three ROI-pair types, with the shuffle null: median r, 95th percentile of r, fraction with FDR < 0.01 and r >= 0.3"""
    U = P.U; I, J = candidate_pairs(U, bins[-1], rows=np.where(U.src < 2)[0]); G = U.gap_um(bins[-1])[I, J]; r = P.r[I, J]; z = P.z[I, J]; q = bh_fdr(P.p_of_z(z))   # hidden merge constituents excluded (they overlap their merged ROI by construction)
    kind = np.where((U.src[I] == 0) & (U.src[J] == 0), 'acc-acc', np.where((U.src[I] == 0) | (U.src[J] == 0), 'acc-rej', 'rej-rej'))
    out = []
    for k in ('acc-acc', 'acc-rej', 'rej-rej'):
        for a, b in zip(bins[:-1], bins[1:]):
            m = (kind == k) & (G >= a) & (G < b)          # first bin (0, 0.01): overlapping or touching footprints (gap 0)
            if m.sum() == 0: continue
            out.append(dict(kind=k, lo=a, hi=b, n=int(m.sum()), med_r=float(np.median(r[m])), p95_r=float(np.quantile(r[m], 0.95)), frac_sig=float(((q[m] < 0.01) & (r[m] >= 0.3)).mean())))
    return out
