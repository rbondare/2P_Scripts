"""Axon grouping layer, step 3: group the ACCEPTED ROIs.

Two levels (design principle 1: only the strict one is ever averaged):
  unit     = fragment merge: footprints within D_FRAG um, residual r >= R_MIN (FDR-significant), and members that pass the homogeneity test (below) -> averaged into one trace
  cluster  = co-active association of units within D_BUNDLE um: annotation only (cluster_id), traces are NEVER averaged

The agglomeration is greedy (like suite2p_analysis/merge_suite2p_rois.agglomerate, re-written here in pure numpy because that module imports suite2p / torch at load time):
repeatedly join the two groups with the highest correlation of their summed (residual) traces as long as the nearest footprints are <= d um apart and r >= r_min; the group trace
is recomputed after every join, so a noisy ROI can join once its group has averaged down. Guard against bridging two different axons: every member of the smaller side must itself
correlate significantly (p < P_MEMBER) and positively with the trace of the other side.

Homogeneity (is the unit mean a faithful stand-in for every member?): leave-one-out residual of every member  e_i = y_i - beta_i * mean(y_others)  on the band-passed dF/F,
beta fitted on one half of the session and the residual judged on the other half (cross-fit). For noisy copies of one signal e_i is noise at the level expected from the member's own
frame-to-frame noise, so  ratio_i = var(e_i) / (kappa (sigma_i^2 + beta_i^2 sigma_others^2)) is ~1; a member with its own events gives ratio >> 1. The threshold is NOT assumed:
ag_validate.calibrate_homogeneity builds it from simulated noisy copies that carry the REAL noise of the session (rejected noise blobs).
"""
import numpy as np
from ag_data import xc, v3_features, fill_rows, DT, GAP_SCAN_UM
import ag_pairs as ap

D_FRAG, D_BUNDLE = 5.0, 30.0            # um (see the gap diagnostic: residual r between accepted ROIs is above the null for gaps up to ~5-20 um)
R_MIN = 0.3                             # residual r of two groups to join (the value of merge_suite2p_rois.CORR_THRESH, there after PC1 removal too)
P_MEMBER = 0.01
R_HALF = 0.15                           # within-unit mean pair r each half of the session must reach (stability)


def _gaps(U, groups):
    G = U.gap_um(GAP_SCAN_UM); n = len(groups); M = np.full((n, n), np.inf)
    for a in range(n):
        for b in range(a + 1, n):
            M[a, b] = M[b, a] = G[np.ix_(groups[a], groups[b])].min()
    return M


def _ztrace(P, members):
    x = P.Zr[members].mean(0); v = P.v; y = np.zeros_like(x, dtype=np.float64); xv = x[v]; s = xv.std()
    if s > 1e-12: y[v] = (xv - xv.mean()) / s
    return y


def _acf(P, x):
    F = np.fft.rfft(x, n=2 * len(x)); return np.fft.irfft(F * np.conj(F))[:ap.MAXLAG_ACF + 1] / P.n


def _pair_z(P, x, y, ax=None, ay=None):
    """r, z, p of two z-scored full-length traces"""
    r = float(x @ y / P.n); ax = _acf(P, x) if ax is None else ax; ay = _acf(P, y) if ay is None else ay
    w = np.r_[1.0, 2.0 * np.ones(ap.MAXLAG_ACF)]; sg = np.sqrt(max(float((ax * w) @ ay), 1.0) / P.n); z = r / sg
    return r, z, float(P.p_of_z(z))


def agglomerate(P, init_groups, d_um, r_min=R_MIN, p_join=1e-3, verbose=False):
    """init_groups: list of lists of ROI rows of the Universe. Returns the merged list of lists (sorted by size) and the merge log."""
    U = P.U; groups = [list(g) for g in init_groups]; n = len(groups)
    tr = [_ztrace(P, g) for g in groups]; ac = [_acf(P, t) for t in tr]
    GM = _gaps(U, groups); alive = [True] * n; blocked = set(); log = []
    def corr_row(a):
        out = {}
        for b in range(len(groups)):
            if b == a or not alive[b] or GM[a, b] > d_um: continue
            r, z, p = _pair_z(P, tr[a], tr[b], ac[a], ac[b]); out[b] = (r, z, p)
        return out
    cand = {}                                                             # (a, b) a < b -> (r, z, p)
    for a in range(n):
        for b, v in corr_row(a).items():
            if a < b: cand[(a, b)] = v
    while True:
        best = None
        for (a, b), (r, z, p) in cand.items():
            if (a, b) in blocked or not (alive[a] and alive[b]): continue
            if r >= r_min and p <= p_join and (best is None or r > best[0]): best = (r, a, b)
        if best is None: break
        r, a, b = best
        # guard: each member of the smaller group must itself correlate (significantly, positively) with the other group's trace
        small, other = (a, b) if len(groups[a]) <= len(groups[b]) else (b, a)
        rr, zz, pp = P.vs_trace(tr[other], rows=groups[small])
        if len(groups[small]) > 1 and not np.all((rr > 0) & (pp <= P_MEMBER)):
            blocked.add((a, b)); continue
        groups[a] = groups[a] + groups[b]; alive[b] = False; tr[a] = _ztrace(P, groups[a]); ac[a] = _acf(P, tr[a])
        GM[a, :] = np.minimum(GM[a, :], GM[b, :]); GM[:, a] = GM[a, :]; GM[a, a] = 0
        for key in [k for k in cand if a in k or b in k]: cand.pop(key)
        blocked = {k for k in blocked if a not in k and b not in k}
        for bb, v in corr_row(a).items(): cand[(min(a, bb), max(a, bb))] = v
        log.append((r, len(groups[a])))
        if verbose: print('   join -> %d ROIs, r %.2f' % (len(groups[a]), r))
    out = [sorted(g) for g, al in zip(groups, alive) if al]; out.sort(key=lambda g: (-len(g), g[0])); return out, log


# ---------------------------------------------------------------------------------------------------------- homogeneity
_KAPPA = None


def kappa():
    """variance left by the 1-30 s band-pass on unit-variance white noise"""
    global _KAPPA
    if _KAPPA is None:
        w = np.random.default_rng(1).standard_normal((1, 40000)); b = xc.bandpass(w, ap.BP_LO, ap.BP_HI)[0]; _KAPPA = float(np.var(b[500:-500]))
    return _KAPPA


def member_ratio(X, hf, valid, halves):
    """X: n x T band-passed traces of the members (NaN -> 0), hf: white-noise sigma of the unfiltered 5 Hz dF/F per member. Returns ratio_i per member (cross-fit over the two halves)."""
    n = len(X); ratio = np.full(n, np.nan); kap = kappa()
    for i in range(n):
        oth = [k for k in range(n) if k != i]; m = X[oth].mean(0); s2m = (hf[oth] ** 2).sum() / len(oth) ** 2; vals = []
        for h in (0, 1):
            fit, test = halves[1 - h], halves[h]
            beta = float(X[i, fit] @ m[fit] / max(m[fit] @ m[fit], 1e-12)); e = X[i, test] - beta * m[test]
            vals.append(e.var() / (kap * (hf[i] ** 2 + beta ** 2 * s2m) + 1e-12))
        ratio[i] = np.mean(vals)
    return ratio


def one_factor_test(P, members):
    """Is the unit's residual correlation structure rank 1 (noisy copies of ONE signal)? Principal-axis one-factor fit to the members' residual correlation matrix; returns the largest
    |z| = |R_ij - lam_i lam_j| / sigma_ij over member pairs (sigma_ij: Bartlett). Needs >= 4 members: a one-factor model has n loadings for n(n-1)/2 correlations, so with 2 or 3 members
    the fit is exact and NOTHING can be tested (NaN) - a pair or triple of fragments cannot reveal from their own correlations whether they are one axon or several.
    Copies of one signal leave only noise; two axons sharing part of their signal leave a block structure (within-axon r > between-axon r)."""
    n = len(members)
    if n < 4: return np.nan
    R = P.r[np.ix_(members, members)].astype(float); S = P.sigma[np.ix_(members, members)].astype(float); h = np.clip(np.abs(R - np.eye(n)).max(1), 0.05, 0.99)
    for _ in range(200):
        Rr = R.copy(); np.fill_diagonal(Rr, h); w, V = np.linalg.eigh(Rr); lam = V[:, -1] * np.sqrt(max(w[-1], 1e-9)); lam = lam * (np.sign(lam.sum()) or 1.0)
        hn = np.clip(lam ** 2, 0.01, 0.999)
        if np.max(np.abs(hn - h)) < 1e-6: break
        h = hn
    E = R - np.outer(lam, lam); iu = np.triu_indices(n, 1); return float(np.max(np.abs(E[iu] / S[iu])))


def split_halves(P):
    idx = np.where(P.v)[0]; mid = idx[len(idx) // 2]; return [idx[idx < mid], idx[idx >= mid]]


def unit_homogeneity(P, members, halves=None):
    """per-member ratio (see module docstring); X = band-passed raw dF/F of the Universe (the domain that would be averaged)"""
    U = P.U; halves = halves or split_halves(P); X = np.where(np.isfinite(P.Y[members]), P.Y[members], 0.0); hf = U.v3()['hf'][members]
    return member_ratio(X, hf, P.v, halves)


# ---------------------------------------------------------------------------------------------------------- the whole accepted-ROI grouping
class Grouping:
    """units (strict), clusters (association) and tiers of the accepted ROIs of one session"""

    def __init__(self, P, d_frag=D_FRAG, d_bundle=D_BUNDLE, r_min=R_MIN, homog_thr=None, verbose=True):
        self.P, U = P, P.U; self.params = dict(d_frag=d_frag, d_bundle=d_bundle, r_min=r_min, level=P.level, npc=P.npc, homog_thr=homog_thr)
        acc = list(np.where(U.src == 0)[0]); halves = split_halves(P); self.halves = halves
        units, self.log_units = agglomerate(P, [[i] for i in acc], d_frag, r_min)
        # association level on the unit traces
        clus, self.log_clusters = agglomerate(P, units, d_bundle, r_min)
        self.units = units; self.clusters_of_units = clus
        self.unit_of = np.full(U.N, -1); self.cluster_of = np.full(U.N, -1)
        for k, g in enumerate(units): self.unit_of[g] = k
        for k, g in enumerate(clus): self.cluster_of[g] = k
        # unit traces (dF/F, plain mean) and the per-unit statistics
        self.unit_dff = np.array([np.nanmean(U.dff[g], 0) for g in units]) if units else np.zeros((0, U.T))
        v3u = v3_features(fill_rows(self.unit_dff)) if len(units) else None
        rows = []
        for k, g in enumerate(units):
            n = len(g); d = dict(unit=k, n=n, members=g, v3_pass=bool(v3u['passes'][k]) if v3u is not None else False)
            if n >= 2:
                C = [P.Zr[g][:, h] @ P.Zr[g][:, h].T / len(h) for h in halves]                       # pair r inside each half (residual domain)
                iu = np.triu_indices(n, 1); d['half_r'] = [float(c[iu].mean()) for c in C]; d['stable'] = bool(min(d['half_r']) >= R_HALF)
                d['ratio'] = unit_homogeneity(P, g, halves); d['ratio_max'] = float(np.nanmax(d['ratio']))      # descriptive: member-specific variance / noise expectation (little power, see ag_validate)
                d['factor_z'] = one_factor_test(P, g)
                d['homog'] = None if (homog_thr is None or n < 4) else bool(d['factor_z'] <= homog_thr(n))   # None = not testable (pair / triple) or no calibrated threshold
            rows.append(d)
        self.unit_info = rows
        # tiers. A = >= 4 members, stable in both halves, group mean passes V3, rank-1 test passed (homogeneity DEMONSTRATED): the only units whose mean is meant to replace the members.
        #        B = multi-member unit that is not A, reasons in d['why']: untestable_n2 / untestable_n3 (homogeneity cannot be checked) | heterogeneous (rank-1 test failed) | unstable | v3_fail
        #        C = singleton passing V3 | D = singleton failing V3
        self.tier = {}
        for d in rows:
            g = d['members']
            if d['n'] >= 2:
                why = []
                if d['n'] < 4: why.append('untestable_n%d' % d['n'])                      # pair / triple: homogeneity cannot be checked from the fragments' own correlations
                if not d['stable']: why.append('unstable')
                if not d['v3_pass']: why.append('v3_fail')
                if d['n'] >= 4 and d['homog'] is False: why.append('heterogeneous')
                if d['n'] >= 4 and d['homog'] is None: why.append('no_threshold')
                d['why'] = why; t = 'A' if not why else 'B'; d['avg_ok'] = (t == 'A')
            else:
                t = 'C' if bool(U.cls[g[0]] == 'sig') else 'D'; d['why'] = []; d['avg_ok'] = False
            for i in g: self.tier[i] = t
        if verbose: self.summary()

    def summary(self):
        U = self.P.U; tiers = list(self.tier.values()); n_units = len(self.units); multi = [d for d in self.unit_info if d['n'] >= 2]
        print('  grouping [%s]: %d accepted ROIs -> %d units (%d with >= 2 members, sizes %s), %d clusters; ROIs in tiers %s' % (
            self.P.level, len(self.tier), n_units, len(multi), sorted([d['n'] for d in multi], reverse=True)[:12], len(self.clusters_of_units),
            {t: tiers.count(t) for t in sorted(set(tiers))}))
