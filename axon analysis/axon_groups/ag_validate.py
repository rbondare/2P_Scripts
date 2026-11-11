"""Axon grouping layer, step 7: validation. Everything here answers one question: what does grouping do to the data?

  null_runs            the whole grouping pipeline on data with NO shared structure: (i) every ROI independently circular-shifted, (ii) independent ROIs + a planted common mode
                       (the real state signal at realistic loadings). Expected: ~0 edges and ~0 groups. Anything else is manufactured by the method.
  calibrate_homogeneity  simulated noisy copies of ONE real signal carrying the REAL noise of the session (rejected noise blobs): the distribution of the member-specific residual ratio
                       (ag_group.member_ratio) under 'all members are the same axon' -> the threshold of the homogeneity gate, per unit size
  sim_hetero           two adjacent axons that share only a fraction of their signal (true correlation rho): how often are they merged, does the gate catch it,
                       how much of each axon's own signal survives in the merged mean?
Simulated fragments are built in the unfiltered grid dF/F domain and pass through exactly the same band-pass / common-mode removal as real ROIs (PairStats.process).
"""
import copy, numpy as np
from ag_data import xc, v3_features, fill_rows, DT
import ag_pairs as ap, ag_group as ag, ag_noise


class SimU:
    """minimal Universe stand-in for simulated fragments: all footprints touching, so only correlation and the guard decide"""

    def __init__(self, U, hf):
        self.N = len(hf); self.T = U.T; self.src = np.zeros(self.N, int); self.cls = np.array(['sig'] * self.N, dtype=object); self.um = U.um; self._hf = hf

    def gap_um(self, d): return np.zeros((self.N, self.N), np.float32)

    def v3(self): return dict(hf=self._hf)


def band_sd(x): return float(np.nanstd(xc.bandpass(np.atleast_2d(x), ap.BP_LO, ap.BP_HI)[0][5:-5]))


class Sim:
    def __init__(self, U, P, seed=0):
        self.U, self.P = U, P; self.rng = np.random.default_rng(seed); D = fill_rows(U.dff); self.D = D
        sig = np.where((U.src == 0) & (U.cls == 'sig'))[0]
        if len(sig) < 3: sig = np.where(U.src == 0)[0]
        self.sigD = D[sig]; self.state = self.sigD.mean(0); self.state_sd = band_sd(self.state)
        pool = np.where(ag_noise.noise_pool(U))[0]
        if len(pool) < 10: raise ValueError('noise pool too small (%d) for simulation' % len(pool))
        self.noiseD = D[pool]; self.T = U.T

    def shifted(self, x):
        return np.roll(x, int(self.rng.integers(ap.SEP, self.T - ap.SEP)))

    def real_signal(self):
        """one real signal-like trace at a random circular shift (events and kinetics of a real axon, decoupled from the state signal and from the other draws)"""
        return self.shifted(self.sigD[int(self.rng.integers(len(self.sigD)))])

    def fragment(self, s, snr=None, state_rel=None):
        """unfiltered fragment of signal s: s + state coupling + real noise blob scaled to the band-passed SNR"""
        snr = float(np.exp(self.rng.uniform(np.log(0.7), np.log(3.0)))) if snr is None else snr
        g = self.rng.uniform(0.5, 2.0) if state_rel is None else state_rel
        nz = self.shifted(self.noiseD[int(self.rng.integers(len(self.noiseD)))]); nz = nz - nz.mean(); sd_s = band_sd(s)
        c = sd_s / (snr * max(band_sd(nz), 1e-9)); return s + g * sd_s / self.state_sd * self.state + c * nz

    def simp(self, frags):
        """(SimP, fragments) for a list of unfiltered fragments"""
        Yraw = np.array(frags); hf = v3_features(Yraw)['hf']; Y, Z = self.P.process(Yraw)
        S = copy.copy(self.P); S.U = SimU(self.U, hf); S.Y = Y; S.Zr = Z; S.r, S.sigma, S.z = S._stats(Z); return S


def _thr_fn(table, q):
    ns = sorted(table); vals = {n: float(np.quantile(table[n], q)) for n in ns}
    def f(n):
        k = min(max(n, ns[0]), ns[-1]); return vals[k]
    f.table = vals; return f


def calibrate_homogeneity(U, P, ns=(4, 5, 6), n_rep=150, q=0.95, seed=0):
    """null distribution of the rank-1 statistic (ag_group.one_factor_test) of simulated homogeneous units of n fragments (noisy copies of one real signal with real noise blobs);
    also the older member-ratio statistic and the merge success. Returns (threshold function of n at quantile q, {n: factor_z array}, {n: ratio array}, merge-success rate per n)"""
    sim = Sim(U, P, seed); table = {}; rtab = {}; success = {}
    for n in ns:
        vals = []; rv = []; ok = 0
        for _ in range(n_rep):
            s = sim.real_signal(); S = sim.simp([sim.fragment(s) for _ in range(n)]); mem = list(range(n))
            vals.append(ag.one_factor_test(S, mem)); rv.append(float(np.nanmax(ag.unit_homogeneity(S, mem))))
            g, _l = ag.agglomerate(S, [[i] for i in mem], 1e9); ok += (len(g) == 1)
        table[n] = np.array(vals); rtab[n] = np.array(rv); success[n] = ok / n_rep
    return _thr_fn(table, q), table, rtab, success


def sim_hetero(U, P, rhos=(0.95, 0.9, 0.8, 0.7, 0.5, 0.3), mixes=((2, 2), (3, 3), (2, 1)), n_rep=40, thr=None, seed=1):
    """two axons A and B (mA, mB fragments) whose signals have true correlation rho (rho = f^2: f c + sqrt(1 - f^2) u_A, f c + sqrt(1 - f^2) u_B).
    Per (mix, rho): fraction of runs in which fragments of A and B end up in ONE group (false merge), fraction of those flagged by the rank-1 test (needs >= 3 members),
    fraction where all fragments of A merge (recovery), R2 of axon A's own signal explained by the merged mean (what averaging loses)."""
    sim = Sim(U, P, seed); rows = []
    for mA, mB in mixes:
        for rho in rhos:
            f = np.sqrt(rho); fm = fg = ntest = rec = 0; r2 = []
            for _ in range(n_rep):
                sdref = band_sd(sim.real_signal()); sc = lambda x: x / max(band_sd(x), 1e-9) * sdref               # equal band-passed SD for the three components
                c, uA, uB = sc(sim.real_signal()), sc(sim.real_signal()), sc(sim.real_signal()); sA = f * c + np.sqrt(1 - rho) * uA; sB = f * c + np.sqrt(1 - rho) * uB
                frags = [sim.fragment(sA) for _ in range(mA)] + [sim.fragment(sB) for _ in range(mB)]
                S = sim.simp(frags); g, _l = ag.agglomerate(S, [[i] for i in range(mA + mB)], 1e9)
                ga = [set(x) for x in g]; A = set(range(mA)); B = set(range(mA, mA + mB)); rec += any(A <= x for x in ga)
                mix = [x for x in ga if (x & A) and (x & B)]
                if mix:
                    fm += 1; mem = sorted(mix[0])
                    if len(mem) >= 4 and thr is not None:
                        ntest += 1; fg += ag.one_factor_test(S, mem) > thr(len(mem))
                    mean = S.Y[mem].mean(0); yA = xc.bandpass(sA[None], ap.BP_LO, ap.BP_HI)[0]; v = S.v; r = np.corrcoef(mean[v], yA[v])[0, 1]; r2.append(r * r)
            rows.append(dict(mix='%d+%d' % (mA, mB), rho=rho, false_merge=fm / n_rep, tested=ntest, flagged=(fg / ntest if ntest else np.nan), recovered_A=rec / n_rep,
                             r2_A_in_merged=float(np.mean(r2)) if r2 else np.nan))
    return rows


def null_runs(U, P, seed=2):
    """(i) every ROI independently circular-shifted; (ii) independent ROIs + the real state signal planted at realistic loadings.
    Returns, per null, the number of FDR-significant candidate edges with r >= R_MIN among accepted pairs within D_BUNDLE, and the number of multi-ROI units / clusters."""
    rng = np.random.default_rng(seed); D = fill_rows(U.dff); T = U.T; out = {}
    shifted = np.array([np.roll(D[i], int(rng.integers(ap.SEP, T - ap.SEP))) for i in range(U.N)])
    sig = np.where((U.src == 0) & (U.cls == 'sig'))[0]; m = D[sig].mean(0); mc = m - m.mean()
    load = np.array([(D[i] - D[i].mean()) @ mc / (mc @ mc) for i in range(U.N)])
    planted = shifted + np.outer(np.clip(load, 0, None), m)
    acc = list(np.where(U.src == 0)[0]); I, J = ap.candidate_pairs(U, ag.D_BUNDLE, rows=np.array(acc))
    for name, X in (('independent circular shifts', shifted), ('independent + planted common mode', planted)):
        P0 = P.with_rows(X); q = ap.bh_fdr(P0.p_of_z(P0.z[I, J])); edges = int(((q < 0.01) & (P0.r[I, J] >= ag.R_MIN)).sum())
        units, _ = ag.agglomerate(P0, [[i] for i in acc], ag.D_FRAG); clus, _ = ag.agglomerate(P0, units, ag.D_BUNDLE)
        out[name] = dict(candidate_pairs=len(I), edges=edges, multi_units=sum(len(u) > 1 for u in units), multi_clusters=sum(len(c) > 1 for c in clus), max_r=float(np.max(P0.r[I, J])))
    # reference: the real data
    q = ap.bh_fdr(P.p_of_z(P.z[I, J])); out['real data'] = dict(candidate_pairs=len(I), edges=int(((q < 0.01) & (P.r[I, J] >= ag.R_MIN)).sum()), max_r=float(np.max(P.r[I, J])))
    return out


def param_sweep(U, levels=('raw', 'pc1', 'pc1beh'), npcs=(1, 2), d_frags=(3.0, 5.0, 10.0), r_mins=(0.3, 0.5), verbose=True):
    """how many accepted ROIs end up in multi-ROI units when the removal level, the number of removed PCs, the fragment distance and r_min change"""
    rows = []; acc = [[i] for i in np.where(U.src == 0)[0]]
    for level in levels:
        for npc in ((1,) if level == 'raw' else npcs):
            try: P = ap.PairStats(U, level, npc=npc, verbose=False)
            except Exception as e:
                if verbose: print('   sweep: level %s failed (%s)' % (level, e))
                continue
            for d in d_frags:
                for r in r_mins:
                    units, _ = ag.agglomerate(P, acc, d, r_min=r); multi = [u for u in units if len(u) > 1]
                    rows.append(dict(level=level, npc=npc, d_frag=d, r_min=r, n_units=len(units), n_multi=len(multi), rois_in_multi=sum(len(u) for u in multi), largest=max(len(u) for u in units)))
    if verbose:
        for r in rows: print('   sweep %-6s npc %d d_frag %4.1f r_min %.1f: %2d units, %2d multi-member (%2d ROIs), largest %d' % (r['level'], r['npc'], r['d_frag'], r['r_min'], r['n_units'], r['n_multi'], r['rois_in_multi'], r['largest']))
    return rows
