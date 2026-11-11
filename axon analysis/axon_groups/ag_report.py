"""Axon grouping layer, step 7 (report card): does grouping help, and does it hide or manufacture anything?  ONE function per question, all on real data.

held_out(U, level, fold): groups are defined on one half of the session (Grouping on a time slice), everything is then judged on the OTHER half:
  A  SNR gain      band-passed SNR of the unit mean vs its members and the ideal sqrt(n), against a control: the mean of n spatially near accepted ROIs that were NOT grouped together
  B  coupling spread   do members of a unit reproducibly differ in their coupling to behaviour? (member-minus-unit coupling in half 1 vs half 2, Spearman over all members)
  C  inflation     PC1 share and mean pair r of the accepted ROIs (ungrouped, ROI level) vs of the unit means in-sample / held-out vs RANDOM merges of the same sizes at the same places.
                   Fewer, cleaner units raise the PC1 share mechanically: the random-merge control is what the grouped value has to be compared with.
"""
import numpy as np
from scipy.stats import spearmanr
from ag_data import xc, v3_features, fill_rows
import ag_pairs as ap, ag_group as ag

BEH_COUPLING = ('face', 'whisk', 'pupil', 'loco')


def snr_band(P, dff_rows):
    """band-passed (1-30 s) SD over the white-noise level expected from the frame-to-frame differences (kappa of the filter): comparable between a member and a unit mean"""
    X = np.atleast_2d(dff_rows); Y = xc.bandpass(X, ap.BP_LO, ap.BP_HI)[:, P.v]; hf = v3_features(fill_rows(X))['hf']
    return np.nanstd(Y, axis=1) / (np.sqrt(ag.kappa()) * np.maximum(hf, 1e-9))


def pc1_share(Z):
    """share of variance of the first PC of z-scored rows (rows = ROIs / units)"""
    if len(Z) < 3: return np.nan
    s = np.linalg.svd(Z - Z.mean(1, keepdims=True), compute_uv=False); return float(s[0] ** 2 / (s ** 2).sum())


def mean_pair_r(Z):
    C = np.corrcoef(Z); iu = np.triu_indices(len(Z), 1); return float(np.mean(C[iu]))


def _z(P, X):
    Y = xc.bandpass(np.atleast_2d(X), ap.BP_LO, ap.BP_HI); return ap.zscore_valid(np.where(np.isfinite(Y), Y, 0.0), P.v)[:, P.v]


def beh_coupling(P, rows_dff):
    """r of band-passed rows with every behaviour (lag 0) on the valid frames: n_rows x n_beh (NaN where the behaviour is missing)"""
    S = P.U.S; Z = _z(P, rows_dff); out = np.full((len(Z), len(BEH_COUPLING)), np.nan)
    for k, nm in enumerate(BEH_COUPLING):
        if nm not in S.beh: continue
        b = xc.bandpass(S.beh[nm], ap.BP_LO, ap.BP_HI)[P.v]
        if not np.isfinite(b).all(): ok = np.isfinite(b)
        else: ok = np.ones(len(b), bool)
        if ok.sum() < 200: continue
        bz = (b[ok] - b[ok].mean()) / b[ok].std(); out[:, k] = Z[:, ok] @ bz / ok.sum()
    return out


def held_out(U, level='pc1', fold=0, homog_thr=None, n_ctrl=100, seed=0, verbose=True):
    rng = np.random.default_rng(seed); T = U.T; mid = T // 2; sl = [(0, mid), (mid, T)]; a, b = sl[fold], sl[1 - fold]
    UA, UB = U.time_slice(*a), U.time_slice(*b)
    PA = ap.PairStats(UA, level, verbose=False); GA = ag.Grouping(PA, homog_thr=homog_thr, verbose=False); PB = ap.PairStats(UB, level, verbose=False)
    acc = np.where(U.src == 0)[0]; gap = U.gap_um(60.0); res = dict(fold=fold, units=[], n_units=len(GA.units))
    snr_roi = dict(zip(acc, snr_band(PB, UB.dff[acc])))
    multi = [d for d in GA.unit_info if d['n'] >= 2]
    beh_roi = dict(zip(acc, beh_coupling(PB, UB.dff[acc]))) if True else {}
    rows = []; spread_h1, spread_h2 = [], []
    for d in multi:
        g = d['members']; n = len(g); mean_dff = np.nanmean(UB.dff[g], 0)
        s_unit = float(snr_band(PB, mean_dff)[0]); s_mem = np.array([snr_roi[i] for i in g]); gain = s_unit / np.median(s_mem)
        # control: n accepted ROIs near the unit (gap <= D_BUNDLE from any member), NOT all from this unit
        near = [i for i in acc if i not in g and (gap[i, g] <= ag.D_BUNDLE).any()]; cg = []; cr = []
        iu = np.triu_indices(n, 1); res_r = float(PB.r[np.ix_(g, g)][iu].mean())                  # residual r among the members in the HELD-OUT half
        if len(near) >= n:
            for _ in range(n_ctrl):
                c = list(rng.choice(near, n, replace=False)); cg.append(float(snr_band(PB, np.nanmean(UB.dff[c], 0))[0] / np.median([snr_roi[i] for i in c]))); cr.append(float(PB.r[np.ix_(c, c)][iu].mean()))
        rows.append(dict(n=n, tier=GA.tier[g[0]], why=d['why'], snr_unit=s_unit, snr_member_median=float(np.median(s_mem)), gain=float(gain), ideal=float(np.sqrt(n)),
                         ctrl_gain_median=float(np.median(cg)) if cg else np.nan, ctrl_gain_p95=float(np.quantile(cg, 0.95)) if cg else np.nan,
                         res_r_heldout=res_r, res_r_ctrl_median=float(np.median(cr)) if cr else np.nan, res_r_ctrl_p95=float(np.quantile(cr, 0.95)) if cr else np.nan, members=[int(U.idx[i]) for i in g]))
    res['units'] = rows
    # B: do members of a unit reproducibly differ in their coupling to behaviour? (member minus unit-mean coupling, half 1 vs half 2)
    U1, U2 = U.time_slice(0, mid), U.time_slice(mid, T); P1 = ap.PairStats(U1, level, verbose=False); P2 = ap.PairStats(U2, level, verbose=False)
    dev1, dev2, snrs = [], [], []
    for d in multi:
        g = d['members']; c1 = beh_coupling(P1, U1.dff[g]); c2 = beh_coupling(P2, U2.dff[g])
        with np.errstate(all='ignore'): dev1.append(c1 - np.nanmean(c1, 0)); dev2.append(c2 - np.nanmean(c2, 0))
        sn = 0.5 * (snr_band(P1, U1.dff[g]) + snr_band(P2, U2.dff[g])); snrs.append(sn - sn.mean())
    res['coupling_spread'] = {}
    if dev1:
        D1, D2 = np.concatenate(dev1), np.concatenate(dev2)
        for k, nm in enumerate(BEH_COUPLING):
            ok = np.isfinite(D1[:, k]) & np.isfinite(D2[:, k])
            if ok.sum() >= 8:
                rho, p = spearmanr(D1[ok, k], D2[ok, k]); sn = np.concatenate(snrs); rs, ps = spearmanr(0.5 * (D1[ok, k] + D2[ok, k]), sn[ok])
                res['coupling_spread'][nm] = dict(rho=float(rho), p=float(p), n=int(ok.sum()), sd_dev=float(np.nanstd(D1[ok, k])), rho_with_snr=float(rs), p_snr=float(ps))
    # C: inflation of shared-signal statistics
    Zroi_B = _z(PB, UB.dff[acc]); Zroi_A = _z(PA, UA.dff[acc])
    unitsA = [d['members'] for d in GA.unit_info]
    Zu_A = _z(PA, np.array([np.nanmean(UA.dff[g], 0) for g in unitsA])); Zu_B = _z(PB, np.array([np.nanmean(UB.dff[g], 0) for g in unitsA]))
    rnd = []
    for _ in range(30):                                                                  # random merges of the same sizes among accepted ROIs, kept near each other
        left = list(acc); mem = []
        for g in sorted(unitsA, key=len, reverse=True):
            if len(g) == 1: mem.append([left.pop(int(rng.integers(len(left))))]) if left else None; continue
            seed = g[0]; pool = [i for i in left if gap[seed, i] <= ag.D_BUNDLE or i == seed]
            if len(pool) < len(g): pool = left
            pick = list(rng.choice(pool, min(len(g), len(pool)), replace=False)); mem.append(pick); left = [i for i in left if i not in pick]
        rnd.append(_z(PB, np.array([np.nanmean(UB.dff[g], 0) for g in mem if len(g)])))
    res['inflation'] = dict(pc1_roi_A=pc1_share(Zroi_A), pc1_roi_B=pc1_share(Zroi_B), pc1_units_in_sample=pc1_share(Zu_A), pc1_units_held_out=pc1_share(Zu_B),
                            pc1_random_merge_held_out=float(np.mean([pc1_share(z) for z in rnd])), r_roi_B=mean_pair_r(Zroi_B), r_units_in_sample=mean_pair_r(Zu_A),
                            r_units_held_out=mean_pair_r(Zu_B), r_random_merge_held_out=float(np.mean([mean_pair_r(z) for z in rnd])), n_roi=len(acc), n_units=len(unitsA))
    if verbose:
        print('  held-out fold %d: %d units from half %d (%d multi-member), judged on the other half' % (fold, len(unitsA), fold, len(multi)))
        for r in rows: print('    n=%d tier %s why %s: SNR gain %.2f (ideal %.2f; near-ROI control median %.2f p95 %.2f) | residual r among members, held-out %.2f (near-ROI control median %.2f p95 %.2f)' % (
            r['n'], r['tier'], ','.join(r['why']) or '-', r['gain'], r['ideal'], r['ctrl_gain_median'], r['ctrl_gain_p95'], r['res_r_heldout'], r['res_r_ctrl_median'], r['res_r_ctrl_p95']))
        print('    PC1 share: ROI level %.2f | units in-sample %.2f -> held-out %.2f | random merges (held-out) %.2f   mean pair r: ROIs %.2f | units in-sample %.2f -> held-out %.2f | random %.2f' % (
            res['inflation']['pc1_roi_B'], res['inflation']['pc1_units_in_sample'], res['inflation']['pc1_units_held_out'], res['inflation']['pc1_random_merge_held_out'],
            res['inflation']['r_roi_B'], res['inflation']['r_units_in_sample'], res['inflation']['r_units_held_out'], res['inflation']['r_random_merge_held_out']))
        for nm, v in res['coupling_spread'].items(): print('    members differ reproducibly in coupling to %s? rho %+.2f (p %.3f, n %d members, SD of member-unit deviation %.2f); deviation vs member SNR rho %+.2f' % (nm, v['rho'], v['p'], v['n'], v['sd_dev'], v['rho_with_snr']))
    return res


def compare_noise(U, P, R, n_draw=200, seed=0, verbose=True):
    """units (>= 2 accepted ROIs) versus the matched noise reference. For a unit of n members the reference is the MEAN OF n matched noise blobs (same n, same averaging, so SNR and n
    are not confounded with signal); n_draw random draws from the brightness-matched caliper plus the deterministic top-20 set. Metrics per trace: band-passed SNR, V3 transient test,
    upward events per minute, |r| with behaviour (mean over face / whisking / pupil / ball, lag 0, band-passed 1-30 s). Behaviours are used only here, to compare, never to select ROIs."""
    rng = np.random.default_rng(seed); draws = R['noise_draws']; units = [(k, g) for k, g in enumerate(R['units']) if len(g) >= 2]; rows = []
    def metrics(X):
        X = np.atleast_2d(X); v3 = v3_features(fill_rows(X)); c = beh_coupling(P, X)
        return dict(snr=snr_band(P, X), v3=v3['passes'].astype(float), ev_per_min=v3['ev_pos'] / (U.T * 0.2 / 60), abs_r_beh=np.nanmean(np.abs(c), axis=1))
    for k, g in units:
        n = len(g); mu = metrics(np.nanmean(U.dff[g], 0)); nz = []
        for d in range(n_draw):
            pick = rng.choice(draws[int(rng.integers(len(draws)))], min(n, draws.shape[1]), replace=False); nz.append(metrics(np.nanmean(U.dff[pick], 0)))
        top = metrics(np.nanmean(U.dff[R['noise_top'][:n]], 0)); row = dict(unit=k, n=n, tier=R['tier'][g[0]])
        for m in ('snr', 'v3', 'ev_per_min', 'abs_r_beh'):
            arr = np.array([z[m][0] for z in nz]); row[m] = dict(unit=float(mu[m][0]), noise_median=float(np.nanmedian(arr)), noise_p95=float(np.nanquantile(arr, 0.95)), frac_noise_ge_unit=float(np.nanmean(arr >= mu[m][0])), top20_mean=float(top[m][0]))
        rows.append(row)
    if verbose:
        for r in rows:
            print('    unit %2d n=%d tier %s | SNR %.1f (noise median %.1f, p95 %.1f) | events/min %.1f (noise %.1f) | |r| beh %.2f (noise %.2f, p95 %.2f; fraction of noise draws >= unit %.2f)' % (
                r['unit'], r['n'], r['tier'], r['snr']['unit'], r['snr']['noise_median'], r['snr']['noise_p95'], r['ev_per_min']['unit'], r['ev_per_min']['noise_median'],
                r['abs_r_beh']['unit'], r['abs_r_beh']['noise_median'], r['abs_r_beh']['noise_p95'], r['abs_r_beh']['frac_noise_ge_unit']))
    return rows
