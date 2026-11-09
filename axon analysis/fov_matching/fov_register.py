"""Step 2-3a: pairwise registration of the mean images of all sessions in microns (masked normalised cross-correlation, several candidate peaks per pair).

For every pair (A, B) two feature images are correlated: AXON scale (DoG 0.8 / 6 um: thin bright axons + puncta) and VESSEL scale (DoG 3 / 15 um: dark vessel shadows,
large-scale brightness). Both are resampled to a COARSE grid (0.5 um/px) over the whole search space; the peak list of the average of the two correlation maps is kept
(top-K, non-maximum suppression 4 um, parabolic sub-pixel), together with the two single-feature values at the peak (agreement check) and a phase-randomised null
(maximum of the same map when B is replaced by a phase-randomised copy of itself: preserves its spatial power spectrum, destroys the arrangement).

Convention: t = (dy, dx) in um is the position of session B's pixel (0, 0) centre in session A's coordinate system, i.e. a point p_B of B is at p_B + t in A.
Result: cache/pairwise.pkl  {(A, B): dict(peaks=[dict(t, r, r_axon, r_vessel, overlap, z)], null_max, ...)} for A before B in SESSIONS.
Run: miniforge3\\envs\\suite2p_venv\\python.exe fov_register.py
"""
import os, sys, time, pickle, itertools, numpy as np
from scipy.signal import fftconvolve
from scipy import ndimage as ndi
from fov_common import *

CG = 0.5                                  # coarse grid for the global search (um/px)
FEATS = {'axon': (0.8, 6.0), 'vessel': (3.0, 15.0)}
MIN_OVERLAP = 0.25                        # fraction of the smaller image's valid area that must overlap
ERODE_UM = 2.0                            # additional border trimming (registration edge artefacts)
TOPK, SEP_UM = 5, 4.0
N_SURR = 6


def features(S, grid=CG):
    """dict feature -> (image on grid, mask on grid), image robust-z-scored inside the mask and zero outside"""
    out = {}
    m = to_grid(S.valid.astype(float), S.py, S.px, order=1, grid=grid) > 0.99
    m = ndi.binary_erosion(m, iterations=max(1, int(round(ERODE_UM / grid))))
    for k, (s1, s2) in FEATS.items():
        a = dog(S.mean, S.valid, s1, s2, S.py, S.px)
        g = to_grid(a, S.py, S.px, order=1, grid=grid)
        x = g[m]; g = (g - np.median(x)) / (np.std(x) + 1e-12)
        out[k] = (np.where(m, g, 0.0), m)
    return out


def _cc(x, y): return fftconvolve(x, y[::-1, ::-1], mode='full')


def masked_ncc(a, ma, b, mb):
    """normalised cross-correlation of a (fixed) and b (moved) over their overlap for every integer shift; returns r map (Ha+Hb-1, Wa+Wb-1) and overlap pixel count map.
    Index (i, j) <-> t = (i - (Hb - 1), j - (Wb - 1)) pixels (position of b's pixel (0, 0) in a)."""
    fa, fb = ma.astype(np.float64), mb.astype(np.float64)
    a = a * ma; b = b * mb
    N = _cc(fa, fb); Sa = _cc(a, fb); Sb = _cc(fa, b); Sab = _cc(a, b); Saa = _cc(a * a, fb); Sbb = _cc(fa, b * b)
    N = np.maximum(np.round(N), 0)
    with np.errstate(all='ignore'):
        num = Sab - Sa * Sb / N
        den = np.sqrt(np.maximum(Saa - Sa ** 2 / N, 0) * np.maximum(Sbb - Sb ** 2 / N, 0))
        r = num / den
    r[~np.isfinite(r)] = -1
    return r, N


def phase_randomise(b, mb, rng):
    F = np.fft.rfft2(b * mb)
    ph = np.exp(1j * rng.uniform(0, 2 * np.pi, F.shape))
    c = np.fft.irfft2(np.abs(F) * ph, s=b.shape)
    c = c * mb; x = c[mb]
    return (c - x.mean()) / (x.std() + 1e-12) * mb


def top_peaks(r, valid, shape_b, grid, k=TOPK, sep_um=SEP_UM):
    rr = np.where(valid, r, -1.0).copy(); out = []
    sep = int(round(sep_um / grid))
    for _ in range(k):
        i = np.unravel_index(np.argmax(rr), rr.shape)
        if rr[i] <= -1: break
        # parabolic sub-pixel refinement
        dy = dx = 0.0
        if 0 < i[0] < rr.shape[0] - 1 and 0 < i[1] < rr.shape[1] - 1:
            c = rr[i]; up, dn, lf, rt = rr[i[0] - 1, i[1]], rr[i[0] + 1, i[1]], rr[i[0], i[1] - 1], rr[i[0], i[1] + 1]
            if min(up, dn, lf, rt) > -1:
                den_y = up - 2 * c + dn; den_x = lf - 2 * c + rt
                dy = 0.5 * (up - dn) / den_y if den_y < 0 else 0.0; dx = 0.5 * (lf - rt) / den_x if den_x < 0 else 0.0
                dy, dx = float(np.clip(dy, -1, 1)), float(np.clip(dx, -1, 1))
        out.append((i, dy, dx, float(rr[i])))
        rr[max(0, i[0] - sep):i[0] + sep + 1, max(0, i[1] - sep):i[1] + sep + 1] = -1
    return out


def register_pair(FA, FB, rng, n_surr=N_SURR, grid=CG):
    """FA, FB: feature dicts of sessions A and B. Returns dict(peaks=[...], null_max=...)."""
    (aA, mA), (aB, mB) = FA['axon'], FB['axon']
    areaA, areaB = mA.sum(), mB.sum(); minN = MIN_OVERLAP * min(areaA, areaB)
    rs = {}; N = None
    for k in FEATS:
        rs[k], N = masked_ncc(FA[k][0], FA[k][1], FB[k][0], FB[k][1])
    valid = N >= minN
    rc = 0.5 * (rs['axon'] + rs['vessel'])
    Hb, Wb = aB.shape
    peaks = []
    for (i, dy, dx, rv) in top_peaks(rc, valid, aB.shape, grid):
        t = ((i[0] - (Hb - 1) + dy) * grid, (i[1] - (Wb - 1) + dx) * grid)
        peaks.append(dict(t=t, r=rv, r_axon=float(rs['axon'][i]), r_vessel=float(rs['vessel'][i]), overlap=float(N[i] / min(areaA, areaB))))
    # null: phase-randomised B (both features), maximum of the combined map under the same overlap constraint
    nm = []
    for s in range(n_surr):
        rr = 0
        for k in FEATS:
            sb = phase_randomise(FB[k][0], FB[k][1], rng)
            r, _ = masked_ncc(FA[k][0], FA[k][1], sb, FB[k][1]); rr = rr + 0.5 * r
        nm.append(float(np.max(np.where(valid, rr, -1))))
    # map statistics for a z-like peak height
    x = rc[valid]
    med, mad = float(np.median(x)), float(np.median(np.abs(x - np.median(x))) * 1.4826 + 1e-9)
    for p in peaks: p['z'] = (p['r'] - med) / mad
    return dict(peaks=peaks, null_max=float(np.max(nm)), null_mean=float(np.mean(nm)), map_med=med, map_mad=mad, n_valid=int(valid.sum()))


def main():
    t0 = time.time(); inv = read_inventory()
    F = {}
    for nm in SESSIONS:
        S = Sess(nm, inv, rois=False); F[nm] = features(S)
        print('features', short(nm), F[nm]['axon'][0].shape, '%.0f s' % (time.time() - t0), flush=True)
    res = {}; rng = np.random.default_rng(0)
    for A, B in itertools.combinations(SESSIONS, 2):
        res[(A, B)] = register_pair(F[A], F[B], rng)
        p = res[(A, B)]['peaks'][0]
        print('%s %s r %.2f (axon %.2f vessel %.2f) null %.2f  t=(%.1f, %.1f) ov %.2f  %.0f s' % (short(A), short(B), p['r'], p['r_axon'], p['r_vessel'], res[(A, B)]['null_max'], p['t'][0], p['t'][1], p['overlap'], time.time() - t0), flush=True)
    with open(os.path.join(CACHE, 'pairwise.pkl'), 'wb') as f: pickle.dump(dict(res=res, sessions=SESSIONS, grid=CG, feats=FEATS), f)
    print('done %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
