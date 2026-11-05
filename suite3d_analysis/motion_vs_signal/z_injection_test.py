"""How small an axial displacement could the 5-plane profile of THIS recording reveal? Injection test on the real movie.

    python z_injection_test.py <diag_dir> [--movie regraw_bin2.npy] [--dz-um 3]

The estimator is the gain + shift model of z_motion_estimate.py (I_t - M ~= dg_t * M + d_t * dM/dplane, band-passed,
closed-form least squares over the whole field). Here a KNOWN displacement or gain is added to the real binned movie
(real noise, real lateral motion and real calcium are all still in it), the estimator is run again and the response is
compared with what was injected:
  * axial shift of delta planes: DeltaI = M(z + delta) - M(z), cubic interpolation of the mean volume along the 5 planes
    (non-linear in delta, so not simply the model that is being fitted); edge planes are extrapolated linearly
  * gain g: DeltaI = g * M  (what a calcium-like, shape-preserving brightening looks like)
Pulses of different amplitude and duration are injected at random times; a pulse counts as detected when the mean estimated
displacement during the pulse (minus the 5 s before it) exceeds the 99th percentile of the same statistic at random
times without injection. Reports the minimum amplitude detected with 80% power for each duration, the estimator's
gain (estimated/injected) and the leak of a gain change into the displacement estimate.
Nothing here validates the physical plane spacing: amplitudes are in PLANES (x dz-um for micrometres if the header's
spacing is right).
"""
import argparse
import os

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter

ap = argparse.ArgumentParser()
ap.add_argument("diag")
ap.add_argument("--movie", default="regraw_bin2.npy")
ap.add_argument("--sigma", type=float, default=6.0)
ap.add_argument("--margin", type=int, default=12)
ap.add_argument("--dz-um", type=float, default=3.0)
ap.add_argument("--fs", type=float, default=5.0017)
a = ap.parse_args()
FS = a.fs

mov = np.load(os.path.join(a.diag, a.movie), mmap_mode="r")
T, nz, ny, nx = mov.shape
M = np.zeros((nz, ny, nx), np.float64)
for i in range(0, T, 250):
    M += np.asarray(mov[i:i + 250], dtype=np.float32).sum(0)
M = (M / T).astype(np.float32)

valid = np.zeros((ny, nx), bool)
valid[a.margin:-a.margin, a.margin:-a.margin] = True
valid &= (M > 0).all(0)


def bandpass(v):
    return v - gaussian_filter(v, (0, a.sigma, a.sigma), mode="reflect")


D = np.gradient(M, axis=0)
BM, BD = bandpass(M) * valid, bandpass(D) * valid
P1, P2 = bandpass(BM), bandpass(BD)
A, B, C = (BM * BM).sum(), (BM * BD).sum(), (BD * BD).sum()
det = A * C - B * B

s1 = np.empty(T)
s2 = np.empty(T)
for i in range(0, T, 200):
    blk = np.asarray(mov[i:i + 200], dtype=np.float32)
    s1[i:i + 200] = np.tensordot(blk, P1, axes=([1, 2, 3], [0, 1, 2]))
    s2[i:i + 200] = np.tensordot(blk, P2, axes=([1, 2, 3], [0, 1, 2]))
m1, m2 = s1.mean(), s2.mean()


def estimate(s1_, s2_):
    x1, x2 = s1_ - m1, s2_ - m2
    return (A * x2 - B * x1) / det, (C * x1 - B * x2) / det          # d, dg


d0, g0 = estimate(s1, s2)
print(f"{T} volumes, Gram |corr(M, D)| = {abs(B) / np.sqrt(A * C):.3f}; real-data d SD {d0.std():.3f} planes, gain SD {g0.std():.4f}")

# projection of an injected shift / gain onto the two regressors, tabulated over a grid of delta
zk = np.arange(nz)
cs = CubicSpline(zk, M, axis=0, bc_type="natural", extrapolate=True)
grid = np.round(np.arange(-1.5, 1.5001, 0.05), 3)
A1 = np.empty(len(grid))
A2 = np.empty(len(grid))
for k, dl in enumerate(grid):
    dM = cs(zk + dl) - M
    A1[k], A2[k] = (dM * P1).sum(), (dM * P2).sum()
G1, G2 = (M * P1).sum(), (M * P2).sum()                              # projection of a unit gain


def inject_shift(delta_t):
    return np.interp(delta_t, grid, A1), np.interp(delta_t, grid, A2)


# linearity of the estimator: response to constant shifts
print("\nresponse to a constant injected shift (mean over all volumes), planes:")
for dl in (0.1, 0.25, 0.5, 1.0):
    i1, i2 = inject_shift(np.full(T, dl))
    d1, _ = estimate(s1 + i1, s2 + i2)
    print(f"  injected {dl:+.2f} -> estimated {np.mean(d1 - d0):+.3f}  ({np.mean(d1 - d0) / dl:.2f}x)")
i1, i2 = G1 * 0.3, G2 * 0.3
d1, g1 = estimate(s1 + i1, s2 + i2)
print(f"injected gain +0.30 -> estimated gain {np.mean(g1 - g0):+.3f}, leak into d {np.mean(d1 - d0):+.4f} planes "
      f"(= {np.mean(d1 - d0) / 0.3:+.3f} planes per unit gain)")

# pulses
rng = np.random.default_rng(0)
amps = [0.1, 0.2, 0.3, 0.5, 1.0]
durs = {"0.2 s (1 vol)": 1, "1 s": 5, "5 s": 25, "20 s": 100}
pre = int(5 * FS)
print("\nfraction of injected axial pulses detected (alpha 1%, null = same statistic at random times):")
print("duration        " + "  ".join(f"{am:4.1f} pl" for am in amps) + "   min amp @80% (planes / um)")
for name, L in durs.items():
    starts_pool = np.arange(pre + 1, T - L - 1)
    # null: statistic at random times in the real data
    st0 = rng.choice(starts_pool, 3000)
    null = np.array([d0[s:s + L].mean() - d0[s - pre:s].mean() for s in st0])
    thr_hi, thr_lo = np.percentile(null, 99.5), np.percentile(null, 0.5)
    row, minamp = [], None
    for am in amps:
        det_ = 0
        n = 300
        for s in rng.choice(starts_pool, n):
            dl = np.zeros(T)
            dl[s:s + L] = am
            i1, i2 = inject_shift(dl)
            d1, _ = estimate(s1 + i1, s2 + i2)
            st = d1[s:s + L].mean() - d1[s - pre:s].mean()
            det_ += st > thr_hi
        frac = det_ / n
        row.append(frac)
        if minamp is None and frac >= 0.8:
            minamp = am
    ma = ("> %.1f" % amps[-1]) if minamp is None else "%.1f / %.1f" % (minamp, minamp * a.dz_um)
    print(f"{name:14s}  " + "  ".join(f"{f:6.2f}" for f in row) + f"   {ma}")
    print(f"                (null 99% band of the statistic: [{thr_lo:+.3f}, {thr_hi:+.3f}] planes)")
