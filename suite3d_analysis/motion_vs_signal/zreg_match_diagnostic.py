"""Steps 1-3 of suite2p-in-depth's volume-based z-registration (Baruchin et al. 2026, Methods "Volume-based
z-registration"), re-implemented on the x/y-registered binned movie, as a feasibility check for this recording:

    python zreg_match_diagnostic.py <diag_dir> [--ref-volumes 400] [--sigma 6] [--margin 12]

  1. reference image per included plane = mean of --ref-volumes consecutive volumes at the temporal midpoint
  2. correlation C[t, p, r] between the frame of plane p in cycle t and reference plane r (band-passed, Pearson)
  3. best-matching acquired plane p*(t, r) = argmax_p C; displacement d(t, r) = p* - r; bounds L_t = min_r d,
     U_t = max_r d; target plane k scored by sum_t [0 <= k + L_t < P] + [0 <= k + U_t < P]; ties -> mean best-match C.

NOT the package (not installed here): no per-frame suite2p rigid registration against each reference (the movie is
already registered in x/y by Suite3D), no blending / reconstruction. It answers: do the five planes differ enough that a
frame can be assigned to a depth, and how much does the assignment fluctuate?
"""
import argparse
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

ap = argparse.ArgumentParser()
ap.add_argument("diag")
ap.add_argument("--ref-volumes", type=int, default=400)
ap.add_argument("--sigma", type=float, default=6.0)
ap.add_argument("--margin", type=int, default=12)
ap.add_argument("--movie", default="reg_bin2.npy")
a = ap.parse_args()

mov = np.load(os.path.join(a.diag, a.movie), mmap_mode="r")
T, P, ny, nx = mov.shape
valid = np.zeros((ny, nx), bool)
valid[a.margin:-a.margin, a.margin:-a.margin] = True


def prep(img):            # band-pass + mask + zero-mean + unit norm
    b = img - gaussian_filter(img, a.sigma, mode="reflect")
    b = b[valid]
    b = b - b.mean()
    return b / (np.linalg.norm(b) + 1e-9)


mid = T // 2
sl = slice(mid - a.ref_volumes // 2, mid + a.ref_volumes // 2)
ref = np.asarray(mov[sl], dtype=np.float32).mean(0)            # (P, y, x)
R = np.stack([prep(ref[r]) for r in range(P)])                 # (P, npix)

C = np.empty((T, P, P), dtype=np.float32)                       # C[t, p, r]
for i in range(0, T, 200):
    blk = np.asarray(mov[i:i + 200], dtype=np.float32)
    for p in range(P):
        F = np.stack([prep(f) for f in blk[:, p]])              # (n, npix)
        C[i:i + 200, p, :] = F @ R.T
print("C[t, p, r]: mean over t")
np.set_printoptions(precision=3, suppress=True, linewidth=140)
print(C.mean(0))

best_p = C.argmax(axis=1)                                        # (T, r): best acquired plane for each reference plane
d = best_p - np.arange(P)[None, :]
L, U = d.min(1), d.max(1)
print("\nplane shift d(t,r) = p*(t,r) - r : fraction of (t,r) pairs for each value")
vals, cnt = np.unique(d, return_counts=True)
for v, c in zip(vals, cnt):
    print(f"   {v:+d}: {c / d.size:.3f}")
print(f"cycles where all five reference matches agree (L = U): {np.mean(L == U):.3f}")
print(f"cycles with U - L >= 2 planes (inconsistent matches): {np.mean(U - L >= 2):.3f}")

score = np.array([np.sum((0 <= k + L) & (k + L < P)) + np.sum((0 <= k + U) & (k + U < P)) for k in range(P)])
meanC = np.array([C[:, :, k].max(1).mean() for k in range(P)])
print("\ntarget-plane score S(k) =", score, " (max = %d)" % (2 * T))
print("mean best-match correlation per reference plane:", meanC)
best_k = int(np.lexsort((meanC, score))[-1])
print("selected target reference plane:", best_k)

# how decisive is the matching? margin between best and 2nd best acquired plane for each reference plane
srt = np.sort(C, axis=1)
margin = (srt[:, -1, :] - srt[:, -2, :])
print("median gap between best and second-best match correlation per reference plane:", np.median(margin, 0))
# noise floor: the same statistic for two independent halves of the data should agree if the assignment is real
half = np.arange(T) % 2 == 0
Ch1, Ch2 = C[half], C[~half]
print("\nodd/even-cycle consistency of the mean correlation profile (should be ~1 for a stable structure):",
      np.corrcoef(Ch1.mean(0).ravel(), Ch2.mean(0).ravel())[0, 1].round(4))

# temporal behaviour of the matches: do they persist (real axial excursions) or flicker (noise)?
cent = (best_p - np.arange(P)[None, :]).mean(1)                 # mean displacement per cycle
ac = [np.corrcoef(cent[:-k], cent[k:])[0, 1] for k in (1, 2, 5, 10, 25, 50)]
print("autocorrelation of the mean plane shift at lags 1,2,5,10,25,50 cycles:", np.round(ac, 3))
np.savez(os.path.join(a.diag, "zreg_match.npz"), C=C, best_p=best_p, L=L, U=U, score=score, meanC=meanC)

fig, ax = plt.subplots(2, 3, figsize=(15, 7.5), constrained_layout=True)
im = ax[0, 0].imshow(C.mean(0), vmin=C.mean(0).min(), vmax=C.mean(0).max(), cmap="viridis")
ax[0, 0].set_xlabel("reference plane r")
ax[0, 0].set_ylabel("acquired plane p")
ax[0, 0].set_title("mean correlation C[p, r]")
plt.colorbar(im, ax=ax[0, 0])
ax[0, 1].hist(d.ravel(), bins=np.arange(-P + 0.5, P, 1))
ax[0, 1].set_xlabel("best-match displacement p* - r (planes)")
ax[0, 1].set_ylabel("count")
tt = np.arange(T) / 5.0017
ax[0, 2].plot(tt, cent, lw=0.5)
ax[0, 2].set_xlabel("time (s)")
ax[0, 2].set_ylabel("mean displacement over reference planes")
ax[1, 0].plot(tt, L, lw=0.5, label="lower bound L")
ax[1, 0].plot(tt, U, lw=0.5, label="upper bound U")
ax[1, 0].legend()
ax[1, 0].set_xlabel("time (s)")
ax[1, 0].set_ylabel("bounds on axial shift (planes)")
ax[1, 1].bar(range(P), score)
ax[1, 1].set_xlabel("candidate target plane k")
ax[1, 1].set_ylabel("recoverability score S(k)")
ax[1, 2].hist(margin.ravel(), bins=60)
ax[1, 2].set_xlabel("best - second-best correlation")
fig.savefig(os.path.join(a.diag, "zreg_match.png"), dpi=90)
print("saved zreg_match.png")
