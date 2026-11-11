"""Axon grouping layer, step 5: the noise reference.

noise_pool(U): rejected ROIs that are plausible "noise blobs" = not hidden merge constituents, not touching / overlapping any accepted ROI (2 px dilation), failing the V3 transient
test, skew < 0.3 (the strict control of paper_figures/fig12_rejected_rois.py). `attached` ROIs (ag_attach, any confidence) are removed by the caller.
match_noise(...): brightness matching of noise blobs to a set of axon members (nearest neighbour without replacement on log F0, log npix, Fneu/F) -> the top-20 reference
plus random matched draws (null).
"""
import os, numpy as np
from scipy.ndimage import binary_dilation

SKEW_STRICT = 0.3


def footprint_geometry(U):
    """overlap / adjacent fractions of every ROI's pixels with the accepted ROIs' mask (and its 2 px dilation); accepted ROIs themselves get NaN"""
    if hasattr(U, '_geom'): return U._geom
    ops = np.load(os.path.join(U.plane_dir, 'ops.npy'), allow_pickle=True).item(); Ly, Lx = np.asarray(ops['meanImg']).shape
    acc = np.zeros((Ly, Lx), bool)
    for i in np.where(U.src == 0)[0]: acc[U.pix[i][:, 0], U.pix[i][:, 1]] = True
    dil = binary_dilation(acc, iterations=2); ov = np.full(U.N, np.nan); ad = np.full(U.N, np.nan)
    for i in np.where(U.src > 0)[0]: ov[i] = acc[U.pix[i][:, 0], U.pix[i][:, 1]].mean(); ad[i] = dil[U.pix[i][:, 0], U.pix[i][:, 1]].mean()
    U._geom = dict(overlap=ov, adjacent=ad); return U._geom


def noise_pool(U, exclude=None):
    """boolean mask over the Universe rows"""
    g = footprint_geometry(U); passes = U.v3()['passes']; sk = U.get_skew()
    m = (U.src == 1) & (~passes) & (sk < SKEW_STRICT) & (g['adjacent'] == 0)
    if exclude is not None: m &= ~np.asarray(exclude, bool)
    return m


def match_noise(U, member_rows, pool_mask, k=20, n_random=200, seed=0):
    """top-k noise blobs matched to `member_rows` (accepted ROIs of the groups) on log F0, log npix, Fneu/F, plus n_random random matched draws.
    Distance = Euclidean in standardised features (SD of the pooled members + pool); one pool ROI per member (without replacement), the k best matches are the reference.
    Returns dict(top=rows, draws=array(n_random, k), smd=standardised mean difference per feature of the top set versus the members)"""
    rng = np.random.default_rng(seed); pool = np.where(pool_mask)[0]; mem = np.asarray(member_rows)
    feat = lambda r: np.c_[np.log(np.maximum(U.F0[r], 1e-3)), np.log(U.npix[r]), np.log(np.maximum(U.F_over_Fneu[r], 1e-3))]
    Fm, Fp = feat(mem), feat(pool); sd = np.concatenate([Fm, Fp]).std(0); sd[sd < 1e-6] = 1.0
    used = set(); pairs = []
    D = (((Fm[:, None, :] - Fp[None, :, :]) / sd) ** 2).sum(2)
    order = np.argsort(D, axis=None)
    taken_m, taken_p = set(), set()
    for f in order:
        a, b = divmod(int(f), D.shape[1])
        if a in taken_m or b in taken_p: continue
        taken_m.add(a); taken_p.add(b); pairs.append((a, b, D[a, b]))
        if len(taken_m) == len(mem) or len(taken_p) == len(pool): break
    pairs.sort(key=lambda t: t[2]); top = pool[[b for _, b, _ in pairs[:k]]]
    # random matched draws: from the pool ROIs inside the matching caliper (nearest-member distance <= median distance of the top-k matches, at least 3 k ROIs)
    dmin = D.min(0); cal = max(np.median([d for _, _, d in pairs[:k]]), np.sort(dmin)[min(3 * k, len(dmin) - 1)])
    ok = pool[dmin <= cal]
    draws = np.array([rng.choice(ok, size=min(k, len(ok)), replace=False) for _ in range(n_random)])
    smd = (feat(top).mean(0) - Fm.mean(0)) / sd
    return dict(top=top, draws=draws, smd=smd, features=['log F0', 'log npix', 'log Fneu/F'])
