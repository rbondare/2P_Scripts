"""Step 3c: sub-micron refinement of the placement + independent landmark check + local-deformation estimate.

1. Fine pose refinement: every placed pair with a consistent coarse peak is re-correlated on the 0.25 um grid (AXON-scale feature only, window +-REFINE_WIN um around the
   coarse prediction, parabolic sub-pixel), then the poses are re-solved by robust weighted least squares (same scheme as fov_pose.refine).
2. Landmark check: the very bright blob (cell-like object with dark halo, ref-frame position ~ (23, 100) um) is located independently in every session that contains it
   (maximum of the 1.5 um-smoothed mean image inside +-LM_WIN of the predicted position); the spread of predicted - found positions is an alignment residual that does
   not come from the registration step itself.
3. Local deformation: for strongly correlated same-day / neighbouring pairs the overlap is cut into blocks and the per-block residual shift relative to the global shift
   is measured (how much a rigid shift misses; sets the tolerance of the ROI matching).
Output: cache/poses_fine.pkl, fov_pose_fine.csv, fov_landmark_check.csv, fov_local_deformation.csv
"""
import os, pickle, csv, numpy as np
from scipy import ndimage as ndi
from fov_common import *
from fov_register import features, masked_ncc, top_peaks
from fov_pose import peaks_between, load

FG = 0.25
REFINE_WIN = 3.0
LM_REF = None            # set from the data: the brightest smoothed-mean location of the reference session
LM_WIN = 8.0


def fine_edge(Fa, Fb, t0):
    """masked NCC at FG um/px between axon-scale features, restricted to |t - t0| <= REFINE_WIN; returns (t, r) or None"""
    (a, ma), (b, mb) = Fa['axon'], Fb['axon']
    r, N = masked_ncc(a, ma, b, mb)
    Hb, Wb = b.shape
    ii, jj = np.meshgrid(np.arange(r.shape[0]), np.arange(r.shape[1]), indexing='ij')
    ty = (ii - (Hb - 1)) * FG; tx = (jj - (Wb - 1)) * FG
    win = (np.hypot(ty - t0[0], tx - t0[1]) <= REFINE_WIN) & (N >= 0.2 * min(ma.sum(), mb.sum()))
    if not win.any(): return None
    rr = np.where(win, r, -1.0); i = np.unravel_index(np.argmax(rr), rr.shape)
    dy = dx = 0.0
    if 0 < i[0] < rr.shape[0] - 1 and 0 < i[1] < rr.shape[1] - 1:
        c = rr[i]; up, dn, lf, rt = rr[i[0] - 1, i[1]], rr[i[0] + 1, i[1]], rr[i[0], i[1] - 1], rr[i[0], i[1] + 1]
        if min(up, dn, lf, rt) > -1:
            dyy = up - 2 * c + dn; dxx = lf - 2 * c + rt
            dy = float(np.clip(0.5 * (up - dn) / dyy, -1, 1)) if dyy < 0 else 0.0; dx = float(np.clip(0.5 * (lf - rt) / dxx, -1, 1)) if dxx < 0 else 0.0
    return np.array([(i[0] - (Hb - 1) + dy) * FG, (i[1] - (Wb - 1) + dx) * FG]), float(rr[i]), float(N[i] / min(ma.sum(), mb.sum()))


def main():
    inv = read_inventory(); res = load()
    pk = pickle.load(open(os.path.join(CACHE, 'poses.pkl'), 'rb')); P = pk['P']; names = [s for s in SESSIONS if s in P]
    Sx = {n: Sess(n, inv, rois=False) for n in names}
    F = {n: features(Sx[n], grid=FG) for n in names}
    # ---- 1. fine edges
    edges = []
    for ia, a in enumerate(names):
        for b in names[ia + 1:]:
            pred = P[b] - P[a]
            cp, info = peaks_between(res, a, b)
            best = min(cp, key=lambda x: np.hypot(*(x[0] - pred)))
            if np.hypot(*(best[0] - pred)) > 2.5 or best[1]['r'] < 0.45: continue
            out = fine_edge(F[a], F[b], pred)
            if out is None: continue
            t, r, ov = out
            if np.hypot(*(t - pred)) > REFINE_WIN: continue
            edges.append((a, b, t, r, ov, best[1]['r'] - info['null_max']))
    print('fine edges:', len(edges), flush=True)
    idx = {s: i for i, s in enumerate(names)}
    Pn = {s: P[s].copy() for s in names}
    for it in range(5):
        A = np.zeros((2 * len(edges) + 2, 2 * len(names))); y = np.zeros(2 * len(edges) + 2); W = np.zeros(2 * len(edges) + 2)
        for k, (a, b, t, r, ov, exc) in enumerate(edges):
            resid = np.hypot(*(t - (Pn[b] - Pn[a]))); hub = 1.0 if resid <= 0.5 else 0.5 / resid
            w = max(r - 0.3, 0.02) ** 2 * max(ov, 0.2) * hub
            for d in range(2):
                A[2 * k + d, 2 * idx[b] + d] = 1; A[2 * k + d, 2 * idx[a] + d] = -1; y[2 * k + d] = t[d]; W[2 * k + d] = w
        r0 = names[0]
        for d in range(2): A[2 * len(edges) + d, 2 * idx[r0] + d] = 1; y[2 * len(edges) + d] = P[r0][d]; W[2 * len(edges) + d] = 1e3
        sw = np.sqrt(W)[:, None]; sol, *_ = np.linalg.lstsq(A * sw, y * sw[:, 0], rcond=None)
        Pn = {s: sol[2 * idx[s]:2 * idx[s] + 2] for s in names}
    ref = pk['ref']; off = Pn[ref].copy(); Pn = {s: p - off for s, p in Pn.items()}
    res_e = [(a, b, float(np.hypot(*(t - (Pn[b] - Pn[a])))), r) for (a, b, t, r, ov, exc) in edges]
    rows = []
    for s in names:
        e = [x for x in res_e if s in (x[0], x[1])]
        rows.append(dict(session=s, short=short(s), y_um=round(Pn[s][0], 2), x_um=round(Pn[s][1], 2), shift_vs_coarse_um=round(float(np.hypot(*(Pn[s] - (P[s] - P[ref])))), 2),
                         n_edges=len(e), rms_resid_um=round(float(np.sqrt(np.mean([x[2] ** 2 for x in e]))), 2) if e else '', max_resid_um=round(max(x[2] for x in e), 2) if e else '',
                         median_r=round(float(np.median([x[3] for x in e])), 2) if e else ''))
    with open(os.path.join(OUT, 'fov_pose_fine.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    for r in rows: print('  %s y=%7.2f x=%7.2f shift_vs_coarse %.2f  edges %2d rms %.2f max %.2f med_r %.2f' % (r['short'], r['y_um'], r['x_um'], r['shift_vs_coarse_um'], r['n_edges'], r['rms_resid_um'], r['max_resid_um'], r['median_r']))
    pickle.dump(dict(P=Pn, ref=ref, edges=edges, resid=res_e), open(os.path.join(CACHE, 'poses_fine.pkl'), 'wb'))

    # ---- 2. landmark: brightest smoothed blob in the reference session
    sm = {}
    for n in names:
        S = Sx[n]; b = np.where(S.valid, S.mean, np.median(S.mean[S.valid]))
        sm[n] = ndi.gaussian_filter(b, (1.5 / S.py, 1.5 / S.px))
    S0 = Sx[ref]; m0 = np.where(S0.valid, sm[ref], -np.inf); i0 = np.unravel_index(np.argmax(m0), m0.shape)
    lm_ref = np.array([i0[0] * S0.py, i0[1] * S0.px]) + Pn[ref]        # reference-frame position (um)
    print('landmark in reference session at (y, x) = (%.1f, %.1f) um' % tuple(lm_ref))
    lrows = []
    for n in names:
        S = Sx[n]; pred = lm_ref - Pn[n]                                  # position in this session's own um coordinates
        py_, px_ = pred[0] / S.py, pred[1] / S.px
        if not (S.yr[0] + 4 < py_ < S.yr[1] - 4 and S.xr[0] + 4 < px_ < S.xr[1] - 4):
            lrows.append(dict(session=n, short=short(n), in_fov='no', found_dy_um='', found_dx_um='', peak_contrast='')); continue
        wy, wx = int(LM_WIN / S.py), int(LM_WIN / S.px)
        y0, y1 = max(0, int(py_) - wy), min(S.shape[0], int(py_) + wy + 1); x0, x1 = max(0, int(px_) - wx), min(S.shape[1], int(px_) + wx + 1)
        sub = sm[n][y0:y1, x0:x1]; j = np.unravel_index(np.argmax(sub), sub.shape)
        found = np.array([(y0 + j[0]) * S.py, (x0 + j[1]) * S.px])
        contrast = (sub.max() - np.median(sm[n][S.valid])) / (np.std(sm[n][S.valid]) + 1e-9)
        lrows.append(dict(session=n, short=short(n), in_fov='yes', found_dy_um=round(found[0] - pred[0], 2), found_dx_um=round(found[1] - pred[1], 2), peak_contrast=round(float(contrast), 1)))
    with open(os.path.join(OUT, 'fov_landmark_check.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(lrows[0].keys())); w.writeheader(); w.writerows(lrows)
    print('landmark check (found - predicted, um; contrast in robust SDs of the smoothed mean image):')
    for r in lrows: print('  %s %s %s %s %s' % (r['short'], r['in_fov'], r['found_dy_um'], r['found_dx_um'], r['peak_contrast']))

    # ---- 3. local deformation (block-wise residual shifts) for strongly correlated pairs
    drows = []
    for (a, b, t, r, ov, exc) in edges:
        if r < 0.8 or ov < 0.5: continue
        (xa, ma), (xb, mb) = F[a]['axon'], F[b]['axon']
        # overlap region in a's grid
        tpix = np.round((Pn[b] - Pn[a]) / FG).astype(int)                 # b origin in a pixels
        ya0, xa0 = max(0, tpix[0]), max(0, tpix[1]); yb0, xb0 = max(0, -tpix[0]), max(0, -tpix[1])
        h = min(xa.shape[0] - ya0, xb.shape[0] - yb0); w_ = min(xa.shape[1] - xa0, xb.shape[1] - xb0)
        if h < 80 or w_ < 160: continue
        A_ = xa[ya0:ya0 + h, xa0:xa0 + w_]; B_ = xb[yb0:yb0 + h, xb0:xb0 + w_]; MA = ma[ya0:ya0 + h, xa0:xa0 + w_] & mb[yb0:yb0 + h, xb0:xb0 + w_]
        sh = []
        for by in range(2):
            for bx in range(3):
                ys, ye = by * h // 2, (by + 1) * h // 2; xs, xe = bx * w_ // 3, (bx + 1) * w_ // 3
                m_ = MA[ys:ye, xs:xe]
                if m_.mean() < 0.5: continue
                try:
                    rb, Nb = masked_ncc(A_[ys:ye, xs:xe], m_, B_[ys:ye, xs:xe], m_)
                    hb, wb = m_.shape; ii, jj = np.meshgrid(np.arange(rb.shape[0]), np.arange(rb.shape[1]), indexing='ij')
                    dy_, dx_ = (ii - (hb - 1)) * FG, (jj - (wb - 1)) * FG               # b shifted by (dy, dx) relative to its nominal place
                    ok = (np.hypot(dy_, dx_) <= 3.0) & (Nb >= 0.5 * m_.size)
                    if not ok.any(): continue
                    q = np.where(ok, rb, -1.0); i = np.unravel_index(np.argmax(q), q.shape)
                    if q[i] < 0.3: continue
                    sh.append(np.array([dy_[i], dx_[i]]))
                except Exception:
                    pass
        if len(sh) >= 3:
            sh = np.array(sh); rel = sh - np.median(sh, 0)
            drows.append(dict(a=short(a), b=short(b), r=round(r, 2), n_blocks=len(sh), block_spread_um=round(float(np.sqrt((rel ** 2).sum(1).mean())), 2), max_block_dev_um=round(float(np.hypot(*rel.T).max()), 2)))
    if drows:
        with open(os.path.join(OUT, 'fov_local_deformation.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(drows[0].keys())); w.writeheader(); w.writerows(drows)
        print('local deformation (block shifts relative to the pair\'s median block shift): median RMS %.2f um, median max %.2f um over %d pairs' % (
            np.median([d['block_spread_um'] for d in drows]), np.median([d['max_block_dev_um'] for d in drows]), len(drows)))


if __name__ == '__main__':
    main()
