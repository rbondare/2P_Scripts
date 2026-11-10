"""Step 4: match the user's accepted suite2p ROIs (iscell == 1) across sessions in the common um frame (poses from cache/poses_fine.pkl).

ROI masks are rasterised on the 0.25 um common grid (reference frame), optionally dilated by DIL_UM (thin axon ROIs are ~1 um wide, residual misalignment 0.3-0.6 um,
see fov_local_deformation.csv). Pair score = IoU of the (dilated) masks. Per session pair: ROIs that lie FULLY inside the other session's valid FOV (margin 0.5 um) are
'evaluable'; Hungarian assignment with IoU >= IOU_MIN; match rate = matched / evaluable (separately for each direction and for the smaller set).
Controls per pair: CHANCE = the same computation after displacing session b (ROIs and its FOV together) by random shifts of 6-25 um (N_NULL draws);
CEILING = same-day same-FOV pairs (identical stage position), where only detection / curation differs.
Outputs: roi_match_pairs.csv, master_axons.csv (+ cache/roi_match.pkl with the per-session cell sets and the cluster table).
Run: miniforge3\\envs\\suite2p_venv\\python.exe fov_roi_match.py [--dil 0.75] [--iou 0.2]
"""
import os, sys, csv, pickle, itertools, argparse, numpy as np
from scipy import sparse
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment
from fov_common import *

FG = 0.25
N_NULL = 40
MARGIN_UM = 0.5


def rasterise(S, P, accepted_only=True, which=None):
    """cell sets of the ROIs on the reference canvas (integer cell indices, canvas origin handled by the caller): list of (i_cells, j_cells) in ABSOLUTE grid coordinates
    (cell index = round(um / FG)) so that no canvas bounds are needed"""
    ids = np.where(S.accepted)[0] if which is None else np.asarray(which)
    out = []
    for i in ids:
        yy, xx, lam = S.roi_pix(i)
        y0, y1, x0, x1 = yy.min(), yy.max(), xx.min(), xx.max()
        loc = np.zeros((y1 - y0 + 1, x1 - x0 + 1), bool); loc[yy - y0, xx - x0] = True
        # cells whose centres fall into the ROI's native pixels
        uy0, uy1 = P[0] + (y0 - 0.5) * S.py, P[0] + (y1 + 0.5) * S.py
        ux0, ux1 = P[1] + (x0 - 0.5) * S.px, P[1] + (x1 + 0.5) * S.px
        ci = np.arange(int(np.ceil(uy0 / FG)), int(np.floor(uy1 / FG)) + 1); cj = np.arange(int(np.ceil(ux0 / FG)), int(np.floor(ux1 / FG)) + 1)
        if len(ci) == 0 or len(cj) == 0:     # ROI smaller than one cell: take the cell of its centroid
            ci = np.array([int(round((P[0] + yy.mean() * S.py) / FG))]); cj = np.array([int(round((P[1] + xx.mean() * S.px) / FG))])
            out.append((int(i), ci, cj, np.ones((1, 1), bool))); continue
        py_ = np.round(((ci * FG) - P[0]) / S.py).astype(int) - y0; px_ = np.round(((cj * FG) - P[1]) / S.px).astype(int) - x0
        ok_y = (py_ >= 0) & (py_ < loc.shape[0]); ok_x = (px_ >= 0) & (px_ < loc.shape[1])
        m = np.zeros((len(ci), len(cj)), bool)
        m[np.ix_(ok_y, ok_x)] = loc[np.ix_(py_[ok_y], px_[ok_x])]
        if not m.any(): ci = np.array([int(round((P[0] + yy.mean() * S.py) / FG))]); cj = np.array([int(round((P[1] + xx.mean() * S.px) / FG))]); m = np.ones((1, 1), bool)
        out.append((int(i), ci, cj, m))
    return out


def dilate(masks, d_cells):
    if d_cells <= 0: return masks
    st = ndi.generate_binary_structure(2, 1); out = []
    for (i, ci, cj, m) in masks:
        mp = np.pad(m, d_cells); md = ndi.binary_dilation(mp, structure=ndi.iterate_structure(st, d_cells))
        ci2 = np.arange(ci[0] - d_cells, ci[0] - d_cells + md.shape[0]); cj2 = np.arange(cj[0] - d_cells, cj[0] - d_cells + md.shape[1])
        out.append((i, ci2, cj2, md))
    return out


class Canvas:
    def __init__(self, ymin, ymax, xmin, xmax):
        self.oy = int(np.floor(ymin / FG)) - 40; self.ox = int(np.floor(xmin / FG)) - 40
        self.H = int(np.ceil(ymax / FG)) + 40 - self.oy; self.W = int(np.ceil(xmax / FG)) + 40 - self.ox

    def sparse(self, masks, shift=(0, 0)):
        rows, cols = [], []
        for k, (i, ci, cj, m) in enumerate(masks):
            a, b = np.nonzero(m)
            yy = ci[a] + shift[0] - self.oy; xx = cj[b] + shift[1] - self.ox
            ok = (yy >= 0) & (yy < self.H) & (xx >= 0) & (xx < self.W)
            rows.append(np.full(ok.sum(), k)); cols.append(yy[ok] * self.W + xx[ok])
        r = np.concatenate(rows) if rows else np.array([], int); c = np.concatenate(cols) if cols else np.array([], int)
        return sparse.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(len(masks), self.H * self.W))


def fov_rect(S, P):
    """valid-FOV rectangle (um, reference frame) shrunk by MARGIN_UM: (ymin, ymax, xmin, xmax)"""
    return (P[0] + (S.yr[0] - 0.5) * S.py + MARGIN_UM, P[0] + (S.yr[1] - 0.5) * S.py - MARGIN_UM, P[1] + (S.xr[0] - 0.5) * S.px + MARGIN_UM, P[1] + (S.xr[1] - 0.5) * S.px - MARGIN_UM)


def in_rect(masks, rect, shift_um=(0, 0)):
    out = np.zeros(len(masks), bool)
    for k, (i, ci, cj, m) in enumerate(masks):
        a, b = np.nonzero(m)
        y = (ci[a] * FG + shift_um[0]); x = (cj[b] * FG + shift_um[1])
        out[k] = (y.min() >= rect[0]) and (y.max() <= rect[1]) and (x.min() >= rect[2]) and (x.max() <= rect[3])
    return out


def match(Ma, Mb, ia, ib, iou_min):
    """Ma, Mb sparse row-binary matrices restricted to evaluable ROIs; returns list of (row_a, row_b, iou)"""
    if Ma.shape[0] == 0 or Mb.shape[0] == 0: return []
    inter = (Ma @ Mb.T).toarray(); na = np.asarray(Ma.sum(1)).ravel(); nb = np.asarray(Mb.sum(1)).ravel()
    iou = inter / (na[:, None] + nb[None, :] - inter + 1e-9)
    if (iou >= iou_min).sum() == 0: return []
    cost = np.where(iou >= iou_min, -iou, 1e3)
    r, c = linear_sum_assignment(cost)
    return [(a, b, float(iou[a, b])) for a, b in zip(r, c) if iou[a, b] >= iou_min]


def pair_stats(sa, sb, shift_cells=(0, 0), shift_um=(0, 0)):
    """evaluable ROI index sets of both sessions for (possibly displaced) session b; returns indices into the masks lists"""
    ea = np.where(in_rect(sa['masks'], sb['rect_shift'](shift_um)))[0]
    eb = np.where(in_rect(sb['masks'], sa['rect'], shift_um=shift_um))[0]
    return ea, eb


def build(P, dil_um, accepted_only=True, names=None):
    inv = read_inventory(); names = names or [s for s in SESSIONS if s in P and s not in EXCLUDE_ROI_SESSIONS]
    sess = {}
    for n in names:
        S = Sess(n, inv)
        masks0 = rasterise(S, P[n], accepted_only)
        # drop ROIs on the bright blob landmark (not an axon)
        keep = []
        for (i, ci, cj, m) in masks0:
            a, b = np.nonzero(m); keep.append(np.hypot(ci[a].mean() * FG - LM[0], cj[b].mean() * FG - LM[1]) > BLOB_R)
        masks0 = [x for x, k in zip(masks0, keep) if k]
        masks = dilate(masks0, int(round(dil_um / FG)))
        rect = fov_rect(S, P[n])
        sess[n] = dict(S=S, masks=masks, masks0=masks0, rect=rect, ids=np.array([m[0] for m in masks]))
        sess[n]['rect_shift'] = (lambda sh, r=rect: (r[0] + sh[0], r[1] + sh[0], r[2] + sh[1], r[3] + sh[1]))
    return sess


def canvas_for(sess):
    ys = []; xs = []
    for n, s in sess.items():
        for (i, ci, cj, m) in s['masks'][:50]: pass
        ys += [s['rect'][0], s['rect'][1]]; xs += [s['rect'][2], s['rect'][3]]
    return Canvas(min(ys) - 30, max(ys) + 30, min(xs) - 30, max(xs) + 30)


def run_pair(a, b, sess, cv, iou_min, shift_um=(0.0, 0.0)):
    sa, sb = sess[a], sess[b]
    shift_cells = (int(round(shift_um[0] / FG)), int(round(shift_um[1] / FG)))
    shb_rect = sb['rect_shift'](shift_um)
    ea = np.where(in_rect(sa['masks'], shb_rect))[0]
    # b's ROIs displaced by shift: in view of a if inside a's rect after displacement
    eb = np.where(in_rect(sb['masks'], sa['rect'], shift_um=shift_um))[0]
    Ma = cv.sparse([sa['masks'][k] for k in ea]); Mb = cv.sparse([sb['masks'][k] for k in eb], shift=shift_cells)
    m = match(Ma, Mb, ea, eb, iou_min)
    return ea, eb, [(int(ea[x]), int(eb[y]), v) for x, y, v in m]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--dil', type=float, default=0.75); ap.add_argument('--iou', type=float, default=0.2)
    ap.add_argument('--tag', default=''); args = ap.parse_args()
    pk = pickle.load(open(os.path.join(CACHE, 'poses_fine.pkl'), 'rb')); P = pk['P']
    sess = build(P, args.dil); cv = canvas_for(sess); names = list(sess)
    inv = read_inventory(); day = {n: inv[n]['tif_first_mtime'][5:10] for n in names}
    same_day = {}
    for g in SAME_DAY.values():
        for a, b in itertools.combinations(g, 2): same_day[(a, b)] = same_day[(b, a)] = True
    rng = np.random.default_rng(1)
    rows = []; allm = {}
    for a, b in itertools.combinations(names, 2):
        ea, eb, mt = run_pair(a, b, sess, cv, args.iou)
        allm[(a, b)] = mt
        # chance: displace b by random shifts 6-25 um
        ch = []
        for _ in range(N_NULL):
            ang = rng.uniform(0, 2 * np.pi); rad = rng.uniform(6, 25); sh = (rad * np.sin(ang), rad * np.cos(ang))
            ea2, eb2, mt2 = run_pair(a, b, sess, cv, args.iou, shift_um=sh)
            n_min = min(len(ea2), len(eb2))
            ch.append(len(mt2) / n_min if n_min else np.nan)
        n_min = min(len(ea), len(eb)); fr = len(mt) / n_min if n_min else np.nan
        rows.append(dict(a=a, b=b, a_short=short(a), b_short=short(b), same_day=int((a, b) in same_day), days_apart=abs((np.datetime64('2026-' + day[a]) - np.datetime64('2026-' + day[b])).astype(int)),
                         n_a_acc=len(sess[a]['masks']), n_b_acc=len(sess[b]['masks']), n_a_inview=len(ea), n_b_inview=len(eb), n_matched=len(mt), frac_of_smaller=round(fr, 3) if n_min else '',
                         frac_of_a=round(len(mt) / len(ea), 3) if len(ea) else '', frac_of_b=round(len(mt) / len(eb), 3) if len(eb) else '',
                         chance_mean=round(float(np.nanmean(ch)), 3), chance_sd=round(float(np.nanstd(ch)), 3), excess=round(fr - float(np.nanmean(ch)), 3) if n_min else ''))
        r = rows[-1]; print('%s %s same_day=%d inview %3d/%3d matched %3d frac %.2f chance %.2f+-%.2f' % (r['a_short'], r['b_short'], r['same_day'], len(ea), len(eb), len(mt), fr if n_min else np.nan, r['chance_mean'], r['chance_sd']), flush=True)
    tag = args.tag
    with open(os.path.join(OUT, 'roi_match_pairs%s.csv' % tag), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    pickle.dump(dict(matches=allm, args=vars(args)), open(os.path.join(CACHE, 'roi_match%s.pkl' % tag), 'wb'))
    sd = [r for r in rows if r['same_day'] and r['frac_of_smaller'] != '']; cd = [r for r in rows if not r['same_day'] and r['frac_of_smaller'] != '' and r['n_a_inview'] >= 5 and r['n_b_inview'] >= 5]
    print('\nSAME-DAY pairs (%d): median frac %.2f (chance %.2f)' % (len(sd), np.median([r['frac_of_smaller'] for r in sd]), np.median([r['chance_mean'] for r in sd])))
    print('CROSS-DAY pairs (%d): median frac %.2f (chance %.2f)' % (len(cd), np.median([r['frac_of_smaller'] for r in cd]), np.median([r['chance_mean'] for r in cd])))


if __name__ == '__main__':
    main()
