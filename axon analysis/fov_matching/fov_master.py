"""Step 4b + 5: master axons (clusters of matched accepted ROIs across sessions) and a detector-independent look at where each master axon is VISIBLE in every session.

Master axons: constrained union-find over all pairwise matches (cache/roi_match.pkl, IoU >= 0.2 on masks dilated by 0.75 um), strongest first; two clusters are only merged if
they contain no ROI of the same session (an axon is one ROI per session), which stops chains of near-misses from merging distinct axons.
Footprint of a master axon = union of its member ROI masks in the reference frame (0.25 um grid).

For every (master axon, session) with the footprint FULLY inside that session's valid FOV (margin 0.5 um) the following is measured on the session's own images resampled to
the reference grid: mean robust-z of (a) the axon-scale DoG of the mean image (tonic structure), (b) suite2p Vcorr (local correlation = activity structure), (c) max_proj minus
mean (peak activity). The reference distribution ('background') for each session and image = the same statistic for master footprints displaced by random 5-20 um shifts that stay
inside the valid FOV. Class per (axon, session):
  accepted                         an accepted ROI of that session is matched to the axon
  rejected ROI present             no accepted ROI, but a suite2p ROI (iscell == 0) with IoU >= 0.2 (dilated)
  structure, no ROI                nothing detected, but the DoG / Vcorr / max-proj contrast is above the 95th percentile of the background
  nothing visible                  contrast at background level ('silent' or absent or out of the focal plane - not separable without a reference channel)
  out of FOV                       footprint not fully inside that session's valid field
Outputs: master_axons.csv, master_presence.csv (long), cache/master.pkl
Run: miniforge3\\envs\\suite2p_venv\\python.exe fov_master.py
"""
import os, sys, csv, pickle, itertools, numpy as np
from scipy import ndimage as ndi
from fov_common import *
from fov_roi_match import build, canvas_for, rasterise, dilate, in_rect, FG
from fov_register import FEATS

N_BG = 400
ACC_IOU = 0.2


def canvas_images(S, P, cv):
    """DoG(axon scale), Vcorr, max_proj - mean, each robust-z-scored inside the valid mask, resampled onto the reference canvas (NaN outside the valid FOV)"""
    m = S.mean
    s1, s2 = FEATS['axon']
    d = dog(m, S.valid, s1, s2, S.py, S.px)
    v = S.vcorr
    mp = (S.max_proj - m) if S.max_proj is not None else None
    out = {}
    ys = (cv.oy + np.arange(cv.H)) * FG - P[0]; xs = (cv.ox + np.arange(cv.W)) * FG - P[1]
    yy, xx = np.meshgrid(ys / S.py, xs / S.px, indexing='ij')
    inside = (yy >= S.yr[0]) & (yy <= S.yr[1] - 1) & (xx >= S.xr[0]) & (xx <= S.xr[1] - 1)
    for k, im in (('dog', d), ('vcorr', v), ('maxp', mp)):
        if im is None: continue
        im = np.array(im, float); ok = S.valid & np.isfinite(im)
        x = im[ok]; med = np.median(x); mad = np.median(np.abs(x - med)) * 1.4826 + 1e-12
        z = np.where(ok, (im - med) / mad, 0.0)
        g = ndi.map_coordinates(z, [yy, xx], order=1, mode='nearest', prefilter=False)
        out[k] = np.where(inside, g, np.nan)
    out['inside'] = inside
    return out


def cluster(sess, matches):
    """constrained union-find; returns list of clusters: dict(session -> row index into sess[session]['masks'])"""
    parent = {}; members = {}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for n, s in sess.items():
        for k in range(len(s['masks'])): key = (n, k); parent[key] = key; members[key] = {n: k}
    edges = []
    for (a, b), mt in matches.items():
        for (ia, ib, v) in mt: edges.append((v, (a, ia), (b, ib)))
    edges.sort(key=lambda e: -e[0]); refused = 0
    for v, x, y in edges:
        rx, ry = find(x), find(y)
        if rx == ry: continue
        if set(members[rx]) & set(members[ry]): refused += 1; continue
        parent[ry] = rx; members[rx].update(members[ry]); del members[ry]
    print('clusters: %d (from %d ROIs); merges refused because of a same-session conflict: %d' % (len(members), sum(len(s['masks']) for s in sess.values()), refused))
    return list(members.values())


def footprint(sess, mem):
    cells = set()
    for n, k in mem.items():
        i, ci, cj, m = sess[n]['masks0'][k]
        a, b = np.nonzero(m); cells.update(zip((ci[a]).tolist(), (cj[b]).tolist()))
    c = np.array(sorted(cells)); return c[:, 0], c[:, 1]


def main():
    pk = pickle.load(open(os.path.join(CACHE, 'poses_fine.pkl'), 'rb')); P = pk['P']
    mt = pickle.load(open(os.path.join(CACHE, 'roi_match.pkl'), 'rb')); matches = mt['matches']; dil = mt['args']['dil']
    sess = build(P, dil); cv = canvas_for(sess); names = list(sess)
    # clusters on the matches computed with the same dilation / threshold
    clusters = cluster(sess, matches)
    clusters.sort(key=lambda mem: (-len(mem), min(mem)))
    inv = read_inventory()
    # images on canvas
    print('resampling images to the canvas ...', flush=True)
    IM = {n: canvas_images(sess[n]['S'], P[n], cv) for n in names}
    # all suite2p ROIs (accepted + rejected) as dilated masks for the 'rejected ROI present' class
    allm = {}
    for n in names:
        S = sess[n]['S']; rej = np.where(~S.accepted)[0]
        allm[n] = dilate(rasterise(S, P[n], which=rej), int(round(dil / FG))) if len(rej) else []
    # footprints
    foot = [footprint(sess, mem) for mem in clusters]
    rect = {n: sess[n]['rect'] for n in names}
    def inside(ci, cj, n, shift=(0, 0)):
        r = rect[n]; y = ci * FG + shift[0]; x = cj * FG + shift[1]
        return y.min() >= r[0] and y.max() <= r[1] and x.min() >= r[2] and x.max() <= r[3]
    def stat(n, key, ci, cj, shift=(0, 0)):
        y = ci + int(round(shift[0] / FG)) - cv.oy; x = cj + int(round(shift[1] / FG)) - cv.ox
        v = IM[n][key][y, x]; return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
    rng = np.random.default_rng(5)
    # background distributions per session and image
    print('background distributions ...', flush=True)
    BG = {n: {k: [] for k in ('dog', 'vcorr', 'maxp')} for n in names}
    for n in names:
        tries = 0
        while len(BG[n]['dog']) < N_BG and tries < 20 * N_BG:
            tries += 1
            f = foot[rng.integers(len(foot))]; ang = rng.uniform(0, 2 * np.pi); rad = rng.uniform(5, 20); sh = (rad * np.sin(ang), rad * np.cos(ang))
            if not inside(f[0], f[1], n, sh): continue
            for k in BG[n]:
                if k in IM[n]: BG[n][k].append(stat(n, k, f[0], f[1], sh))
        for k in BG[n]: BG[n][k] = np.array(BG[n][k])
    thr = {n: {k: (np.nanpercentile(v, 95) if len(v) > 20 else np.nan) for k, v in BG[n].items()} for n in names}
    # per (master, session)
    rows = []; mrows = []
    for q, (mem, (ci, cj)) in enumerate(zip(clusters, foot)):
        cy, cx = ci.mean() * FG, cj.mean() * FG
        n_in = 0; n_acc = len(mem); cl = {}
        for n in names:
            if not inside(ci, cj, n):
                rows.append(dict(master=q, session=n, short=short(n), cls='out of FOV')); continue
            n_in += 1
            z = {k: stat(n, k, ci, cj) for k in ('dog', 'vcorr', 'maxp') if k in IM[n]}
            if n in mem:
                cls = 'accepted'; roi = int(sess[n]['ids'][mem[n]]); iou = np.nan
            else:
                roi = -1; iou = 0.0
                # rejected suite2p ROI overlapping the footprint
                if allm[n]:
                    inter = []
                    fset = set(zip(ci.tolist(), cj.tolist()))
                    best = (0.0, -1)
                    for (i2, ci2, cj2, m2) in allm[n]:
                        a, b = np.nonzero(m2); cells = set(zip((ci2[a]).tolist(), (cj2[b]).tolist()))
                        if not cells: continue
                        u = len(cells & fset)
                        if u == 0: continue
                        v = u / (len(cells) + len(fset) - u)
                        if v > best[0]: best = (v, i2)
                    iou, roi = best
                if iou >= ACC_IOU: cls = 'rejected ROI present'
                else:
                    roi = -1
                    above = [z[k] > thr[n][k] for k in z if np.isfinite(thr[n][k])]
                    cls = 'structure, no ROI' if any(above) else 'nothing visible'
            rows.append(dict(master=q, session=n, short=short(n), cls=cls, roi_id=roi, iou_rejected=iou, **{'z_' + k: round(v, 2) for k, v in z.items()}))
            cl[n] = cls
        mrows.append(dict(master=q, n_sessions_accepted=n_acc, n_sessions_in_view=n_in, y_um=round(cy, 2), x_um=round(cx, 2), n_cells=len(ci), area_um2=round(len(ci) * FG * FG, 2),
                          sessions=';'.join(short(n) for n in names if n in mem)))
    with open(os.path.join(OUT, 'master_axons.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(mrows[0].keys())); w.writeheader(); w.writerows(mrows)
    keys = ['master', 'session', 'short', 'cls', 'roi_id', 'iou_rejected', 'z_dog', 'z_vcorr', 'z_maxp']
    with open(os.path.join(OUT, 'master_presence.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows: w.writerow({k: r.get(k, '') for k in keys})
    pickle.dump(dict(clusters=[{n: int(sess[n]['ids'][k]) for n, k in mem.items()} for mem in clusters], rows=rows, mrows=mrows, thr=thr, names=names, BG={n: {k: v for k, v in d.items()} for n, d in BG.items()}),
                open(os.path.join(CACHE, 'master.pkl'), 'wb'))
    # per-session summary
    print('\nper session: master axons in view / accepted / rejected ROI present / structure, no ROI / nothing visible')
    for n in names:
        r = [x for x in rows if x['session'] == n and x['cls'] != 'out of FOV']
        c = {k: sum(x['cls'] == k for x in r) for k in ('accepted', 'rejected ROI present', 'structure, no ROI', 'nothing visible')}
        print('  %s in view %4d  acc %4d  rej %4d  struct %4d  nothing %4d' % (short(n), len(r), c['accepted'], c['rejected ROI present'], c['structure, no ROI'], c['nothing visible']))


if __name__ == '__main__':
    main()
