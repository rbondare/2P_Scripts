"""Registration QC figures from cache/poses.pkl:
 Q1_fov_footprints   outline of every session's valid FOV in the reference frame (um), coloured by day
 Q2_registered_gallery  every session's mean image on the common frame (same crop), with its accepted ROIs as outlines
 Q3_pair_overlays    red / green overlay of a few informative pairs (within day, adjacent days, the weakest placed sessions)
Run after fov_register.py + fov_pose.py:  miniforge3\\envs\\suite2p_venv\\python.exe fov_qc.py
"""
import os, pickle, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from fov_common import *

DAYCOL = {'1021': '#1b9e77', '1028': '#66a61e', '1022': '#7570b3', '1023': '#d95f02', '1027': '#e7298a', '1029': '#e6ab02', '1030': '#a6761d', '1104': '#1f78b4', '1105': '#33a02c', '1106': '#e31a1c'}


def day_of(inv, nm): return inv[nm]['tif_first_mtime'][5:7] + inv[nm]['tif_first_mtime'][8:10]


def main(grid=0.5):
    P = pickle.load(open(os.path.join(CACHE, 'poses.pkl'), 'rb'))['P']; inv = read_inventory()
    names = [s for s in SESSIONS if s in P]
    S = {n: Sess(n, inv) for n in names}
    # common canvas
    ext = []
    for n in names:
        s = S[n]; y0, y1 = s.yr; x0, x1 = s.xr
        ext.append((P[n][0] + y0 * s.py, P[n][0] + (y1 - 1) * s.py, P[n][1] + x0 * s.px, P[n][1] + (x1 - 1) * s.px))
    ext = np.array(ext); ymin, ymax, xmin, xmax = ext[:, 0].min() - 3, ext[:, 1].max() + 3, ext[:, 2].min() - 3, ext[:, 3].max() + 3
    # Q1 footprints
    fig, ax = plt.subplots(figsize=(7, 7 * (ymax - ymin) / (xmax - xmin) + 0.6))
    for n, e in zip(names, ext):
        d = day_of(inv, n)
        ax.add_patch(Rectangle((e[2], e[0]), e[3] - e[2], e[1] - e[0], fill=False, ec=DAYCOL.get(d, 'k'), lw=1.1, alpha=0.9))
        ax.text(e[2] + 1, e[0] + 3, short(n), fontsize=5, color=DAYCOL.get(d, 'k'))
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymax, ymin); ax.set_aspect('equal'); ax.set_xlabel('x (um, reference = 11-05 1520)'); ax.set_ylabel('y (um)')
    ax.set_title('Valid field of view of every session in the common frame')
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Q1_fov_footprints.png'), dpi=200); plt.close(fig)
    # Q2 gallery on the common canvas
    ncol = 4; nrow = int(np.ceil(len(names) / ncol))
    H = int(np.ceil((ymax - ymin) / grid)); W = int(np.ceil((xmax - xmin) / grid))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.2 * nrow * H / W + 0.3 * nrow))
    axs = np.atleast_2d(axs)
    for k, n in enumerate(names):
        a = axs[k // ncol, k % ncol]; s = S[n]
        img = np.full((H, W), np.nan)
        g = to_grid(np.where(s.valid, s.mean, np.nan), s.py, s.px, order=1, grid=grid, cval=np.nan)
        gm = to_grid(s.valid.astype(float), s.py, s.px, order=1, grid=grid) > 0.99
        oy = int(round((P[n][0] - ymin) / grid)); ox = int(round((P[n][1] - xmin) / grid))
        hh, ww = g.shape; y1, x1 = min(H, oy + hh), min(W, ox + ww)
        sub = np.where(gm, g, np.nan)[:y1 - oy, :x1 - ox]; img[oy:y1, ox:x1] = sub
        lo, hi = np.nanpercentile(img, [1, 99.7]); a.imshow(img, cmap='gray', vmin=lo, vmax=hi, extent=(xmin, xmax, ymax, ymin), interpolation='nearest')
        acc = np.where(s.accepted)[0]
        for i in acc:
            yy, xx, lam = s.roi_pix(i)
            a.plot(P[n][1] + xx.mean() * s.px, P[n][0] + yy.mean() * s.py, '.', ms=1.6, color='#ff3b30', alpha=0.8)
        a.set_title('%s  (%d acc. ROIs)' % (short(n), s.accepted.sum()), fontsize=7); a.set_xticks([]); a.set_yticks([])
    for k in range(len(names), nrow * ncol): axs[k // ncol, k % ncol].axis('off')
    fig.suptitle('Mean images on the common frame (same crop, um); red dots = accepted ROI centroids', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Q2_registered_gallery.png'), dpi=170); plt.close(fig)
    # Q3 pair overlays
    pairs = [('AnimalRB24_261105_1423', 'AnimalRB24_261105_1520'), ('AnimalRB24_261105_1423', 'AnimalRB24_261106_1317'), ('AnimalRB24_261106_1317', 'AnimalRB24_261106_1604'),
             ('AnimalRB24_261105_1423', 'AnimalRB24_261023_1538'), ('AnimalRB24_261105_1423', 'AnimalRB24_261027_1456'), ('AnimalRB24_261105_1423', 'AnimalRB24_261030_1744'),
             ('AnimalRB24_261105_1423', 'AnimalRB24_261021_1630'), ('AnimalRB24_261105_1423', 'AnimalRB24_261022_2026'), ('AnimalRB24_261105_1423', 'AnimalRB24_261021_1826')]
    fig, axs = plt.subplots(3, 3, figsize=(13, 9)); axs = axs.ravel()
    for a, (A, B) in zip(axs, pairs):
        if A not in P or B not in P: a.axis('off'); continue
        canvas = []
        for n in (A, B):
            s = S[n]; g = to_grid(np.where(s.valid, s.mean, np.nan), s.py, s.px, order=1, grid=grid, cval=np.nan)
            gm = to_grid(s.valid.astype(float), s.py, s.px, order=1, grid=grid) > 0.99
            f = np.full((H, W), np.nan); oy = int(round((P[n][0] - ymin) / grid)); ox = int(round((P[n][1] - xmin) / grid))
            y1, x1 = min(H, oy + g.shape[0]), min(W, ox + g.shape[1]); f[oy:y1, ox:x1] = np.where(gm, g, np.nan)[:y1 - oy, :x1 - ox]
            # axon-scale band-pass for display
            fz = np.where(np.isfinite(f), f, np.nanmedian(f)); bp = ndi.gaussian_filter(fz, 1.5) - ndi.gaussian_filter(fz, 10); bp[~np.isfinite(f)] = np.nan
            lo, hi = np.nanpercentile(bp, [2, 99.5]); canvas.append(np.clip((bp - lo) / (hi - lo), 0, 1))
        rgb = np.zeros((H, W, 3)); rgb[..., 0] = np.nan_to_num(canvas[0]); rgb[..., 1] = np.nan_to_num(canvas[1]); rgb[..., 2] = 0.5 * np.nan_to_num(canvas[1])
        a.imshow(rgb, extent=(xmin, xmax, ymax, ymin)); a.set_title('red %s / green %s' % (short(A), short(B)), fontsize=8); a.set_xticks([]); a.set_yticks([])
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Q3_pair_overlays.png'), dpi=150); plt.close(fig)
    print('written Q1-Q3 in', OUT)


if __name__ == '__main__':
    main()
