"""Figures of the RB24 cross-session FOV / axon-location analysis -> plots/fov_matching/F1..F6 (PNG 200 dpi + PDF).
 F1 registration (footprints, landmark residual, alignment evidence)      F2 registered mean-image gallery with accepted ROIs
 F3 ROI matching vs chance / same-day ceiling                             F4 presence matrix of the confirmed master axons (class per session)
 F5 what differs between sessions / days                                  F6 example axons across all sessions (band-pass image, fixed colour scale)
Run after fov_analysis.py:  miniforge3\\envs\\suite2p_venv\\python.exe fov_figures.py [F1 F2 ...]
"""
import os, sys, csv, pickle, collections, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap
from scipy import ndimage as ndi
from scipy import stats
from fov_common import *
from fov_roi_match import Canvas, rasterise, fov_rect, FG
from fov_master import canvas_images

plt.rcParams.update({'font.size': 7, 'axes.titlesize': 8, 'axes.titleweight': 'bold', 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.dpi': 200})
INV = read_inventory()
ORDER = sorted(SESSIONS, key=lambda s: INV[s]['tif_first_mtime'])
RORDER = [s for s in ORDER if s not in EXCLUDE_ROI_SESSIONS]       # sessions used in the ROI-based analyses
DAY = {s: INV[s]['tif_first_mtime'][5:10] for s in SESSIONS}
DAYS = sorted(set(DAY.values()))
DCOL = dict(zip(DAYS, ['#1b9e77', '#7570b3', '#d95f02', '#e7298a', '#e6ab02', '#a6761d', '#66a61e', '#1f78b4', '#33a02c', '#e31a1c']))
CLS = ['accepted', 'rejected ROI present', 'structure, no ROI', 'nothing visible', 'out of FOV']
CCOL = {'accepted': '#1F5FBF', 'rejected ROI present': '#D98A00', 'structure, no ROI': '#7FC8C0', 'nothing visible': '#D1D5DB', 'out of FOV': '#FFFFFF'}
LM = (24.0, 97.3)          # landmark (bright blob) in the reference frame, from fov_refine.py


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=200, bbox_inches='tight'); fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight'); plt.close(fig); print('written', name, flush=True)


def load_all():
    P = pickle.load(open(os.path.join(CACHE, 'poses_fine.pkl'), 'rb'))
    M = pickle.load(open(os.path.join(CACHE, 'master.pkl'), 'rb')); A = pickle.load(open(os.path.join(CACHE, 'analysis.pkl'), 'rb'))
    return P, M, A


def make_canvas(P):
    ys, xs = [], []
    for n in SESSIONS:
        S = Sess(n, INV, rois=False); r = fov_rect(S, P['P'][n]); ys += [r[0], r[1]]; xs += [r[2], r[3]]
    return Canvas(min(ys) - 30, max(ys) + 30, min(xs) - 30, max(xs) + 30)


def mean_on_canvas(S, Pn, cv, y0, y1, x0, x1):
    """raw mean image on the reference grid inside the crop [y0,y1]x[x0,x1] um (NaN outside the valid FOV)"""
    ys = np.arange(y0, y1, FG) - Pn[0]; xs = np.arange(x0, x1, FG) - Pn[1]
    yy, xx = np.meshgrid(ys / S.py, xs / S.px, indexing='ij')
    b = np.where(S.valid, S.mean, np.median(S.mean[S.valid]))
    g = ndi.map_coordinates(b, [yy, xx], order=1, mode='nearest')
    inside = (yy >= S.yr[0]) & (yy <= S.yr[1] - 1) & (xx >= S.xr[0]) & (xx <= S.xr[1] - 1)
    return np.where(inside, g, np.nan)


def outline(masks, shape, oy, ox):
    ov = np.zeros(shape, bool)
    for (i, ci, cj, m) in masks:
        e = m & ~ndi.binary_erosion(m, iterations=1) if m.sum() > 4 else m
        a, b = np.nonzero(e); yy = ci[a] - oy; xx = cj[b] - ox; ok = (yy >= 0) & (yy < shape[0]) & (xx >= 0) & (xx < shape[1]); ov[yy[ok], xx[ok]] = True
    return ov


# ---------------------------------------------------------------- F1
def fig1(P, M, A):
    pf = {r['session']: r for r in csv.DictReader(open(os.path.join(OUT, 'fov_pose_fine.csv')))}
    lm = {r['session']: r for r in csv.DictReader(open(os.path.join(OUT, 'fov_landmark_check.csv')))}
    fig = plt.figure(figsize=(12, 4.6)); gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], wspace=0.35)
    ax = fig.add_subplot(gs[0])
    for n in ORDER:
        S = Sess(n, INV, rois=False); r = fov_rect(S, P['P'][n]); d = DAY[n]
        ax.add_patch(Rectangle((r[2], r[0]), r[3] - r[2], r[1] - r[0], fill=False, ec=DCOL[d], lw=1.0, alpha=0.9))
    for d in DAYS: ax.plot([], [], color=DCOL[d], label=d)
    ax.plot(LM[1], LM[0], 'k*', ms=8, label='landmark (blob)')
    ax.set_xlim(-18, 195); ax.set_ylim(150, -65); ax.set_aspect('equal'); ax.set_xlabel('x (um, reference = 11-05 1520)'); ax.set_ylabel('y (um)')
    ax.legend(ncol=2, fontsize=6, loc='lower left', frameon=False); ax.set_title('a  Valid field of every session in the common frame')
    ax = fig.add_subplot(gs[1]); xs = np.arange(len(ORDER))
    for k, n in enumerate(ORDER):
        r = lm[n]
        if r['in_fov'] == 'yes':
            v = np.hypot(float(r['found_dy_um']), float(r['found_dx_um'])); ax.bar(k, v, color=DCOL[DAY[n]])
        else: ax.text(k, 0.1, 'not in\nFOV', ha='center', fontsize=5, rotation=90)
    ax.axhline(0.5, color='k', lw=0.5, ls=':'); ax.set_xticks(xs); ax.set_xticklabels([short(n) for n in ORDER], rotation=90, fontsize=6)
    ax.set_ylabel('landmark position error after alignment (um)'); ax.set_title('b  Independent landmark check')
    ax = fig.add_subplot(gs[2])
    for k, n in enumerate(ORDER):
        r = pf[n]; ax.scatter(k, float(r['median_r']), s=10 + 8 * float(r['n_edges']), color=DCOL[DAY[n]], alpha=0.85)
    ax.set_xticks(xs); ax.set_xticklabels([short(n) for n in ORDER], rotation=90, fontsize=6); ax.set_ylabel('median image correlation (NCC) over its edges')
    ax.set_title('c  Alignment evidence (dot size = number of edges)'); ax.set_ylim(0, 0.9)
    save(fig, 'F1_registration')


# ---------------------------------------------------------------- F2
def fig2(P, M, A, cv):
    y0, y1, x0, x1 = -22.0, 58.0, -14.0, 116.0
    ncol = 4; nrow = 5; fig, axs = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 4.4 * nrow * (y1 - y0) / (x1 - x0) + 0.5))
    ylm, xlm = LM
    for k, n in enumerate(ORDER):
        a = axs[k // ncol, k % ncol]; S = Sess(n, INV); Pn = P['P'][n]
        im = mean_on_canvas(S, Pn, cv, y0, y1, x0, x1)
        yy, xx = np.meshgrid(np.arange(y0, y1, FG), np.arange(x0, x1, FG), indexing='ij'); blob = np.hypot(yy - ylm, xx - xlm) < 9
        v = im[np.isfinite(im) & ~blob]
        if len(v) > 100:
            lo, hi = np.percentile(v, 1), np.percentile(v, 99.5); disp = np.clip((im - lo) / (hi - lo + 1e-9), 0, 1) ** 0.7
        else: disp = np.zeros_like(im)
        a.imshow(np.where(np.isfinite(im), disp, 0.0), cmap='gray', extent=(x0, x1, y1, y0), interpolation='bilinear', vmin=0, vmax=1)
        if n in EXCLUDE_ROI_SESSIONS:
            a.set_title('%s  (ROIs not used: horizontal bands)' % short(n), fontsize=7, color=DCOL[DAY[n]])
        else:
            masks = rasterise(S, Pn); oy, ox = int(round(y0 / FG)), int(round(x0 / FG))
            ov = outline(masks, im.shape, oy, ox)
            rgba = np.zeros(im.shape + (4,)); rgba[ov] = (1.0, 0.25, 0.2, 0.9); a.imshow(rgba, extent=(x0, x1, y1, y0), interpolation='nearest')
            a.set_title('%s  %d accepted ROIs' % (short(n), S.accepted.sum()), fontsize=7, color=DCOL[DAY[n]])
        a.set_xticks([]); a.set_yticks([])
    for k in range(len(ORDER), nrow * ncol): axs[k // ncol, k % ncol].axis('off')
    fig.suptitle('Registered mean images of all 20 recordings, same window of the common frame (%d x %d um); each panel scaled on its own 1-99.5 percentile (blob excluded); red = accepted ROIs' % (x1 - x0, y1 - y0), fontsize=8)
    fig.tight_layout(); save(fig, 'F2_registered_gallery')


# ---------------------------------------------------------------- F3
def fig3(P, M, A):
    rows = list(csv.DictReader(open(os.path.join(OUT, 'roi_match_pairs.csv'))))
    idx = {s: i for i, s in enumerate(RORDER)}; Mx = np.full((len(RORDER), len(RORDER)), np.nan)
    cat = collections.defaultdict(list); chance = collections.defaultdict(list)
    for r in rows:
        a, b = r['a'], r['b']
        if r['frac_of_smaller'] == '' or min(int(r['n_a_inview']), int(r['n_b_inview'])) < 5: continue
        f = float(r['frac_of_smaller']); Mx[idx[a], idx[b]] = Mx[idx[b], idx[a]] = f
        nov = lambda s: DAY[s] >= '11-04'
        gap = int(r['days_apart'])
        if r['same_day'] == '1': c = 'same day'
        elif nov(a) and nov(b): c = 'Nov, %d day%s apart' % (gap, '' if gap == 1 else 's')
        elif (not nov(a)) and (not nov(b)): c = 'Oct, different days'
        else: c = 'Oct vs Nov (>= 6 days)'
        cat[c].append(f); chance[c].append(float(r['chance_mean']))
    fig = plt.figure(figsize=(11.5, 4.8)); gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1], wspace=0.3)
    ax = fig.add_subplot(gs[0]); im = ax.imshow(Mx, cmap='viridis', vmin=0, vmax=0.9)
    ax.set_xticks(range(len(RORDER))); ax.set_xticklabels([short(s) for s in RORDER], rotation=90, fontsize=6); ax.set_yticks(range(len(RORDER))); ax.set_yticklabels([short(s) for s in RORDER], fontsize=6)
    plt.colorbar(im, ax=ax, fraction=0.04, label='fraction of the smaller evaluable ROI set matched (IoU >= 0.2)'); ax.set_title('a  Matched accepted ROIs per session pair (blank: < 5 ROIs in both fields)')
    ax = fig.add_subplot(gs[1]); order = ['same day', 'Nov, 1 day apart', 'Nov, 2 days apart', 'Oct, different days', 'Oct vs Nov (>= 6 days)']
    rng = np.random.default_rng(0)
    for k, c in enumerate(order):
        v = np.array(cat.get(c, [])); ch = np.array(chance.get(c, []))
        if len(v) == 0: continue
        ax.scatter(k + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color='#1F5FBF', alpha=0.7); ax.hlines(np.mean(v), k - 0.3, k + 0.3, color='k', lw=1.4)
        ax.scatter(k + 0.35 + rng.uniform(-0.05, 0.05, len(ch)), ch, s=8, color='#8A8F98', alpha=0.6)
        ax.text(k, -0.07, 'n=%d' % len(v), ha='center', fontsize=6)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=20, ha='right'); ax.set_ylim(-0.1, 1.0); ax.set_ylabel('matched fraction')
    ax.set_title('b  By relation (blue pairs, black mean; grey = random displacement control)')
    save(fig, 'F3_roi_matching')


# ---------------------------------------------------------------- F4
def fig4(P, M, A):
    cl, rows = M['clusters'], M['rows']; pres = collections.defaultdict(dict)
    for r in rows: pres[r['master']][r['session']] = r['cls']
    conf = [i for i, c in enumerate(cl) if len(c) >= 2]
    nov = [s for s in RORDER if DAY[s] >= '11-04']
    conf.sort(key=lambda i: (-sum(pres[i][s] == 'accepted' for s in nov), -len(cl[i]), i))
    code = {c: k for k, c in enumerate(CLS)}
    Z = np.array([[code[pres[i][s]] for s in RORDER] for i in conf])
    fig = plt.figure(figsize=(8.5, 10)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 14], width_ratios=[20, 1.6], hspace=0.03, wspace=0.03)
    ax = fig.add_subplot(gs[1, 0]); ax.imshow(Z, aspect='auto', cmap=ListedColormap([CCOL[c] for c in CLS]), vmin=-0.5, vmax=len(CLS) - 0.5, interpolation='nearest')
    ax.set_xticks(range(len(RORDER))); ax.set_xticklabels([short(s) for s in RORDER], rotation=90, fontsize=6); ax.set_ylabel('%d confirmed master axons (accepted in >= 2 sessions)' % len(conf)); ax.set_yticks([])
    ax.spines['left'].set_visible(False); ax.spines['bottom'].set_visible(False)
    for k in range(1, len(RORDER)):
        if DAY[RORDER[k]] != DAY[RORDER[k - 1]]: ax.axvline(k - 0.5, color='k', lw=0.6)
    at = fig.add_subplot(gs[0, 0]); at.bar(range(len(RORDER)), [sum(pres[i][s] == 'accepted' for i in conf) for s in RORDER], color=[DCOL[DAY[s]] for s in RORDER]); at.set_xlim(-0.5, len(RORDER) - 0.5)
    at.set_xticks([]); at.set_ylabel('accepted'); at.set_title('Where each confirmed axon was accepted / detected / visible', fontsize=8)
    ar = fig.add_subplot(gs[1, 1]); ar.barh(range(len(conf)), [len(cl[i]) for i in conf], color='#444'); ar.set_ylim(len(conf) - 0.5, -0.5); ar.set_yticks([]); ar.set_xlabel('n sessions', fontsize=6)
    lg = fig.add_subplot(gs[0, 1]); lg.axis('off')
    for k, c in enumerate(CLS[:4]): lg.add_patch(Rectangle((0, 1 - 0.25 * (k + 1)), 0.2, 0.2, color=CCOL[c])); lg.text(0.25, 1 - 0.25 * (k + 1) + 0.03, c, fontsize=5)
    save(fig, 'F4_presence_matrix')


# ---------------------------------------------------------------- F5
def fig5(P, M, A):
    sm = {r['session']: r for r in csv.DictReader(open(os.path.join(OUT, 'fov_session_summary.csv')))}
    fig = plt.figure(figsize=(14, 8.2)); gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.32)
    ax = fig.add_subplot(gs[0, 0]); bot = np.zeros(len(RORDER))
    for c, key in (('accepted', 'recall'), ('rejected ROI present', 'frac_rejected_present'), ('structure, no ROI', 'frac_structure_no_roi'), ('nothing visible', 'frac_nothing')):
        v = np.array([float(sm[s][key]) for s in RORDER]); ax.bar(range(len(RORDER)), v, bottom=bot, color=CCOL[c], label=c, edgecolor='none'); bot += v
    ax.set_xticks(range(len(RORDER))); ax.set_xticklabels([short(s) + ' (%s)' % sm[s]['n_confirmed_in_view'] for s in RORDER], rotation=90, fontsize=5.5); ax.set_ylabel('fraction of confirmed axons in view')
    ax.legend(fontsize=5.5, ncol=4, loc='lower center', bbox_to_anchor=(0.5, 1.0), frameon=False, columnspacing=0.8, handlelength=1.0); ax.set_title('a  What happened to the axons confirmed by the OTHER sessions (n in brackets)', pad=16, fontsize=7)
    ax = fig.add_subplot(gs[0, 1])
    for k, s in enumerate(RORDER):
        for key, col, off in (('z_dog_accepted', '#1F5FBF', -0.2), ('z_dog_missing', '#D98A00', 0.0), ('z_dog_nothing', '#8A8F98', 0.2)):
            try: v = float(sm[s][key])
            except Exception: continue
            if np.isfinite(v): ax.scatter(k + off, v, s=12, color=col)
        ax.scatter(k, float(sm[s]['bg95_dog']), marker='_', s=60, color='k')
    ax.scatter([], [], color='#1F5FBF', label='accepted'); ax.scatter([], [], color='#D98A00', label='not accepted'); ax.scatter([], [], color='#8A8F98', label='nothing visible'); ax.scatter([], [], marker='_', color='k', label='background 95th pct')
    ax.legend(fontsize=5.5, frameon=False, loc='upper left'); ax.set_xticks(range(len(RORDER))); ax.set_xticklabels([short(s) for s in RORDER], rotation=90, fontsize=5.5)
    ax.set_ylabel('axon-scale contrast of the mean image at the footprint (robust z)'); ax.set_title('b  Tonic visibility of the confirmed axons')
    ax = fig.add_subplot(gs[0, 2]); k = np.arange(1, len(A['f_obs']) + 1)
    ax.bar(k - 0.2, A['f_obs'], 0.4, color='#1F5FBF', label='observed'); ax.bar(k + 0.2, A['f_mod'], 0.4, color='#C9CED6', yerr=[np.maximum(A['f_mod'] - A['f_lo'], 0), np.maximum(A['f_hi'] - A['f_mod'], 0)], error_kw=dict(lw=0.6), label='random subsets, same rates')
    ax.set_xlabel('in how many of the %d November sessions the axon was accepted' % len(A['nov'])); ax.set_ylabel('axons (in view of all)'); ax.legend(fontsize=6, frameon=False); ax.set_title('c  Stable core + sporadic axons')
    ax = fig.add_subplot(gs[1, 0]); C = A['C']; names = A['names']; idx = {s: i for i, s in enumerate(names)}
    day = {s: np.datetime64('2026-' + DAY[s]) for s in names}; pts = collections.defaultdict(list); nov = A['nov']
    for s in nov:
        for t in nov:
            if s != t and np.isfinite(C[idx[s], idx[t]]): pts[abs((day[s] - day[t]).astype(int))].append(C[idx[s], idx[t]])
    rng = np.random.default_rng(1)
    for g in sorted(pts):
        ax.scatter(g + rng.uniform(-0.12, 0.12, len(pts[g])), pts[g], s=9, color='#1F5FBF', alpha=0.6); ax.hlines(np.mean(pts[g]), g - 0.25, g + 0.25, color='k', lw=1.4)
    ax.axhline(0.13, color='#8A8F98', ls=':', lw=0.8); ax.text(2.35, 0.14, 'chance', color='#8A8F98', fontsize=6, ha='right')
    ax.set_xticks(sorted(pts)); ax.set_xlabel('days between the two sessions (11-04..11-06 block)'); ax.set_ylabel('P(accepted in t | accepted in s, in view)'); ax.set_title('d  Carry-over of accepted axons'); ax.set_ylim(0, 1)
    ax = fig.add_subplot(gs[1, 1])
    for s in RORDER:
        r = sm[s]
        if DAY[s] < '11-04' or not np.isfinite(float(r['recall'])): continue
        ax.scatter(float(r['minutes']), float(r['recall']), s=22, color=DCOL[DAY[s]]); ax.annotate(short(s)[6:], (float(r['minutes']), float(r['recall'])), fontsize=5.5, xytext=(2, 2), textcoords='offset points')
    ax.set_xscale('log'); ax.set_xlabel('recording length (min, log)'); ax.set_ylabel('recall of confirmed axons'); ax.set_title('e  Recall vs recording length (Nov block)'); ax.set_ylim(0, 0.65)
    ax = fig.add_subplot(gs[1, 2]); d6 = [s for s in RORDER if DAY[s] == '11-06']; ref = 'AnimalRB24_261106_1317'
    rr = {}
    for (a, b, t, r, ov, exc) in P['edges']:
        if a == ref: rr[b] = r
        if b == ref: rr[a] = r
    tm = [int(INV[s]['tif_first_mtime'][11:13]) * 60 + int(INV[s]['tif_first_mtime'][14:16]) for s in d6]
    ax.plot(tm, [rr.get(s, 1.0) for s in d6], 'o-', color='#33a02c', label='image similarity to 13:17 (NCC)')
    ax.plot(tm, [float(sm[s]['recall']) for s in d6], 's-', color='#1F5FBF', label='recall of confirmed axons')
    ax.plot(tm, [float(sm[s]['frac_nothing']) for s in d6], '^-', color='#8A8F98', label='fraction nothing visible')
    for x, s, dy in zip(tm, d6, (0.02, 0.02, 0.09, 0.02)): ax.annotate(short(s)[6:], (x, dy), fontsize=5.5, ha='center')
    ax.set_xticks(tm); ax.set_xticklabels(['%d:%02d' % (t // 60, t % 60) for t in tm], fontsize=6)
    ax.set_xlabel('start time on 11-06'); ax.set_ylim(0, 1.05); ax.legend(fontsize=5.5, frameon=False, loc='upper right', bbox_to_anchor=(1.0, 0.98)); ax.set_title('f  Within-day change on 11-06 (same stage position)')
    save(fig, 'F5_what_differs')


# ---------------------------------------------------------------- F6
def fig6(P, M, A, cv):
    cl, rows = M['clusters'], M['rows']; pres = collections.defaultdict(dict)
    for r in rows: pres[r['master']][r['session']] = r['cls']
    common = A['common']; seen = A['seen']
    nov = [s for s in RORDER if DAY[s] >= '11-04']
    pick = []
    core = [i for i, k in zip(common, seen) if k >= 8]; mid = [i for i, k in zip(common, seen) if 4 <= k <= 6]; one = [i for i, k in zip(common, seen) if k == 1]
    rng = np.random.default_rng(3)
    pick = list(rng.choice(core, min(6, len(core)), replace=False)) + list(rng.choice(mid, min(4, len(mid)), replace=False)) + list(rng.choice(one, min(4, len(one)), replace=False))
    labels = ['seen in >=8'] * min(6, len(core)) + ['seen in 4-6'] * min(4, len(mid)) + ['seen once'] * min(4, len(one))
    S_all = {n: Sess(n, INV, rois=False) for n in RORDER}
    IM = {n: canvas_images(S_all[n], P['P'][n], cv) for n in RORDER}
    # footprints from the cluster members: re-rasterise the member ROIs (accepted ones)
    ROI = {n: Sess(n, INV) for n in RORDER}
    ras = {n: {m[0]: m for m in rasterise(ROI[n], P['P'][n])} for n in RORDER}
    mrow = M['mrows']
    W = 12.0; hw = int(W / FG)
    fig, axs = plt.subplots(len(pick), len(RORDER), figsize=(1.05 * len(RORDER), 1.05 * len(pick) + 0.6))
    for r_, (q, lab) in enumerate(zip(pick, labels)):
        mem = cl[q]; cells = []
        for n, rid in mem.items():
            i, ci, cj, m = ras[n][rid]; a, b = np.nonzero(m); cells += list(zip(ci[a].tolist(), cj[b].tolist()))
        cells = np.array(sorted(set(cells))); cy, cx = int(round(cells[:, 0].mean())), int(round(cells[:, 1].mean()))
        fp = np.zeros((2 * hw, 2 * hw), bool); yy = cells[:, 0] - (cy - hw); xx = cells[:, 1] - (cx - hw); ok = (yy >= 0) & (yy < 2 * hw) & (xx >= 0) & (xx < 2 * hw); fp[yy[ok], xx[ok]] = True
        edge = fp & ~ndi.binary_erosion(fp)
        for c_, n in enumerate(RORDER):
            a = axs[r_, c_]; a.set_xticks([]); a.set_yticks([])
            cls = pres[q][n]
            if cls == 'out of FOV': a.set_facecolor('#f4f4f4'); continue
            img = IM[n]['dog'][cy - cv.oy - hw:cy - cv.oy + hw, cx - cv.ox - hw:cx - cv.ox + hw]
            a.imshow(img, cmap='magma', vmin=-1, vmax=6, interpolation='bilinear')
            rgba = np.zeros(fp.shape + (4,)); col = (0.2, 0.9, 1.0, 1.0) if cls == 'accepted' else (1, 1, 1, 0.55)
            rgba[edge] = col; a.imshow(rgba, interpolation='nearest')
            for sp in a.spines.values(): sp.set_visible(True); sp.set_linewidth(1.6); sp.set_color(CCOL[cls] if cls != 'nothing visible' else '#BBBBBB')
            if r_ == 0: a.set_title(short(n), fontsize=5, rotation=90, color=DCOL[DAY[n]])
        axs[r_, 0].set_ylabel(lab, fontsize=5.5, rotation=0, ha='right', va='center')
    fig.suptitle('Example master axons across all sessions (%d x %d um windows; axon-scale band-pass image of each session in robust z, fixed scale -1..6; outline = footprint, cyan = accepted ROI; frame colour = class as in F4)' % (W, W), fontsize=7)
    fig.subplots_adjust(wspace=0.04, hspace=0.04); save(fig, 'F6_example_axons')


if __name__ == '__main__':
    which = sys.argv[1:] or ['F1', 'F2', 'F3', 'F4', 'F5', 'F6']
    P, M, A = load_all(); cv = make_canvas(P)
    if 'F1' in which: fig1(P, M, A)
    if 'F2' in which: fig2(P, M, A, cv)
    if 'F3' in which: fig3(P, M, A)
    if 'F4' in which: fig4(P, M, A)
    if 'F5' in which: fig5(P, M, A)
    if 'F6' in which: fig6(P, M, A, cv)
