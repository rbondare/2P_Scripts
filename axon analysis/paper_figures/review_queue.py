"""ROI review queue: for every aggregated session with rejected suite2p ROIs, a ranked thumbnail sheet of the rejected ROIs that look most like axons / boutons, and a decision sheet to fill in.
usage: python review_queue.py [session ...]          (default: every session that has a cache/rejroi_<session>.pkl table; build the tables with rejected_rois.py first)
output (plots/summary_figures/rejected_roi_review/):  review_<session>.pdf   labels_review_<session>.csv   review_queue_summary.csv   README.txt

Queue = rejected ROIs (iscell 0) that are not hidden constituents of a merged ROI (suite2p inmerge), do not overlap an accepted ROI, and are activity-like: they pass the axon transient test used for the
accepted ROIs OR have skewness > 0.5 (transient-rich trace).  Ranked by an activity score = mean percentile rank (within the session's rejected ROIs) of skewness, slow variance (transient test), local contrast
and F / Fneu.  At most MAX_Q ROIs per session.  The thumbnails show only morphology and the dF/F trace (no behaviour), so the decision cannot be driven by the coupling that is being tested.
Rejected ROIs with a median neuropil-corrected fluorescence <= 5 counts (dim) are not part of the tables and therefore not in the queue.
Fill in the 'decision' column of labels_review_<session>.csv with accept / reject / unsure (blank = no decision), then run apply_curation.py.
"""
import os, sys, csv, glob, pickle, textwrap, numpy as np
from fig_common import *
import xs_common as xc
from scipy.stats import rankdata
from matplotlib.backends.backend_pdf import PdfPages

REV = os.path.join(OUT, 'rejected_roi_review'); MAX_Q = 160; PER_PAGE = 32; C_CAND, C_ACC = '#FFD400', '#00E5FF'


def percentile_rank(x):
    x = np.asarray(x, float); ok = np.isfinite(x); r = np.full(len(x), 0.5)
    if ok.sum() > 1: r[ok] = (rankdata(x[ok]) - 1) / (ok.sum() - 1)
    return r


def queue_of(d):
    g = d['grp']; hid = d.get('inmerge', np.zeros(len(g))) > 0; R = (g == 1) & ~hid; ids = np.where(R)[0]
    score = np.mean([percentile_rank(d['skew'][ids]), percentile_rank(np.log10(np.maximum(d['slow_ratio'][ids], 1e-3))), percentile_rank(d['contrast'][ids]), percentile_rank(np.log2(np.maximum(d['F_over_Fneu'][ids], 1e-3)))], axis=0)
    cand = (d['passes'][ids] | (d['skew'][ids] > 0.5)) & (d['overlap'][ids] == 0)
    order = ids[cand][np.argsort(-score[cand])][:MAX_Q]; sc = dict(zip(ids, score))
    return order, sc, dict(n_rej=int(R.sum()), n_hidden=int(((g == 1) & hid).sum()), n_cand=int(cand.sum()), n_queued=len(order))


def build(name):
    d = pickle.load(open(os.path.join(CACHE, 'rejroi_%s.pkl' % name), 'rb')); order, sc, info = queue_of(d); os.makedirs(REV, exist_ok=True)
    p = os.path.join(xc.DATA_DIR, name, 'suite2p', 'plane0'); S = xc.Session(name, with_rejected=False, verbose=False)
    stat = np.load(os.path.join(p, 'stat.npy'), allow_pickle=True); ops = np.load(os.path.join(p, 'ops.npy'), allow_pickle=True).item(); mimg = np.asarray(ops['meanImg'], float); Ly, Lx = mimg.shape
    ic = np.load(os.path.join(p, 'iscell.npy'))[:, 0] > 0; accmask = np.zeros((Ly, Lx), bool)
    for k in np.where(ic)[0]: accmask[stat[k]['ypix'], stat[k]['xpix']] = True
    F = np.load(os.path.join(p, 'F.npy'), mmap_mode='r'); Fn = np.load(os.path.join(p, 'Fneu.npy'), mmap_mode='r'); sidx = d['idx'][order]
    Fg = xc.bin_rows(S.tca, np.asarray(F[sidx]).astype(np.float32), S.edges); Fng = xc.bin_rows(S.tca, np.asarray(Fn[sidx]).astype(np.float32), S.edges); dff = xc.best_dff_grid(Fg, Fng)[0]
    # decision sheet (never overwritten once it contains a decision or a comment)
    f_lab = os.path.join(REV, 'labels_review_%s.csv' % name); filled = False
    if os.path.exists(f_lab):
        rows_ = list(csv.reader(open(f_lab, newline=''))); filled = any((len(r) > 6 and r[6].strip()) or (len(r) > 7 and r[7].strip()) for r in rows_[1:])
    if filled and not os.environ.get('FORCE_SHEETS'): print('   %s: decision sheet already contains decisions: kept as it is' % name)
    else:
        with open(f_lab, 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['rank', 'suite2p_roi_index', 'activity_score', 'skewness', 'passes_transient_test', 'touches_accepted_roi', 'decision (accept / reject / unsure)', 'comment'])
            for r_, i in enumerate(order, 1): w.writerow([r_, int(d['idx'][i]), '%.2f' % sc[i], '%.2f' % d['skew'][i], int(d['passes'][i]), int(d['adjacent'][i] > 0), '', ''])
    with PdfPages(os.path.join(REV, 'review_%s.pdf' % name)) as pdf:
        for pg in range(0, max(len(order), 1), PER_PAGE):
            fig = plt.figure(figsize=(W2, 250 * MM))
            fig.text(0.01, 0.992, textwrap.fill('%s - rejected ROIs ranked by how axon-like they look (rank 1 = most), %d-%d of %d.  Left: mean image (yellow = this ROI, cyan = accepted ROIs); right: session dF/F (top) and 90 s around its largest transient (bottom).  Decide in labels_review_%s.csv: accept / reject / unsure.  '
                                                '#rank  roi index   T = passes the transient test, adj = touches an accepted ROI.' % (short(name), pg + 1, min(pg + PER_PAGE, len(order)), len(order), name), 150), fontsize=6.2, fontweight='bold', va='top')
            gs = fig.add_gridspec(8, 4, left=0.01, right=0.99, top=0.955, bottom=0.01, hspace=0.30, wspace=0.07)
            for q, i in enumerate(order[pg:pg + PER_PAGE]):
                r_, c_ = divmod(q, 4); sub = gs[r_, c_].subgridspec(2, 2, width_ratios=[1, 1.45], height_ratios=[1, 1], wspace=0.05, hspace=0.12); a = fig.add_subplot(sub[:, 0]); b1 = fig.add_subplot(sub[0, 1]); b2 = fig.add_subplot(sub[1, 1])
                s = stat[d['idx'][i]]; y, x = s['ypix'], s['xpix']; cy, cx = int(np.median(y)), int(np.median(x)); y0, y1, x0, x1 = max(cy - 25, 0), min(cy + 26, Ly), max(cx - 25, 0), min(cx + 26, Lx)
                crop = mimg[y0:y1, x0:x1]; lo, hi = np.percentile(crop, [1, 99.7]); a.imshow(crop, cmap='gray', vmin=lo, vmax=hi, interpolation='nearest')
                m = np.zeros(crop.shape); m[np.clip(y - y0, 0, crop.shape[0] - 1), np.clip(x - x0, 0, crop.shape[1] - 1)] = 1; a.contour(m, levels=[0.5], colors=C_CAND, linewidths=0.7)
                if accmask[y0:y1, x0:x1].any(): a.contour(accmask[y0:y1, x0:x1].astype(float), levels=[0.5], colors=C_ACC, linewidths=0.35)
                a.axis('off'); a.text(0.02, 0.98, '#%d' % (pg + q + 1), transform=a.transAxes, color='white', fontsize=6.5, fontweight='bold', va='top', bbox=dict(fc='black', ec='none', alpha=0.6, pad=0.8))
                a.text(0.02, 0.02, 'roi %d%s%s' % (d['idx'][i], '  T' if d['passes'][i] else '', '  adj' if d['adjacent'][i] > 0 else ''), transform=a.transAxes, color='white', fontsize=4.8, va='bottom', bbox=dict(fc='black', ec='none', alpha=0.55, pad=0.6))
                tr = xc.nsmooth(dff[pg + q], 1.0 / xc.DT); b1.plot(S.t / 60, tr, color=INK, lw=0.35); b1.set_xlim(S.t[0] / 60, S.t[-1] / 60); b1.set_yticks([]); b1.tick_params(labelsize=4.2, length=1.2)
                k = int(np.nanargmax(np.where(np.isfinite(tr), tr, -np.inf))); j0 = max(k - int(45 / xc.DT), 0); j1 = min(k + int(45 / xc.DT), len(tr)); b2.plot(S.t[j0:j1] - S.t[j0], tr[j0:j1], color=INK, lw=0.5); b2.set_yticks([]); b2.tick_params(labelsize=4.2, length=1.2); b2.set_xlim(0, S.t[j1 - 1] - S.t[j0])
                for ax in (b1, b2):
                    for sp in ('left', 'top', 'right'): ax.spines[sp].set_visible(False)
                if r_ < 7: b1.set_xticklabels([]); b2.set_xticklabels([])
            pdf.savefig(fig)
            if pg == 0: fig.savefig(os.path.join(REV, 'preview_%s.png' % name), dpi=100, facecolor='white')
            plt.close(fig)
    return dict(session=name, **info)


if __name__ == '__main__':
    names = sys.argv[1:] or [os.path.basename(f)[7:-4] for f in sorted(glob.glob(os.path.join(CACHE, 'rejroi_*.pkl')))]
    for n in names:
        r = build(n); print('%-26s rejected %3d (hidden merge constituents %2d) -> %3d candidates, %3d in the queue' % (n, r['n_rej'], r['n_hidden'], r['n_cand'], r['n_queued']), flush=True)
    os.makedirs(REV, exist_ok=True); rows = []                  # summary of ALL sessions (cheap: from the cached tables), whatever subset was rebuilt
    for f in sorted(glob.glob(os.path.join(CACHE, 'rejroi_*.pkl'))):
        d = pickle.load(open(f, 'rb')); rows.append(dict(session=d['name'], **queue_of(d)[2]))
    with open(os.path.join(REV, 'review_queue_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['session', 'n_rej', 'n_hidden', 'n_cand', 'n_queued']); w.writeheader(); [w.writerow(r) for r in rows]
    open(os.path.join(REV, 'README.txt'), 'w').write(__doc__ + '\nAfter filling in the decisions: python apply_curation.py [--apply] [session ...]   (dry run without --apply)\n')
