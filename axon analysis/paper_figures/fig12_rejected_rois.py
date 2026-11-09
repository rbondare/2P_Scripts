"""Figure 12 - what are the 'rejected' suite2p ROIs? (are they unpicked boutons / axon fragments, or background?)
a-c  feature distributions: accepted signal-like ROIs | rejected, activity-like (pass the axon transient test) | rejected, strict control (fail the test, low skew, not touching accepted ROIs)
d    coupling to whisking against activity-likeness (rejected ROIs ranked by skew, 5 bins)     e  fraction of ROIs coupled positively / negatively to whisking in each group, per session
f    coupling against distance to the nearest accepted ROI      g  neuropil test: r(whisking) with dF/F from F - 0.5 Fneu (as used everywhere), from F alone, from Fneu alone
h    r(whisking) against F / Fneu     i  session table
Also writes the blind labelling galleries (plots/summary_figures/rejected_roi_gallery/): thumbnails of accepted + rejected ROIs in random order (mean image + ROI outline + dF/F trace, no behaviour), a labels sheet to fill in and a key.
usage: python fig12_rejected_rois.py (after python rejected_rois.py) -> Fig12_rejected_rois.pdf/.png, gallery, ledger entries
"""
import os, pickle, csv, numpy as np
from fig_common import *
import xs_common as xc
from scipy.stats import spearmanr

C_ACT = '#0F8B8D'; C_CTRL = '#5B6270'; SKEW_STRICT = 0.3
GAL = os.path.join(OUT, 'rejected_roi_gallery'); GAL_SESS = ['AnimalRB24_261021_1630', 'AnimalRB22_261027_1138', 'AnimalRB23_261023_1642']


def load_all():
    T = {}
    for n in UNFLAGGED:
        f = os.path.join(CACHE, 'rejroi_%s.pkl' % n)
        if os.path.exists(f): T[n] = pickle.load(open(f, 'rb'))
    return T


def masks(d):
    g = d['grp']; hid = d.get('inmerge', np.zeros(len(g))) > 0; R = (g == 1) & ~hid; A = (g == 0) & d['sig_acc']; act = R & d['passes']; strict = R & ~d['passes'] & (d['skew'] < SKEW_STRICT) & (d['adjacent'] == 0)
    return dict(acc=A, rej=R, act=act, strict=strict)


def posneg(d, m, key='whisk', suffix=''):
    r = d['r_%s%s' % (key, suffix)]; p = d['p_%s%s' % (key, suffix)]; ok = m & np.isfinite(r)
    return (100 * np.mean((r[ok] > 0) & (p[ok] < 0.01)), 100 * np.mean((r[ok] < 0) & (p[ok] < 0.01)), int(ok.sum())) if ok.sum() else (np.nan, np.nan, 0)


def violin(ax, data, x0, color, w=0.38, lims=None):
    data = np.asarray(data, float); data = data[np.isfinite(data)]
    if len(data) < 5: return
    lo, hi = np.percentile(data, [1, 99])
    if lims: lo, hi = max(lo, lims[0]), min(hi, lims[1]); data = data[(data >= lo) & (data <= hi)]
    kde_x = np.linspace(lo, hi, 80)
    from scipy.stats import gaussian_kde
    k = gaussian_kde(data)(kde_x); k = k / k.max() * w
    ax.fill_betweenx(kde_x, x0 - k, x0 + k, color=color, alpha=0.55, lw=0); q = np.percentile(data, [25, 50, 75]); ax.plot([x0 - w * 0.6, x0 + w * 0.6], [q[1]] * 2, color=INK, lw=1.2); ax.plot([x0, x0], [q[0], q[2]], color=INK, lw=0.8)


def main():
    T = load_all(); M = {n: masks(T[n]) for n in T}; names = list(T)
    fig = plt.figure(figsize=(W2, 226 * MM)); gs = fig.add_gridspec(3, 3, left=0.075, right=0.985, top=0.935, bottom=0.085, hspace=0.50, wspace=0.34)
    rng = np.random.default_rng(3)
    # ------------------------------------------------ a-c  distributions
    feats = [('skew', 'skewness of dF/F', lambda v: v), ('F_over_Fneu', 'median F / median Fneu (log2)', lambda v: np.log2(np.maximum(v, 1e-3))), ('contrast', 'local contrast (mean image, ROI / ring)', lambda v: v)]
    for j, (k, lab, tf) in enumerate(feats):
        a = fig.add_subplot(gs[0, j])
        for x0, key, col in ((0, 'acc', C_SIG), (1, 'act', C_ACT), (2, 'strict', C_CTRL)):
            v = np.concatenate([tf(T[n][k][M[n][key]]) for n in names]); violin(a, v, x0, col, lims={'skew': (-1, 3), 'F_over_Fneu': (-1, 4), 'contrast': (0.5, 3.2)}[k])
        a.set_xticks([0, 1, 2]); a.set_xticklabels(['accepted\nsignal-like', 'rejected,\nactivity-like', 'rejected,\nstrict control'], fontsize=4.9); a.set_ylabel(lab, fontsize=5.4); a.set_xlim(-0.6, 2.6)
        if j == 0: plabel(a, 'a', dx=-26, dy=8)
        else: plabel(a, 'bc'[j - 1], dx=-26, dy=8)
        if k == 'contrast': a.set_ylim(0.5, 3.2)
        if k == 'skew': a.set_ylim(-1, 3)
    n_act = {n: M[n]['act'].sum() / M[n]['rej'].sum() for n in names}; h_act = pooled_unflagged(n_act, n_boot=4000)
    fig.text(0.075, 0.997, wrap12('13 unflagged sessions: accepted signal-like %d ROIs; rejected %d ROIs, of which %.0f%% [%.0f, %.0f] (mean over sessions) pass the axon transient test (activity-like) and %d form the strict control (fail the test, skew < %.1f, not touching accepted ROIs)' % (
        sum(M[n]['acc'].sum() for n in names), sum(M[n]['rej'].sum() for n in names), 100 * h_act['est'], 100 * h_act['lo'], 100 * h_act['hi'], sum(M[n]['strict'].sum() for n in names), SKEW_STRICT), 170), fontsize=5.2, color=MUTED, va='top')
    ledger_add('Fig12', 'rej_activity_like', 'fraction of rejected ROIs that pass the axon transient criterion (5 Hz grid version)', h_act['est'], h_act['lo'], h_act['hi'], h_act['n_sessions'], h_act['n_animals'], 'ambiguous', note='possible unpicked boutons / axon fragments in the control; the grid version of the test is lenient (accepted signal-like ROIs pass it 40-90%)', unit='fraction of rejected ROIs')
    # ------------------------------------------------ d  sensitivity curve (rejected ROIs ranked by skew)
    ad = fig.add_subplot(gs[1, 0]); NB = 5; P_ = []; N_ = []
    for n in names:
        d = T[n]; m = M[n]['rej'] & np.isfinite(d['r_whisk']); idx = np.where(m)[0]
        if len(idx) < 25: continue
        o = idx[np.argsort(d['skew'][idx])]; bins = np.array_split(o, NB); pos = [100 * np.mean((d['r_whisk'][b] > 0) & (d['p_whisk'][b] < 0.01)) for b in bins]; neg = [100 * np.mean((d['r_whisk'][b] < 0) & (d['p_whisk'][b] < 0.01)) for b in bins]
        P_.append(pos); N_.append(neg); ad.plot(range(NB), pos, color=FAINT, lw=0.6); ad.plot(range(NB), neg, color=FAINT, lw=0.6, ls='--')
    ad.plot(range(NB), np.mean(P_, 0), color=INK, lw=1.6, label='significantly positive'); ad.plot(range(NB), np.mean(N_, 0), color=INK, lw=1.6, ls='--', label='significantly negative')
    pa = np.mean([posneg(T[n], M[n]['acc'])[0] for n in names if M[n]['acc'].sum() >= 5]); na = np.mean([posneg(T[n], M[n]['acc'])[1] for n in names if M[n]['acc'].sum() >= 5])
    ad.axhline(pa, color=C_SIG, lw=0.9, ls='-', alpha=0.8); ad.axhline(na, color=C_SIG, lw=0.9, ls='--', alpha=0.8); ad.text(NB - 0.95, pa + 2, 'accepted axons', color=C_SIG, fontsize=4.8, ha='right')
    ad.set_xticks(range(NB)); ad.set_xticklabels(['lowest\nskew', '2', '3', '4', 'highest\nskew'], fontsize=4.9); ad.set_xlabel('rejected ROIs ranked by skewness (per session)', fontsize=5.4); ad.set_ylabel('% of ROIs coupled to whisking (p < 0.01)', fontsize=5.4); ad.set_ylim(0, 100)
    ad.legend(frameon=False, fontsize=4.6, loc='upper center', handlelength=1.6, borderpad=0.1, bbox_to_anchor=(0.5, 1.0)); plabel(ad, 'd', dx=-26, dy=8)
    # ------------------------------------------------ e  groups, per session
    ae = fig.add_subplot(gs[1, 1]); groups = [('acc', 'accepted\nsignal-like', C_SIG), ('rej', 'rejected,\nall', C_REJ), ('act', 'rejected,\nactivity-like', C_ACT), ('strict', 'rejected,\nstrict control', C_CTRL)]; xt = []; xl = []
    for i, (key, lab, col) in enumerate(groups):
        for q, (kind, ls) in enumerate((('pos', 0), ('neg', 1))):
            x0 = i * 2.6 + q * 0.95; vals = {}
            for n in names:
                pn = posneg(T[n], M[n][key])
                if pn[2] >= 8: vals[n] = pn[q]
            for n, v in vals.items():
                if q == 0: sdot(ae, x0 + rng.uniform(-0.2, 0.2), v, n, color=col, ms=3.3)
                else: ae.plot(x0 + rng.uniform(-0.2, 0.2), v, MARK[animal(n)], mfc='white', mec=col, mew=0.8, ms=3.3, zorder=3)
            h = hboot(list(vals.values()), list(vals.keys()), n_boot=3000)
            if h['n_sessions']: ae.plot([x0 - 0.38, x0 + 0.38], [h['est']] * 2, color=INK, lw=1.4, zorder=5); ae.plot([x0, x0], [h['lo'], h['hi']], color=INK, lw=0.8, zorder=5)
            if key == 'strict': ledger_add('Fig12', 'strict_%s' % kind, 'strict control ROIs significantly %s coupled to whisking' % ('positively' if q == 0 else 'negatively'), h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'ambiguous', note='strict control = fail transient test, skew < %.1f, not touching accepted ROIs; accepted signal-like: %s' % (SKEW_STRICT, '%.0f%% positive, %.0f%% negative' % (pa, na)), unit='% of ROIs')
        xt += [i * 2.6 + 0.47]; xl.append(lab)
    ae.set_xticks(xt); ae.set_xticklabels(xl, fontsize=4.7); ae.set_ylim(0, 106); ae.set_yticks([0, 20, 40, 60, 80, 100]); ae.set_xlim(-0.8, 3 * 2.6 + 1.9); ae.set_ylabel('% of ROIs coupled to whisking (p < 0.01)', fontsize=5.4); plabel(ae, 'e', dx=-26, dy=14)
    ae.set_title('filled = positive, open = negative; dot = session (>= 8 ROIs)', fontsize=4.8, loc='left', fontweight='normal', color=MUTED, pad=3)
    # ------------------------------------------------ f  distance to accepted ROIs
    af = fig.add_subplot(gs[1, 2]); edges = [('touching\n(<= 2 px)', None), ('3-30 px', (0, 30)), ('30-80 px', (30, 80)), ('> 80 px', (80, 1e9))]; res = []
    allR = {k: np.concatenate([T[n][k][M[n]['rej']] for n in names]) for k in ('adjacent', 'nn_dist', 'r_whisk', 'p_whisk')}
    for i, (lab, rg) in enumerate(edges):
        m = (allR['adjacent'] > 0) if rg is None else ((allR['adjacent'] == 0) & (allR['nn_dist'] >= rg[0]) & (allR['nn_dist'] < rg[1])); m &= np.isfinite(allR['r_whisk'])
        pos = 100 * np.mean((allR['r_whisk'][m] > 0) & (allR['p_whisk'][m] < 0.01)); neg = 100 * np.mean((allR['r_whisk'][m] < 0) & (allR['p_whisk'][m] < 0.01)); res.append((pos, neg, int(m.sum())))
        af.bar(i - 0.2, pos, 0.38, color=INK, alpha=0.85); af.bar(i + 0.2, neg, 0.38, color='white', ec=INK, lw=0.8); af.text(i, 101, 'n=%d' % m.sum(), ha='center', va='bottom', fontsize=4.5, color=MUTED)
    af.set_xticks(range(4)); af.set_xticklabels([e[0] for e in edges], fontsize=4.9); af.set_ylim(0, 100); af.set_xlabel('distance to the nearest accepted ROI', fontsize=5.4); af.set_ylabel('% of rejected ROIs coupled', fontsize=5.4); plabel(af, 'f', dx=-26, dy=8)
    af.set_title('filled = positive, open = negative (pooled ROIs)', fontsize=4.8, loc='left', fontweight='normal', color=MUTED, pad=9)
    # ------------------------------------------------ g  neuropil test
    ag = fig.add_subplot(gs[2, 0]); bins = np.linspace(-1, 1, 41); ctr = 0.5 * (bins[1:] + bins[:-1]); txt = []; hl = []
    for suffix, lab, ls, col in (('', 'F − 0.5 Fneu (as used)', '-', C_CTRL), ('_rawF', 'F alone', '--', C_CTRL), ('_Fneu', 'Fneu alone', ':', INK)):
        r = np.concatenate([T[n]['r_whisk' + suffix][M[n]['rej']] for n in names]); p = np.concatenate([T[n]['p_whisk' + suffix][M[n]['rej']] for n in names]); ok = np.isfinite(r)
        cnt, _ = np.histogram(r[ok], bins=bins); hl.append(ag.plot(ctr, 100 * cnt / cnt.sum(), color=col, ls=ls, lw=1.3 if suffix != '_Fneu' else 1.6)[0]); txt.append((lab, 100 * np.mean((r[ok] > 0) & (p[ok] < 0.01)), 100 * np.mean((r[ok] < 0) & (p[ok] < 0.01))))
        if suffix == '_rawF':
            ledger_add('Fig12', 'neuropil_test_rawF_neg', 'rejected ROIs significantly NEGATIVELY coupled to whisking with dF/F from F alone (no neuropil subtraction)', np.nan, np.nan, np.nan, len(names), 3, 'ambiguous', note='%.0f%% negative / %.0f%% positive with F alone vs %.0f%% / %.0f%% with F - 0.5 Fneu: the sign split is not a subtraction artefact' % (100 * np.mean((r[ok] < 0) & (p[ok] < 0.01)), 100 * np.mean((r[ok] > 0) & (p[ok] < 0.01)), txt[0][2], txt[0][1]), unit='% of ROIs')
    ag.axvline(0, color=MUTED, lw=0.4, ls=':'); ag.set_xlabel('r with whisking (band-passed)', fontsize=5.4); ag.set_ylabel('% of rejected ROIs', fontsize=5.4); ag.legend(hl, ['%s: %.0f%% pos / %.0f%% neg' % (lab, po, ne) for lab, po, ne in txt], frameon=False, fontsize=4.6, loc='upper left', handlelength=1.8, borderpad=0.1); plabel(ag, 'g', dx=-26, dy=8); ag.set_ylim(0, 24)
    # ------------------------------------------------ h  r vs F / Fneu
    ah = fig.add_subplot(gs[2, 1]); xr = np.concatenate([np.log2(np.maximum(T[n]['F_over_Fneu'][M[n]['rej']], 1e-3)) for n in names]); yr = np.concatenate([T[n]['r_whisk'][M[n]['rej']] for n in names]); ok = np.isfinite(yr) & np.isfinite(xr)
    ah.plot(xr[ok], yr[ok], 'o', ms=1.4, mfc=C_REJ, mec='none', alpha=0.25, rasterized=True)
    xa = np.concatenate([np.log2(np.maximum(T[n]['F_over_Fneu'][M[n]['acc']], 1e-3)) for n in names]); ya = np.concatenate([T[n]['r_whisk'][M[n]['acc']] for n in names]); oka = np.isfinite(ya)
    ah.plot(xa[oka], ya[oka], 'o', ms=2.4, mfc=C_SIG, mec='white', mew=0.2, alpha=0.9); rho = spearmanr(xr[ok], yr[ok])[0]
    ah.axhline(0, color=MUTED, lw=0.4); ah.set_xlabel('median F / median Fneu (log2; 0 = ROI as bright as its neuropil)', fontsize=5.0); ah.set_ylabel('r with whisking', fontsize=5.4); ah.set_xlim(-1.5, 4.5); ah.set_ylim(-1, 1); plabel(ah, 'h', dx=-26, dy=8)
    ah.text(0.02, 0.03, 'rejected (grey): Spearman %+.2f; accepted signal-like (blue)' % rho, transform=ah.transAxes, fontsize=4.6, va='bottom')
    # ------------------------------------------------ i  table
    ai = fig.add_subplot(gs[2, 2]); ai.axis('off'); plabel(ai, 'i', dx=-6, dy=8)
    xs_ = (0.0, 0.54, 0.72, 0.90, 1.05)
    for x_, h_ in zip(xs_, ('session', 'accepted\nsig.-like', 'rejected', 'activity-\nlike', 'strict\ncontrol')): ai.text(x_, 1.02, h_, fontsize=4.8, fontweight='bold', va='top', ha='left' if x_ == 0 else 'right')
    for i, n in enumerate(names):
        for x_, v_ in zip(xs_, (short(n), '%d' % M[n]['acc'].sum(), '%d' % M[n]['rej'].sum(), '%.0f%%' % (100 * n_act[n]), '%d' % M[n]['strict'].sum())): ai.text(x_, 0.90 - i * 0.065, v_, fontsize=4.8, va='top', ha='left' if x_ == 0 else 'right')
    fig.text(0.075, 0.004, wrap12('All ROIs of the session\'s suite2p plane0 output with the same F / Fneu and dF/F recipe; rejected = iscell 0, not a hidden constituent of a merged ROI (suite2p inmerge > 0), median Fc > 5 counts. Activity-like = passes the transient test used for the accepted axons (slow variance >= 4 x white-noise expectation, '
                                  'upward events >= 2 x downward events) applied on the 5 Hz grid (approximate: accepted ROIs pass 40-90%). Coupling = r (lag 0, band-passed 1-30 s) with whisking and its circular-shift p < 0.01 (shift >= 120 s). Panel g: the same ROIs (first 220 rejected ROIs per session) '
                                  'with dF/F recomputed from F alone or from the neuropil trace alone. 13 unflagged sessions, 3 mice.'), fontsize=4.6, color=MUTED)
    savefig(fig, 'Fig12_rejected_rois')
    return T


def wrap12(s, w=212):
    import textwrap
    return textwrap.fill(s.replace('\n', ' '), w)


# ------------------------------------------------------------------ blind labelling gallery
def gallery(name, T, n_act=20, n_rest=20, n_acc=10, per_page=32):
    from matplotlib.backends.backend_pdf import PdfPages
    if os.path.exists(os.path.join(GAL, 'labels_%s.csv' % name)) and not os.environ.get('FORCE_GALLERY'): print('gallery', name, 'exists: not overwritten (a labels sheet may already be filled in; FORCE_GALLERY=1 to rebuild)'); return
    d = T[name]; M = masks(d); rng = np.random.default_rng(11); os.makedirs(GAL, exist_ok=True)
    pick = [(i, 'rejected, activity-like') for i in rng.choice(np.where(M['act'])[0], min(n_act, M['act'].sum()), replace=False)] + [(i, 'rejected, not activity-like') for i in rng.choice(np.where(M['rej'] & ~M['act'])[0], n_rest, replace=False)] + \
           [(i, 'accepted') for i in rng.choice(np.where(d['grp'] == 0)[0], min(n_acc, (d['grp'] == 0).sum()), replace=False)]
    order = rng.permutation(len(pick)); pick = [pick[i] for i in order]
    S = xc.Session(name, with_rejected=False, verbose=False); p = os.path.join(xc.DATA_DIR, name, 'suite2p', 'plane0'); stat = np.load(os.path.join(p, 'stat.npy'), allow_pickle=True); ops = np.load(os.path.join(p, 'ops.npy'), allow_pickle=True).item(); mimg = np.asarray(ops['meanImg'], float)
    F = np.load(os.path.join(p, 'F.npy'), mmap_mode='r'); Fn = np.load(os.path.join(p, 'Fneu.npy'), mmap_mode='r'); ids = np.array([d['idx'][i] for i, _ in pick])
    Fg = xc.bin_rows(S.tca, np.asarray(F[ids]).astype(np.float32), S.edges); Fng = xc.bin_rows(S.tca, np.asarray(Fn[ids]).astype(np.float32), S.edges); dff = xc.best_dff_grid(Fg, Fng)[0]; Ly, Lx = mimg.shape
    with open(os.path.join(GAL, 'labels_%s.csv' % name), 'w', newline='') as fl, open(os.path.join(GAL, 'KEY_%s_do_not_open_before_labelling.csv' % name), 'w', newline='') as fk:
        wl = csv.writer(fl); wk = csv.writer(fk); wl.writerow(['thumb', 'label (bouton / axon segment / neuropil-background / noise-artifact / unclear)', 'comment']); wk.writerow(['thumb', 'suite2p roi index', 'group', 'skew', 'passes transient test', 'adjacent to accepted'])
        for t_, (i, grp) in enumerate(pick, 1): wl.writerow([t_, '', '']); wk.writerow([t_, int(d['idx'][i]), grp, '%.2f' % d['skew'][i], bool(d['passes'][i]), bool(d['adjacent'][i] > 0) if d['grp'][i] == 1 else ''])
    with PdfPages(os.path.join(GAL, 'gallery_%s.pdf' % name)) as pdf:
        for pg in range(0, len(pick), per_page):
            fig = plt.figure(figsize=(W2, 250 * MM)); fig.text(0.01, 0.992, '%s - blind labelling sheet, thumbnails %d-%d (mean image 51 x 51 px with ROI outline | session dF/F trace; no behaviour shown)' % (short(name), pg + 1, min(pg + per_page, len(pick))), fontsize=7, fontweight='bold', va='top')
            gs = fig.add_gridspec(8, 4, left=0.01, right=0.99, top=0.965, bottom=0.01, hspace=0.18, wspace=0.06)
            for q, (i, grp) in enumerate(pick[pg:pg + per_page]):
                r_, c_ = divmod(q, 4); sub = gs[r_, c_].subgridspec(1, 2, width_ratios=[1, 1.35], wspace=0.04); a = fig.add_subplot(sub[0]); b = fig.add_subplot(sub[1]); s = stat[d['idx'][i]]; y, x = s['ypix'], s['xpix']
                cy, cx = int(np.median(y)), int(np.median(x)); y0, y1, x0, x1 = max(cy - 25, 0), min(cy + 26, Ly), max(cx - 25, 0), min(cx + 26, Lx); crop = mimg[y0:y1, x0:x1]; lo, hi = np.percentile(crop, [1, 99.7])
                a.imshow(crop, cmap='gray', vmin=lo, vmax=hi, interpolation='nearest'); m = np.zeros(crop.shape); m[np.clip(y - y0, 0, crop.shape[0] - 1), np.clip(x - x0, 0, crop.shape[1] - 1)] = 1; a.contour(m, levels=[0.5], colors='#FFD400', linewidths=0.6); a.axis('off')
                a.text(0.03, 0.97, str(pg + q + 1), transform=a.transAxes, color='white', fontsize=7, fontweight='bold', va='top', bbox=dict(fc='black', ec='none', alpha=0.6, pad=1))
                tr = xc.nsmooth(dff[[j for j, (ii, _) in enumerate(pick) if ii == i][0]], 1.0 / xc.DT); b.plot(S.t / 60, tr, color=INK, lw=0.4); b.set_xlim(S.t[0] / 60, S.t[-1] / 60); b.set_yticks([]); b.tick_params(labelsize=4.5, length=1.5)
                for sp in ('left', 'top', 'right'): b.spines[sp].set_visible(False)
                if r_ < 7: b.set_xticklabels([])
            for pdf_fig in (fig,): pdf.savefig(pdf_fig); plt.close(pdf_fig)
    print('gallery', name, len(pick), 'thumbnails')


if __name__ == '__main__':
    T = main()
    for n in GAL_SESS:
        if n in T: gallery(n, T)
