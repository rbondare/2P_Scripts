"""Figure 4 - heterogeneity across ROIs and spatial structure.
a  FOV coloured by single-ROI coupling to face motion (1-plane 1630, 5-plane 1651)   b  ROI x ROI correlation (1651), ordered by PC1 loading   c  example ROIs from the two ends of the coupling distribution
d  variance explained per principal component   e  single-ROI coupling by ROI class   f  what PC2 carries   g  do nearby ROIs share more signal?
usage: python fig4_heterogeneity.py -> plots/summary_figures/Fig4_heterogeneity.pdf/.png ; ledger entries
"""
import os, pickle, numpy as np
from fig_common import *
import xs_common as xc
from scipy.ndimage import binary_erosion
from matplotlib.colors import TwoSlopeNorm

Z = r'Z:/Group Members/Rima/DATA_2P/'
EX1, EX5 = 'AnimalRB24_261021_1630', 'AnimalRB24_261104_1651'


def vf(r, k=6): v = np.asarray(r['var_frac'], float)[:k]; return np.concatenate([v, np.full(k - len(v), np.nan)])


def load_roi(n):
    f = os.path.join(CACHE, 'roi_%s.pkl' % n); return pickle.load(open(f, 'rb')) if os.path.exists(f) else None


def fov_value_map(ax, name, R, cmap, norm):
    p = Z + name + '/suite2p/plane0/'; ops = np.load(p + 'ops.npy', allow_pickle=True).item(); stat = np.load(p + 'stat.npy', allow_pickle=True); ic = np.load(p + 'iscell.npy')[:, 0] > 0
    img = ops.get('max_proj'); yr, xr = ops['yrange'], ops['xrange']; sub = img if img is not None else ops['meanImg']; sub = sub.astype(float)
    lo, hi = np.percentile(sub, [5, 99.6]); g = np.clip((sub - lo) / (hi - lo), 0, 1) ** 0.7; rgb = np.dstack([g, g, g]) * 0.5
    acc = np.where(ic)[0]; val = np.full(len(acc), np.nan)
    for grp in ('sig', 'noise'):
        if ('idx_' + grp) in R: val[R['idx_' + grp]] = R['r_' + grp]['face']
    for k, i in enumerate(acc):
        m = np.zeros((ops['Ly'], ops['Lx']), bool); m[stat[i]['ypix'], stat[i]['xpix']] = True; m = m[yr[0]:yr[1], xr[0]:xr[1]]
        if not np.isfinite(val[k]): continue
        c = np.array(cmap(norm(val[k]))[:3]); edge = m & ~binary_erosion(m); rgb[m] = 0.35 * rgb[m] + 0.65 * c; rgb[edge] = c
    ax.imshow(rgb, extent=(xr[0], xr[1], yr[1], yr[0]), interpolation='nearest', aspect='equal'); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    return val


def main():
    st = load_state(); R = {n: load_roi(n) for n in ALL}
    fig = plt.figure(figsize=(W2, 150 * MM))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.78, 1.0], hspace=0.30, left=0.075, right=0.985, top=0.945, bottom=0.095)
    top = gs[0].subgridspec(1, 4, wspace=0.14, width_ratios=[1.35, 1.35, 0.8, 1.35])
    cmap = matplotlib.colormaps['RdBu_r']; norm = TwoSlopeNorm(vmin=-0.9, vcenter=0.0, vmax=0.9)
    a1 = fig.add_subplot(top[0]); a2 = fig.add_subplot(top[1])
    v1 = fov_value_map(a1, EX1, R[EX1], cmap, norm); v5 = fov_value_map(a2, EX5, R[EX5], cmap, norm)
    a1.set_title('%s\n(1 plane)' % short(EX1), fontsize=5.8, loc='left', fontweight='normal'); a2.set_title('%s\n(5 planes summed)' % short(EX5), fontsize=5.8, loc='left', fontweight='normal'); plabel(a1, 'a', dx=-14, dy=12)
    cax = a2.inset_axes([0.15, -0.14, 0.7, 0.05]); cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal'); cb.set_label('single-ROI r with face motion (0 s lag)', fontsize=5.4, labelpad=1); cb.ax.tick_params(labelsize=5, length=1.5)
    # b  correlation matrix
    ab = fig.add_subplot(top[2]); r5 = R[EX5]; ordr = np.argsort(r5['loadings'][:, 0])[::-1]; Cm = r5['C'][np.ix_(ordr, ordr)]
    im = ab.imshow(Cm, cmap='RdBu_r', vmin=-1, vmax=1, interpolation='nearest'); ab.set_xticks([]); ab.set_yticks([]); ab.set_title('ROI \u00d7 ROI r (%d signal-like ROIs,\nordered by PC1 loading)' % len(ordr), fontsize=5.8, loc='left', fontweight='normal'); plabel(ab, 'b', dx=-6, dy=12)
    # c  example ROIs: best and worst coupled signal-like ROIs of 1651
    S = xc.Session(EX5, with_rejected=False, verbose=False); sig = np.where(S.sig)[0]; rf = r5['r_sig']['face']; o = np.argsort(rf); hi = sig[o[-3:][::-1]]; lo = sig[o[:3]]
    ac = fig.add_subplot(top[3]); t = S.t - S.t[0]; sel = (t >= 180) & (t < 360); tt = t[sel] - t[sel][0]; off = 0.0
    face = xc.nsmooth(S.beh['face'], 1.0 / xc.DT)[sel]
    def unit(y): lo_, hi_ = np.nanpercentile(y, [1, 99]); return (y - lo_) / max(hi_ - lo_, 1e-6)
    for grp, ids, col in (('high', hi, C_SIG), ('low', lo, C_REJ)):
        for k in ids:
            y = unit(xc.nsmooth(S.dff[k], 1.0 / xc.DT)[sel]); ac.plot(tt, y + off, color=col, lw=0.7); ac.text(-4, off + 0.45, 'r = %.2f' % rf[list(sig).index(k)], fontsize=4.8, ha='right', va='center', color=col); off -= 1.25
        off -= 0.3
    ac.plot(tt, unit(face) + off, color=C_BEH['face'], lw=0.6); ac.text(-4, off + 0.45, 'face motion', fontsize=4.8, ha='right', va='center'); ac.set_xlim(-34, tt[-1]); ac.set_ylim(off - 1.0, 1.3); ac.axis('off')
    scalebar(ac, tt[-1] - 30, off - 0.55, dx=30, lx='30 s'); ac.set_title('best (blue) / worst (grey) coupled ROIs\n%s, 3 min, own scaling' % short(EX5), fontsize=5.8, loc='left', fontweight='normal'); plabel(ac, 'c', dx=-6, dy=12)

    # ------------------------------------------------ row 2
    bot = gs[1].subgridspec(1, 4, wspace=0.50, width_ratios=[1.0, 1.15, 1.0, 1.25]); ad = fig.add_subplot(bot[0]); ae = fig.add_subplot(bot[1]); af = fig.add_subplot(bot[2]); ag = fig.add_subplot(bot[3])
    # d  scree
    for n in ALL:
        r = R[n]
        if not r or 'var_frac' not in r: continue
        fl = is_flagged(n); ad.plot(np.arange(1, 7), vf(r), color=FAINT if fl else C_SIG, lw=0.5, alpha=1 if fl else 0.5, ls='--' if fl else '-', marker=MARK[animal(n)], ms=2.4, mfc='white' if fl else C_SIG, mec=MUTED if fl else C_SIG, mew=0.4)
        c = st[n]['coh'].get('rej')
        if c and not fl: ad.plot(1.15, np.nanmean([x['pc1'] for x in c if x]), 'o', color=C_REJ, ms=2.6, mew=0)
    ad.plot(np.arange(1, 7), np.nanmean([vf(R[n]) for n in UNFLAGGED if R[n] and 'var_frac' in R[n]], 0), color=C_SIG, lw=1.8); ad.plot([], [], 'o', color=C_REJ, ms=2.6, label='PC1 of random rejected sets'); ad.legend(frameon=False, fontsize=4.9, loc='upper right', handletextpad=0.1)
    ad.set_xlabel('principal component'); ad.set_ylabel('fraction of variance (signal-like ROIs)'); ad.set_ylim(0, 0.85); plabel(ad, 'd', dx=-22, dy=6)
    pc1 = {n: R[n]['var_frac'][0] for n in UNFLAGGED if R[n] and 'var_frac' in R[n]}; h = pooled_unflagged(pc1, n_boot=4000); rj = {n: np.nanmean([x['pc1'] for x in st[n]['coh']['rej'] if x]) for n in UNFLAGGED if st[n]['coh'].get('rej')}; hr = pooled_unflagged(rj, n_boot=4000)
    ledger_add('Fig4', 'pc1_share', 'PC1 share of variance among signal-like ROIs (band-passed 1-30 s)', h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'supported', note='random rejected-ROI sets of the same size: %.2f' % hr['est'], unit='fraction')
    # e  single-ROI r by class (state analysis arrays, r at the population lag)
    rng = np.random.default_rng(5); pos = {'sig': 0, 'noise': 1, 'rej': 2}; cols = {'sig': C_SIG, 'noise': C_NOISE, 'rej': C_REJ}
    for grp in ('sig', 'noise', 'rej'):
        allr = []
        for n in UNFLAGGED:
            r = st[n]['roi'].get('face', {}).get(grp)
            if not r: continue
            allr.append(np.asarray(r['r'])); med = np.median(r['r']); sdot(ae, pos[grp] + rng.uniform(-0.28, 0.28), med, n, color=cols[grp], hollow_color=cols[grp], ms=3.6, zorder=4)
        allr = np.concatenate(allr); kde_x = np.linspace(-0.8, 1.0, 120); from scipy.stats import gaussian_kde
        kd = gaussian_kde(allr, bw_method=0.25)(kde_x); kd = kd / kd.max() * 0.42; ae.fill_betweenx(kde_x, pos[grp] - kd, pos[grp] + kd, color=cols[grp], alpha=0.25, lw=0, zorder=1)
        ae.text(pos[grp], 1.05, '%d ROIs' % len(allr), ha='center', fontsize=4.9, color=MUTED, va='bottom')
    ae.axhline(0, color=MUTED, lw=0.4); ae.set_xticks([0, 1, 2]); ae.set_xticklabels(['signal-like', 'noise-like', 'rejected'], fontsize=5.4); ae.set_ylim(-0.7, 1.05); ae.set_ylabel('single-ROI r with face motion'); plabel(ae, 'e', dx=-22, dy=6)
    ae.text(0.0, 1.08, 'violins: all ROIs of the 13 unflagged sessions; dots: session medians', transform=ae.transAxes, fontsize=4.8, color=MUTED, va='bottom')
    # f  what PC2 carries
    for n in ALL:
        r = R[n]
        if not r or 'pc_r' not in r or 'face' not in r['pc_r']: continue
        a, b = r['pc_r']['face']; sdot(af, abs(a), abs(b), n, color=C_SIG, hollow_color=MUTED, ms=3.8)
    af.plot([0, 1], [0, 1], color=MUTED, lw=0.4, ls=':'); af.set_xlim(0, 1); af.set_ylim(0, 1); af.set_xlabel('|r| of PC1 with face motion'); af.set_ylabel('|r| of PC2 with face motion'); plabel(af, 'f', dx=-22, dy=6)
    p1 = {n: abs(R[n]['pc_r']['face'][0]) for n in UNFLAGGED if R[n] and 'pc_r' in R[n] and 'face' in R[n]['pc_r']}; p2 = {n: abs(R[n]['pc_r']['face'][1]) for n in p1}
    af.text(0.03, 0.96, 'PC2 / PC1 variance:\nmedian %.2f' % np.median([R[n]['pc2_over_pc1'] for n in UNFLAGGED if R[n] and 'pc2_over_pc1' in R[n]]), transform=af.transAxes, fontsize=5.0, color=MUTED, ha='left', va='top')
    h2 = pooled_unflagged(p2, n_boot=3000); ledger_add('Fig4', 'pc2_face', '|r| of PC2 of the signal-like ROIs with face motion', h2['est'], h2['lo'], h2['hi'], h2['n_sessions'], h2['n_animals'], 'ambiguous', note='PC1: %.2f' % pooled_unflagged(p1, n_boot=3000)['est'], unit='|r|')
    # g  spatial structure
    vals = {}
    for n in ALL:
        r = R[n]; sp = r.get('space') if r else None
        if not sp: continue
        x = rng.uniform(-0.18, 0.18); sdot(ag, x, sp['rho'], n, color=C_SIG, hollow_color=MUTED, ms=3.8)
        ag.plot([x - 0.04, x + 0.04], [-sp['null95']] * 2, color=C_REJ, lw=0.4) if False else None
        if not is_flagged(n): vals[n] = sp['rho']
    ns_ = np.mean([R[n]['space']['null95'] for n in ALL if R[n] and R[n].get('space')]); ag.axhspan(-ns_, ns_, color=C_REJ, alpha=0.3, lw=0, zorder=0)
    h3 = pooled_unflagged(vals, n_boot=4000); ag.plot([-0.3, 0.3], [h3['est']] * 2, color=INK, lw=1.4, zorder=5); ag.plot([0, 0], [h3['lo'], h3['hi']], color=INK, lw=0.9, zorder=5)
    k_neg = sum(1 for n in UNFLAGGED if R[n] and R[n].get('space') and R[n]['space']['p'] < 0.01 and R[n]['space']['rho'] < 0)
    ag.axhline(0, color=MUTED, lw=0.4); ag.set_xlim(-0.6, 0.6); ag.set_xticks([]); ag.set_ylabel('Spearman \u03c1: pairwise r vs distance'); ag.set_ylim(-0.9, 0.3); plabel(ag, 'g', dx=-22, dy=6)
    ag.text(0.5, 1.03, '%d of %d unflagged sessions: \u03c1 < 0, p < 0.01' % (k_neg, len(vals)), transform=ag.transAxes, fontsize=5.0, color=MUTED, ha='center', va='bottom'); ag.text(0.97, 0.03, 'grey: mean 95% range of\nROI-position permutations', transform=ag.transAxes, fontsize=4.8, color=MUTED, ha='right', va='bottom')
    ledger_add('Fig4', 'space_rho', 'Spearman rho between pairwise ROI correlation and ROI distance', h3['est'], h3['lo'], h3['hi'], h3['n_sessions'], h3['n_animals'], 'ambiguous', note='%d/%d sessions p<0.01; cannot separate nearby fragments of one axon / local neuropil from real topography' % (k_neg, len(vals)), unit='rho')
    fig.text(0.075, 0.006, 'Signal-like ROIs: V3 axon criterion; band-passed 1-30 s best dF/F on a 5 Hz grid; sessions with < 8 signal-like ROIs have no spatial test. Marker shape = animal, filled = unflagged, hollow = flagged (not pooled);\n'
             'black bars = mean of per-animal means ± 95% hierarchical bootstrap CI. Nearby ROIs may be fragments of one axon or share local neuropil / z-motion: panel g shows a spatial relationship, not a topographic map.', fontsize=4.9, color=MUTED)
    savefig(fig, 'Fig4_heterogeneity')


if __name__ == '__main__':
    main()
