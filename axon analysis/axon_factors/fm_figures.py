"""python fm_figures.py [--domain band] [session ...]   -> axon analysis/plots/axon_factors/<session>_factors_<domain>.png
 a  held-out R2 versus number of factors (leave-block-out and leave-ROI-out)     b  what each factor correlates with (behaviours, nuisance signals)
 c  loading maps (varimax) of the first three factors: each accepted ROI coloured by its loading    d  factor time courses (first 12 min) with face motion
Style: reference palette of the other figure code (fixed order, light surface); diverging map for signed loadings (two hues, neutral midpoint)."""
import os, sys, pickle, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import cm as mcm
import fm_fit as fm
import ag_data as ad
from ag_data import xc

SURFACE, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8984'; SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
plt.rcParams.update({'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE, 'text.color': INK, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.edgecolor': INK2, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlesize': 9, 'axes.titleweight': 'regular', 'axes.titlelocation': 'left', 'font.size': 8,
                     'legend.frameon': False, 'axes.grid': True, 'grid.color': '#e6e5e1', 'grid.linewidth': 0.5})
OUT = os.path.join(ad.REPO, 'axon analysis', 'plots', 'axon_factors')


def fig_session(name, domain='band'):
    with open(os.path.join(fm.CACHE, '%s_fm_%s.pkl' % (name, domain)), 'rb') as fh: o = pickle.load(fh)
    U = ad.load_universe(name, verbose=False); acc = np.where(U.src == 0)[0]; k = o['kstar']; os.makedirs(OUT, exist_ok=True)
    ops = np.load(os.path.join(U.plane_dir, 'ops.npy'), allow_pickle=True).item(); img = np.asarray(ops['meanImg'], float); lo, hi = np.percentile(img, [1, 99.5]); base = np.clip((img - lo) / (hi - lo), 0, 1) ** 0.7
    fig = plt.figure(figsize=(11, 10.5)); gs = fig.add_gridspec(3, 3, height_ratios=[1, 0.62, 0.9], hspace=0.5, wspace=0.28)
    ax = fig.add_subplot(gs[0, 0]); ks = np.arange(len(o['med_R2_sig'])); ax.plot(ks, o['med_R2_sig'], '-o', color=SERIES[0], ms=3.5, lw=1.2, label='signal-like ROIs'); ax.plot(ks, o['med_R2_all'], '-o', color=SERIES[1], ms=3.5, lw=1.2, label='all accepted ROIs')
    ax.axvline(k, color=MUTED, ls=':', lw=0.8); ax.text(k + 0.15, 0.6, 'k* = %d' % k, fontsize=7, color=INK2); ax.set_xlabel('number of shared factors'); ax.set_ylabel('held-out R2 (median ROI)'); ax.set_title('a  how many shared signals?'); ax.legend(fontsize=7, loc='center right')
    ax = fig.add_subplot(gs[0, 1:]); cols = list(o['corr_beh']) + list(o['corr_nuis']); lab = {'mean neuropil': 'neuropil', 'rejected-ROI PC1': 'rejected PC1', 'rigid shift': 'rigid shift'}; M = np.array([[(o['corr_beh'] if c in o['corr_beh'] else o['corr_nuis'])[c][j] for c in cols] + [o['ridge_r2_beh'][j], o['share'][j]] for j in range(k)])
    im = ax.imshow(M[:, :len(cols)], cmap='coolwarm', vmin=-1, vmax=1, aspect='auto'); ax.set_xticks(range(len(cols) + 2)); ax.set_xticklabels([lab.get(c, c) for c in cols] + ['behaviour ridge R2', 'variance share'], fontsize=7, rotation=20); ax.set_yticks(range(k)); ax.set_yticklabels(['factor %d' % (j + 1) for j in range(k)]); ax.grid(False)
    for j in range(k):
        for c in range(len(cols) + 2): ax.text(c, j, '%+.2f' % M[j, c] if c < len(cols) else '%.2f' % M[j, c], ha='center', va='center', fontsize=7, color=INK)
    ax.imshow(np.full((k, 2), 0.5), extent=(len(cols) - 0.5, len(cols) + 1.5, k - 0.5, -0.5), cmap='Greys', vmin=0, vmax=3, aspect='auto', alpha=0.15); ax.set_xlim(-0.5, len(cols) + 1.5)
    ax.set_title('b  correlation of each factor with behaviour and nuisance signals; behaviour ridge R2 (held-out); variance share', fontsize=8)
    for j in range(min(k, 3)):
        a = fig.add_subplot(gs[1, j]); a.grid(False); rgb = np.dstack([base] * 3) * 0.5 + 0.08; vmax = np.abs(o['Lv'][:, j]).max()
        for i, r in enumerate(acc):
            c = np.array(mcm.coolwarm((o['Lv'][i, j] / vmax + 1) / 2)[:3]); y, x = U.pix[r][:, 0], U.pix[r][:, 1]; rgb[y, x] = c
        a.imshow(rgb, interpolation='nearest'); a.set_xticks([]); a.set_yticks([]); a.set_title('c  loadings of factor %d (varimax), +-%.2f' % (j + 1, vmax))
    ax = fig.add_subplot(gs[2, :]); pos = np.where(fm.get_matrix(U, domain)[1])[0]; t = (U.t[pos] - U.t[0]) / 60; m = t < 12
    for j in range(min(k, 3)): ax.plot(t[m], o['Fv'][j][m] / o['Fv'][j].std() + 6 * (min(k, 3) - 1 - j), color=SERIES[j], lw=0.8, label='factor %d' % (j + 1))
    if 'face' in U.S.beh:
        f = xc.bandpass(U.S.beh['face'], 1, 30)[pos] if domain == 'band' else U.S.beh['face'][pos]; f = (f - np.nanmean(f)) / np.nanstd(f); ax.plot(t[m], f[m] + 6 * min(k, 3), color=MUTED, lw=0.8, label='face motion')
    ax.set_xlabel('time (min)'); ax.set_yticks([]); ax.set_title('d  factor time courses (z, offset), first 12 min'); ax.legend(fontsize=7, ncol=4, loc='upper right')
    fig.suptitle('%s: factor model of %d accepted ROIs (%s domain)' % (name, o['n_roi'], domain), x=0.01, ha='left', fontsize=10); fig.savefig(os.path.join(OUT, '%s_factors_%s.png' % (name, domain)), dpi=150, bbox_inches='tight'); plt.close(fig)


if __name__ == '__main__':
    args = sys.argv[1:]; domain = 'band'
    if '--domain' in args: domain = args[args.index('--domain') + 1]; args = [a for a in args if a not in ('--domain', domain)]
    for nm in (args or ['AnimalRB24_261021_1630', 'AnimalRB22_261027_1138', 'AnimalRB23_261023_1642', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520']): fig_session(nm, domain)
