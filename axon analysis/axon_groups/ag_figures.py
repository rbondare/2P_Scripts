"""Axon grouping layer, figures. python ag_figures.py <session> [...]   (needs the outputs of ag_run.py; --validate outputs for the report card)
  <session>_map.png     FOV: units (colour = unit, solid = tier A averaged unit, dashed = tier B association only), attached rejected ROIs (fill strength = confidence), top-20 noise blobs, orphans
  <session>_units.png   member traces and unit mean of the multi-member units
  <session>_report.png  report card: residual r vs gap, held-out SNR gain, simulation operating curve, inflation of PC1 share
Style: reference data-viz palette of suite2p_analysis/merge_suite2p_rois.py (fixed order, light surface), thin marks, direct labels, grid recessive.
"""
import os, sys, pickle, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import ag_data as ad, ag_pairs as ap, ag_group as ag

SURFACE, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a8984'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
plt.rcParams.update({'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE, 'text.color': INK, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.edgecolor': INK2, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlesize': 9, 'axes.titleweight': 'regular', 'axes.titlelocation': 'left',
                     'font.size': 8, 'legend.frameon': False, 'axes.grid': True, 'grid.color': '#e6e5e1', 'grid.linewidth': 0.5})
OUT = ad.OUT


def _load(name):
    with open(os.path.join(ad.CACHE, '%s_groups.pkl' % name), 'rb') as fh: return pickle.load(fh)


def _val(name, level='pc1'):
    f = os.path.join(ad.CACHE, '%s_validation_%s.pkl' % (name, level))
    if not os.path.exists(f): return None
    with open(f, 'rb') as fh: return pickle.load(fh)


def fig_map(name):
    U = ad.load_universe(name); R = _load(name); ops = np.load(os.path.join(U.plane_dir, 'ops.npy'), allow_pickle=True).item(); img = np.asarray(ops['meanImg'], float)
    lo, hi = np.percentile(img, [1, 99.5]); base = np.clip((img - lo) / (hi - lo), 0, 1) ** 0.7; rgb = np.dstack([base] * 3) * 0.55 + 0.1
    multi = [k for k, g in enumerate(R['units']) if len(g) >= 2]; col = {k: np.array(matplotlib.colors.to_rgb(SERIES[n % len(SERIES)])) for n, k in enumerate(multi)}
    fig, ax = plt.subplots(figsize=(11, 11 * img.shape[0] / img.shape[1] + 1.2)); ax.grid(False)
    def paint(rows, color, a):
        for i in rows: y, x = U.pix[i][:, 0], U.pix[i][:, 1]; rgb[y, x] = (1 - a) * rgb[y, x] + a * color
    # accepted singletons: neutral; units: colour
    for i in np.where(U.src == 0)[0]:
        if R['unit_of'][i] not in col: paint([i], np.array([0.85, 0.85, 0.82]), 0.5)
    for k in multi: paint(R['units'][k], col[k], 0.75)
    conf_a = {'high': 0.7, 'medium': 0.45, 'low': 0.2}
    for j, c in R['attach']['conf'].items():
        if c in conf_a and R['attach']['best'].get(j) and R['attach']['best'][j]['unit'] in col and U.src[j] == 1: paint([j], col[R['attach']['best'][j]['unit']], conf_a[c])
    ax.imshow(rgb, interpolation='nearest'); ax.set_xticks([]); ax.set_yticks([])
    def outline(rows, color, ls, lw=0.8):
        for i in rows:
            m = np.zeros(img.shape, bool); m[U.pix[i][:, 0], U.pix[i][:, 1]] = True; ax.contour(m, levels=[0.5], colors=[color], linestyles=ls, linewidths=lw)
    for k in multi:
        solid = R['unit_info'][k].get('avg_ok', False); outline(R['units'][k], tuple(col[k]), 'solid' if solid else 'dashed')
    outline(list(R['noise_top']), '#ffffff', 'solid', 0.9)
    orph = [j for j, o in R['attach']['orphan'].items() if o]; outline(orph, '#eda100', 'solid', 0.8)
    n_conf = {c: sum(1 for j, v in R['attach']['conf'].items() if U.src[j] == 1 and v == c) for c in ('high', 'medium', 'low')}
    handles = [Patch(color=SERIES[0], label='unit with >= 2 accepted ROIs (colour = unit; solid outline = tier A, averaged; dashed = tier B, association only)'),
               Patch(color='#d9d9d1', label='accepted ROI without a partner'), Patch(color=SERIES[0], alpha=0.5, label='attached rejected ROI (fill strength = confidence: high %d / medium %d / low %d)' % (n_conf['high'], n_conf['medium'], n_conf['low'])),
               Line2D([0], [0], color='white', lw=1.2, label='top-20 brightness-matched noise blobs'), Line2D([0], [0], color='#eda100', lw=1.2, label='orphan: active rejected ROI, attached to nothing (review candidates)')]
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.01), fontsize=7, ncol=1)
    ax.set_title('%s: %d accepted ROIs -> %d units (%d multi-member), level %s' % (name, (U.src == 0).sum(), len(R['units']), len(multi), R['level']), loc='left')
    fig.savefig(os.path.join(OUT, name + '_map.png'), dpi=170, bbox_inches='tight'); plt.close(fig)


def fig_units(name, max_units=6, win_s=240):
    U = ad.load_universe(name); R = _load(name); multi = [k for k, g in enumerate(R['units']) if len(g) >= 2][:max_units]
    if not multi: return
    fig, axs = plt.subplots(len(multi), 1, figsize=(10, 1.5 * len(multi) + 0.8))
    axs = np.atleast_1d(axs); t = U.t - U.t[0]
    for ax, k in zip(axs, multi):
        g = R['units'][k]; D = U.dff[g]; mean = np.nanmean(D, 0); sd = np.nanstd(mean)
        # window with the largest unit activity
        w = int(win_s / ad.DT); s = np.nanvar(np.lib.stride_tricks.sliding_window_view(np.nan_to_num(mean), w)[::25], axis=1); a = int(np.argmax(s)) * 25
        for r in D: ax.plot(t[a:a + w] - t[a], r[a:a + w], color=MUTED, lw=0.5, alpha=0.8)
        ax.plot(t[a:a + w] - t[a], mean[a:a + w], color=SERIES[0], lw=1.2); d = R['unit_info'][k]; ax.set_xlim(0, win_s)
        ax.text(0.995, 0.92, 'unit %d: %d ROIs, tier %s%s, half-r %s' % (k, len(g), R['tier'][g[0]], (' (' + ','.join(d['why']) + ')') if d['why'] else '', '/'.join('%.2f' % x for x in d['half_r'])), transform=ax.transAxes, ha='right', va='top', fontsize=7, color=INK2)
        ax.set_ylabel('dF/F')
    axs[-1].set_xlabel('time within the window (s); each unit shows its own most active window'); fig.suptitle('%s: members (grey) and unit mean (blue), most active %d s window of each multi-member unit' % (name, win_s), x=0.01, ha='left', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, name + '_units.png'), dpi=150); plt.close(fig)


def fig_report(name, level='pc1'):
    V = _val(name, level)
    if V is None: print('no validation output for', name); return
    U = ad.load_universe(name); P = ap.PairStats(U, level, verbose=False); diag = ap.gap_diagnostic(P)
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.6))
    # (a) residual r vs gap
    ax = axs[0, 0]; kinds = (('acc-acc', SERIES[0], 'accepted - accepted'), ('acc-rej', SERIES[1], 'accepted - rejected'), ('rej-rej', SERIES[2], 'rejected - rejected'))
    for kd, c, lab in kinds:
        rows = [d for d in diag if d['kind'] == kd]; x = np.arange(len(rows))
        if not rows: continue
        ax.plot([d['hi'] for d in rows], [d['med_r'] for d in rows], '-o', color=c, lw=1.2, ms=3.5, label=lab + ' (median)'); ax.plot([d['hi'] for d in rows], [d['p95_r'] for d in rows], '--', color=c, lw=0.8)
    ax.axhline(0, color=INK2, lw=0.6); ax.axhline(ag.R_MIN, color=MUTED, lw=0.6, ls=':'); ax.text(ax.get_xlim()[1], ag.R_MIN, ' r_min', va='bottom', ha='right', fontsize=7, color=MUTED)
    ax.set_xscale('symlog', linthresh=2); ax.set_xlabel('footprint gap, upper edge of bin (um)'); ax.set_ylabel('residual r (after %s removal)' % level); ax.set_title('a  correlation vs distance (dashed = 95th percentile)'); ax.legend(fontsize=7)
    # (b) held-out SNR gain
    ax = axs[0, 1]; rows = [r for ho in V['held_out'] for r in ho['units']]
    if rows:
        x = np.arange(len(rows)); ax.bar(x, [r['gain'] for r in rows], color=SERIES[0], width=0.55, label='unit mean vs members (held-out half)')
        ax.plot(x, [r['ctrl_gain_median'] for r in rows], 'o', color=SERIES[1], ms=4, label='control: n near ungrouped ROIs'); ax.plot(x, [r['ideal'] for r in rows], '_', color=INK, ms=10, mew=1.2, label='ideal sqrt(n)')
        ax.set_xticks(x); ax.set_xticklabels(['n=%d\n%s' % (r['n'], r['tier']) for r in rows], fontsize=6)
    ax.set_ylabel('SNR gain'); ax.set_title('b  what averaging buys (groups from one half, SNR on the other)'); ax.legend(fontsize=7)
    # (c) simulation operating curve
    ax = axs[1, 0]; H = V['hetero']; mixes = sorted(set(h['mix'] for h in H))
    for n, mx in enumerate(mixes):
        rr = [h for h in H if h['mix'] == mx]; ax.plot([h['rho'] for h in rr], [h['false_merge'] for h in rr], '-o', color=SERIES[n], lw=1.2, ms=3.5, label='%s fragments: merged' % mx)
        ax.plot([h['rho'] for h in rr], [h['flagged'] for h in rr], '--s', color=SERIES[n], lw=0.9, ms=3, mfc='none', label='%s: flagged heterogeneous' % mx)
    ax.set_xlabel('true correlation of the two axons\' signals (simulated, real noise)'); ax.set_ylabel('fraction of runs'); ax.set_ylim(-0.03, 1.05); ax.invert_xaxis(); ax.set_title('c  do two DIFFERENT axons get merged, and caught?'); ax.legend(fontsize=6.5, ncol=1)
    # (d) inflation
    ax = axs[1, 1]; I = [ho['inflation'] for ho in V['held_out']]
    keys = [('pc1_roi_B', 'ROIs (no grouping)'), ('pc1_units_in_sample', 'units, in-sample'), ('pc1_units_held_out', 'units, held-out'), ('pc1_random_merge_held_out', 'random merges, held-out')]
    vals = [np.nanmean([i[k] for i in I]) for k, _ in keys]; ax.bar(np.arange(4), vals, color=[MUTED, SERIES[0], SERIES[0], SERIES[1]], width=0.55)
    for n, v in enumerate(vals): ax.text(n, v, '%.2f' % v, ha='center', va='bottom', fontsize=7)
    ax.set_xticks(np.arange(4)); ax.set_xticklabels([l for _, l in keys], fontsize=7); ax.set_ylabel('PC1 share of variance'); ax.set_title('d  does grouping manufacture shared signal? (accepted ROIs, band-passed)')
    n0 = V['null']; txt = '\n'.join('%s: %d edges / %d pairs' % (k, v['edges'], v['candidate_pairs']) for k, v in n0.items())
    fig.text(0.01, 0.005, 'null runs (edges with FDR < 0.01 and r >= %.1f among accepted pairs within %d um):  ' % (ag.R_MIN, ag.D_BUNDLE) + ' | '.join('%s: %d of %d' % (k, v['edges'], v['candidate_pairs']) for k, v in n0.items()), fontsize=7, color=INK2)
    fig.suptitle('%s  grouping report card (level %s)' % (name, level), x=0.01, ha='left', fontsize=10); fig.tight_layout(rect=(0, 0.02, 1, 0.97)); fig.savefig(os.path.join(OUT, name + '_report.png'), dpi=150); plt.close(fig)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for nm in [a for a in sys.argv[1:] if not a.startswith('--')]:
        fig_map(nm); fig_units(nm); fig_report(nm)
