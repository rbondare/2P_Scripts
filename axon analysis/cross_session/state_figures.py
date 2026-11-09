"""Figures of cross-session analysis 1 (state signal). usage: python state_figures.py   (reads state_part*.pkl / state_signal_results.pkl in plots/cross_session, writes S*.png there)

Colour = what the thing is, never its rank: signal-like ROIs blue, noise-like ROIs (accepted by the experimenter on purpose, no visible axon signal) amber, rejected suite2p ROIs grey;
animals RB24 / RB23 / RB22 in the session labels' colours (blue / orange / teal). Every figure has the same session order (animal, date).
"""
import os, glob, pickle, sys, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
FIG = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session'
plt.rcParams.update({'font.family': 'Arial' if any('Arial' in f.name for f in matplotlib.font_manager.fontManager.ttflist) else 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.75, 'axes.spines.top': False, 'axes.spines.right': False,
                     'xtick.direction': 'out', 'ytick.direction': 'out', 'figure.dpi': 110, 'savefig.dpi': 200, 'axes.titlesize': 8.5, 'axes.titleweight': 'normal'})
C_SIG, C_NOISE, C_REJ, INK, MUTED = '#1F5FBF', '#D98A00', '#8A8F98', '#222222', '#6B7280'
C_ANIMAL = {'RB24': '#1F5FBF', 'RB23': '#D1495B', 'RB22': '#1B998B'}
BEH_LAB = {'pupil': 'pupil area', 'dpupil': 'pupil dilation rate', 'face': 'face motion', 'whisk': 'whisking', 'eyelid': 'eyelid opening', 'loco': 'locomotion (ball)', 'shift': 'brain shift'}


def load():
    res = {}
    for f in sorted(glob.glob(os.path.join(FIG, 'state_part*.pkl'))): res.update(pickle.load(open(f, 'rb')))
    if not res and os.path.exists(os.path.join(FIG, 'state_signal_results.pkl')): res = pickle.load(open(os.path.join(FIG, 'state_signal_results.pkl'), 'rb'))
    return res


def short(nm): a, d, t = nm.replace('Animal', '').split('_'); return '%s  %s-%s  %s' % (a, d[2:4], d[4:6], t)
def animal(nm): return nm.replace('Animal', '')[:4]
def order(res): return sorted(res, key=lambda n: (animal(n)[::-1], n))   # RB22, RB23, RB24 blocks; by date inside


def save(fig, name): fig.savefig(os.path.join(FIG, name), bbox_inches='tight', facecolor='white'); plt.close(fig); print('saved', name)


def ylabels(ax, names, fs=7):
    ax.set_yticks(range(len(names))); ax.set_yticklabels([short(n) for n in names], fontsize=fs)
    for t, n in zip(ax.get_yticklabels(), names): t.set_color(C_ANIMAL[animal(n)])
    ax.set_ylim(len(names) - 0.4, -0.6)


def fig_overview(res, use, excl):
    fig, axs = plt.subplots(1, 6, figsize=(13, 0.34 * len(use) + 1.6), sharey=True, gridspec_kw=dict(width_ratios=[1.3, 0.8, 1, 1, 1, 1], wspace=0.12))
    y = np.arange(len(use))
    a = axs[0]; ylabels(a, use)
    ns = np.array([res[n]['n_sig'] for n in use]); nn = np.array([res[n]['n_noise'] for n in use]); npl = np.array([res[n]['n_neuropil'] for n in use])
    a.barh(y, ns, color=C_SIG, height=0.7, label='signal-like'); a.barh(y, nn, left=ns, color=C_NOISE, height=0.7, label='noise-like'); a.barh(y, npl, left=ns + nn, color=C_REJ, height=0.7, label='Fc<=0 (neuropil-limited)')
    a.set_xlabel('accepted ROIs'); a.legend(frameon=False, fontsize=6.5, loc='lower right'); a.set_title('ROIs')
    axs[1].barh(y, [res[n]['dur_min'] for n in use], color=MUTED, height=0.7); axs[1].set_xlabel('minutes'); axs[1].set_title('recording')
    for ax, key, ttl, lim in ((axs[2], 'pc1', 'PC1 share of variance', (0, 1)), (axs[3], 'same_sign', 'same-sign fraction (PC1)', (0.4, 1.02)), (axs[4], 'r_med', 'median pairwise r', (-0.1, 0.8))):
        for i, n in enumerate(use):
            c = res[n]['coh']; v = c['sig'][key] if c['sig'] else np.nan
            rj = [x[key] for x in c['rej'] if x]
            if rj: ax.plot([min(rj), max(rj)], [i, i], color=C_REJ, lw=3, alpha=0.5, solid_capstyle='round'); ax.plot(np.mean(rj), i, 'o', mfc='white', mec=C_REJ, ms=4)
            ax.plot(v, i, 'o', color=C_SIG, ms=5)
            if c['noise']: ax.plot(c['noise'][key], i, 's', color=C_NOISE, ms=3.5)
        ax.set_xlim(lim); ax.set_title(ttl); ax.axvline(0.5 if key == 'same_sign' else 0, color=MUTED, lw=0.5, ls=':')
    axs[4].plot([], [], 'o', color=C_SIG, label='signal-like ROIs'); axs[4].plot([], [], 's', color=C_NOISE, ms=3.5, label='noise-like ROIs'); axs[4].plot([], [], 'o', mfc='white', mec=C_REJ, label='random sets of rejected ROIs\n(same size; range of 10 sets)')
    axs[4].legend(frameon=False, fontsize=6.5, loc='lower right')
    ax = axs[5]; ax.axis('off')
    ax.text(0, 0.98, 'how to read', fontsize=8, fontweight='bold', va='top', transform=ax.transAxes)
    ax.text(0, 0.9, 'PC1 share: how much of the ROI-to-ROI variance\nis ONE shared signal. Same-sign: fraction of ROIs with\nthe sign of the PC1 loading. Rejected sets have the\nsame number of ROIs; if the signal were only optics /\nneuropil they would look alike.\n\nExcluded: ' + ('; '.join('%s (%s)' % (short(n), res[n]['why']) for n in excl) if excl else 'none'),
            fontsize=6.8, va='top', transform=ax.transAxes, color=INK, wrap=True)
    fig.suptitle('1 | Is there ONE shared signal across the axons, and does it exceed what random non-accepted ROIs show?  (band-passed 1-30 s, best dF/F)', x=0.01, ha='left', fontsize=9.5)
    save(fig, 'S1_session_overview_coherence.png')


def fig_heatmaps(res, use):
    behs = ['pupil', 'face', 'whisk', 'eyelid', 'loco', 'shift']
    def mat(kind):
        M = np.full((len(use), len(behs)), np.nan); S = np.zeros_like(M, bool); L = np.full_like(M, np.nan)
        for i, n in enumerate(use):
            for j, b in enumerate(behs):
                d = res[n]['xc'].get(b)
                if d is None: continue
                if kind == 'sig': M[i, j] = d['peak_r']; S[i, j] = d['p'] < 0.01; L[i, j] = d['peak_lag_s']
                elif kind == 'rej' and d['rej_peaks']: M[i, j] = np.mean([p[1] for p in d['rej_peaks']])
                elif kind == 'noise' and d['r_noise'] is not None: k = np.argmax(np.abs(d['r_noise'])); M[i, j] = d['r_noise'][k]
        return M, S, L
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 0.33 * len(use) + 1.8), sharey=True, gridspec_kw=dict(wspace=0.08))
    for ax, kind, ttl in ((axs[0], 'sig', 'signal-like ROIs (population mean)'), (axs[1], 'noise', 'noise-like ROIs left in by the experimenter'), (axs[2], 'rej', 'random sets of rejected ROIs')):
        M, S, L = mat(kind)
        im = ax.imshow(np.ma.masked_invalid(M), cmap='RdBu_r', vmin=-0.8, vmax=0.8, aspect='auto')
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if np.isfinite(M[i, j]): ax.text(j, i, ('%.2f' % M[i, j]).replace('0.', '.') + ('*' if S[i, j] else ''), ha='center', va='center', fontsize=6, color='white' if abs(M[i, j]) > 0.45 else INK)
        ax.set_xticks(range(len(behs))); ax.set_xticklabels([BEH_LAB[b] for b in behs], rotation=45, ha='right', fontsize=7); ax.set_title(ttl, fontsize=8)
        ax.set_facecolor('#EEEEEE')
    ylabels(axs[0], use)
    cb = fig.colorbar(im, ax=axs, fraction=0.012, pad=0.01); cb.set_label('peak cross-correlation r (|lag| <= 10 s)', fontsize=7)
    fig.suptitle('2 | Coupling of the population signal to behaviour (peak r within +-10 s; * = beyond all lags >= 120 s away, p<0.01; grey = not available / ball implausible)', x=0.01, ha='left', fontsize=9.3)
    save(fig, 'S2_coupling_heatmaps.png')


def fig_leadlag(res, use):
    behs = ['pupil', 'face', 'whisk', 'loco']
    fig, axs = plt.subplots(2, len(behs), figsize=(3.9 * len(behs), 6.6), gridspec_kw=dict(height_ratios=[1.0, 1.25], hspace=0.38, wspace=0.85))
    for j, b in enumerate(behs):
        ax = axs[0, j]; curves = []
        for n in use:
            d = res[n]['xc'].get(b)
            if d is None: continue
            ax.plot(d['lags_s'], d['r'], color=C_ANIMAL[animal(n)], lw=0.7, alpha=0.55); curves.append(d['r']); lag_ref = d['lags_s']
        if curves: ax.plot(lag_ref, np.median(np.stack(curves), 0), color=INK, lw=1.8)
        ax.axvline(0, color=MUTED, lw=0.6, ls=':'); ax.axhline(0, color=MUTED, lw=0.6, ls=':'); ax.set_xlim(-10, 10); ax.set_title('%s  (n = %d sessions)' % (BEH_LAB[b], len(curves)), fontsize=8)
        ax.set_xlabel('lag (s)   > 0: behaviour later than dF/F'); ax.set_ylabel('cross-correlation r') if j == 0 else None
        ax = axs[1, j]; rows = [(n, res[n]['xc'][b]) for n in use if b in res[n]['xc'] and res[n]['xc'][b]['p'] < 0.05]
        rows.sort(key=lambda t: t[1]['peak_lag_s'])
        for i, (n, d) in enumerate(rows):
            ax.plot(d['lag_ci'], [i, i], color=C_ANIMAL[animal(n)], lw=1.2, alpha=0.8); ax.plot(d['peak_lag_s'], i, 'o', color=C_ANIMAL[animal(n)], ms=4)
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([short(n) for n, _ in rows], fontsize=5.8); ax.set_ylim(len(rows) - 0.4, -0.8)
        for t, (n, _) in zip(ax.get_yticklabels(), rows): t.set_color(C_ANIMAL[animal(n)])
        ax.axvline(0, color=MUTED, lw=0.6, ls=':'); ax.set_xlim(-10, 10); ax.set_xlabel('peak lag (s), 95% CI (block bootstrap)')
        if rows: ax.set_title('median peak lag %+.1f s (sessions with p<0.05)' % np.median([d['peak_lag_s'] for _, d in rows]), fontsize=7.3)
    fig.suptitle('3 | Lead / lag between the population signal and behaviour.  Top: r(lag) per session (colour = animal), black = median.  Bottom: peak lag with CI.\n'
                 'Face motion, whisking and locomotion are synchronous with the population signal (lag ~ 0); the pupil follows by ~1 s, which is the slow pupil kinetics, not proof that the signal leads. '
                 '(Pupil dilation RATE is left out: its curve is the derivative of the pupil curve, its peak only marks the steepest rise.)', x=0.01, ha='left', fontsize=9, y=1.0)
    save(fig, 'S3_lead_lag.png')


def fig_sign(res, use):
    behs = ['pupil', 'face', 'whisk', 'loco']
    fig, axs = plt.subplots(1, len(behs), figsize=(3.0 * len(behs), 0.33 * len(use) + 1.6), sharey=True, gridspec_kw=dict(wspace=0.1))
    for ax, b in zip(axs, behs):
        for i, n in enumerate(use):
            r = res[n]['roi'].get(b)
            if not r: continue
            vals = {k: (r[k]['pos'] - r[k]['neg']) / max(r[k]['n'], 1) for k in r}
            xs = [v for v in vals.values()]
            if len(xs) > 1: ax.plot([min(xs), max(xs)], [i, i], color='#CCCCCC', lw=1.5, zorder=1)
            if 'rej' in vals: ax.plot(vals['rej'], i, 'o', mfc='white', mec=C_REJ, ms=5, zorder=2)
            if 'noise' in vals: ax.plot(vals['noise'], i, 's', color=C_NOISE, ms=4, zorder=3)
            if 'sig' in vals: ax.plot(vals['sig'], i, 'o', color=C_SIG, ms=5, zorder=4)
        ax.axvline(0, color=MUTED, lw=0.6, ls=':'); ax.set_xlim(-1, 1); ax.set_title(BEH_LAB[b]); ax.set_xlabel('(significant + minus significant -) / ROIs')
    ylabels(axs[0], use)
    axs[-1].plot([], [], 'o', color=C_SIG, label='signal-like'); axs[-1].plot([], [], 's', color=C_NOISE, ms=4, label='noise-like'); axs[-1].plot([], [], 'o', mfc='white', mec=C_REJ, label='rejected (<=200 ROIs)'); axs[-1].legend(frameon=False, fontsize=6.5, loc='lower left')
    fig.suptitle('4 | Sign test per ROI at the population lag: signal-like ROIs should all go the SAME way (+1 or -1); mixed-sign coupling (0) is the signature of motion / neuropil artefacts', x=0.01, ha='left', fontsize=9.2)
    save(fig, 'S4_sign_control.png')


def fig_variance(res, use):
    use2 = [n for n in use if 'ridge' in res[n]]
    fig, axs = plt.subplots(1, 2, figsize=(11, 0.34 * len(use2) + 1.7), sharey=True, gridspec_kw=dict(width_ratios=[1.0, 1.15], wspace=0.08))
    ax = axs[0]
    for i, n in enumerate(use2):
        g = res[n]['ridge']; rj = g.get('rej', [])
        if len(rj): ax.plot([np.nanmin(rj), np.nanmax(rj)], [i, i], color=C_REJ, lw=3, alpha=0.5, solid_capstyle='round'); ax.plot(np.nanmean(rj), i, 'o', mfc='white', mec=C_REJ, ms=4)
        if 'noise' in g: ax.plot(g['noise'], i, 's', color=C_NOISE, ms=4)
        ax.plot(g['sig']['full'], i, 'o', color=C_SIG, ms=5.5)
    ax.axvline(0, color=MUTED, lw=0.6, ls=':'); ax.set_xlabel('held-out R2 (10 contiguous folds)'); ax.set_xlim(-0.2, 0.85); ax.set_title('behaviour explains the population signal')
    ax.plot([], [], 'o', color=C_SIG, label='signal-like mean'); ax.plot([], [], 's', color=C_NOISE, ms=4, label='noise-like mean'); ax.plot([], [], 'o', mfc='white', mec=C_REJ, label='random sets of rejected ROIs (range)'); ax.legend(frameon=False, fontsize=6.5, loc='lower right')
    ylabels(ax, use2)
    ax = axs[1]; fams = ['pupil', 'face', 'whisk', 'eyelid', 'loco']; cols = {'pupil': '#C0392B', 'face': '#E08E0B', 'whisk': '#2E8B57', 'eyelid': '#7D5BA6', 'loco': '#2C5AA0'}
    for i, n in enumerate(use2):
        u = res[n]['ridge']['sig']['unique']
        for k, f in enumerate(fams):
            if f in u and np.isfinite(u[f]): ax.barh(i + (k - 2) * 0.15, u[f], height=0.14, color=cols[f], label=BEH_LAB[f] if i == 0 else None)
    ax.axvline(0, color=MUTED, lw=0.6); ax.set_xlabel('unique R2 (full model minus model without the family)'); ax.set_title('what each behaviour adds on top of the others'); ax.legend(frameon=False, fontsize=6.5, loc='lower right')
    fig.suptitle('5 | How much of the shared signal does behaviour explain?  Blue = axon ROIs; if amber / grey (non-axon) are as high, the explained part is common-mode, not axon calcium', x=0.01, ha='left', fontsize=9.2)
    save(fig, 'S5_variance_explained.png')


def fig_states(res, use):
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 0.34 * len(use) + 1.5), sharey=True, gridspec_kw=dict(wspace=0.1))
    for ax, key, ttl in ((axs[0], 'face_active_minus_quiet', 'active minus quiet (face-motion terciles)'), (axs[1], 'run_minus_still', 'running minus still (ball, where plausible)')):
        for i, n in enumerate(use):
            s = res[n]['states'].get(key)
            if s is None: continue
            est, ci = s; ax.plot(ci, [i, i], color=C_SIG, lw=1.5); ax.plot(est, i, 'o', color=C_SIG, ms=5)
            rj = res[n]['states'].get(key + '_rej')
            if rj: ax.plot(rj, [i] * len(rj), '.', color=C_REJ, ms=4, alpha=0.7)
        ax.axvline(0, color=MUTED, lw=0.6, ls=':'); ax.set_title(ttl); ax.set_xlabel('mean population dF/F difference (dF/F units)')
    ylabels(axs[0], use)
    fig.suptitle('7 | State dependence of the population signal (blue = signal-like, block-bootstrap 95% CI; grey dots = random rejected sets, face-motion contrast only)', x=0.01, ha='left', fontsize=9.2)
    save(fig, 'S7_states.png')


def fig_examples(res, use, k=4):
    cand = sorted([n for n in use if 'face' in res[n]['xc']], key=lambda n: -abs(res[n]['xc']['face']['peak_r']))[:k]
    fig, axs = plt.subplots(len(cand), 1, figsize=(11, 2.5 * len(cand)), sharex=False, gridspec_kw=dict(hspace=0.75))
    axs = np.atleast_1d(axs)
    for ax, n in zip(axs, cand):
        r = res[n]; t = r['t']; P = r['P']; bg = r['beh_grid']
        T0 = t[0] + 0.25 * (t[-1] - t[0]); sel = (t >= T0) & (t < T0 + 300)
        z = lambda v: (v - np.nanmean(v)) / np.nanstd(v)
        ax.plot(t[sel] - T0, z(P)[sel] * 1.0 + 6, color=C_SIG, lw=0.9)
        for off, b, c, lab in ((3, 'face', '#E08E0B', 'face motion'), (0, 'whisk', '#2E8B57', 'whisking'), (-3, 'pupil', '#C0392B', 'pupil')):
            if b in bg: ax.plot(t[sel] - T0, z(bg[b])[sel] + off, color=c, lw=0.7)
            ax.text(301, off, lab, fontsize=6.5, color=c, va='center')
        if 'loco' in bg and r['loco_ok']: ax.plot(t[sel] - T0, z(bg['loco'])[sel] - 6, color='#2C5AA0', lw=0.7); ax.text(301, -6, 'locomotion', fontsize=6.5, color='#2C5AA0', va='center')
        ax.text(301, 6, 'axon population', fontsize=6.5, color=C_SIG, va='center'); ax.set_yticks([]); ax.set_xlim(0, 300); ax.spines['left'].set_visible(False)
        ax.set_title('%s | face r %+.2f, pupil r %+.2f | 5 min excerpt, z-scored' % (short(n), r['xc']['face']['peak_r'], r['xc'].get('pupil', {}).get('peak_r', np.nan)), fontsize=8, loc='left', color=C_ANIMAL[animal(n)])
    axs[-1].set_xlabel('time (s)')
    fig.suptitle('6 | Examples: the sessions with the strongest face-motion coupling (population mean of the signal-like ROIs vs behaviour)', x=0.01, ha='left', fontsize=9.3, y=1.0)
    save(fig, 'S6_examples.png')


if __name__ == '__main__':
    res = load(); names = order(res)
    use = [n for n in names if res[n].get('usable')]; excl = [n for n in names if not res[n].get('usable')]
    print('usable:', len(use), '| excluded:', [(short(n), res[n].get('why')) for n in excl])
    fig_overview(res, use, excl); fig_heatmaps(res, use); fig_leadlag(res, use); fig_sign(res, use); fig_variance(res, use); fig_examples(res, use); fig_states(res, use)
