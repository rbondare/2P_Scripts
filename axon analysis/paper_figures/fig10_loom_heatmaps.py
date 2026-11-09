"""Figure 10 - looming discs: heatmaps, best examples, averages and the question 'does any session respond clearly?'
top (4 examples: best clean expanding, largest expanding, same session receding, largest receding): single trials | ROIs x time (trial average) | axon mean +- SEM (raw, movement removed) | pupil / face / whisking
bottom: all sessions x time (trial-averaged axon mean) for expanding and receding discs, with the response z against the shifted-onset null (raw and movement removed)
usage: python fig10_loom_heatmaps.py -> plots/summary_figures/Fig10_loom_heatmaps.pdf/.png ; ledger entries
Data: cross_session/stim_response_results.pkl (population mean, behaviours, movement-removed residual, nulls) + loom_roi.py (ROI-level, same alignment / high-pass / baseline).
"""
import numpy as np
from fig_common import *
from loom_roi import build, STIMS
from matplotlib.colors import TwoSlopeNorm

C_MV = '#8FB3E8'; ZTHR = 2.33; AW = (0, 4)
EXAMPLES = [('AnimalRB24_261023_1538', 'loom_expand', 'clean example'), ('AnimalRB24_261021_1630', 'loom_expand', 'largest raw response'),
            ('AnimalRB24_261023_1538', 'loom_recede', 'same session, receding'), ('AnimalRB24_261106_1604', 'loom_recede', 'largest receding')]
TTL = {'loom_expand': 'expanding', 'loom_recede': 'receding'}


def sem(M): return np.nanstd(M, 0) / np.sqrt(np.isfinite(M).sum(0).clip(1))


def main():
    R = load_pickle('stim_response_results.pkl'); fig = plt.figure(figsize=(W2, 244 * MM))
    gst = fig.add_gridspec(1, 1, left=0.105, right=0.985, top=0.935, bottom=0.415)
    top = gst[0].subgridspec(4, 4, height_ratios=[1.0, 1.35, 1.25, 1.0], hspace=0.14, wspace=0.24); norm = TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.5)
    for j, (nm, key, why) in enumerate(EXAMPLES):
        d = R[nm]['keys'][key]; L = build(nm)['keys'][key]; rel = d['rel']; A = L['A']; fl = is_flagged(nm); n_tr = d['n']
        a0 = fig.add_subplot(top[0, j]); a1 = fig.add_subplot(top[1, j], sharex=a0); a2 = fig.add_subplot(top[2, j], sharex=a0); a3 = fig.add_subplot(top[3, j], sharex=a0)
        # single trials (axon mean)
        a0.imshow(d['P'], aspect='auto', cmap='RdBu_r', norm=norm, extent=(rel[0], rel[-1], n_tr + 0.5, 0.5), interpolation='nearest')
        blk = L['blk']
        for b in np.where(np.diff(blk) != 0)[0]: a0.axhline(b + 1.5, color='white', lw=0.8)
        a0.set_yticks([1, n_tr]); a0.set_ylabel('trial', fontsize=5.4, labelpad=1); a0.tick_params(labelbottom=False)
        a0.set_title('%s, %s%s\n%s; %d trials, %d ROIs' % (TTL[key], short(nm), ' †' if fl else '', why, n_tr, A.shape[0]), fontsize=5.0, loc='left', fontweight='bold', pad=3, color=MUTED if fl else INK)
        if j == 0: plabel(a0, 'a', dx=-34, dy=14)
        # ROIs x time, trial average, sorted by response
        M = np.nanmean(A, 1); resp = np.nanmean(M[:, (rel >= AW[0]) & (rel < AW[1])], 1); M = M[np.argsort(resp)[::-1]]
        im = a1.imshow(M, aspect='auto', cmap='RdBu_r', norm=norm, extent=(rel[0], rel[-1], M.shape[0] + 0.5, 0.5), interpolation='nearest')
        a1.set_yticks([1, M.shape[0]]); a1.set_ylabel('ROI (sorted)', fontsize=5.4, labelpad=1); a1.tick_params(labelbottom=False)
        a1.text(0.985, 0.04, '%d / %d ROIs > 0' % ((resp > 0).sum(), len(resp)), transform=a1.transAxes, ha='right', va='bottom', fontsize=4.8, color=INK, bbox=dict(fc='white', ec='none', alpha=0.75, pad=0.8))
        if j == 0: plabel(a1, 'b', dx=-34, dy=4)
        for a in (a0, a1): a.axvline(0, color=INK, lw=0.5, ls='--'); a.spines['left'].set_visible(False)
        # mean traces
        m = np.nanmean(d['P'], 0); s = sem(d['P']); a2.fill_between(rel, m - s, m + s, color=C_SIG, alpha=0.25, lw=0); a2.plot(rel, m, color=C_SIG, lw=1.4, label='raw')
        if 'res_mv' in d: a2.plot(rel, np.nanmean(d['res_mv'], 0), color=C_MV, lw=1.1, label='movement removed')
        a2.axvspan(*AW, color=C_REJ, alpha=0.15, lw=0); a2.axvline(0, color=INK, lw=0.5, ls='--'); a2.axhline(0, color=MUTED, lw=0.4); a2.tick_params(labelbottom=False)
        a2.set_ylim(-0.3 if j else -0.3, 0.62 if j != 1 else 0.7); a2.set_ylabel('axon mean dF/F' if j == 0 else '', fontsize=5.4)
        zr, zm = d['z']['P'], d['z']['res_mv']; a2.text(0.02, 0.97, 'z = %.1f raw, %.1f movement removed' % (zr, zm), transform=a2.transAxes, va='top', fontsize=4.8, color=INK)
        if j == 0: plabel(a2, 'c', dx=-34, dy=4); a2.legend(frameon=False, fontsize=4.6, loc='lower right', handlelength=1.0, borderpad=0.1, labelspacing=0.15)
        if j > 0: a2.tick_params(labelleft=False)
        # behaviour
        for b, col in (('pupil', C_BEH['pupil']), ('face', C_BEH['face']), ('whisk', C_BEH['whisk'])):
            if b in d['beh']:
                Bm = np.nanmean(d['beh'][b], 0); a3.plot(rel, Bm, color=col, lw=1.0, label={'pupil': 'pupil', 'face': 'face motion', 'whisk': 'whisking'}[b])
        a3.axvline(0, color=INK, lw=0.5, ls='--'); a3.axhline(0, color=MUTED, lw=0.4); a3.axvspan(*AW, color=C_REJ, alpha=0.15, lw=0); a3.set_xlabel('time from loom onset (s)', fontsize=5.6)
        a3.set_ylim(-1.2, 2.2); a3.set_ylabel('behaviour (z)' if j == 0 else '', fontsize=5.4)
        if j == 0: plabel(a3, 'd', dx=-34, dy=4); a3.legend(frameon=False, fontsize=4.6, loc='upper left', handlelength=1.0, borderpad=0.1, labelspacing=0.15, ncol=1)
        if j > 0: a3.tick_params(labelleft=False)
    cax = fig.add_axes([0.105, 0.978, 0.12, 0.007]); cb = fig.colorbar(im, cax=cax, orientation='horizontal', ticks=[-0.5, 0, 0.5]); cb.ax.tick_params(labelsize=4.8, length=1.2)
    fig.text(0.235, 0.9815, 'dF/F (high-passed, baseline-subtracted)      † = flagged session (not pooled)', fontsize=5.0, color=MUTED, ha='left', va='center')

    # ------------------------------------------------ e..h  all sessions
    bot = fig.add_gridspec(1, 5, left=0.155, right=0.985, top=0.335, bottom=0.080, width_ratios=[1.0, 0.45, 0.66, 1.0, 0.45], wspace=0.06); rng = np.random.default_rng(2); panel = {}
    for q, key in enumerate(STIMS):
        sess = [n for n in ALL if key in R[n].get('keys', {})]; sess.sort(key=lambda n: -R[n]['keys'][key]['z']['P'])
        rel = R[sess[0]]['keys'][key]['rel']; rows = [np.nanmean(R[n]['keys'][key]['P'], 0) for n in sess]
        unfl = [r for n, r in zip(sess, rows) if not is_flagged(n)]; rows.append(np.nanmean(unfl, 0)); labels = ['%s%s (%d)' % (short(n), ' †' if is_flagged(n) else '', R[n]['keys'][key]['n']) for n in sess] + ['mean, %d unflagged' % len(unfl)]
        ah = fig.add_subplot(bot[3 * q]); ad = fig.add_subplot(bot[3 * q + 1], sharey=ah); panel[key] = (ah, ad)
        Mx = np.array(rows); ah.imshow(Mx, aspect='auto', cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-0.4, vcenter=0, vmax=0.4), extent=(rel[0], rel[-1], len(rows) - 0.5, -0.5), interpolation='nearest')
        ah.axhline(len(rows) - 1.5, color='white', lw=1.5); ah.axvline(0, color=INK, lw=0.5, ls='--'); ah.set_yticks(range(len(rows))); ah.set_yticklabels(labels, fontsize=4.8); ah.set_xlabel('time from loom onset (s)', fontsize=5.6)
        for tl, n in zip(ah.get_yticklabels(), sess + [None]):
            if n is not None and is_flagged(n): tl.set_color(MUTED)
        ah.set_title('%s loom: axon mean per session' % TTL[key], fontsize=5.6, loc='left', fontweight='bold', pad=3); plabel(ah, 'e' if q == 0 else 'f', dx=-70, dy=6)
        for i, n in enumerate(sess):
            d = R[n]['keys'][key]
            for ser, col in (('P', C_SIG), ('res_mv', C_MV)):
                z = d['z'][ser]; mk = MARK[animal(n)]; fl = is_flagged(n)
                ad.plot(z, i + (-0.12 if ser == 'P' else 0.12), mk, mfc='white' if fl else col, mec=MUTED if fl else 'white', mew=0.8 if fl else 0.3, ms=4.0)
            ad.plot([d['z']['P'], d['z']['res_mv']], [i - 0.12, i + 0.12], color=FAINT, lw=0.6, zorder=0)
        ad.axvline(ZTHR, color=INK, lw=0.5, ls='--'); ad.axvline(0, color=MUTED, lw=0.4); ad.set_xlim(-2.5, 5.0); ad.set_xlabel('response z\nvs shifted-onset null', fontsize=5.4); ad.tick_params(labelleft=False); ad.set_ylim(len(rows) - 0.5, -0.5)
        ad.text(ZTHR + 0.1, len(rows) - 0.9, 'z = 2.33', fontsize=4.4, color=MUTED, ha='left', va='bottom')
        if q == 0: ad.plot([], [], 'o', color=C_SIG, ms=3.6, label='raw'); ad.plot([], [], 'o', color=C_MV, ms=3.6, label='movement removed'); ad.legend(frameon=False, fontsize=4.4, loc='center right', handletextpad=0.1, borderpad=0.1, labelspacing=0.2)
        unf = [n for n in sess if not is_flagged(n)]; k_clear = sum(1 for n in unf if R[n]['keys'][key]['z']['P'] > ZTHR and R[n]['keys'][key]['z']['res_mv'] > ZTHR and np.nanmean(R[n]['keys'][key]['amp']['P']) > 0); k_raw = sum(1 for n in unf if R[n]['keys'][key]['z']['P'] > ZTHR)
        ledger_add('Fig10', key + '_clear_sessions', 'unflagged sessions with a clear %s-loom response (z > 2.33 raw AND after movement removal)' % TTL[key], k_clear / max(len(unf), 1), np.nan, np.nan, len(unf), len({animal(n) for n in unf}), 'ambiguous' if k_clear else 'not seen', note='%d/%d unflagged sessions pass raw only; largest raw z of all sessions: %s%s' % (k_raw, len(unf), short(sess[0]), ' (flagged)' if is_flagged(sess[0]) else ''), unit='fraction of sessions')
    fig.text(0.105, 0.004, 'Top: population = mean best dF/F of the signal-like ROIs, high-passed > 30 s, baseline -2..0 s subtracted; dashed line = loom onset; grey = response window 0-4 s; white lines in the trial raster separate blocks.\n'
             'b: each ROI averaged over the trials, sorted by its 0-4 s response. Movement removed = face + whisking + ball ridge model fitted outside all stimulus windows and subtracted (pupil not used: it reflects the stimulus).\n'
             'e, f: rows ordered by raw z (trials in brackets); z = response against onsets circularly shifted inside their own block (1000 draws). Examples: largest responses (columns 1, 2, 4) and a matched receding disc (3).\n'
             'Pupil, face and whisking are high-passed z. Marker shape = animal (circle RB24, square RB23, triangle RB22), hollow / † = flagged (not pooled).', fontsize=4.6, color=MUTED)
    savefig(fig, 'Fig10_loom_heatmaps')


if __name__ == '__main__':
    main()
