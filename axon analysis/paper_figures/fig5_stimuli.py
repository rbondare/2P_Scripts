"""Figure 5 - visual stimulus responses of the axon mean.
a  single-trial rasters (RB24 10-23 1538: the one unflagged session with all four visual stimuli, not the strongest session)   b  trial averages (axon mean raw / movement removed, pupil, face)
c  response size per session vs shifted-onset null   d  expanding minus receding   e  trial-wise loom response vs evoked pupil change
usage: python fig5_stimuli.py -> plots/summary_figures/Fig5_stimuli.pdf/.png ; ledger entries
Data: cross_session/stim_response_results.pkl (population mean P high-passed > 30 s; 'movement removed' = ridge model of face + whisking + ball fitted OUTSIDE the stimulus windows and subtracted).
"""
import numpy as np
from fig_common import *
from scipy.stats import spearmanr
from matplotlib.colors import TwoSlopeNorm

STIM = [('flash', 'full-field /\nlocal flash'), ('bar', 'moving bar'), ('loom_expand', 'expanding\nloom'), ('loom_recede', 'receding\nloom'), ('square', 'square\n(task / passive)')]
AWIN = {'flash': (0, 2), 'bar': (0.5, 4), 'loom_expand': (0, 4), 'loom_recede': (0, 4), 'square': (0, 3)}
EX = 'AnimalRB24_261023_1538'
C_MV = '#8FB3E8'


def main():
    R = load_pickle('stim_response_results.pkl'); fig = plt.figure(figsize=(W2, 218 * MM))
    gs = fig.add_gridspec(3, 1, height_ratios=[34, 62, 62], hspace=0.40, left=0.085, right=0.985, top=0.955, bottom=0.075)
    # ------------------------------------------------ a  rasters
    ga = gs[0].subgridspec(1, 4, wspace=0.28); r = R[EX]['keys']; vmax = 0
    for k in ('flash', 'bar', 'loom_expand', 'loom_recede'): vmax = max(vmax, np.nanpercentile(np.abs(r[k]['P']), 98))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    for j, k in enumerate(('flash', 'bar', 'loom_expand', 'loom_recede')):
        a = fig.add_subplot(ga[j]); d = r[k]; rel = d['rel']; im = a.imshow(d['P'], aspect='auto', cmap='RdBu_r', norm=norm, extent=(rel[0], rel[-1], d['n'] + 0.5, 0.5), interpolation='nearest')
        a.axvline(0, color=INK, lw=0.6, ls='--'); a.set_xlabel('time from onset (s)', fontsize=5.8); a.set_title('%s (n = %d)' % (dict(STIM)[k].replace('\n', ' '), d['n']), fontsize=5.8, loc='left', fontweight='normal')
        a.set_ylabel('trial' if j == 0 else ''); a.set_yticks([1, d['n']])
        if j == 0: plabel(a, 'a', dx=-24, dy=8)
    cax = fig.add_axes([0.80, 0.982, 0.11, 0.007]); cb = fig.colorbar(im, cax=cax, orientation='horizontal'); cb.ax.tick_params(labelsize=4.5, length=1.2); fig.text(0.795, 0.9855, 'axon mean dF/F, high-passed', fontsize=4.8, color=MUTED, ha='right', va='center')
    fig.text(0.085, 0.995, '%s: single trials' % short(EX), fontsize=5.8, color=MUTED, va='top')

    # ------------------------------------------------ b  trial averages
    gb = gs[1].subgridspec(2, 5, hspace=0.14, wspace=0.34, height_ratios=[1.5, 1]); axb = []
    for j, (k, lab) in enumerate(STIM):
        a1 = fig.add_subplot(gb[0, j]); a2 = fig.add_subplot(gb[1, j], sharex=a1); axb.append(a1); n_s = 0; M = []; Mv = []; Bp = []; Bf = []; rel = None
        for n in ALL:
            if k not in R[n].get('keys', {}): continue
            d = R[n]['keys'][k]; rel = d['rel']; m = np.nanmean(d['P'], 0); fl = is_flagged(n)
            a1.plot(rel, m, color=FAINT if fl else C_SIG, lw=0.5, alpha=1 if fl else 0.45, ls='--' if fl else '-', zorder=2)
            if not fl:
                n_s += 1; M.append(m)
                if 'res_mv' in d: Mv.append(np.nanmean(d['res_mv'], 0))
                if 'pupil' in d['beh']: Bp.append(np.nanmean(d['beh']['pupil'], 0))
                if 'face' in d['beh']: Bf.append(np.nanmean(d['beh']['face'], 0))
        if M:
            a1.plot(rel, np.nanmean(M, 0), color=C_SIG, lw=1.7, zorder=4, label='raw')
            if Mv: a1.plot(rel, np.nanmean(Mv, 0), color=C_MV, lw=1.2, ls='-', zorder=3, label='movement removed')
        w = AWIN[k]; a1.axvspan(w[0], w[1], color=C_REJ, alpha=0.15, lw=0, zorder=0)
        a1.axvline(0, color=INK, lw=0.5, ls=':'); a1.axhline(0, color=MUTED, lw=0.4); a1.set_title('%s\n(%d unflagged sessions)' % (lab.replace('\n', ' '), n_s), fontsize=5.8, loc='left', fontweight='bold'); a1.tick_params(labelbottom=False)
        if j == 0: a1.set_ylabel('axon mean dF/F (high-passed)'); plabel(a1, 'b', dx=-24, dy=14); a1.legend(frameon=False, fontsize=4.8, loc='upper right', handlelength=1.0, borderpad=0.1)
        if Bp: a2.plot(rel, np.nanmean(Bp, 0), color=C_BEH['pupil'], lw=1.0, label='pupil')
        if Bf: a2.plot(rel, np.nanmean(Bf, 0), color=C_BEH['face'], lw=1.0, label='face')
        a2.axvline(0, color=INK, lw=0.5, ls=':'); a2.axhline(0, color=MUTED, lw=0.4); a2.set_xlabel('time from onset (s)', fontsize=5.8); a2.axvspan(w[0], w[1], color=C_REJ, alpha=0.15, lw=0, zorder=0)
        if j == 0: a2.set_ylabel('behaviour (high-passed z)'); a2.legend(frameon=False, fontsize=4.8, loc='upper left', handlelength=1.0, borderpad=0.1)
    # y limits shared between the axon panels
    lim = max(abs(np.array(a.get_ylim())).max() for a in axb); [a.set_ylim(-0.3 * lim, lim) for a in axb[:0]]

    # ------------------------------------------------ c  response size per session vs null
    gc = gs[2].subgridspec(1, 3, wspace=0.42, width_ratios=[1.5, 0.9, 1.1]); ac = fig.add_subplot(gc[0]); ad = fig.add_subplot(gc[1]); ae = fig.add_subplot(gc[2]); rng = np.random.default_rng(4)
    for j, (k, lab) in enumerate(STIM):
        for q, (ser, col) in enumerate((('P', C_SIG), ('res_mv', C_MV))):
            vals = {}; nul = []; ksig = 0
            for n in ALL:
                d = R[n].get('keys', {}).get(k)
                if not d or ser not in d['amp']: continue
                A = float(np.nanmean(d['amp'][ser])); x = j + (q - 0.5) * 0.4 + rng.uniform(-0.07, 0.07); sdot(ac, x, A, n, color=col, hollow_color=MUTED if q == 0 else MUTED, ms=3.4)
                if not is_flagged(n): vals[n] = A; nul.append(d['null'][ser]); ksig += int(d['z'][ser] > 2.33 and A > 0)
            h = pooled_unflagged(vals, n_boot=4000)
            if h['n_sessions']:
                L = min(len(x) for x in nul); pn = np.mean([x[:L] for x in nul], 0); lo_n, hi_n = np.percentile(pn, [2.5, 97.5]); p = float((1 + (np.abs(pn) >= abs(np.mean(list(vals.values())))).sum()) / (1 + L)); x0 = j + (q - 0.5) * 0.4
                ac.fill_between([x0 - 0.2, x0 + 0.2], lo_n, hi_n, color=C_REJ, alpha=0.35, lw=0, zorder=0); ac.plot([x0 - 0.17, x0 + 0.17], [h['est']] * 2, color=INK, lw=1.4, zorder=5); ac.plot([x0, x0], [h['lo'], h['hi']], color=INK, lw=0.9, zorder=5)
                if q == 0: ac.text(j, 0.205, '%d/%d' % (ksig, h['n_sessions']), ha='center', va='bottom', fontsize=4.8, color=MUTED)
                verdict = 'supported' if (h['lo'] > 0 and p < 0.01 and ksig >= 0.5 * h['n_sessions']) else 'ambiguous' if (h['lo'] > 0 or p < 0.05) else 'not seen'
                ledger_add('Fig5', '%s_%s' % (k, ser), 'axon response to %s (%s)' % (lab.replace('\n', ' '), 'raw' if ser == 'P' else 'movement removed'), h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], verdict, note='pooled shifted-onset p = %.3f; sessions with z > 2.3: %d/%d' % (p, ksig, h['n_sessions']), unit='dF/F (high-passed)')
    ac.axhline(0, color=MUTED, lw=0.4); ac.set_xticks(range(len(STIM))); ac.set_xticklabels([s[1] for s in STIM], fontsize=5.2); ac.set_xlim(-0.6, len(STIM) - 0.4); ac.set_ylim(-0.12, 0.28); ac.set_ylabel('response (dF/F)'); plabel(ac, 'c', dx=-24, dy=8)
    ac.text(0.0, 1.03, 'left of pair: raw (blue); right: movement removed (pale);\ngrey = null of the pooled mean; n/N = unflagged sessions with z > 2.3', transform=ac.transAxes, fontsize=4.6, color=MUTED, va='bottom')
    # d  expanding - receding within sessions
    for q, (ser, col, lab) in enumerate((('P', C_SIG, 'raw'), ('res_mv', C_MV, 'movement\nremoved'))):
        vals = {}
        for n in ALL:
            k = R[n].get('keys', {})
            if 'loom_expand' in k and 'loom_recede' in k and ser in k['loom_expand']['amp'] and ser in k['loom_recede']['amp']:
                v = float(np.nanmean(k['loom_expand']['amp'][ser]) - np.nanmean(k['loom_recede']['amp'][ser])); sdot(ad, q + rng.uniform(-0.15, 0.15), v, n, color=col, hollow_color=MUTED, ms=3.6)
                if not is_flagged(n): vals[n] = v
        h = pooled_unflagged(vals, n_boot=4000)
        if h['n_sessions']:
            ad.plot([q - 0.3, q + 0.3], [h['est']] * 2, color=INK, lw=1.4, zorder=5); ad.plot([q, q], [h['lo'], h['hi']], color=INK, lw=0.9, zorder=5)
            ledger_add('Fig5', 'expand_minus_recede_' + ser, 'expanding minus receding loom response (%s)' % ('raw' if ser == 'P' else 'movement removed'), h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'supported' if h['lo'] > 0 else 'not seen' if h['lo'] < 0 < h['hi'] else 'ambiguous', note='sessions with both stimuli', unit='dF/F')
    ad.axhline(0, color=MUTED, lw=0.4); ad.set_xticks([0, 1]); ad.set_xticklabels(['raw', 'movement\nremoved'], fontsize=5.4); ad.set_xlim(-0.6, 1.6); ad.set_ylabel('expanding − receding (dF/F)'); plabel(ad, 'd', dx=-24, dy=8)
    ad.text(0.0, 1.03, 'sessions with\nboth stimuli', transform=ad.transAxes, fontsize=4.6, color=MUTED, va='bottom')
    # e  trial-wise loom amplitude vs evoked pupil change
    xs, ys, per = [], [], {}
    for n in ALL:
        d = R[n].get('keys', {}).get('loom_expand')
        if not d or 'b_pupil' not in d['amp']: continue
        a = d['amp']['P']; b = d['amp']['b_pupil']; ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 5: continue
        fl = is_flagged(n); mk = MARK[animal(n)]
        ae.plot(b[ok], a[ok], mk, mfc='white' if fl else C_SIG, mec=MUTED if fl else 'white', mew=0.5 if not fl else 0.7, ms=3.0, alpha=0.9, zorder=2)
        per[n] = spearmanr(a[ok], b[ok])[0]
        if not fl: xs += list(b[ok]); ys += list(a[ok])
    xs = np.array(xs); ys = np.array(ys); rho, p = spearmanr(xs, ys); m, c0 = np.polyfit(xs, ys, 1); xx = np.linspace(xs.min(), xs.max(), 10); ae.plot(xx, m * xx + c0, color=INK, lw=1.0, zorder=5)
    hr = pooled_unflagged(per, n_boot=3000); ae.axhline(0, color=MUTED, lw=0.4); ae.axvline(0, color=MUTED, lw=0.4); ae.set_xlabel('evoked pupil change (z)'); ae.set_ylabel('axon response per trial (dF/F)'); plabel(ae, 'e', dx=-24, dy=8)
    ae.text(0.03, 1.03, 'unflagged: ρ = %+.2f (%d trials);\nwithin session %+.2f [%.2f, %.2f]' % (rho, len(xs), hr['est'], hr['lo'], hr['hi']), transform=ae.transAxes, fontsize=4.8, va='bottom')
    ledger_add('Fig5', 'loom_vs_pupil', 'trial-wise expanding-loom axon response vs evoked pupil change (within-session Spearman)', hr['est'], hr['lo'], hr['hi'], hr['n_sessions'], hr['n_animals'], 'ambiguous', note='confound: trials with larger arousal change have larger responses in every measured channel; pooled trials rho %.2f (n = %d)' % (rho, len(xs)), unit='rho')
    fig.text(0.085, 0.004, 'Response = mean(axon, window after onset) - mean(baseline) per trial, population mean of the signal-like ROIs high-passed > 30 s; windows: flash 0-2 s, bar 0.5-4 s, loom 0-4 s, square 0-3 s.\n'
             'Null: onsets circularly shifted inside their own block (1000 draws). Movement removed = ridge model of face motion + whisking + ball (lags -1..+3 s) fitted outside all stimulus windows and subtracted;\n'
             'pupil is deliberately not used (flashes constrict the pupil: a pupil-adjusted residual would be spurious). Marker shape = animal, filled = unflagged, hollow / dashed = flagged (not pooled);\n'
             'black bars = mean of per-animal means with 95% hierarchical bootstrap CI (at most 3 mice, most stimuli 2-3).', fontsize=4.7, color=MUTED)
    savefig(fig, 'Fig5_stimuli')


if __name__ == '__main__':
    main()
