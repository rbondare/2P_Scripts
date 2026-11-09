"""Figures + table of cross-session analysis 2 (stimulus-evoked responses and behaviour-adjusted responses). usage: python stim_figures.py
Reads plots/cross_session/stim_response_results.pkl, writes L*.png and stim_summary.csv there.
z = (session mean amplitude - mean of the shifted-onset null) / SD of the null: comparable over sessions with different dF/F scales. Filled marker = p < 0.05 against the null.
"""
import os, sys, pickle, csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stim_response as sr
from state_figures import C_SIG, C_NOISE, C_REJ, INK, MUTED, C_ANIMAL, short, animal, ylabels
FIG = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session'
plt.rcParams.update({'font.family': 'Arial' if any('Arial' in f.name for f in matplotlib.font_manager.fontManager.ttflist) else 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.75, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 110, 'savefig.dpi': 200, 'axes.titlesize': 8.5, 'axes.titleweight': 'normal'})
BEHC = {'pupil': '#C0392B', 'face': '#E08E0B', 'whisk': '#2E8B57', 'loco': '#2C5AA0'}
BEHL = {'pupil': 'pupil', 'face': 'face motion', 'whisk': 'whisking', 'loco': 'locomotion'}


def p_of(d, s):
    n = d['null'][s]; m = np.nanmean(d['amp'][s]); return (1 + (np.abs(n) >= abs(m)).sum()) / (1 + len(n))


def save(fig, name): fig.savefig(os.path.join(FIG, name), bbox_inches='tight', facecolor='white'); plt.close(fig); print('saved', name)


def sess_with(res, key): return sorted([n for n in res if key in res[n].get('keys', {})], key=lambda n: (animal(n)[::-1], n))


def fig_traces(res, keys):
    rows = ['P'] + ['pupil', 'face', 'whisk', 'loco']
    fig, axs = plt.subplots(len(rows), len(keys), figsize=(2.9 * len(keys), 8.2), sharex='col', gridspec_kw=dict(hspace=0.12, wspace=0.28))
    for j, key in enumerate(keys):
        ss = sess_with(res, key)
        for i, rw in enumerate(rows):
            ax = axs[i, j]; M = []
            for n in ss:
                d = res[n]['keys'][key]; C = d['P'] if rw == 'P' else d['beh'].get(rw)
                if C is None: continue
                mu = np.nanmean(C, 0); M.append(mu); ax.plot(d['rel'], mu, color=C_ANIMAL[animal(n)], lw=0.7, alpha=0.6)
            if M:
                L = min(len(m) for m in M); ax.plot(d['rel'][:L], np.nanmean(np.stack([m[:L] for m in M]), 0), color=INK, lw=1.8)
            ax.axvline(0, color=MUTED, lw=0.6, ls=':'); ax.axhline(0, color=MUTED, lw=0.6, ls=':'); aw = sr.AWIN[key]; ax.axvspan(aw[0], aw[1], color='#EAF1FB', zorder=0)
            if i == 0: ax.set_title('%s\n%d sessions, %d trials' % (sr.LAB[key], len(ss), sum(res[n]['keys'][key]['n'] for n in ss)), fontsize=8)
            if j == 0: ax.set_ylabel('population dF/F' if rw == 'P' else BEHL[rw] + ' (SD)')
            if i == len(rows) - 1: ax.set_xlabel('time from onset (s)')
    fig.suptitle('L1 | Trial-averaged population dF/F and behaviour per stimulus (thin = session, colour = animal; black = mean over sessions; blue band = response window)', x=0.01, ha='left', fontsize=9.3, y=0.94)
    save(fig, 'L1_trial_averages.png')


def fig_amplitudes(res, keys):
    series = [('P', 'dF/F, raw', C_SIG), ('res_mv', 'dF/F after removing\nface + whisking + locomotion', '#0B7A75'), ('res_all', 'dF/F after removing\nthese + pupil (secondary)', '#8A5BB5')]
    beh_rows = ['pupil', 'face', 'whisk', 'loco']
    nr = len(series) + len(beh_rows)
    fig, axs = plt.subplots(nr, len(keys), figsize=(2.9 * len(keys), 1.35 * nr + 0.8), sharey='row', gridspec_kw=dict(hspace=0.5, wspace=0.12))
    for j, key in enumerate(keys):
        ss = sess_with(res, key)
        for i in range(nr):
            ax = axs[i, j]; s = series[i][0] if i < len(series) else 'b_' + beh_rows[i - len(series)]
            col = series[i][2] if i < len(series) else BEHC[beh_rows[i - len(series)]]
            xs = []; zs = []
            for k, n in enumerate(ss):
                d = res[n]['keys'][key]
                if s not in d['z']: continue
                sig = p_of(d, s) < 0.05
                ax.plot(k, d['z'][s], 'o', mfc=col if sig else 'white', mec=col, ms=5); xs.append(k); zs.append(d['z'][s])
            ax.axhline(0, color=MUTED, lw=0.6, ls=':'); ax.axhspan(-2, 2, color='#EEEEEE', zorder=0)
            if zs: ax.axhline(np.mean(zs), color=col, lw=1.4, xmin=0.02, xmax=0.98)
            ax.set_xticks(range(len(ss))); ax.set_xticklabels([short(n)[:4] + ' ' + short(n)[6:] if i == nr - 1 else '' for n in ss], rotation=90, fontsize=5.5)
            if i == nr - 1:
                for t, n in zip(ax.get_xticklabels(), ss): t.set_color(C_ANIMAL[animal(n)])
            if j == 0: ax.set_ylabel((series[i][1] if i < len(series) else 'evoked ' + BEHL[beh_rows[i - len(series)]]) + '\nz vs null', fontsize=6.8)
            if i == 0: ax.set_title('%s' % sr.LAB[key], fontsize=8)
    fig.suptitle('L2 | Evoked amplitude per session as z against onsets shifted inside the block (grey band = +-2; filled = p<0.05; line = mean over sessions)', x=0.01, ha='left', fontsize=9.3, y=0.94)
    save(fig, 'L2_amplitudes_raw_and_behaviour_adjusted.png')


def fig_expand_recede(res):
    both = [n for n in sess_with(res, 'loom_expand') if 'loom_recede' in res[n]['keys']]
    fig, axs = plt.subplots(1, 4, figsize=(12.5, 3.6), gridspec_kw=dict(wspace=0.4))
    for ax, (s, ttl) in zip(axs, (('P', 'dF/F raw'), ('res_mv', 'dF/F, movement removed'), ('b_pupil', 'evoked pupil'), ('b_loco', 'evoked locomotion'))):
        for n in both:
            de, dr = res[n]['keys']['loom_expand'], res[n]['keys']['loom_recede']
            if s not in de['z'] or s not in dr['z']: continue
            ax.plot([0, 1], [de['z'][s], dr['z'][s]], '-o', color=C_ANIMAL[animal(n)], ms=4, lw=1, alpha=0.9)
        ax.axhspan(-2, 2, color='#EEEEEE', zorder=0); ax.axhline(0, color=MUTED, lw=0.6, ls=':'); ax.set_xticks([0, 1]); ax.set_xticklabels(['expanding', 'receding']); ax.set_xlim(-0.3, 1.3); ax.set_ylabel('z vs null'); ax.set_title(ttl)
    fig.suptitle('L3 | Expanding vs receding disc within each session that has both (%d sessions; colour = animal)' % len(both), x=0.01, ha='left', fontsize=9.3, y=1.02)
    save(fig, 'L3_expanding_vs_receding.png')


def fig_trials(res, key='loom_expand'):
    ss = sess_with(res, key); fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.4), gridspec_kw=dict(wspace=0.32))
    for ax, b in zip(axs, ('pupil', 'face', 'loco')):
        X, Y = [], []
        for n in ss:
            d = res[n]['keys'][key]; a = d['amp']['P']; ab = d['amp'].get('b_' + b)
            if ab is None: continue
            ok = np.isfinite(a) & np.isfinite(ab); za = (a[ok] - np.mean(a[ok])) / (np.std(a[ok]) + 1e-9); zb = (ab[ok] - np.mean(ab[ok])) / (np.std(ab[ok]) + 1e-9)
            ax.plot(zb, za, 'o', color=C_ANIMAL[animal(n)], ms=3, alpha=0.6); X += list(zb); Y += list(za)
        if len(X) > 5:
            from scipy.stats import spearmanr; rho, p = spearmanr(X, Y); ax.set_title('%s: Spearman rho %+.2f (p = %.3g, n = %d trials)' % (BEHL[b], rho, p, len(X)), fontsize=7.6)
        ax.axhline(0, color=MUTED, lw=0.5, ls=':'); ax.axvline(0, color=MUTED, lw=0.5, ls=':'); ax.set_xlabel('evoked %s (z within session)' % BEHL[b]); ax.set_ylabel('dF/F response (z within session)')
    fig.suptitle('L4 | Trial by trial: does a bigger %s response go with a bigger behaviour change? (z within session, all sessions pooled)' % sr.LAB[key], x=0.01, ha='left', fontsize=9.3, y=1.03)
    save(fig, 'L4_trial_relations_%s.png' % key)


def table(res):
    rows = []
    for key in sr.KEYS:
        for n in sess_with(res, key):
            d = res[n]['keys'][key]; r = dict(session=n, stimulus=key, n_trials=d['n'], n_signal_like=res[n]['n_sig'])
            for s in ('P', 'res_mv', 'res_all', 'Phat_mv') + tuple('b_' + b for b in sr.BEH_ORDER):
                if s in d['amp']: r[s + '_amp'] = round(float(np.nanmean(d['amp'][s])), 4); r[s + '_z'] = round(d['z'][s], 2); r[s + '_p'] = round(p_of(d, s), 4)
            rows.append(r)
    cols = sorted({k for r in rows for k in r}, key=lambda k: (k not in ('session', 'stimulus', 'n_trials', 'n_signal_like'), k))
    with open(os.path.join(FIG, 'stim_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); [w.writerow(r) for r in rows]
    print('saved stim_summary.csv', len(rows), 'rows')
    # pooled
    print('\npooled over sessions (mean of session means; null = mean of independent shifts):')
    for key in sr.KEYS:
        for s in ('P', 'res_mv', 'res_all', 'b_pupil', 'b_face', 'b_loco'):
            p = sr.pooled(res, key, s)
            if p: print('  %-12s %-7s %2d sessions mean %+.3f  null %+.3f..%+.3f  p=%.3f' % (key, s, p['n_sessions'], p['mean'], *p['null95'], p['p']))


if __name__ == '__main__':
    res = pickle.load(open(os.path.join(FIG, 'stim_response_results.pkl'), 'rb'))
    keys = [k for k in ('flash', 'loom_expand', 'loom_recede', 'bar', 'square') if len(sess_with(res, k)) >= 2]
    fig_traces(res, keys); fig_amplitudes(res, keys); fig_expand_recede(res); fig_trials(res, 'loom_expand'); table(res)
