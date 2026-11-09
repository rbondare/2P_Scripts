"""Figures of cross-session analysis 4 (task engagement / context). usage: python task_figures.py  (reads task_results.pkl, writes T*.png + task_summary.csv in plots/cross_session)

Sessions are told apart by marker (circle, square, triangle, diamond, down-triangle; the flagged RB23 session is a grey cross and is never pooled).
Colour = what the thing is: axons blue (signal-like mean), axons with the concurrent movement removed light blue / hollow, behaviour dark ink, licks brown, hit trials blue, miss trials orange.
Pooled statistic = Stouffer combination of the per-session z (each against its own circular-shift / permutation null).
"""
import os, csv, pickle, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
FIG = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session'
plt.rcParams.update({'font.family': 'Arial' if any('Arial' in f.name for f in matplotlib.font_manager.fontManager.ttflist) else 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.75, 'axes.spines.top': False, 'axes.spines.right': False,
                     'xtick.direction': 'out', 'ytick.direction': 'out', 'figure.dpi': 110, 'savefig.dpi': 200, 'axes.titlesize': 8.5, 'axes.titleweight': 'normal'})
C_AX, C_RES, INK, MUTED, C_LICK, C_HIT, C_MISS, C_FLAG = '#1F5FBF', '#8FB3E8', '#222222', '#6B7280', '#8C5A2B', '#1F5FBF', '#E07B1F', '#9CA3AF'
R = pickle.load(open(os.path.join(FIG, 'task_results.pkl'), 'rb'))
NAMES = [n for n in R if 'epochs' in R[n]]
POOL = [n for n in NAMES if not R[n].get('flag')]
MK = {n: m for n, m in zip(NAMES, ['o', 's', '^', 'D', 'v', 'X'])}
def short(nm): a, d, t = nm.replace('Animal', '').split('_'); return '%s %s-%s %s' % (a, d[2:4], d[4:6], t)
def col_of(n, base): return C_FLAG if R[n].get('flag') else base
def save(fig, name): fig.savefig(os.path.join(FIG, name), bbox_inches='tight', facecolor='white'); plt.close(fig); print('saved', name)
def stouffer(z): z = np.array([x for x in z if np.isfinite(x)]); return (float(z.sum() / np.sqrt(len(z))), len(z)) if len(z) else (np.nan, 0)
def pval(z): from math import erf, sqrt; return 1 - erf(abs(z) / sqrt(2))
SUMMARY = []


# ------------------------------------------------------------------ T1 timeline
def fig_timeline():
    n = len(NAMES); fig, axs = plt.subplots(5, n, figsize=(2.75 * n, 9.2), gridspec_kw=dict(hspace=0.28, wspace=0.3, height_ratios=[1, 1, 1, 1, 1]))
    for j, nm in enumerate(NAMES):
        r = R[nm]; tl = r['timeline']; tm = (tl['t'] - tl['t'][0] + 30) / 60; tr = r['trials']; t0 = tl['t'][0] - 30
        sub0 = (np.atleast_1d(tr['sub_t0']) - t0) / 60; sub1 = (np.atleast_1d(tr['sub_t1']) - t0) / 60
        for a in axs[:, j]:
            for k in range(len(sub0)): a.axvspan(sub0[k], sub1[k], color='#F3F4F6', lw=0)
        a = axs[0, j]; on = (tr['on'] - t0) / 60; hit = tr['hit'] > 0
        a.plot(on[hit], np.ones(hit.sum()), '|', color=C_HIT, ms=7, mew=1.0); a.plot(on[~hit], np.zeros((~hit).sum()), '|', color=C_MISS, ms=7, mew=1.0)
        hr = np.array([np.mean(hit[max(0, i - 4):i + 5]) for i in range(len(hit))]); a.plot(on, hr * 0.8 + 0.1, color=INK, lw=0.9)
        a.set_yticks([0, 1]); a.set_yticklabels(['miss', 'hit'], fontsize=6.5); a.set_title(short(nm) + ('\n(not pooled)' if r.get('flag') else '\n%d signal-like ROIs' % r['n_sig']), fontsize=7.5, color=C_FLAG if r.get('flag') else INK)
        axs[1, j].plot(tm, tl['lick_rate'], color=C_LICK, lw=0.9); axs[2, j].plot(tm, tl['P_sd'], color=C_AX, lw=0.9); axs[3, j].plot(tm, tl['P'], color=C_AX, lw=0.9)
        for bn, c in (('pupil', INK), ('face', MUTED)):
            if bn in tl: axs[4, j].plot(tm, tl[bn], color=c, lw=0.9, label=bn if j == 0 else None)
        for i in range(5):
            axs[i, j].set_xlim(0, tm[-1] + 0.5)
            if i < 4: axs[i, j].set_xticklabels([])
        axs[4, j].set_xlabel('minutes', fontsize=7)
    for i, lab in enumerate(['trials (ticks) and\n10-trial hit rate (line)', 'licks / min', 'axon fluctuation SD\n(1-30 s, dF/F)', 'axon mean dF/F', 'pupil (dark), face motion (grey), z']): axs[i, 0].set_ylabel(lab, fontsize=6.8)
    axs[4, 0].legend(frameon=False, fontsize=6, loc='upper left')
    fig.suptitle('T1 | The lick task over the session: engagement (hits, licking) next to the axon signal and arousal  (grey bands = subsessions, white gaps = ~100 s breaks without trials)', x=0.01, ha='left', fontsize=9.5)
    save(fig, 'T1_task_timeline.png')


# ------------------------------------------------------------------ T2 subsessions
def fig_subsessions():
    items = [('hit', 'hit rate (lick in the square)', 100), ('lat', 'median first-lick latency (s)', 1), ('timeouts', 'time-outs per trial (licks in the wait)', 1), ('lick_rate', 'licks / min (whole subsession)', 1),
             ('P_sd', 'axon fluctuation SD (dF/F)', 1), ('P_mean', 'axon mean dF/F', 1), ('P_rate', 'axon events / min', 1), ('pupil', 'pupil area (z)', 1)]
    fig, axs = plt.subplots(2, 4, figsize=(12.5, 6.0), gridspec_kw=dict(hspace=0.5, wspace=0.38)); axs = axs.ravel()
    for a, (k, ttl, sc) in zip(axs, items):
        for nm in NAMES:
            S_ = R[nm]['subs']; x = [s['sub'] for s in S_]; y = [s.get(k, np.nan) * sc for s in S_]
            col = col_of(nm, C_LICK if k in ('lick_rate',) else (INK if k in ('hit', 'lat', 'timeouts', 'pupil') else C_AX))
            a.plot(x, y, '-' + MK[nm], color=col, lw=0.9, ms=4.5, alpha=0.9)
        a.set_title(ttl, fontsize=7.8); a.set_xticks([1, 2, 3, 4]); a.set_xlabel('subsession', fontsize=7)
    for nm in NAMES: axs[0].plot([], [], '-' + MK[nm], color=INK if not R[nm].get('flag') else C_FLAG, ms=4.5, lw=0.9, label=short(nm))
    axs[0].legend(frameon=False, fontsize=5.8, loc='lower left')
    fig.suptitle('T2 | Engagement over the session (top row: behaviour of the animal) and the axon signal per subsession (bottom row)', x=0.01, ha='left', fontsize=9.5)
    save(fig, 'T2_task_subsessions.png')


# ------------------------------------------------------------------ T3 epochs
def fig_epochs():
    order = ['wait (hit)', 'square (hit)', 'post (hit)', 'wait (miss)', 'square (miss)', 'post (miss)', 'break']
    series = [('P', 'axons (signal-like mean)', C_AX), ('Pres', 'axons after removing face + whisk + ball', C_RES), ('face', 'face motion', INK), ('lick', 'lick rate', C_LICK)]
    fig, axs = plt.subplots(1, 4, figsize=(13.5, 4.0), sharey=False, gridspec_kw=dict(wspace=0.28))
    for a, (sk, ttl, col) in zip(axs, series):
        for nm in NAMES:
            E = R[nm]['epochs']
            if sk not in E['series']: continue
            si = E['series'].index(sk)
            for ei, e in enumerate(order):
                if e not in E['names']: continue
                i = E['names'].index(e); z = E['z'][i, si]; v = E['obs'][i, si]
                a.plot(ei + (list(NAMES).index(nm) - 2.5) * 0.09, v, MK[nm], color=col_of(nm, col), mfc=col_of(nm, col) if abs(z) > 3 else 'white', ms=5, mew=0.9, alpha=0.95)
        for x in (2.5,): a.axvline(x, color=MUTED, lw=0.4, ls=':'); a.axvline(5.5, color=MUTED, lw=0.4, ls=':')
        a.axhline(0, color=MUTED, lw=0.5); a.set_xticks(range(len(order))); a.set_xticklabels([o.replace(' (', '\n(') for o in order], fontsize=6.3, rotation=0); a.set_title(ttl, fontsize=8, color=INK)
        a.set_ylabel('mean in the epoch minus mean over all epochs (session SD)', fontsize=6.3)
    for nm in NAMES: axs[0].plot([], [], MK[nm], color=INK if not R[nm].get('flag') else C_FLAG, ms=4.5, label=short(nm))
    axs[3].legend(frameon=False, fontsize=6, loc='upper left', bbox_to_anchor=(1.02, 1.0))
    zs_ = {}
    for sk in ('P', 'Pres', 'face', 'lick'):
        for e in order:
            zz = [R[n]['epochs']['z'][R[n]['epochs']['names'].index(e), R[n]['epochs']['series'].index(sk)] for n in POOL if e in R[n]['epochs']['names'] and sk in R[n]['epochs']['series']]
            mm = [R[n]['epochs']['obs'][R[n]['epochs']['names'].index(e), R[n]['epochs']['series'].index(sk)] for n in POOL if e in R[n]['epochs']['names'] and sk in R[n]['epochs']['series']]
            z, k = stouffer(zz); SUMMARY.append(dict(analysis='epoch', series=sk, item=e, mean_effect_SD=round(float(np.mean(mm)), 3) if mm else np.nan, pooled_z=round(z, 2), n_sessions=k))
    fig.suptitle('T3 | Trial epochs and breaks: filled marker = beyond chance (|z| > 3 vs the whole epoch pattern shifted in time), hollow = not.  Pupil is left out: it dilates in the breaks because the screen differs there', x=0.01, ha='left', fontsize=9)
    save(fig, 'T3_task_epochs.png')


# ------------------------------------------------------------------ T4 engagement bins and T5 pre-trial predictors share a dot-plot
def dotplot(ax, rows, getter, title, xlabel):
    for i, (key, lab, col) in enumerate(rows):
        for nm in NAMES:
            d = getter(nm, key)
            if d is None or not np.isfinite(d['rho']): continue
            yy = i + (list(NAMES).index(nm) - 2.5) * 0.11; sig = abs((d['rho'] - d['null_mean']) / max(d['null_sd'], 1e-9)) > 2
            ax.plot([d['null_mean'] - 1.96 * d['null_sd'], d['null_mean'] + 1.96 * d['null_sd']], [yy, yy], color='#D1D5DB', lw=2.0, solid_capstyle='butt', zorder=1)
            ax.plot(d['rho'], yy, MK[nm], color=col_of(nm, col), mfc=col_of(nm, col) if sig else 'white', ms=4.8, mew=0.9, zorder=3)
    ax.axvline(0, color=INK, lw=0.5); ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[1] for r in rows], fontsize=7.2); ax.set_ylim(len(rows) - 0.4, -0.6); ax.set_title(title, fontsize=8); ax.set_xlabel(xlabel, fontsize=7)
    for t, r in zip(ax.get_yticklabels(), rows): t.set_color(r[2] if r[2] != C_RES else C_AX)


def fig_bins():
    rows = [('P_mean', 'axon mean dF/F', C_AX), ('P_sd', 'axon fluctuation SD', C_AX), ('P_rate', 'axon event rate', C_AX), ('Pres_mean', 'axon mean, movement removed', C_RES), ('Pres_sd', 'axon SD, movement removed', C_RES), ('pupil', 'pupil', INK), ('face', 'face motion', INK), ('whisk', 'whisking', INK)]
    fig, ax = plt.subplots(1, 1, figsize=(6.6, 4.6))
    dotplot(ax, rows, lambda nm, k: R[nm]['bin_corr'].get(k), 'T4 | Engagement over the session: lick rate per minute vs axon / behaviour per minute', 'Spearman rho (linear time trend removed from both; break minutes left out)\ngrey bar = 95% of circular shifts of the lick-rate series; filled = beyond it')
    for nm in NAMES: ax.plot([], [], MK[nm], color=INK if not R[nm].get('flag') else C_FLAG, ms=4.5, label='%s (%d bins)' % (short(nm), R[nm]['bin_corr']['P_mean']['n']))
    ax.legend(frameon=False, fontsize=5.8, loc='lower left', bbox_to_anchor=(1.0, 0.0))
    for k, lab, c in rows:
        z, nn = stouffer([(R[n]['bin_corr'][k]['rho'] - R[n]['bin_corr'][k]['null_mean']) / max(R[n]['bin_corr'][k]['null_sd'], 1e-9) for n in POOL if k in R[n]['bin_corr']]); SUMMARY.append(dict(analysis='engagement bins (lick rate vs)', series=k, item='rho', mean_effect_SD=round(float(np.nanmean([R[n]['bin_corr'][k]['rho'] for n in POOL if k in R[n]['bin_corr']])), 3), pooled_z=round(z, 2), n_sessions=nn))
    save(fig, 'T4_task_engagement_bins.png')


def fig_pre():
    rows = [('P_hp', 'axon (30 s high-passed)', C_AX), ('P_res', 'axon, movement removed', C_RES), ('pupil', 'pupil', INK), ('face', 'face motion', INK), ('whisk', 'whisking', INK), ('loco', 'locomotion', INK)]
    outs = [('hit', 'hit (lick in the square)'), ('log latency', 'log latency of the first lick (hits)'), ('licks before square', 'licks in the 5 s before the square')]
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.9), sharey=True, gridspec_kw=dict(wspace=0.12))
    for a, (ok, ttl) in zip(axs, outs):
        dotplot(a, rows, lambda nm, k, ok=ok: R[nm]['pre_corr'].get(ok, {}).get(k), ttl, 'Spearman rho, state in [-6, -1] s before the square')
        for k, lab, c in rows:
            z, nn = stouffer([(R[n]['pre_corr'][ok][k]['rho'] - R[n]['pre_corr'][ok][k]['null_mean']) / max(R[n]['pre_corr'][ok][k]['null_sd'], 1e-9) for n in POOL if k in R[n]['pre_corr'][ok]]); SUMMARY.append(dict(analysis='pre-trial state vs ' + ok, series=k, item='rho', mean_effect_SD=round(float(np.nanmean([R[n]['pre_corr'][ok][k]['rho'] for n in POOL if k in R[n]['pre_corr'][ok]])), 3), pooled_z=round(z, 2), n_sessions=nn))
    for nm in NAMES: axs[0].plot([], [], MK[nm], color=INK if not R[nm].get('flag') else C_FLAG, ms=4.5, label=short(nm))
    axs[0].legend(frameon=False, fontsize=5.6, loc='lower left')
    fig.suptitle('T5 | Does the state just before the square predict how the animal performs?  (grey bar = 95% of the outcome series shifted by >= 10 trials; filled = beyond it)', x=0.01, ha='left', fontsize=9.5)
    save(fig, 'T5_task_pre_trial_state.png')


# ------------------------------------------------------------------ T6 hit vs miss
def fig_hitmiss():
    n = len(NAMES); fig = plt.figure(figsize=(2.75 * n, 8.4)); gs = fig.add_gridspec(3, n, hspace=0.5, wspace=0.32, height_ratios=[1, 1, 1.15])
    for j, nm in enumerate(NAMES):
        r = R[nm]; rel = r['trace_rel']
        for i, (kh, km, lab) in enumerate((('trace_hit', 'trace_miss', 'axon (30 s high-passed)'), ('trace_res_hit', 'trace_res_miss', 'axon, movement removed'))):
            a = fig.add_subplot(gs[i, j])
            if kh in r:
                for key, col, tag in ((kh, C_HIT, 'hit'), (km, C_MISS, 'miss')):
                    M = r[key]; m = np.nanmean(M, 0); se = np.nanstd(M, 0) / np.sqrt(np.maximum((np.isfinite(M).sum(0)), 1))
                    a.fill_between(rel, m - se, m + se, color=col, alpha=0.2, lw=0); a.plot(rel, m, color=col, lw=1.2, label='%s (%d)' % (tag, len(M)))
            a.axvline(0, color=MUTED, lw=0.5, ls=':'); a.axvspan(0, 3, color='#F3F4F6', lw=0, zorder=0); a.axhline(0, color=MUTED, lw=0.4)
            if i == 0: a.set_title(short(nm) + ('\n(not pooled)' if r.get('flag') else ''), fontsize=7.5, color=C_FLAG if r.get('flag') else INK)
            if j == 0: a.set_ylabel(lab + '\ndF/F vs -3..-0.5 s', fontsize=6.6)
            a.set_xlabel('s from the square onset', fontsize=6.8); a.legend(frameon=False, fontsize=6, loc='upper left')
    for c, (ki, ttl) in enumerate((('response', 'hit - miss: response 0-3 s'), ('pre level', 'hit - miss: level before the square'))):
        a = fig.add_subplot(gs[2, 0:3] if c == 0 else gs[2, 3:n]); rows_ = []
        for i, (pk, lab, col) in enumerate((('P_hp', 'axon (30 s high-passed)', C_AX), ('P_res', 'axon, movement removed', C_RES))):
            for nm in NAMES:
                d = R[nm]['hit_miss'].get('%s %s' % (pk, ki))
                if not d: continue
                yy = i + (list(NAMES).index(nm) - 2.5) * 0.11; sg = abs(d['z']) > 2
                a.plot([d['diff'] - 1.96 * d['null_sd'] * 0 - 1.96 * d['null_sd'], d['diff'] + 1.96 * d['null_sd']], [yy, yy], color='#D1D5DB', lw=2.0, solid_capstyle='butt', zorder=1)
                a.plot(d['diff'], yy, MK[nm], color=col_of(nm, col), mfc=col_of(nm, col) if sg else 'white', ms=5, mew=0.9, zorder=3)
            zz, kk = stouffer([R[n]['hit_miss']['%s %s' % (pk, ki)]['z'] for n in POOL if '%s %s' % (pk, ki) in R[n]['hit_miss']])
            SUMMARY.append(dict(analysis='hit vs miss ' + ki, series=pk, item='hit - miss (dF/F)', mean_effect_SD=round(float(np.mean([R[n]['hit_miss']['%s %s' % (pk, ki)]['diff'] for n in POOL])), 3), pooled_z=round(zz, 2), n_sessions=kk))
        a.axvline(0, color=INK, lw=0.5); a.set_yticks([0, 1]); a.set_yticklabels(['axon (30 s high-passed)', 'axon, movement removed'] if c == 0 else ['', '']); a.set_ylim(1.5, -0.5); a.set_title(ttl, fontsize=8); a.set_xlabel('mean dF/F difference (grey = 95% of label permutations; filled |z| > 2)', fontsize=6.8)
    fig.suptitle('T6 | Trials with a lick (hit) vs without (miss): the axon trace aligned to the square, and the hit - miss difference with and without the concurrent movement removed', x=0.01, y=0.95, ha='left', fontsize=9.5)
    save(fig, 'T6_task_hit_vs_miss.png')


if __name__ == '__main__':
    fig_timeline(); fig_subsessions(); fig_epochs(); fig_bins(); fig_pre(); fig_hitmiss()
    with open(os.path.join(FIG, 'task_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['analysis', 'series', 'item', 'mean_effect_SD', 'pooled_z', 'n_sessions']); w.writeheader(); [w.writerow(r) for r in SUMMARY]
    print('pooled (Stouffer z over the %d non-flagged sessions); |z| > 2 listed:' % len(POOL))
    for r in SUMMARY:
        if np.isfinite(r['pooled_z']) and abs(r['pooled_z']) > 2: print('  %-32s %-8s %-18s effect %+.3f  z %+.2f (n=%d)' % (r['analysis'], r['series'], r['item'], r['mean_effect_SD'], r['pooled_z'], r['n_sessions']))
