"""Figures of cross-session analysis 3 (grey before / after the looming blocks). usage: python grey_figures.py  (reads grey_results.pkl, writes G*.png in plots/cross_session)

Colour = what the thing is: signal-like axons blue, noise-like ROIs amber, rejected suite2p ROIs grey, behaviour dark ink; session lines in the animal's colour (RB24 blue, RB23 red, RB22 teal).
"""
import os, pickle, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
FIG = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session'
plt.rcParams.update({'font.family': 'Arial' if any('Arial' in f.name for f in matplotlib.font_manager.fontManager.ttflist) else 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.75, 'axes.spines.top': False, 'axes.spines.right': False,
                     'xtick.direction': 'out', 'ytick.direction': 'out', 'figure.dpi': 110, 'savefig.dpi': 200, 'axes.titlesize': 8.5, 'axes.titleweight': 'normal'})
C_SIG, C_NOISE, C_REJ, INK, MUTED = '#1F5FBF', '#D98A00', '#8A8F98', '#222222', '#6B7280'
C_ANIMAL = {'RB24': '#1F5FBF', 'RB23': '#D1495B', 'RB22': '#1B998B'}
D = pickle.load(open(os.path.join(FIG, 'grey_results.pkl'), 'rb')); R = D['results']; CON = D['contrasts']
USE = [n for n in R if 'windows' in R[n]]


def short(nm): a, d, t = nm.replace('Animal', '').split('_'); return '%s  %s-%s  %s' % (a, d[2:4], d[4:6], t)
def animal(nm): return nm.replace('Animal', '')[:4]
USE = sorted(USE, key=lambda n: (animal(n)[::-1], n))
def save(fig, name): fig.savefig(os.path.join(FIG, name), bbox_inches='tight', facecolor='white'); plt.close(fig); print('saved', name)


def metric_val(w, pop, m):
    if not w['valid']: return np.nan
    if pop in ('sig', 'noise', 'rej'): d = w.get(pop); return np.nan if d is None else d[m]
    if pop == 'res': d = w.get('res'); return np.nan if d is None else d[m]
    if pop == 'beh': return w['beh'].get(m, np.nan)
    if pop == 'coupling': return w['coupling'].get(m, np.nan)


# ------------------------------------------------------------------ G1 time course around the looming blocks
def fig_timecourse():
    panels = [('sig', 'signal-like axons (mean best dF/F)', C_SIG), ('noise', 'noise-like ROIs (accepted on purpose)', C_NOISE), ('rej', 'rejected suite2p ROIs (random sets)', C_REJ),
              ('b_face', 'face motion', INK), ('b_whisk', 'whisking', INK), ('b_pupil', 'pupil area', INK), ('b_eyelid', 'eyelid opening', INK)]
    nb = 12; edges = np.linspace(0, 120, nb + 1); ctr = 0.5 * (edges[1:] + edges[:-1])
    fig, axs = plt.subplots(2, 4, figsize=(13, 5.6), gridspec_kw=dict(wspace=0.32, hspace=0.5)); axs = axs.ravel()
    for a, (key, ttl, col) in zip(axs, panels):
        allm = []
        for nm in USE:
            pk = R[nm]['peri']
            if key not in pk['pre']: continue
            rel = pk['rel']; pre = pk['pre'][key]; post = pk['post'][key]
            def bins(M): return np.array([[np.nanmean(M[j][(rel >= edges[i]) & (rel < edges[i + 1])]) if np.isfinite(M[j][(rel >= edges[i]) & (rel < edges[i + 1])]).any() else np.nan for i in range(nb)] for j in range(len(M))])
            bp, bq = bins(pre), bins(post); allv = np.concatenate([bp.ravel(), bq.ravel()]); mu, sd = np.nanmean(allv), np.nanstd(allv)
            if not np.isfinite(sd) or sd == 0: continue
            mp, mq = (np.nanmean(bp, 0) - mu) / sd, (np.nanmean(bq, 0) - mu) / sd
            a.plot(ctr - 120, mp, color=C_ANIMAL[animal(nm)], lw=0.7, alpha=0.45); a.plot(ctr + 14, mq, color=C_ANIMAL[animal(nm)], lw=0.7, alpha=0.45)
            allm.append(np.concatenate([mp, mq]))
        if allm:
            m = np.nanmean(np.array(allm), 0); se = np.nanstd(np.array(allm), 0) / np.sqrt(len(allm))
            for sl, off in ((slice(0, nb), -120), (slice(nb, 2 * nb), 14)):
                a.fill_between(ctr + off, (m - se)[sl], (m + se)[sl], color=col, alpha=0.18, lw=0); a.plot(ctr + off, m[sl], color=col, lw=2.0)
        a.axvspan(0, 14, color='#E5E7EB', lw=0); a.text(7, 0.02, 'loom\nblock', transform=a.get_xaxis_transform(), ha='center', va='bottom', fontsize=6, color=MUTED)
        a.axhline(0, color=MUTED, lw=0.4, ls=':'); a.set_title(ttl, fontsize=8, color=INK)
        a.set_xticks([-120, -60, 0, 14, 74, 134]); a.set_xticklabels(['-120', '-60', '', '', '+60', '+120']); a.set_xlabel('s before the block start | after the block end', fontsize=6.5); a.set_ylabel('z (session grey)', fontsize=7)
    ax = axs[-1]; ax.axis('off')
    ax.text(0, 1, 'How to read', fontweight='bold', va='top', transform=ax.transAxes, fontsize=8)
    ax.text(0, 0.9, 'Mean over the loom blocks of a session (thin line,\nanimal colour) and over sessions (thick, +- s.e.m.)\nin 10 s bins; 120 s of grey before and after.\nValues are z-scored per session over all its\ngrey bins so sessions are comparable.\n\nA response to the LOOM would be a step between\nthe left and the right half; slow time-in-session\ndrift is a tilt across both halves.\nThe grey of one block and the next block are\none continuous 270 s stretch when looms follow\neach other.\n\nsessions: ' + ', '.join(short(n) for n in USE), fontsize=6.6, va='top', transform=ax.transAxes, color=INK)
    fig.suptitle('G1 | Spontaneous grey before and after the looming blocks: axons and behaviour  (%d sessions, %d loom blocks)' % (len(USE), sum(len(R[n]['peri_blocks']) for n in USE)), x=0.01, ha='left', fontsize=9.5)
    save(fig, 'G1_grey_timecourse.png')


# ------------------------------------------------------------------ G2 paired pre / post per block
def fig_paired():
    items = [('sig', 'mean', 'signal-like: mean dF/F', C_SIG), ('sig', 'sd', 'signal-like: fluctuation SD (1-30 s)', C_SIG), ('sig', 'rate', 'signal-like: event rate', C_SIG), ('sig', 'active', 'signal-like: time above 1 SD', C_SIG),
             ('res', 'mean', 'signal-like after removing\nface + whisk + ball (mean)', C_SIG), ('beh', 'face', 'face motion', INK), ('beh', 'whisk', 'whisking', INK), ('beh', 'pupil', 'pupil area', INK),
             ('beh', 'eyelid', 'eyelid opening', INK), ('coupling', 'r_face', 'axon - face-motion r in the window', INK)]
    fig, axs = plt.subplots(2, 5, figsize=(13, 5.8), gridspec_kw=dict(wspace=0.45, hspace=0.55)); axs = axs.ravel()
    for a, (pop, m, ttl, col) in zip(axs, items):
        sesm = []
        for nm in USE:
            W = R[nm]['windows']; vals = np.array([metric_val(w, pop, m) for w in W]); sc, mu = np.nanstd(vals), np.nanmean(vals)
            if not np.isfinite(sc) or sc == 0: continue
            nbk = max(w['block'] for w in W) + 1; P = {(w['block'], w['kind']): (metric_val(w, pop, m) - mu) / sc for w in W}
            dif = []
            for j in range(nbk):
                x0, x1 = P.get((j, 'pre'), np.nan), P.get((j, 'post'), np.nan)
                if np.isfinite(x0) and np.isfinite(x1):
                    a.plot([0, 1], [x0, x1], color=C_ANIMAL[animal(nm)], lw=0.6, alpha=0.5); a.plot([0, 1], [x0, x1], 'o', color=C_ANIMAL[animal(nm)], ms=2.3, alpha=0.6, mew=0); dif.append((x0, x1))
            if dif:
                d = np.array(dif).mean(0); a.plot([0, 1], d, color=C_ANIMAL[animal(nm)], lw=1.8, marker='o', ms=5, mec='white', mew=0.6); sesm.append(d)
        c = CON.get((pop, m), {}).get('loom')
        txt = '' if not c else 'mean %+.2f SD  [%+.2f, %+.2f]\nsign-flip p = %.3f, %d sessions' % (c['mean'], c['ci'][0], c['ci'][1], c['p'], c['n'])
        a.set_title(ttl, fontsize=7.5, color=col); a.text(0.5, 1.0, txt, transform=a.transAxes, ha='center', va='bottom', fontsize=5.8, color=MUTED) if False else a.set_xlabel(txt, fontsize=5.8, color=INK)
        a.set_xticks([0, 1]); a.set_xticklabels(['before\nloom', 'after\nloom']); a.set_xlim(-0.3, 1.3); a.axhline(0, color=MUTED, lw=0.4, ls=':')
        if m in ('mean',) or True: a.set_ylabel('z within session', fontsize=7)
    fig.suptitle('G2 | The same block: grey before vs grey after (thin = block, thick = session mean; z-scored within session over its windows)', x=0.01, ha='left', fontsize=9.5)
    save(fig, 'G2_grey_pre_post_paired.png')


# ------------------------------------------------------------------ G3 forest: across a loom block vs across plain grey vs difference
def fig_forest():
    rows = [('sig', 'mean', 'signal-like: mean dF/F', C_SIG), ('sig', 'sd', 'signal-like: fluctuation SD', C_SIG), ('sig', 'rate', 'signal-like: event rate', C_SIG), ('sig', 'active', 'signal-like: time above 1 SD', C_SIG),
            ('res', 'mean', 'signal-like, behaviour removed: mean', C_SIG), ('res', 'sd', 'signal-like, behaviour removed: SD', C_SIG),
            ('noise', 'sd', 'noise-like ROIs: fluctuation SD', C_NOISE), ('rej', 'sd', 'rejected ROIs: fluctuation SD', C_REJ), ('rej', 'rate', 'rejected ROIs: event rate', C_REJ),
            ('beh', 'face', 'face motion', INK), ('beh', 'whisk', 'whisking', INK), ('beh', 'pupil', 'pupil area', INK), ('beh', 'eyelid', 'eyelid opening', INK),
            ('coupling', 'r_face', 'axon - face r', INK), ('coupling', 'slope_face', 'axon per face-motion unit (slope)', INK)]
    cols = [('loom', 'across a loom block\n(after - before)'), ('grey', 'across plain grey\n(next block\'s before - this after)'), ('diff', 'difference\n(loom - plain grey)')]
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 5.8), sharey=True, gridspec_kw=dict(wspace=0.3)); y = np.arange(len(rows))
    for a, (ck, ctl) in zip(axs, cols):
        for i, (pop, m, lab, col) in enumerate(rows):
            c = CON.get((pop, m), {}).get(ck)
            if not c: continue
            ps = list(c['per_session'].values())
            a.plot(ps, [i] * len(ps), '.', color=col, ms=4, alpha=0.5)
            a.plot([c['ci'][0], c['ci'][1]], [i, i], color=col, lw=2.2, solid_capstyle='round'); a.plot(c['mean'], i, 'o', color=col, ms=6, mec='white', mew=0.7)
            a.text(1.01, i, 'p=%.2f' % c['p'], transform=a.get_yaxis_transform(), fontsize=6, va='center', color=INK if c['p'] < 0.05 else MUTED)
        a.axvline(0, color=INK, lw=0.6); a.set_title(ctl, fontsize=8); a.set_xlabel('standardised change (SD of the metric across windows)', fontsize=7)
    axs[0].set_yticks(y); axs[0].set_yticklabels([r[2] for r in rows], fontsize=7.5); axs[0].set_ylim(len(rows) - 0.4, -0.6)
    for t, r in zip(axs[0].get_yticklabels(), rows): t.set_color(r[3])
    fig.suptitle('G3 | Does a loom block change the grey that follows it beyond elapsed time?  (dots = sessions, bar = bootstrap 95% CI over sessions, p = exact sign-flip over sessions)', x=0.01, ha='left', fontsize=9)
    save(fig, 'G3_grey_contrasts.png')


# ------------------------------------------------------------------ G4 epochs: grey vs spontaneous block vs loom vs other stimuli
def fig_epochs():
    eps = ['spont block', 'grey (before / after loom)', 'loom blocks', 'other stimulus blocks']; short_e = ['spont.\nblock', 'grey\n(loom pre / post)', 'loom\nblocks', 'other\nstimuli']
    items = [('sd', 'signal-like fluctuation SD (dF/F)'), ('r_face', 'axon - face-motion r'), ('slope_face', 'axon per face-motion unit'), ('mean_face', 'face motion (z, session)')]
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.6), gridspec_kw=dict(wspace=0.38))
    MK = dict(zip(USE, ['o', 's', '^', 'D', 'v', 'P']))
    for a, (k, ttl) in zip(axs, items):
        for nm in USE:
            E = R[nm].get('epochs', {}); xs, ys = [], []
            for i, e in enumerate(eps):
                if e in E and k in E[e]: xs.append(i); ys.append(E[e][k])
            if len(xs) >= 2: a.plot(xs, ys, '-' + MK[nm], color=C_ANIMAL[animal(nm)], lw=0.9, ms=4.5, alpha=0.8)
            elif xs: a.plot(xs, ys, MK[nm], color=C_ANIMAL[animal(nm)], ms=4.5)
        a.set_xticks(range(4)); a.set_xticklabels(short_e, fontsize=6.5); a.set_title(ttl, fontsize=8); a.set_xlim(-0.4, 3.4)
        if k in ('r_face',): a.axhline(0, color=MUTED, lw=0.4, ls=':')
    for nm in USE: axs[0].plot([], [], '-' + MK[nm], color=C_ANIMAL[animal(nm)], ms=4, label=short(nm))
    axs[0].legend(frameon=False, fontsize=5.8, loc='upper left', bbox_to_anchor=(0, -0.28), ncol=2)
    fig.suptitle('G4 | The same signal in four contexts of one session: spontaneous block, grey before / after the loom, loom blocks, other stimulus blocks', x=0.01, ha='left', fontsize=9.5)
    save(fig, 'G4_grey_vs_other_epochs.png')


# ------------------------------------------------------------------ G5 timeline of the windows in session time
def fig_timeline():
    n = len(USE); fig, axs = plt.subplots(3, n, figsize=(2.9 * n, 6.4), sharex='col', gridspec_kw=dict(hspace=0.15, wspace=0.34))
    for j, nm in enumerate(USE):
        r = R[nm]; W = [w for w in r['windows'] if w['valid']]
        for i, (pop, m, ttl, col) in enumerate((('sig', 'sd', 'axon SD', C_SIG), ('sig', 'mean', 'axon mean dF/F', C_SIG), ('beh', 'face', 'face motion z', INK))):
            a = axs[i, j]
            for b in r['blocks']:
                if b.get('pre_ok') and b.get('post_ok'): a.axvspan(b['t0'] / 60, b['t1'] / 60, color='#E5E7EB', lw=0)
            for w in W:
                v = metric_val(w, pop, m)
                if not np.isfinite(v): continue
                tc = 0.5 * (w['t0'] + w['t1']) / 60; a.plot([w['t0'] / 60, w['t1'] / 60], [v, v], color=col, lw=2.2, solid_capstyle='butt', alpha=0.95 if w['kind'] == 'post' else 0.5)
            a.set_ylabel(ttl if j == 0 else '', fontsize=7)
            if i == 0: a.set_title(short(nm), fontsize=7.5, color=C_ANIMAL[animal(nm)])
            if i == 2: a.set_xlabel('session time (min)', fontsize=7)
    fig.suptitle('G5 | Each 100 s grey window in session time (bar = window; dark = after a loom block, light = before; shaded = loom blocks)', x=0.01, ha='left', fontsize=9.5)
    save(fig, 'G5_grey_windows_timeline.png')


if __name__ == '__main__':
    fig_timecourse(); fig_paired(); fig_forest(); fig_epochs(); fig_timeline()
