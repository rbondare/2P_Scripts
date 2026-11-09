"""Figures + table of the behaviour-matched analysis (matched_behaviour.py, matched_validate.py, matched_pool.py).

usage: python matched_figures.py            -> OUT/M1_method_validation, M2_spout_vs_no_spout, M3_other_cases, M4_robustness, M5_roi (pdf + png), matched_results.csv, matched_pooled.csv
Conventions as in paper_figures: axon blue, behaviour hues, shape = animal, filled = unflagged (pooled), hollow = flagged (drawn, never pooled);
estimator markers: naive event-minus-control = open grey circle, regression (no matching) = open blue diamond, matched = filled, behaviour-only benchmark = black bar.
"""
import os, sys, pickle, csv, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import matched_pool as mp
import matched_behaviour as mb
from fig_common import (plt, W1, W2, MM, C_SIG, C_NOISE, C_REJ, C_BEH, BEH_LAB, INK, MUTED, FAINT, MARK, animal, short, is_flagged, FLAG, plabel, clean, fmt_p, hboot)
from matplotlib.lines import Line2D
DT = mb.DT
C_CTRL = '#9AA7BD'                    # matched control windows (axon): same hue family as the axons, lighter
C_NAIVE, C_REG, C_BENCH = '#6B7280', '#5B8DEF', '#111111'
OUTD = mb.OUT


def save(fig, name):
    for ext in ('pdf', 'png'): fig.savefig(os.path.join(OUTD, name + '.' + ext), facecolor='white')
    plt.close(fig); print('saved', name, flush=True)


testable = mp.testable


def mk(ax, x, y, s, filled=True, ms=4.6, color=C_SIG, alpha=1.0, zorder=4):
    m = MARK[animal(s)]
    if filled and not is_flagged(s): ax.plot(x, y, m, mfc=color, mec='white', mew=0.3, ms=ms, alpha=alpha, zorder=zorder)
    else: ax.plot(x, y, m, mfc='white', mec=color, mew=0.9, ms=ms, alpha=alpha, zorder=zorder)


def forest(ax, res, key, sessions, title=None, xlabel='dF/F', show_pool=True, labels=True, xlim=None, norm=False, fs=5.2):
    rows = []
    for s in sessions:
        c = mp.get(res, s, key)
        if c is not None: rows.append((s, c))
    if not rows: ax.axis('off'); ax.text(0.5, 0.5, 'no session', ha='center', va='center', transform=ax.transAxes, color=MUTED); return
    n = len(rows); ypos = {s: n - i + (1.4 if show_pool else 0) - 1 for i, (s, _) in enumerate(rows)}
    ax.axvline(0, color=FAINT, lw=0.7, zorder=0)
    for s, c in rows:
        y = ypos[s]; o = c['obs']; st = c['stats']['matched']; sc = res[s]['sd_P'] if norm else 1.0; ok = testable(c); a = 1.0 if ok else 0.35
        ax.plot([(st['nmean'] - 1.96 * st['nsd']) / sc, (st['nmean'] + 1.96 * st['nsd']) / sc], [y, y], color=FAINT, lw=3.2, solid_capstyle='butt', zorder=1, alpha=0.9)
        ss = c.get('stats_sim', {}).get('matched')
        if ss: ax.plot([ss['nmean'] / sc] * 2, [y - 0.32, y + 0.32], color=C_BENCH, lw=1.1, zorder=3, alpha=a)
        ax.plot(o['naive'] / sc, y + 0.0, 'o', mfc='white', mec=C_NAIVE, mew=0.7, ms=3.0, zorder=2, alpha=a)
        ax.plot(o['reg'] / sc, y, 'D', mfc='white', mec=C_REG, mew=0.7, ms=2.6, zorder=2, alpha=a)
        mk(ax, o['matched'] / sc, y, s, ms=4.6, alpha=a)
        if labels: ax.text(0.995, y, '%d/%d' % (o['n_matched'], o['n_focal']), transform=ax.get_yaxis_transform(), fontsize=fs - 0.4, color=MUTED if ok else '#C0392B', va='center', ha='right')
    ax.set_yticks([ypos[s] for s, _ in rows]); ax.set_yticklabels([short(s) for s, _ in rows], fontsize=fs)
    for t, (s, _) in zip(ax.get_yticklabels(), rows): t.set_color(MUTED if is_flagged(s) else INK)
    if show_pool:
        P = mp.pooled(res, key, sessions, norm=norm)
        if P is not None and len(P['sessions']) >= 2:
            y = 0.2; ax.errorbar(P['est'], y, xerr=[[P['est'] - P['lo']], [P['hi'] - P['est']]], fmt='D', color=C_SIG, mfc=C_SIG, ms=4.2, lw=1.1, capsize=0, zorder=5)
            ax.text(0.995, y, 'n=%d' % len(P['sessions']), transform=ax.get_yaxis_transform(), fontsize=fs - 0.4, color=INK, va='center', ha='right'); ax.set_yticks(list(ax.get_yticks()) + [y]); ax.set_yticklabels([t.get_text() for t in ax.get_yticklabels()][:-1] + ['pooled'], fontsize=fs)
    ax.set_ylim(-0.7, n + (1.6 if show_pool else 0) - 0.3); ax.set_xlabel(xlabel)
    if title: ax.set_title(title, loc='left')
    if xlim: ax.set_xlim(xlim)
    elif labels: x0, x1 = ax.get_xlim(); ax.set_xlim(x0, x1 + 0.28 * (x1 - x0))


def legend_forest(fig, x=0.01, y=0.005, shift=True):
    h = [Line2D([], [], marker='o', ls='', mfc='white', mec=C_NAIVE, ms=3.2, label='naive: event - all controls'), Line2D([], [], marker='D', ls='', mfc='white', mec=C_REG, ms=3, label='regression, no matching'),
         Line2D([], [], marker='o', ls='', mfc=C_SIG, mec='white', ms=4.6, label='matched (filled = unflagged, hollow = flagged)'), Line2D([], [], color=C_BENCH, lw=1.1, label='behaviour-only benchmark'),
         Line2D([], [], color=FAINT, lw=3.2, label='shift null 95%')]
    if not shift: h = h[:-1]
    fig.legend(handles=h, loc='lower left', bbox_to_anchor=(x, y), ncol=5, frameon=False, fontsize=5.4, handletextpad=0.3, columnspacing=1.2)


# ====================================================================================================== pooled aligned traces
def pooled_traces(res, key, sessions, series):
    F, R, D = [], [], []; tau = None
    for s in sessions:
        c = mp.get(res, s, key)
        if c is None or is_flagged(s) or not testable(c) or 'traces' not in c or series not in c['traces']: continue
        t = c['traces'][series]; F.append(t['focal']); R.append(t['ref']); D.append(t['diff_mean']); tau = c['tau']
    if not F: return None
    F, R, D = map(lambda a: np.stack(a, 0), (F, R, D)); sem = lambda a: np.nanstd(a, 0) / np.sqrt(np.isfinite(a).sum(0))
    return dict(tau=tau, F=np.nanmean(F, 0), R=np.nanmean(R, 0), D=np.nanmean(D, 0), Fs=sem(F), Rs=sem(R), Ds=sem(D), n=len(F))


def trace_panel(ax, T, color, spec=None, ylabel=None, ctrl=True, lw=0.9):
    if T is None: ax.axis('off'); return
    t = T['tau']
    if spec is not None: ax.axvspan(spec[0], spec[1], color='#F3F4F6', zorder=0)
    ax.fill_between(t, T['F'] - T['Fs'], T['F'] + T['Fs'], color=color, alpha=0.18, lw=0); ax.plot(t, T['F'], color=color, lw=lw)
    if ctrl: ax.fill_between(t, T['R'] - T['Rs'], T['R'] + T['Rs'], color=C_CTRL, alpha=0.25, lw=0); ax.plot(t, T['R'], color=C_CTRL, lw=lw, ls='-')
    ax.axhline(0, color=FAINT, lw=0.5, zorder=0); ax.axvline(0, color=MUTED, lw=0.5, ls=':', zorder=0)
    if ylabel: ax.set_ylabel(ylabel)


# ====================================================================================================== figure M1
def fig_m1(res, val):
    from matched_behaviour import Rec
    fig = plt.figure(figsize=(W2, 150 * MM)); gs = fig.add_gridspec(2, 3, left=0.07, right=0.985, top=0.95, bottom=0.085, wspace=0.42, hspace=0.5)
    name, key = 'AnimalRB24_261030_1744', 'post'; c = res[name]['contrasts'][key]; R = Rec(name); spec = c['spec']; a, b = int(round(spec[0] / DT)), int(round(spec[1] / DT))
    # (a) example segment
    inner = gs[0, 0:2].subgridspec(4, 1, hspace=0.12, height_ratios=[1.6, 1, 1, 1]); fb = c['focal_bins']; rb = np.concatenate(c['ref_bins']); t = R.S.t - R.S.t[0]
    t0 = 300.0; i0 = int(t0 / DT); i1 = i0 + int(200 / DT); axs = [fig.add_subplot(inner[i]) for i in range(4)]
    series = [('P', R.P, C_SIG, 'axon\ndF/F')] + [(bn, R.beh[bn], C_BEH[bn], BEH_LAB[bn].split(' (')[0].replace('locomotion', 'ball')) for bn in ('face', 'whisk') if bn in R.beh] + [('pupil', R.beh['pupil'], C_BEH['pupil'], 'pupil')]
    for ax, (nm, x, col, lab) in zip(axs, series):
        ax.plot(t[i0:i1], x[i0:i1], color=col, lw=0.6); ax.set_ylabel(lab, fontsize=5.5)
        for f in fb:
            if i0 <= f < i1: ax.axvspan(t[f + a], t[min(f + b, len(t) - 1)], color='#F59E0B', alpha=0.30, lw=0)
        for r in rb:
            if i0 <= r < i1: ax.axvspan(t[r + a], t[min(r + b, len(t) - 1)], color='#64748B', alpha=0.22, lw=0)
        clean(ax, bottom=(ax is axs[-1]))
        if ax is not axs[-1]: ax.spines['bottom'].set_visible(False); ax.tick_params(axis='x', bottom=False, labelbottom=False)
    axs[-1].set_xlabel('time in session (s)'); axs[0].set_title('Example RB24 10-30, 3-10 s after square onset: event windows orange, matched controls grey (darker = overlapping)', loc='left', fontsize=5.6)
    plabel(axs[0], 'a')
    # (b) feature space
    ax = fig.add_subplot(gs[0, 2]); F, y, ok = R.table(spec); pool = c['pool_bins']; ix = {n: i for i, n in enumerate(R.names)}; i1_, i2_ = ix['face'], ix['whisk']
    ax.plot(F[pool, i1_], F[pool, i2_], '.', color=FAINT, ms=2.2, mew=0, label='control pool')
    ax.plot(F[fb, i1_], F[fb, i2_], 'o', mfc='#F59E0B', mec='none', ms=3.0, alpha=0.8, label='event windows')
    for f, r in zip(fb, c['ref_bins']):
        for rr in r: ax.plot([F[f, i1_], F[rr, i1_]], [F[f, i2_], F[rr, i2_]], color='#64748B', lw=0.4, alpha=0.7)
    ax.set_xlabel('face motion (window mean, z)'); ax.set_ylabel('whisking (window mean, z)'); ax.legend(frameon=False, fontsize=5.2, loc='upper left', handletextpad=0.2); plabel(ax, 'b', dx=-26)
    # (c) balance
    ax = fig.add_subplot(gs[1, 0]); xs, ys = [], []
    for s, r in res.items():
        for k, cc in r['contrasts'].items():
            o = cc.get('obs')
            if not o or o.get('n_matched', 0) < 5 or 'smd_before' not in o: continue
            xs += list(np.abs(o['smd_before'])); ys += list(np.abs(o['smd_after']))
    ax.plot(xs, ys, 'o', mfc=C_SIG, mec='none', ms=2.4, alpha=0.5); lim = max(max(xs), 1.0) * 1.05; ax.plot([0, lim], [0, lim], color=FAINT, lw=0.6); ax.axhline(0.1, color=MUTED, lw=0.5, ls=':')
    ax.set_xlabel('|imbalance| before matching (SD)'); ax.set_ylabel('after matching (SD)'); ax.set_xscale('symlog', linthresh=0.3); ax.set_yscale('symlog', linthresh=0.3); ax.set_title('Balance: every feature, every session x contrast', loc='left', fontsize=6.2); plabel(ax, 'c')
    # (d) simulation: bias
    rows = [r for s in val.values() for r in s['rows']]
    ax = fig.add_subplot(gs[1, 1]); scen = ['lin', 'nonlin', 'noisy05', 'noisy10']; lab = ['linear', 'non-linear', 'noisy proxy\n0.5 SD', 'noisy proxy\n1 SD']; w = 0.2
    for j, (est, col) in enumerate((('naive', C_NAIVE), ('reg', C_REG), ('matched', C_SIG))):
        for i, sc in enumerate(scen):
            v = [(r[est] - r['truth']) / r['sd_win'] for r in rows if r['scen'] == sc]
            if v: ax.plot(np.full(len(v), i + (j - 1) * w) + np.random.default_rng(1).uniform(-0.04, 0.04, len(v)), v, 'o', mfc=col, mec='none', ms=2.4, alpha=0.55); ax.plot([i + (j - 1) * w - 0.08, i + (j - 1) * w + 0.08], [np.mean(v)] * 2, color=col, lw=1.4)
    ax.axhline(0, color=INK, lw=0.5); ax.set_xticks(range(len(scen))); ax.set_xticklabels(lab, fontsize=5.2); ax.set_ylabel('bias (SD of a window mean)'); ax.set_title('Axon = behaviour only: bias', loc='left', fontsize=6.2)
    ax.legend(handles=[Line2D([], [], color=C_NAIVE, lw=1.4, label='naive'), Line2D([], [], color=C_REG, lw=1.4, label='regression'), Line2D([], [], color=C_SIG, lw=1.4, label='matched')], frameon=False, fontsize=5.2, loc='upper left'); plabel(ax, 'd')
    # (e) recovery of an injected event response
    ax = fig.add_subplot(gs[1, 2]); lim = 0
    for est, col, mkr in (('naive', C_NAIVE, 'o'), ('reg', C_REG, 'D'), ('matched', C_SIG, 's')):
        xs = [r['truth'] / r['sd_win'] for r in rows if r['scen'] in ('resp025', 'resp050', 'resp100')]; ys = [r[est] / r['sd_win'] for r in rows if r['scen'] in ('resp025', 'resp050', 'resp100')]
        ax.plot(xs, ys, mkr, mfc=col, mec='none', ms=2.8, alpha=0.7); lim = max(lim, max(xs) if xs else 0)
    ax.plot([0, lim * 1.05], [0, lim * 1.05], color=INK, lw=0.5); ax.set_xlabel('injected response (SD of a window mean)'); ax.set_ylabel('estimated'); ax.set_title('Known response added: recovery', loc='left', fontsize=6.2); plabel(ax, 'e')
    save(fig, 'M1_method_validation')


# ====================================================================================================== figure M2
def fig_m2(res):
    fig = plt.figure(figsize=(W2, 215 * MM)); gs = fig.add_gridspec(5, 4, left=0.115, right=0.985, top=0.955, bottom=0.07, wspace=0.42, hspace=0.75, height_ratios=[1.0, 0.75, 1.5, 1.5, 1.25])
    groups = (('task', 'spout + reward', mp.TASK, 0), ('passive', 'no spout', mp.PASSIVE, 2))
    for gk, gname, sess, c0 in groups:
        for j, key in enumerate(('cue', 'post')):
            spec = mb.SPEC[key]; col = c0 + j
            T = pooled_traces(res, key, sess, 'P'); ax = fig.add_subplot(gs[0, col]); axb = fig.add_subplot(gs[1, col])
            ax.set_title('%s | %s%s' % (gname, 'cue 0-3 s' if key == 'cue' else '3-10 s', '' if T is None else '  (n=%d)' % T['n']), loc='left', fontsize=5.8)
            if col == 0: plabel(ax, 'a', dx=-26); ax.set_ylabel('axon dF/F'); axb.set_ylabel('behaviour (z)')
            if T is None:
                for a_ in (ax, axb): a_.axis('off')
                ax.text(0.5, 0.5, 'not testable:\nno behaviour-matched\ncontrol windows', ha='center', va='center', transform=ax.transAxes, fontsize=5.6, color='#C0392B'); continue
            trace_panel(ax, T, C_SIG, spec)
            for bn in ('face', 'whisk'):
                Tb = pooled_traces(res, key, sess, bn)
                if Tb is not None: axb.plot(Tb['tau'], Tb['F'], color=C_BEH[bn], lw=0.9); axb.plot(Tb['tau'], Tb['R'], color=C_BEH[bn], lw=0.9, ls=':')
            axb.axvspan(spec[0], spec[1], color='#F3F4F6', zorder=0); axb.axhline(0, color=FAINT, lw=0.5); axb.axvline(0, color=MUTED, lw=0.5, ls=':'); axb.set_xlabel('time from square onset (s)', fontsize=5.5)
            if col == 0: axb.text(0.02, 0.97, 'face motion', color=C_BEH['face'], transform=axb.transAxes, fontsize=5.0, va='top'); axb.text(0.02, 0.80, 'whisking', color=C_BEH['whisk'], transform=axb.transAxes, fontsize=5.0, va='top')
    fig.text(0.115, 0.915 - 0.0, '', fontsize=5)
    # forests: rows = task / passive, columns = cue / post (each 2 grid columns wide)
    for r, (gk, gname, sess, c0) in enumerate(groups):
        sub = gs[2 + r, :].subgridspec(1, 2, wspace=0.62)
        for k, key in enumerate(('cue', 'post')):
            ax = fig.add_subplot(sub[0, k]); forest(ax, res, key, sess, title='%s | %s' % (gname, mb.LABEL[key]), xlim=(-0.14, 0.14), labels=True, fs=5.6)
            if r == 0 and k == 0: plabel(ax, 'b', dx=-26)
    # pooled comparison
    ax = fig.add_subplot(gs[4, 0:2]); xt, xl = [], []; x = 0
    for key in ('cue', 'post'):
        for gname, sess, col in (('task', mp.TASK, C_SIG), ('passive', mp.PASSIVE, C_REJ)):
            P = mp.pooled(res, key, sess)
            if P is None: ax.text(x, 0, 'not\ntestable', ha='center', va='center', fontsize=5.0, color='#C0392B'); xt.append(x); xl.append('%s\n%s' % ('cue' if key == 'cue' else '3-10 s', gname)); x += 1; continue
            ax.errorbar(x, P['est'], yerr=[[P['est'] - P['lo']], [P['hi'] - P['est']]], fmt='D', color=col, ms=4.6, lw=1.1, capsize=0)
            for s_, v in zip(P['sessions'], P['vals']): mk(ax, x + 0.2, v, s_, ms=3.4, color=col)
            if 'bench' in P: ax.plot([x - 0.25, x + 0.25], [P['bench']] * 2, color=C_BENCH, lw=1.1)
            xt.append(x); xl.append('%s\n%s' % ('cue' if key == 'cue' else '3-10 s', gname)); x += 1
        x += 0.6
    ax.axhline(0, color=FAINT, lw=0.7); ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=5.4); ax.set_ylabel('matched effect (dF/F)'); ax.set_title('Pooled (unflagged); bar = behaviour-only benchmark', loc='left', fontsize=5.8); plabel(ax, 'c', dx=-26)
    ax = fig.add_subplot(gs[4, 2:4]); ax.axis('off'); lines = []
    for key in ('cue', 'post'):
        d = mp.group_diff(res, key, mp.TASK, mp.PASSIVE, field='matched')
        lines.append('%s: task - passive = %+.3f dF/F (%d vs %d sessions,\n   exact permutation %s, smallest possible p %.3f)' % (mb.LABEL[key], d['diff'], d['n_a'], d['n_b'], fmt_p(d['p']), d['min_p']) if d else '%s: task - passive not testable (passive has too few matched windows)' % mb.LABEL[key])
    for key in ('cue', 'post'):
        for gname, sess in (('task', mp.TASK), ('passive', mp.PASSIVE)):
            P = mp.pooled(res, key, sess)
            if P: lines.append('%s %s: %+.3f [%+.3f, %+.3f]; %d session(s), %d animal(s); p vs shift null %s,\n   p vs behaviour-only benchmark %s; %.0f%% of windows matched' % (gname, key, P['est'], P['lo'], P['hi'], len(P['sessions']), P['n_animals'], fmt_p(P.get('p_shift')).replace('p ', ''), fmt_p(P.get('p_sim')).replace('p ', ''), 100 * P['coverage']))
    ax.text(0, 1, '\n'.join(lines), va='top', fontsize=5.0, color=INK, linespacing=1.5, transform=ax.transAxes)
    legend_forest(fig); save(fig, 'M2_spout_vs_no_spout')


# ====================================================================================================== figure M3
def fig_m3(res):
    items = [('flash', mp.VISUAL), ('loom_expand', mp.VISUAL), ('loom_recede', mp.VISUAL), ('bar', mp.VISUAL), ('expand_vs_recede', mp.VISUAL), ('rew_vs_spont', mp.TASK), ('hit_vs_miss', mp.TASK)]
    fig = plt.figure(figsize=(W2, 175 * MM)); gs = fig.add_gridspec(2, 4, left=0.105, right=0.985, top=0.95, bottom=0.09, wspace=1.25, hspace=0.5)
    for i, (key, sess) in enumerate(items):
        ax = fig.add_subplot(gs[i // 4, i % 4]); forest(ax, res, key, sess, title=mb.LABEL[key], xlim=None, labels=True, fs=5.0)
        if i == 0: plabel(ax, 'a', dx=-34)
    ax = fig.add_subplot(gs[1, 3]); ax.axis('off'); lines = ['Pooled (unflagged, testable)', 'matched [95% CI]', 'n sessions; p shift / bench', '']
    short_lab = {'flash': 'flash', 'loom_expand': 'exp. loom', 'loom_recede': 'rec. loom', 'bar': 'bar', 'expand_vs_recede': 'exp - rec', 'rew_vs_spont': 'rew - spont', 'hit_vs_miss': 'hit - miss'}
    for key, sess in items:
        P = mp.pooled(res, key, sess)
        lines.append('%s %+.3f [%+.3f, %+.3f]\n   n=%d; shift %s, bench %s' % (short_lab[key], P['est'], P['lo'], P['hi'], len(P['sessions']), fmt_p(P.get('p_shift')).replace('p = ', 'p ').replace('p ', '', 1) if False else fmt_p(P.get('p_shift')), fmt_p(P.get('p_sim'))) if P else '%s: not testable' % short_lab[key])
    ax.text(-0.15, 1, '\n'.join(lines), va='top', fontsize=4.8, linespacing=1.35, transform=ax.transAxes)
    legend_forest(fig); save(fig, 'M3_other_cases')


# ====================================================================================================== overview (M0) + pooled table
OVERVIEW = [('cue', mp.TASK, 'square cue 0-3 s | spout + reward'), ('post', mp.TASK, 'square 3-10 s | spout + reward'), ('cue', mp.PASSIVE, 'square cue 0-3 s | no spout'), ('post', mp.PASSIVE, 'square 3-10 s | no spout'),
            ('rew_vs_spont', mp.TASK, 'rewarded - spontaneous lick bout'), ('hit_vs_miss', mp.TASK, 'hit - miss trial'),
            ('flash', mp.VISUAL, 'flash'), ('bar', mp.VISUAL, 'moving bar'), ('loom_expand', mp.VISUAL, 'expanding loom'), ('loom_recede', mp.VISUAL, 'receding loom'), ('expand_vs_recede', mp.VISUAL, 'expanding - receding')]


def pooled_rows(res, norm=False):
    rows = []
    for key, sess, lab in OVERVIEW:
        P = mp.pooled(res, key, sess, norm=norm)
        n_all = sum(1 for s in sess if mp.get(res, s, key) is not None and not is_flagged(s))
        rows.append(dict(contrast=key, label=lab, group='task' if sess is mp.TASK else 'passive' if sess is mp.PASSIVE else 'visual', n_unflagged_with_result=n_all, P=P))
    return rows


def write_pooled_csv(res):
    out = []
    for norm in (False, True):
        for r in pooled_rows(res, norm):
            P = r['P']; d = dict(units='SD of the session population trace' if norm else 'dF/F', contrast=r['contrast'], label=r['label'], group=r['group'], n_unflagged_with_result=r['n_unflagged_with_result'])
            if P is None: d.update(n_testable=0)
            else: d.update(n_testable=len(P['sessions']), n_animals=P['n_animals'], naive=P['naive'], reg=P['reg'], matched=P['est'], ci_lo=P['lo'], ci_hi=P['hi'], sessions_positive=P['k_pos'], bench=P.get('bench'), excess=P.get('excess'),
                           excess_lo=P.get('excess_lo'), excess_hi=P.get('excess_hi'), p_shift=P.get('p_shift'), p_behaviour_only=P.get('p_sim'), coverage=P['coverage'])
            out.append(d)
    keys = []
    for d in out:
        for k in d:
            if k not in keys: keys.append(k)
    with open(os.path.join(OUTD, 'matched_pooled.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for d in out: w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in d.items()})


def fig_m0(res):
    rows = pooled_rows(res); n = len(rows)
    fig, ax = plt.subplots(figsize=(W2, 105 * MM), gridspec_kw=dict(left=0.245, right=0.60, top=0.93, bottom=0.14))
    for i, r in enumerate(rows):
        y = n - i; P = r['P']
        if i in (4, 6): ax.axhline(y + 0.5, color=FAINT, lw=0.5)
        if P is None:
            ax.text(0.0, y, 'not testable (too few behaviour matches)', fontsize=5.4, color='#C0392B', va='center', ha='center'); continue
        ax.plot(P['naive'], y + 0.17, 'o', mfc='white', mec=C_NAIVE, mew=0.8, ms=3.8); ax.plot(P['reg'], y + 0.0, 'D', mfc='white', mec=C_REG, mew=0.8, ms=3.2)
        ax.errorbar(P['est'], y - 0.17, xerr=[[P['est'] - P['lo']], [P['hi'] - P['est']]], fmt='o', color=C_SIG, ms=4.4, lw=1.1, capsize=0)
        if 'bench' in P: ax.plot([P['bench']] * 2, [y - 0.4, y + 0.4], color=C_BENCH, lw=1.2)
        for s, v in zip(P['sessions'], P['vals']): mk(ax, v, y - 0.17, s, ms=2.6, alpha=0.55, zorder=3)
    ax.axvline(0, color=FAINT, lw=0.7); ax.set_yticks(range(n, 0, -1)); ax.set_yticklabels([r['label'] for r in rows], fontsize=5.8); ax.set_ylim(0.4, n + 0.7); ax.set_xlabel('event-window axon dF/F minus control windows')
    ax.set_title('What is left of the axon response when the controls have the same behaviour', loc='left', fontsize=6.6)
    for i, r in enumerate(rows):
        y = n - i; P = r['P']
        t = ('%d/%d sess, %d animal%s, %.0f%% matched\n   p shift %s | p benchmark %s' % (len(P['sessions']), r['n_unflagged_with_result'], P['n_animals'], '' if P['n_animals'] == 1 else 's', 100 * P['coverage'],
                                                                              fmt_p(P.get('p_shift')).replace('p ', ''), fmt_p(P.get('p_sim')).replace('p ', ''))) if P else 'n=0/%d sessions' % r['n_unflagged_with_result']
        ax.text(1.03, y, t, transform=ax.get_yaxis_transform(), fontsize=4.6, va='center', color=MUTED)
    legend_forest(fig, x=0.01, y=0.0, shift=False)
    save(fig, 'M0_overview')


# ====================================================================================================== figure M4 (robustness)
VAR_ORDER = [('default', 'default'), ('cal025', 'caliper 0.25'), ('cal100', 'caliper 1.0'), ('k1', '1 control'), ('k4', '4 controls'), ('fsd', '+ within-window SD'), ('fdyn', '+ within-window slope'),
             ('fhist', '+ 10 s history'), ('pupilin', 'pupil inside window'), ('pupilnone', 'no pupil'), ('facewhisk', 'face + whisk only'), ('noloco', 'no ball'), ('nohp', 'no 30 s high-pass'), ('hp10', '10 s high-pass'),
             ('nolagfit', 'no behaviour lag'), ('ext1', 'window -0/+2 s'), ('step5', 'controls >= 5 s apart')]
M4_ITEMS = [('cue', mp.TASK, 'square cue - task'), ('post', mp.TASK, 'square post - task'), ('cue', mp.PASSIVE, 'square cue - passive'), ('bar', mp.VISUAL, 'moving bar'),
            ('loom_expand', mp.VISUAL, 'expanding loom'), ('flash', mp.VISUAL, 'flash'), ('rew_vs_spont', mp.TASK, 'rewarded - spontaneous lick bout'), ('loom_recede', mp.VISUAL, 'receding loom')]


def fig_m4():
    caches = {v: mp.load('cache' if v == 'default' else 'cache_' + v) for v, _ in VAR_ORDER if os.path.isdir(os.path.join(OUTD, 'cache' if v == 'default' else 'cache_' + v))}
    fig = plt.figure(figsize=(W2, 200 * MM)); gs = fig.add_gridspec(len(M4_ITEMS), 2, left=0.085, right=0.985, top=0.96, bottom=0.17, wspace=0.12, hspace=0.35, width_ratios=[1, 1])
    names = [v for v, _ in VAR_ORDER if v in caches]; xs = np.arange(len(names))
    for i, (key, sess, title) in enumerate(M4_ITEMS):
        for j, (field, lab) in enumerate((('matched', 'matched effect (dF/F)'), ('excess', 'minus behaviour-only benchmark'))):
            ax = fig.add_subplot(gs[i, j])
            for x, v in zip(xs, names):
                res = caches[v]; P = mp.pooled(res, key, sess, n_boot=2000)
                if P is None: continue
                if field == 'matched': est, lo, hi = P['est'], P['lo'], P['hi']
                else:
                    if 'excess' not in P: continue
                    est, lo, hi = P['excess'], P['excess_lo'], P['excess_hi']
                ax.errorbar(x, est, yerr=[[est - lo], [hi - est]], fmt='o', color=C_SIG if v == 'default' else INK, ms=3.4, lw=0.9, capsize=0)
                ax.text(x, ax.get_ylim()[1], '', fontsize=4)
            ax.axhline(0, color=FAINT, lw=0.6); ax.set_xlim(-0.6, len(names) - 0.4)
            if i == 0: ax.set_title(lab, loc='left', fontsize=6)
            if j == 0: ax.set_ylabel(title, fontsize=5.4)
            if i < len(M4_ITEMS) - 1: ax.set_xticklabels([]); ax.tick_params(axis='x', labelbottom=False)
            else: ax.set_xticks(xs); ax.set_xticklabels([dict(VAR_ORDER)[v] for v in names], rotation=60, ha='right', fontsize=5.2)
    fig.text(0.01, 0.015, 'Each point: pooled over the unflagged sessions of the group (animal-clustered 95% bootstrap CI); one setting changed at a time relative to the default (blue).', fontsize=5.4, color=MUTED)
    save(fig, 'M4_robustness')


# ====================================================================================================== figure M5 (single ROIs)
def fig_m5(res):
    keys = [('cue', mp.TASK + mp.PASSIVE), ('post', mp.TASK + mp.PASSIVE), ('loom_expand', mp.VISUAL), ('flash', mp.VISUAL), ('rew_vs_spont', mp.TASK)]
    fig, axs = plt.subplots(1, len(keys), figsize=(W2, 78 * MM), gridspec_kw=dict(left=0.075, right=0.99, bottom=0.32, top=0.9, wspace=0.4))
    for ax, (key, sess) in zip(axs, keys):
        rows = [(s, mp.get(res, s, key)) for s in sess]; rows = [(s, c) for s, c in rows if c is not None and testable(c) and 'roi_effect' in c]; x = 0; xt = []
        for s, c in rows:
            e = c['roi_effect'] / res[s]['sd_P']
            ax.plot(np.full(len(e), x) + np.random.default_rng(0).uniform(-0.2, 0.2, len(e)), e, '.', color=C_SIG if not is_flagged(s) else C_REJ, ms=2.2, alpha=0.6)
            ax.plot([x - 0.3, x + 0.3], [np.mean(e)] * 2, color=INK, lw=1.2); xt.append(x); x += 1
        ax.axhline(0, color=FAINT, lw=0.6); ax.set_xticks(xt); ax.set_xticklabels([short(s) for s, _ in rows], rotation=70, ha='right', fontsize=4.6); ax.set_title(mb.LABEL[key], loc='left', fontsize=6)
        if key == 'cue': ax.set_ylabel('per-ROI matched effect\n(SD of the session population trace)')
    save(fig, 'M5_single_roi')


def main():
    res = mp.load('cache'); print(len(res), 'sessions')
    mp.write_csv(res, os.path.join(OUTD, 'matched_results.csv')); write_pooled_csv(res)
    vp = os.path.join(OUTD, 'validation.pkl'); val = pickle.load(open(vp, 'rb')) if os.path.exists(vp) else None
    if val: fig_m1(res, val)
    fig_m0(res); fig_m2(res); fig_m3(res); fig_m4(); fig_m5(res)


if __name__ == '__main__':
    main()
