"""Figure 11 - the same square with a spout + reward (lick task) vs without (passive), and engaged vs not engaged trials, with whisking / pupil / face / ball regressed out.
a raw excerpts: task (RB24 10-30 1744) and passive (RB24 11-04 1651): axon mean, the part predicted by the behaviours, the residual, behaviours, square onsets / licks
b square-aligned averages: axon mean raw, axon mean with behaviours regressed out (face + whisking + pupil + ball), behaviours; task sessions | passive sessions
c, d joint-model contrasts per session: cue - pre (0-3 vs -3..0 s) and post - pre (3-10 vs -3..0 s), raw and with the behaviours regressed out, task vs passive
e engagement inside the task sessions: hit-trial cycles minus miss-trial cycles, and hit rate of the surrounding 10 trials, raw / regressed, with the per-session shuffle null
f how much of the axon signal the behaviours explain, task vs passive
usage: python fig11_task_engagement.py (after python context_state.py) -> plots/summary_figures/Fig11_task_engagement.pdf/.png ; ledger entries
"""
import os, pickle, itertools, numpy as np
from fig_common import *
import xs_common as xc
from context_state import SESS, WIN

C_MV = '#8FB3E8'; C_RES = '#0B2F73'; C_FIT = '#9AA3B2'
TASK = [n for n in SESS if n in TASK_SESSIONS]; PASS = [n for n in SESS if n not in TASK_SESSIONS]
EX_T, EX_P = 'AnimalRB24_261030_1744', 'AnimalRB24_261104_1651'


def wrap(s, w=212):
    import textwrap
    return textwrap.fill(s.replace('\n', ' '), w)


def load(n): return pickle.load(open(os.path.join(CACHE, 'ctx_%s.pkl' % n), 'rb'))


def perm_diff(a, b):
    """exact two-sided permutation test of mean(a) - mean(b) over sessions"""
    v = np.array(list(a) + list(b), float); na = len(a); obs = abs(np.mean(a) - np.mean(b)); cnt = 0; tot = 0
    for idx in itertools.combinations(range(len(v)), na):
        m = np.zeros(len(v), bool); m[list(idx)] = True; tot += 1; cnt += abs(v[m].mean() - v[~m].mean()) >= obs - 1e-12
    return cnt / tot


def sm(x, s=0.5): return xc.nsmooth(np.asarray(x, float), s / xc.DT)


def raw_panel(fig, spec, name, letter, ttl, d):
    t = np.asarray(d['t'], float); ex = d['ex']; t0, t1 = 120.0, 420.0; sl = (t >= t0) & (t < t1); tt = t[sl] - t0
    gs = spec.subgridspec(6, 1, height_ratios=[0.9, 2.0, 1.3, 1.0, 1.0, 1.0], hspace=0.14)
    a0 = fig.add_subplot(gs[0]); ons = np.asarray(ex['ons'], float)
    for k, o in enumerate(ons):
        if t0 <= o < t1:
            a0.axvspan(o - t0, o - t0 + 3, color=C_REJ, alpha=0.35, lw=0)
            if d['task']: a0.plot(o - t0, 1.35, 'v', color=INK, mfc=INK if d['hit'][k] == 1 else 'white', ms=3.0, mew=0.6, clip_on=False)
    if d['task']: lk = d['licks']; lk = lk[(lk >= t0) & (lk < t1)] - t0; a0.vlines(lk, 0, 1, color=INK, lw=0.3)
    else: a0.text(0.5, 0.5, 'no spout: no licks, no reward', transform=a0.transAxes, ha='center', va='center', fontsize=5.0, color=MUTED)
    a0.set_ylim(0, 1.5); a0.set_yticks([]); a0.spines['left'].set_visible(False); a0.spines['bottom'].set_visible(False); a0.tick_params(labelbottom=False, bottom=False)
    a0.set_ylabel('licks' if d['task'] else 'spout', fontsize=5.2, rotation=0, ha='right', va='center', labelpad=3); plabel(a0, letter, dx=-32, dy=12)
    a0.set_title(ttl, fontsize=5.8, loc='left', fontweight='normal', pad=6)
    a1 = fig.add_subplot(gs[1], sharex=a0); a1.plot(tt, sm(ex['P'])[sl], color=C_SIG, lw=1.0, label='axon mean'); a1.plot(tt, sm(ex['fitB'])[sl], color=C_FIT, lw=0.9, label='predicted from behaviours')
    a1.set_yticks([]); a1.spines['left'].set_visible(False); a1.spines['bottom'].set_visible(False); a1.tick_params(labelbottom=False, bottom=False); a1.set_ylabel('axon\ndF/F', fontsize=5.2, rotation=0, ha='right', va='center', labelpad=3)
    a1.legend(frameon=False, fontsize=4.6, loc='upper right', ncol=2, handlelength=1.0, columnspacing=0.8, borderpad=0.1, bbox_to_anchor=(1.0, 1.18))
    a2 = fig.add_subplot(gs[2], sharex=a0); a2.plot(tt, sm(ex['res'])[sl], color=C_RES, lw=0.9); a2.axhline(0, color=MUTED, lw=0.3); a2.set_yticks([]); a2.spines['left'].set_visible(False); a2.spines['bottom'].set_visible(False); a2.tick_params(labelbottom=False, bottom=False)
    a2.set_ylabel('residual', fontsize=5.2, rotation=0, ha='right', va='center', labelpad=3)
    for k, (b, lab) in enumerate((('whisk', 'whisking'), ('face', 'face motion'), ('pupil', 'pupil'))):
        a = fig.add_subplot(gs[3 + k], sharex=a0)
        if b in ex['behs']: a.plot(tt, sm(ex['behs'][b], 0.3)[sl], color=C_BEH[b], lw=0.6)
        a.set_yticks([]); a.spines['left'].set_visible(False); a.set_ylabel(lab, fontsize=5.2, rotation=0, ha='right', va='center', labelpad=3)
        if k < 2: a.spines['bottom'].set_visible(False); a.tick_params(labelbottom=False, bottom=False)
    a.set_xlabel('time (s)'); a.set_xlim(0, t1 - t0)


def group_stat(ax, x0, vals, names, col, ms=3.8, rng=None, mean_only=False):
    """session dots + pooled bar for one condition; returns hboot result of the unflagged sessions"""
    v = {n: vals[n] for n in names if np.isfinite(vals.get(n, np.nan))}
    for n, y in v.items(): sdot(ax, x0 + rng.uniform(-0.17, 0.17), y, n, color=col, hollow_color=MUTED, ms=ms)
    h = hboot([v[n] for n in names if n in v and not is_flagged(n)], [n for n in names if n in v and not is_flagged(n)], n_boot=6000)
    if h['n_sessions']: ax.plot([x0 - 0.32, x0 + 0.32], [h['est']] * 2, color=INK, lw=1.4, zorder=5); ax.plot([x0, x0], [h['lo'], h['hi']], color=INK, lw=0.9, zorder=5)
    return h


def main():
    D = {n: load(n) for n in SESS}; fig = plt.figure(figsize=(W2, 246 * MM)); rng = np.random.default_rng(8)
    gs = fig.add_gridspec(3, 1, height_ratios=[52, 78, 58], hspace=0.30, left=0.105, right=0.985, top=0.967, bottom=0.105)
    top = gs[0].subgridspec(1, 2, wspace=0.20)
    raw_panel(fig, top[0], EX_T, 'a', '%s: lick task (spout + reward)' % short(EX_T), D[EX_T])
    raw_panel(fig, top[1], EX_P, 'b', '%s: passive square (no spout)' % short(EX_P), D[EX_P])

    # ------------------------------------------------ c  square-aligned averages
    mid = gs[1].subgridspec(3, 2, hspace=0.16, wspace=0.14, height_ratios=[1.3, 1.3, 1.0]); lims = {}
    for j, (grp, nm_g) in enumerate(((TASK, 'lick task (spout + reward)'), (PASS, 'passive (no spout)'))):
        for i, (key, lab, col) in enumerate((('P_tr', 'axon mean dF/F', C_SIG), ('res_tr', 'behaviours\nregressed out', C_RES))):
            a = fig.add_subplot(mid[i, j], sharey=lims.get(i)); lims[i] = lims.get(i, a); rel = D[grp[0]]['rel']; M = []
            for n in grp:
                m = np.nanmean(D[n][key], 0); fl = is_flagged(n); a.plot(rel, m, color=FAINT if fl else col, lw=0.5, alpha=1 if fl else 0.45, ls='--' if fl else '-', zorder=2)
                if not fl: M.append(m)
            a.plot(rel, np.nanmean(M, 0), color=col, lw=1.7, zorder=4); a.axvspan(0, 3, color=C_REJ, alpha=0.25, lw=0, zorder=0); a.axvspan(3, 10, color=C_REJ, alpha=0.10, lw=0, zorder=0); a.axvline(0, color=INK, lw=0.5, ls=':'); a.axhline(0, color=MUTED, lw=0.4)
            a.tick_params(labelbottom=False)
            if j == 0: a.set_ylabel(lab, fontsize=5.4)
            else: a.tick_params(labelleft=False)
            if i == 0: a.set_title('%s: %d unflagged sessions (thin) + %d flagged (dashed)' % (nm_g, len(M), len(grp) - len(M)), fontsize=5.6, loc='left', fontweight='bold')
            if i == 0 and j == 0: plabel(a, 'c', dx=-34, dy=8)
        a = fig.add_subplot(mid[2, j], sharey=lims.get(2)); lims[2] = lims.get(2, a)
        for b, lab in (('pupil', 'pupil'), ('face', 'face motion'), ('whisk', 'whisking')):
            ms = [np.nanmean(D[n]['beh_tr'][b], 0) for n in grp if b in D[n]['beh_tr'] and not is_flagged(n)]
            if ms: ms = np.array(ms); a.plot(rel, np.nanmean(ms, 0), color=C_BEH[b], lw=1.2, label=lab if j == 0 else None)
        a.axvspan(0, 3, color=C_REJ, alpha=0.25, lw=0, zorder=0); a.axvspan(3, 10, color=C_REJ, alpha=0.10, lw=0, zorder=0); a.axvline(0, color=INK, lw=0.5, ls=':'); a.axhline(0, color=MUTED, lw=0.4); a.set_xlabel('time from square onset (s)', fontsize=5.6)
        if j == 0: a.set_ylabel('behaviour (z)', fontsize=5.4); a.legend(frameon=False, fontsize=4.6, loc='upper right', handlelength=1.0, borderpad=0.1, labelspacing=0.15, ncol=3, columnspacing=0.8)
        else: a.tick_params(labelleft=False)

    # ------------------------------------------------ d..g  statistics
    bot = gs[2].subgridspec(1, 4, wspace=0.50, width_ratios=[1.0, 1.0, 1.15, 0.8]); ac = fig.add_subplot(bot[0]); ad = fig.add_subplot(bot[1]); ae = fig.add_subplot(bot[2]); af = fig.add_subplot(bot[3]); res = {}
    for a, ci, ttl, letter, ledkey in ((ac, 0, 'cue (0-3 s) − pre (−3..0 s)', 'd', 'cue'), (ad, 1, 'post (3-10 s) − pre', 'e', 'post')):
        pos = {('task', 'raw'): 0, ('task', 'all'): 1.1, ('pass', 'raw'): 3.0, ('pass', 'all'): 4.1}; H = {}
        for (g, mdl), x0 in pos.items():
            names = TASK if g == 'task' else PASS; vals = {n: D[n]['epochs'][mdl]['contrast'][ci] for n in names}
            H[(g, mdl)] = group_stat(a, x0, vals, names, C_SIG if mdl == 'raw' else C_RES, rng=rng); ksig = sum(1 for n in names if not is_flagged(n) and abs(D[n]['epochs'][mdl]['c_z'][ci]) > 2.33 and D[n]['epochs'][mdl]['contrast'][ci] > 0)
            h = H[(g, mdl)]; a.text(x0, 1.0, '%d/%d' % (ksig, h['n_sessions']), transform=a.get_xaxis_transform(), ha='center', va='bottom', fontsize=4.6, color=MUTED)
            ledger_add('Fig11', 'sq_%s_%s_%s' % (ledkey, g, mdl), 'axon response to the square, %s, %s sessions (%s)' % ('cue − pre' if ci == 0 else 'post − pre', 'lick-task' if g == 'task' else 'passive', 'raw' if mdl == 'raw' else 'whisking, face, pupil, ball regressed'), h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'],
                       'ambiguous' if (h['lo'] > 0 and ksig >= 1) else 'not seen', note='%d/%d unflagged sessions with z > 2.33 (shift null); %s' % (ksig, h['n_sessions'], 'one animal' if h['n_animals'] < 2 else '%d animals' % h['n_animals']), unit='dF/F (high-passed)')
        res[ledkey] = H
        # task - passive (unflagged sessions, raw and regressed)
        for mdl, y0 in (('raw', 0.93), ('all', 0.84)):
            va = [D[n]['epochs'][mdl]['contrast'][ci] for n in TASK if not is_flagged(n)]; vb = [D[n]['epochs'][mdl]['contrast'][ci] for n in PASS if not is_flagged(n)]
            a.text(0.03, y0, '%s: task − passive %+.3f, perm. p = %.2f' % ('raw' if mdl == 'raw' else 'regressed', np.mean(va) - np.mean(vb), perm_diff(va, vb)), transform=a.transAxes, fontsize=4.5, va='top')
            ledger_add('Fig11', 'sq_%s_diff_%s' % (ledkey, mdl), 'square response %s: lick-task minus passive sessions (%s)' % ('cue − pre' if ci == 0 else 'post − pre', 'raw' if mdl == 'raw' else 'behaviours regressed'), np.mean(va) - np.mean(vb), np.nan, np.nan, len(va) + len(vb), 2, 'not seen' if perm_diff(va, vb) > 0.05 else 'ambiguous',
                       note='exact permutation over sessions (%d vs %d): p = %.2f; groups differ in animal and training day' % (len(va), len(vb), perm_diff(va, vb)), unit='dF/F')
        a.axhline(0, color=MUTED, lw=0.4); a.set_xticks(list(pos.values())); a.set_xticklabels(['task\nraw', 'task\nreg.', 'passive\nraw', 'passive\nreg.'], fontsize=4.8); a.set_xlim(-0.6, 4.7); a.set_ylim(-0.06, 0.2 if ci == 1 else 0.1)
        a.set_ylabel('dF/F (high-passed)', fontsize=5.4); plabel(a, letter, dx=-24, dy=12); a.set_title(ttl, fontsize=5.6, loc='left', fontweight='bold', pad=11)
    # e  engagement
    ae.set_title('engagement (task sessions)', fontsize=5.6, loc='left', fontweight='bold', pad=11); plabel(ae, 'f', dx=-24, dy=12); cols = [('H', 'raw', 0), ('H', 'all', 1), ('E', 'raw', 2.8), ('E', 'all', 3.8)]; H = {}
    for lab, mdl, x0 in cols:
        vals = {}
        for n in TASK:
            r = D[n]['engage'][mdl].get(lab)
            if not r: continue
            vals[n] = r['beta']; ae.plot([x0 + 0.0, x0 + 0.0], [-1.96 * r['null_sd'] + r['null_mean'], 1.96 * r['null_sd'] + r['null_mean']], color=C_REJ, lw=0, alpha=0)
        names = [n for n in TASK if n in vals]; xs_ = {n: x0 + (i - (len(names) - 1) / 2) * 0.17 for i, n in enumerate(names)}
        for n in names:
            r = D[n]['engage'][mdl][lab]; ae.plot([xs_[n]] * 2, [r['null_mean'] - 1.96 * r['null_sd'], r['null_mean'] + 1.96 * r['null_sd']], color=C_REJ, lw=1.5, alpha=0.7, solid_capstyle='butt', zorder=1); sdot(ae, xs_[n], vals[n], n, color=C_SIG if mdl == 'raw' else C_RES, hollow_color=MUTED, ms=3.6)
        un = [n for n in names if not is_flagged(n)]; h = hboot([vals[n] for n in un], un, n_boot=6000); H[(lab, mdl)] = h
        ae.plot([x0 - 0.35, x0 + 0.35], [h['est']] * 2, color=INK, lw=1.4, zorder=5)
        ksig = sum(1 for n in un if abs(D[n]['engage'][mdl][lab]['z']) > 2.33)
        ledger_add('Fig11', 'engage_%s_%s' % (lab, mdl), 'axon activity in hit-trial vs miss-trial cycles (%s)' % ('raw' if mdl == 'raw' else 'behaviours regressed') if lab == 'H' else 'axon activity vs hit rate of the surrounding 10 trials (%s)' % ('raw' if mdl == 'raw' else 'behaviours regressed'), h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'],
                   'ambiguous' if (h['lo'] > 0 and mdl == 'raw') else 'not seen', note='%d/%d unflagged task sessions with |z| > 2.33 vs shifted hit sequence; one animal (RB24)' % (ksig, h['n_sessions']), unit='dF/F per hit-vs-miss' if lab == 'H' else 'dF/F per SD of hit rate')
    ae.axhline(0, color=MUTED, lw=0.4); ae.set_xticks([0, 1, 2.8, 3.8]); ae.set_xticklabels(['hit−miss\nraw', 'hit−miss\nreg.', 'hit rate\nraw', 'hit rate\nreg.'], fontsize=4.8); ae.set_xlim(-0.6, 4.4); ae.set_ylabel('dF/F (high-passed)', fontsize=5.4)
    ae.set_ylim(-0.065, 0.085); ae.text(0.03, 0.02, 'grey bar = 95% of shifted hit\nsequences, one per session', transform=ae.transAxes, fontsize=4.5, color=MUTED, va='bottom')
    # f  variance explained by behaviours
    for g, names, x0 in (('task', TASK, 0), ('pass', PASS, 1.4)):
        vals = {n: D[n]['epochs']['all']['r2_B'] for n in names}; h = group_stat(af, x0, vals, names, C_BEH['whisk'] if False else INK, rng=rng, ms=3.8)
        ledger_add('Fig11', 'r2B_' + g, 'share of axon variance explained by whisking, face, pupil, ball (in-sample, lags -1..+3 s), %s sessions' % ('lick-task' if g == 'task' else 'passive'), h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'supported', note='1317 / 1529 (11-06, flagged) explain 4% / 2%: little behavioural range', unit='R2')
    af.set_xticks([0, 1.4]); af.set_xticklabels(['task', 'passive'], fontsize=5.2); af.set_xlim(-0.6, 2.0); af.set_ylim(0, 0.8); af.set_ylabel('variance explained by behaviours (R²)', fontsize=5.4); af.set_title('behaviour share', fontsize=5.6, loc='left', fontweight='bold', pad=11); plabel(af, 'g', dx=-24, dy=12)
    fig.text(0.105, 0.004, wrap('a, b: ▼ = square onset (task: filled hit, open miss), grey = 3 s after onset; axon mean and its prediction from pupil + face + whisking + ball, residual = difference. Signal: mean best dF/F of the signal-like ROIs, high-passed > 30 s.\n'
             'c-g: behaviours (pupil, face, whisking, ball where plausible) at lags -1..+3 s enter one ridge model together with the epoch regressors pre (-3..0 s), cue (0-3 s) and post (3-10 s) after square onset (reference: rest of the cycle);\n'
             'raw = epochs only, regress. = with behaviours. Null: epoch pattern circularly shifted by >= 60 s (500 draws); n/N = unflagged sessions with z > 2.33. c: residual = axon mean minus the behaviour part of the model. Task: 1752, 1744, 1423,\n'
             '1520 (RB24) + flagged 1317, RB23 1704; passive: 1456, 1651 (RB24), RB23 1423 + flagged 1529. f: hit-trial cycle minus miss-trial cycle and hit rate of the surrounding 10 trials (per SD), same model, null = hit sequence shifted by >= 5 trials.\n'
             'Black = mean of per-animal means, 95% hierarchical bootstrap CI (task: one animal); task and passive differ in animal and training day. In 1456, 1752 and RB23 1423 / 1704 a pupil / face transient follows the square at ~11 s (end of grey).'), fontsize=4.6, color=MUTED)
    savefig(fig, 'Fig11_task_engagement')


if __name__ == '__main__':
    main()
