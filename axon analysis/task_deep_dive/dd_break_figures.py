"""Tables + figures of the grey-break analysis (dd_break.py, dd_break_fir.py).   usage: python dd_break_figures.py      (needs the pickles in OUTDIR/cache)

Figures (PNG + PDF in plots/task_deep_dive):
  B1_break_state       what the break is: example, pupil / axon / lick rate aligned to the break start and end, lick rate inside the break, pupil step per break
  B2_break_axons       axon activity in the break: level change at onset / offset (vs random task alignments), amplitude and behaviour coupling inside the break vs task stretches
  B3_break_licks       break licks as a lick class: bout-triggered averages, FIR kernels (break-wide terms removed), behaviour-matched contrasts vs the other classes
Tables: break_pooled_contrasts.csv (+ .txt), break_fir.csv, break_lta.csv, (break_table.csv, break_aligned.csv, break_state.csv are written by dd_break.py)
Pooling convention of the task deep dive: unit = session, all unflagged sessions are ONE animal (RB24), pooled = mean of the session estimates, CI = bootstrap over sessions (n = 4, wide by construction),
Stouffer z of the sessions' own shift-null z; flagged sessions (1317, RB23 1704) are drawn hollow and never pooled.
"""
import os, sys, csv, pickle, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'paper_figures'))
import dd_common as dd
import dd_break as db
import fig_common as fc
import matplotlib.pyplot as plt
from xs_common import DT
import xs_common as xc

CACHE = os.path.join(dd.OUTDIR, 'cache')
SESS = db.SESS; UNFL = db.UNFLAGGED; FLAG = set(dd.FLAGGED)
SCOL = dict(zip(['AnimalRB24_261109_1506', 'AnimalRB24_261105_1423', 'AnimalRB24_261030_1744', 'AnimalRB24_261029_1752', 'AnimalRB24_261106_1317', 'AnimalRB23_261029_1704'],
                ['#E11D48', '#1F5FBF', '#D98A00', '#7E57C2', '#6B7280', '#0891B2']))
CLS = ['wait_early', 'ant_to', 'ant_sq', 'win', 'post_a', 'post_b', 'brk']
CLS_LAB = {'wait_early': 'wait\n(early)', 'ant_to': 'antic.\n-> time-out', 'ant_sq': 'antic.\n-> square', 'win': 'reward\nwindow', 'post_a': 'post\n3-6 s', 'post_b': 'post\n> 6 s', 'brk': 'grey\nbreak', 'tail': 'end\ntail'}
CLS_COL = {'wait_early': '#6B7280', 'ant_to': '#C2185B', 'ant_sq': '#E57AA3', 'win': '#1F5FBF', 'post_a': '#2E8B57', 'post_b': '#8C6D31', 'brk': '#000000', 'tail': '#7C3AED'}
rng = np.random.default_rng(5)


def save(fig, name):
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(dd.OUTDIR, name + '.' + ext), bbox_inches='tight', dpi=300)
    plt.close(fig); print('saved', name)
def label(ax, s): ax.text(-0.16, 1.1, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
def ld(prefix, nm):
    f = os.path.join(CACHE, '%s_%s.pkl' % (prefix, nm)); return pickle.load(open(f, 'rb')) if os.path.exists(f) else None
def sess_lab(nm): return fc.short(nm)[5:]
def boot_mean(x, n=5000): x = np.asarray(x, float); return np.percentile([rng.choice(x, len(x)).mean() for _ in range(n)], [2.5, 97.5])


# ============================================================================================================= tables
def fir_stat(r, o, a, b=None):
    B = r['out'][o]; est = B['sum'][a] - (B['sum'][b] if b else 0); bt = B['boot'][a] - (B['boot'][b] if b else 0); nu = B['null'][a] - (B['null'][b] if b else 0)
    z = (est - nu.mean()) / (nu.std() + 1e-12); return est, np.percentile(bt, [2.5, 97.5]), z


def tables():
    # ---- matched contrasts, pooled
    keys = ['brk_vs_wait_early', 'pre_brk_vs_wait_early', 'brkLate_vs_wait_early', 'pre_brkLate_vs_wait_early', 'brkEarly_vs_wait_early', 'brk_vs_ant_to', 'pre_brk_vs_ant_to', 'brk_vs_ant_sq', 'pre_brk_vs_ant_sq', 'brk_vs_rew', 'pre_brk_vs_rew', 'brk_vs_post_a', 'pre_brk_vs_post_a', 'brk_vs_post_b', 'pre_brk_vs_post_b',
            'brk_state_10s', 'tail_vs_wait_early', 'tail_vs_post_b', 'tail_state_10s']
    rows, txt = [], []
    for o in ('P', 'control'):
        for k in keys:
            est, zs_, det = [], [], []
            for n in SESS:
                m = ld('matchedbrk', n); c = m['contrasts'].get(o, {}).get(k) if m else None
                if not c or not c.get('obs') or not np.isfinite(c['obs'].get('matched', np.nan)): continue
                ob = c['obs']
                if n in FLAG or ob['n_matched'] < 8: continue
                est.append(ob['matched']); zs_.append(c['stats']['matched']['z']); det.append('%s %+.3f %d/%d' % (sess_lab(n)[:11], ob['matched'], ob['n_matched'], ob['n_focal']))
            if len(est) < 2: rows.append([o, k, len(est), '', '', '', '', '', '; '.join(det)]); continue
            est = np.array(est); zz = np.array(zs_); ok = np.isfinite(zz); ci = boot_mean(est); st = zz[ok].sum() / np.sqrt(ok.sum())
            rows.append([o, k, len(est), round(est.mean(), 4), round(ci[0], 4), round(ci[1], 4), int((est > 0).sum()), round(st, 2), '; '.join(det)])
    with open(os.path.join(dd.OUTDIR, 'break_pooled_contrasts.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['outcome', 'contrast', 'n_sessions', 'pooled_matched_dFF', 'ci_lo', 'ci_hi', 'n_positive', 'stouffer_z', 'sessions (estimate matched/focal)']); w.writerows(rows)
    with open(os.path.join(dd.OUTDIR, 'break_pooled_contrasts.txt'), 'w') as fh:
        for r in rows: fh.write('%-8s %-24s n %s pooled %s [%s, %s] pos %s Stouffer %s | %s\n' % tuple(r))
    for r in rows:
        if r[0] == 'P': print('%-24s n %s pooled %6s [%6s, %6s] pos %s z %5s | %s' % (r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]))
    # ---- FIR (plain and with break-wide terms)
    rows = []
    for var, prefix in (('plain', 'firbrk'), ('boxcar', 'firbrk2')):
        for n in SESS:
            r = ld(prefix, n)
            if not r: continue
            for o in ('P', 'control'):
                if o not in r['out']: continue
                for c in ('wait_early', 'ant_to', 'ant_sq', 'win', 'post_a', 'post_b', 'brk', 'tail'):
                    if r['n_licks'].get(c, 0) < 20: continue
                    e, ci, z = fir_stat(r, o, 'lick_' + c); rows.append([var, n, o, 'lick_' + c, r['n_licks'][c], round(e, 5), round(ci[0], 5), round(ci[1], 5), round(z, 2)])
                for c in ('wait_early', 'ant_to', 'win', 'post_b'):
                    if r['n_licks'].get(c, 0) < 20 or r['n_licks'].get('brk', 0) < 20: continue
                    e, ci, z = fir_stat(r, o, 'lick_brk', 'lick_' + c); rows.append([var, n, o, 'lick_brk - lick_' + c, r['n_licks']['brk'], round(e, 5), round(ci[0], 5), round(ci[1], 5), round(z, 2)])
                for nm in ('brk_on', 'brk_off'):
                    if nm in r['out'][o]['sum']:
                        e, ci, z = fir_stat(r, o, nm); rows.append([var, n, o, nm, 0, round(e, 5), round(ci[0], 5), round(ci[1], 5), round(z, 2)])
    with open(os.path.join(dd.OUTDIR, 'break_fir.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['variant', 'session', 'outcome', 'term', 'n_licks', 'mean_kernel_0_2s_dFF_per_lick', 'ci_lo', 'ci_hi', 'shift_z']); w.writerows(rows)
    return rows


def lta_table():
    res = pickle.load(open(os.path.join(CACHE, 'break_lta.pkl'), 'rb')); tau = db.TAU; rows = []
    pre = (tau >= -3) & (tau <= -1); post = (tau >= 0) & (tau <= 2)
    for s in SESS:
        for c, d in res[s].items():
            for k in ('P', 'ctl', 'jaw'):
                if k not in d: continue
                A = d[k]; per = np.nanmean(A[:, post], 1) - np.nanmean(A[:, pre], 1); m = np.nanmean(per)
                rows.append([s, c, k, d['n'], round(float(m), 4), round(float(np.nanstd(per) / np.sqrt(np.isfinite(per).sum())), 4)])
    with open(os.path.join(dd.OUTDIR, 'break_lta.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['session', 'class', 'series', 'n_bouts', 'mean_0_2s_minus_pre_-3_-1s', 'sem']); w.writerows(rows)
    return res, rows


# ============================================================================================================= B1
def fig_b1():
    s = 'AnimalRB24_261109_1506'; S, T = db.load(s); x = T.brk[1]
    a, b = x['t0'] - 45, x['t1'] + 40; m = (S.t >= a) & (S.t <= b); tt = S.t[m] - x['t0']
    sd = db.sess_series(S, T); lk = T.lick[(T.lick >= a) & (T.lick <= b)] - x['t0']
    tau = np.arange(-60, 110, DT); ser_all = {}
    for s_ in UNFL:
        S_, T_ = db.load(s_); sd_ = db.sess_series(S_, T_)
        for kind in ('onset', 'offset'):
            for nm in ('pupil', 'P', 'ctl'):
                if nm not in sd_: continue
                for x_ in [y for y in T_.brk if y['kind'] == 'brk']:
                    t0 = x_['t0'] if kind == 'onset' else x_['t1']; i = np.round((t0 + tau - S_.edges[0]) / DT - 0.5).astype(int); o = np.full(len(tau), np.nan); ok = (i >= 0) & (i < S_.nb); o[ok] = sd_[nm][i[ok]]
                    o = o - np.nanmean(o[(tau >= -30) & (tau <= -5)]); ser_all.setdefault((kind, nm), []).append(o)
    # ---- rebuild with a simpler layout (clean grid: a | b-onset | b-offset on the top, c | d below)
    fig = plt.figure(figsize=(fc.W2, 4.3)); gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.42, left=0.075, right=0.985, top=0.94, bottom=0.1, width_ratios=[1.5, 1, 1])
    ex = [fig.add_subplot(gs[i, 0]) for i in range(3)]
    ex[0].plot(tt, sd['P'][m], color=fc.C_SIG, lw=0.7); ex[0].set_ylabel('signal-like\n(dF/F)')
    ex[1].plot(tt, S.beh['pupil'][m], color=fc.C_BEH['pupil'], lw=0.9); ex[1].set_ylabel('pupil area')
    ex[2].plot(tt, sd['jaw'][m], color=fc.INK, lw=0.5); ex[2].set_ylabel('jaw motion (z)'); ex[2].plot(lk, np.full(len(lk), -3.3), '|', color='#8C6D31', ms=4, mew=0.7); ex[2].set_ylim(-3.6, 5)
    for ax in ex:
        ax.axvspan(0, x['t1'] - x['t0'], color='0.92', zorder=0, lw=0); ax.margins(x=0)
    ex[0].set_xticklabels([]); ex[1].set_xticklabels([]); ex[2].set_xlabel('time from the break start (s); grey = break, ticks = licks')
    ex[0].set_title('Example: 11-09 15:06, break 2 (one with an onset rise; the average, b, has none)', loc='left', fontsize=6.3); label(ex[0], 'a')
    for ci, kind in enumerate(("onset", "offset")):
        for ri, (nm, col, lab) in enumerate((('pupil', fc.C_BEH['pupil'], 'pupil (z)'), ('P', fc.C_SIG, 'signal-like (dF/F)'), ('ctl', fc.C_NOISE, 'noise-like (dF/F)'))):
            ax = fig.add_subplot(gs[ri, 1 + ci]); A = np.array(ser_all[(kind, nm)]); mu = np.nanmean(A, 0); se = np.nanstd(A, 0) / np.sqrt(np.isfinite(A).sum(0)); sm = xc.nsmooth(mu, 1.0 / DT) if nm != 'pupil' else mu
            ax.plot(tau, sm, color=col, lw=0.9); ax.fill_between(tau, sm - se, sm + se, color=col, alpha=0.2, lw=0); ax.axvline(0, color='k', lw=0.4); ax.axhline(0, color='0.7', lw=0.4); ax.set_xlim(-45, 60 if kind == 'onset' else 45)
            if nm != 'pupil': ax.set_ylim(-0.2, 0.2)
            if ri == 0: ax.set_title('aligned to the break %s' % kind, loc='left', fontsize=6.3)
            if ci == 0: ax.set_ylabel(lab, fontsize=5.8)
            if ri < 2: ax.set_xticklabels([])
            else: ax.set_xlabel('s from the break %s' % kind)
            if ri == 0 and ci == 0: label(ax, 'b')
    fig.text(0.52, 0.015, 'b: mean +- SEM of the unflagged breaks (n = 10, 4 sessions), each referenced to its own -30..-5 s; slow-baseline dF/F (Fc / running median +-300 s)', fontsize=5.4, color=fc.MUTED)
    save(fig, 'B1_break_state')

    # ---- second figure part: lick rate in the break and pupil step per break  (B1 lower panels as their own small figure, appended to B1 name)
    tab = list(csv.DictReader(open(os.path.join(dd.OUTDIR, 'break_table.csv'))))
    fig, axs = plt.subplots(1, 3, figsize=(fc.W2, 2.2), gridspec_kw=dict(wspace=0.45, left=0.07, right=0.99, top=0.85, bottom=0.2))
    ax = axs[0]; label(ax, 'c'); edges_t = [0, 20, 40, 60]
    for s_ in SESS:
        rr = [r for r in tab if r['session'] == s_ and r['kind'] == 'brk']
        if not rr: continue
        v = np.array([[float(r[k]) if r[k] not in ('', 'nan') else np.nan for k in ('rate_0_20s', 'rate_20_40s', 'rate_40_60s', 'rate_60_end')] for r in rr]); mu = np.nanmean(v, 0)
        ax.plot([10, 30, 50, 70], mu, '-o', ms=2.5, lw=0.8, color=SCOL[s_], mfc=SCOL[s_] if s_ not in FLAG else 'white', label=sess_lab(s_)); ax.axhline(float(rr[0]['wait_rate_hz']), color=SCOL[s_], lw=0.4, ls=':')
    ax.set_xlabel('s since the break start'); ax.set_ylabel('licks / s'); ax.set_title('Licking in the break (dotted = same session, task waits)', loc='left', fontsize=6.0); ax.legend(frameon=False, fontsize=4.6, ncol=2, loc='upper right')
    ax = axs[1]; label(ax, 'd')
    for i, k in enumerate(('pupil_ratio', 'face_ratio', 'whisk_ratio', 'jaw_ratio')):
        for r in tab:
            if r['kind'] != 'brk' or r[k] in ('', 'nan'): continue
            ax.plot(i + rng.uniform(-0.12, 0.12), float(r[k]), 'o', ms=2.6, color=SCOL[r['session']], mfc=SCOL[r['session']] if r['session'] not in FLAG else 'white', mew=0.6)
    ax.axhline(1, color='k', lw=0.4); ax.set_xticks(range(4)); ax.set_xticklabels(['pupil\narea', 'face\nmotion', 'whisking', 'jaw\nmotion']); ax.set_ylabel('break level / task wait level'); ax.set_title('Pupil changes, movement does not (per break)', loc='left', fontsize=6.0)
    ax = axs[2]; label(ax, 'e'); bw = []
    for r in tab:
        if r['kind'] != 'brk': continue
        ax.plot([0, 1, 2], [float(r['pupil_pre']), float(r['pupil_in']), float(r['pupil_post']) if r['pupil_post'] not in ('', 'nan') else np.nan], '-o', ms=2.2, lw=0.6, color=SCOL[r['session']], alpha=0.8, mfc='white' if r['session'] in FLAG else SCOL[r['session']])
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(['before\n(-12..-2 s)', 'in the break\n(10 s .. end-5 s)', 'after\n(+8..+25 s)']); ax.set_ylabel('pupil area (a.u.)'); ax.set_title('Pupil falls back when the next subsession starts', loc='left', fontsize=6.0)
    save(fig, 'B1b_break_licking_pupil')


# ============================================================================================================= B2
def fig_b2():
    al = pickle.load(open(os.path.join(CACHE, 'break_aligned.pkl'), 'rb')); st = pickle.load(open(os.path.join(CACHE, 'break_state.pkl'), 'rb'))
    rows = list(csv.DictReader(open(os.path.join(dd.OUTDIR, 'break_aligned.csv')))); srows = list(csv.DictReader(open(os.path.join(dd.OUTDIR, 'break_state.csv'))))
    fig = plt.figure(figsize=(fc.W2, 3.0)); gs = fig.add_gridspec(1, 3, wspace=0.5, left=0.07, right=0.99, top=0.84, bottom=0.2, width_ratios=[1.2, 1, 1])
    ax = fig.add_subplot(gs[0]); label(ax, 'a'); xs = []
    spec = [('onset', 'early', '0-10 s'), ('onset', 'mid', '10-40 s'), ('onset', 'late', '40 s-end'), ('offset', 'post', '0-20 s')]
    for ni, (nm, col) in enumerate((('P', fc.C_SIG), ('ctl', fc.C_NOISE))):
        for i, (kind, w, lab) in enumerate(spec):
            r = [q for q in rows if q['alignment'] == kind and q['series'] == nm and q['window'] == w and q['group'] == 'pooled']
            if not r: continue
            r = r[0]; xp = i + (ni - 0.5) * 0.3; est = float(r['change_vs_pre']); ax.plot([xp, xp], [float(r['ci_lo']), float(r['ci_hi'])], color=col, lw=1); ax.plot(xp, est, 'o', color=col, ms=3.5)
            if ni == 0: ax.fill_between([i - 0.4, i + 0.4], float(r['null_mean']) - 1.96 * float(r['null_sd']), float(r['null_mean']) + 1.96 * float(r['null_sd']), color='0.85', lw=0, zorder=0)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(4)); ax.set_xticklabels(['onset\n0-10 s', 'onset\n10-40 s', 'onset\n40 s-end', 'offset\n0-20 s'], fontsize=5.6); ax.set_ylabel('change vs the 30 s before (dF/F)')
    ax.set_title('Level change at the break (blue signal-like, amber noise-like;\ngrey band = 95% of random task alignments)', loc='left', fontsize=5.8)
    # b amplitude ratio / coupling
    ax = fig.add_subplot(gs[1]); label(ax, 'b'); names = [('sd_P', 'population'), ('sd_ctl', 'noise-like'), ('sd_face', 'face'), ('sd_whisk', 'whisking'), ('sd_jaw', 'jaw'), ('sd_lick', 'lick rate')]
    for i, (k, lab) in enumerate(names):
        r = [q for q in srows if q['statistic'] == k and q['group'] == 'pooled']
        if not r: continue
        r = r[0]; col = fc.C_SIG if k == 'sd_P' else (fc.C_NOISE if k == 'sd_ctl' else fc.INK); ax.plot([i, i], [float(r['ci_lo']), float(r['ci_hi'])], color=col, lw=1); ax.plot(i, float(r['in_break']), 'o', color=col, ms=3.5)
    ax.axhline(1, color='k', lw=0.4); ax.set_xticks(range(len(names))); ax.set_xticklabels([n[1] for n in names], rotation=40, ha='right', fontsize=5.4); ax.set_ylabel('SD in the break / SD in the task'); ax.set_title('Fluctuation amplitude (30 s high-passed)', loc='left', fontsize=6.0)
    ax = fig.add_subplot(gs[2]); label(ax, 'c'); names = [('r_face', 'face'), ('r_whisk', 'whisking'), ('r_jaw', 'jaw'), ('r_dpupil', 'pupil dilation\nrate'), ('r_lick', 'lick rate')]
    for i, (k, lab) in enumerate(names):
        r = [q for q in srows if q['statistic'] == k and q['group'] == 'pooled']
        if not r: continue
        r = r[0]; tm = float(r['task_mean']); sdn = float(r['task_sd_of_pooled_mean']); ax.fill_between([i - 0.35, i + 0.35], tm - 1.96 * sdn, tm + 1.96 * sdn, color='0.85', lw=0, zorder=0); ax.plot([i, i], [float(r['ci_lo']), float(r['ci_hi'])], color=fc.C_SIG, lw=1); ax.plot(i, float(r['in_break']), 'o', color=fc.C_SIG, ms=3.5)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(len(names))); ax.set_xticklabels([n[1] for n in names], rotation=40, ha='right', fontsize=5.4); ax.set_ylabel('r(population, behaviour)'); ax.set_title('Coupling inside the break (grey = task windows)', loc='left', fontsize=6.0)
    save(fig, 'B2_break_axons')


# ============================================================================================================= B3
def fig_b3():
    res, ltarows = lta_table(); tau = db.TAU
    fig = plt.figure(figsize=(fc.W2, 5.4)); gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.5, left=0.07, right=0.99, top=0.93, bottom=0.12, height_ratios=[1, 1.05])
    # a: LTA
    ax = fig.add_subplot(gs[0, 0]); label(ax, 'a'); ax2 = fig.add_subplot(gs[1, 0]); label(ax2, 'b')
    for c in ('wait_early', 'ant_to', 'win', 'post_b', 'brk'):
        cl = 'rew' if c == 'win' else c
        for k, a_ in (('P', ax), ('jaw', ax2)):
            M = []
            for s in UNFL:
                if cl in res[s] and k in res[s][cl]:
                    A = res[s][cl][k]; M.append(np.nanmean(A, 0) - np.nanmean(np.nanmean(A[:, (tau >= -3) & (tau <= -1)], 1)))
            if not M: continue
            M = np.array(M); mu = np.nanmean(M, 0); a_.plot(tau, xc.nsmooth(mu, 1.0 / DT) if k == 'P' else mu, color=CLS_COL[c], lw=1.5 if c == 'brk' else 0.9, label={'win': 'reward bout', 'brk': 'break', 'ant_to': 'antic. -> time-out', 'wait_early': 'wait (early)', 'post_b': 'post > 6 s'}[c])
    for a_, t in ((ax, 'Population around bout onsets'), (ax2, 'Jaw motion around bout onsets (z)')):
        a_.axvline(0, color='k', lw=0.4); a_.axhline(0, color='0.7', lw=0.4); a_.set_xlabel('s from the bout onset'); a_.set_title(t, loc='left', fontsize=6.0)
    ax.set_ylabel('dF/F vs -3..-1 s'); ax2.set_ylabel('jaw (z) vs -3..-1 s'); ax.legend(frameon=False, fontsize=4.8, ncol=2, loc='upper right')
    # c: FIR kernels (per lick), boxcar variant, by session
    ax = fig.add_subplot(gs[0, 1:]); label(ax, 'c'); cl_l = ['wait_early', 'ant_to', 'ant_sq', 'win', 'post_a', 'post_b', 'brk']
    for var, prefix, dx, mk in (('break-wide terms removed', 'firbrk2', 0.12, 'o'), ('plain', 'firbrk', -0.12, 's')):
        for ci, c in enumerate(cl_l):
            est = []
            for si, s in enumerate(SESS):
                r = ld(prefix, s)
                if not r or r['n_licks'].get(c, 0) < 20: continue
                e, cil, z = fir_stat(r, 'P', 'lick_' + c); est.append((s, e))
                xp = ci + dx + (si - 2.5) * 0.035
                ax.plot(xp, e, mk, ms=2.6, color=SCOL[s], mfc=(SCOL[s] if (s not in FLAG and prefix == 'firbrk2') else 'white'), mew=0.6, alpha=0.9)
            ee = [e for s, e in est if s not in FLAG]
            if ee: ax.plot([ci + dx - 0.12, ci + dx + 0.12], [np.mean(ee)] * 2, color='k', lw=1.1 if prefix == 'firbrk2' else 0.5)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(len(cl_l))); ax.set_xticklabels([CLS_LAB[c] for c in cl_l], fontsize=5.4); ax.set_ylabel('mean kernel 0-2 s (dF/F per lick)')
    ax.set_title('Encoding model, response per lick by class: filled circles = break-wide terms removed (right of each tick), squares = plain model (left);\nsession colours as in B1, hollow = flagged, bar = mean of the unflagged sessions', loc='left', fontsize=5.4)
    # d: matched contrasts
    ax = fig.add_subplot(gs[1, 1:]); label(ax, 'd'); keys = [('brk_vs_wait_early', 'break vs\nwait'), ('pre_brk_vs_wait_early', 'before: break\nvs wait'), ('brk_vs_ant_to', 'break vs\nantic. -> TO'), ('brk_vs_ant_sq', 'break vs\nantic. -> sq'), ('brk_vs_rew', 'break vs\nreward'), ('brk_vs_post_b', 'break vs\npost > 6 s'),
                                                                             ('brk_state_10s', 'break 10 s\nwindows vs task'), ('tail_vs_wait_early', 'tail vs\nwait')]
    for ki, (k, lab) in enumerate(keys):
        ee = []
        for si, s in enumerate(SESS):
            m = ld('matchedbrk', s); c = m['contrasts']['P'].get(k) if m else None
            if not c or not c.get('obs') or not np.isfinite(c['obs'].get('matched', np.nan)): continue
            o = c['obs']; cov = o['n_matched'] / max(o['n_focal'], 1)
            if o['n_matched'] < 8 and s not in FLAG and False: continue
            xp = ki + (si - 2.5) * 0.11; z = c['stats']['matched']['z']; sg = abs(z) > 2.5
            ax.plot(xp, o['matched'], 'o', ms=2.0 + 4.0 * cov, color=SCOL[s], mfc=SCOL[s] if (s not in FLAG) else 'white', mew=0.7 if sg is False else 1.3, alpha=0.9)
            if s not in FLAG and o['n_matched'] >= 8: ee.append(o['matched'])
        if len(ee) >= 2: ax.plot([ki - 0.4, ki + 0.4], [np.mean(ee)] * 2, color='k', lw=1.2)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(len(keys))); ax.set_xticklabels([k[1] for k in keys], fontsize=5.0); ax.set_ylabel('matched difference (dF/F)')
    ax.set_title('Behaviour-matched bout contrasts (0-3 s after the bout onset; marker size = share of bouts with a partner;\nbar = mean of the unflagged sessions with >= 8 matched bouts)', loc='left', fontsize=5.4)
    save(fig, 'B3_break_licks')


if __name__ == '__main__':
    tables()
    fig_b1(); fig_b2(); fig_b3()
