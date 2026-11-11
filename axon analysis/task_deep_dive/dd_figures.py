"""Figures + tables of the lick-task deep dive (task_deep_dive).  usage: python dd_figures.py   (needs the pickles of dd_fir / dd_matched / dd_hitmiss / dd_extras in OUTDIR/cache)

Figures (PNG + PDF in plots/task_deep_dive):
  D1_task_and_lick_classes   what the task does in RB24_261109_1506: wait segments / time-outs, rig-count check of the segments, lick rate vs time to the wait end, lick rate around the square
  D2_lick_class_kernels      FIR kernels per lick class (population + ROI groups + noise-like control), class differences, ROI sign structure, stability over the subsessions
  D3_matched_contrasts       behaviour-matched class contrasts of 1506 and the replication in the other lick sessions
  D4_hit_vs_miss             hit vs miss: time course, windows raw vs matched, pre-trial prediction, ROI scan against the noise-like ROIs, other sessions
Conventions as in paper_figures: blue = signal-like axons, amber = noise-like ROIs (control), session = unit, flagged sessions hollow and never pooled.
"""
import os, sys, pickle, csv, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'paper_figures'))
import dd_common as dd
import fig_common as fc
import matplotlib.pyplot as plt
from xs_common import DT

PRIMARY = 'AnimalRB24_261109_1506'
CACHE = os.path.join(dd.OUTDIR, 'cache')
CLS_COL = {'wait_early': '#6B7280', 'ant_to': '#C2185B', 'ant_sq': '#E57AA3', 'win': '#1F5FBF', 'post_a': '#2E8B57', 'post_b': '#8C6D31'}
CLS_LAB = {'wait_early': 'wait licks (early)', 'ant_to': 'anticipatory -> time-out', 'ant_sq': 'anticipatory -> square', 'win': 'reward window (0.1-3 s)', 'post_a': 'post 3-6 s', 'post_b': 'post > 6 s'}
OUT_COL = {'P': fc.C_SIG, 'g1': '#4C78A8', 'g2': '#2E8B57', 'g3': '#9467BD', 'g4': '#8C564B', 'control': fc.C_NOISE}
OUT_LAB = {'P': 'all signal-like', 'g1': 'group 1', 'g2': 'group 2', 'g3': 'group 3', 'g4': 'group 4', 'control': 'noise-like (control)'}
SESS_ORDER = [PRIMARY, 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261030_1744', 'AnimalRB24_261029_1752', 'AnimalRB24_261106_1317', 'AnimalRB23_261029_1704']
FLAGGED = set(dd.FLAGGED)
SESS_COL = dict(zip(SESS_ORDER, ['#E11D48', '#1F5FBF', '#2E8B57', '#D98A00', '#7E57C2', '#6B7280', '#0891B2']))


def sh(nm): return fc.short(nm)
def load(prefix, nm):
    f = os.path.join(CACHE, '%s_%s.pkl' % (prefix, nm)); return pickle.load(open(f, 'rb')) if os.path.exists(f) else None
def save(fig, name):
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(dd.OUTDIR, name + '.' + ext), bbox_inches='tight', dpi=300)
    plt.close(fig); print('saved', name)
def label(ax, s): ax.text(-0.18, 1.08, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')


# ================================================================================================= D1
def fig_d1():
    T = dd.Task(PRIMARY); s = T.seg
    fig = plt.figure(figsize=(fc.W2, 4.8)); gs = fig.add_gridspec(2, 3, hspace=0.85, wspace=0.55, left=0.07, right=0.985, top=0.90, bottom=0.10)
    # a: one example trial with its licks, coloured by class
    ax = fig.add_subplot(gs[0, :2]); label(ax, 'a')
    t_ex = int(np.argmax((T.n_to >= 5) & T.hit)); t0 = T.tr0[t_ex]
    segs = np.where(s['trial'] == t_ex)[0]
    for j in segs:
        ax.axvspan(s['t0'][j] - t0, s['t1'][j] - t0, ymin=0.0, ymax=0.14, color='#E5E7EB' if j % 2 == 0 else '#CBD5E1', lw=0)
        ax.text((s['t0'][j] + s['t1'][j]) / 2 - t0, 0.07, ('wait %.0f s\n%s' % (s['X'][j], 'time-out' if s['to'][j] else 'square')), ha='center', va='center', fontsize=5, transform=ax.get_xaxis_transform())
    ax.axvspan(T.on[t_ex] - t0, T.st_end[t_ex] - t0, ymin=0, ymax=0.14, color='#06B6D4', lw=0, alpha=0.8)
    ax.axvspan(T.st_end[t_ex] - t0, T.tr1[t_ex] - t0, ymin=0, ymax=0.14, color='#F3F4F6', lw=0)
    ax.text((T.on[t_ex] + T.st_end[t_ex]) / 2 - t0, 0.07, 'square', ha='center', va='center', fontsize=5, transform=ax.get_xaxis_transform(), color='white')
    ax.text((T.st_end[t_ex] + T.tr1[t_ex]) / 2 - t0, 0.07, 'grey 8.5 s', ha='center', va='center', fontsize=5, transform=ax.get_xaxis_transform())
    m = (T.lick_tr == t_ex + 1); lt = T.lick[m]; cl, _ = T.classify(lt); cl = np.array(['win' if c in ('win_other', 'win_early') else c for c in cl], dtype=object)
    for c in CLS_COL:
        k = cl == c
        if k.any(): ax.vlines(lt[k] - t0, 0.2, 0.9, color=CLS_COL[c], lw=0.8, transform=ax.get_xaxis_transform(), label=CLS_LAB[c])
    if np.isfinite(T.rew_t[t_ex]): ax.plot([T.rew_t[t_ex] - t0], [0.97], marker='v', color='k', ms=4, transform=ax.get_xaxis_transform(), clip_on=False); ax.text(T.rew_t[t_ex] - t0 + 0.6, 0.97, 'pump', fontsize=5.5, va='center', transform=ax.get_xaxis_transform())
    ax.set_yticks([]); ax.set_xlabel('time in trial %d (s)' % (t_ex + 1)); ax.spines['left'].set_visible(False)
    ax.set_title('One trial of the session: %d silent time-outs, then the square (licks coloured by class)' % int(T.n_to[t_ex]), loc='left', fontsize=6.5)
    ax.legend(ncol=6, frameon=False, fontsize=5, loc='upper center', bbox_to_anchor=(0.5, -0.22), columnspacing=1.0, handlelength=1.2)
    # b: rig count check
    ax = fig.add_subplot(gs[0, 2]); label(ax, 'b')
    ok = np.isfinite(s['rig']); X_ = s['X']; y = s['rig']
    ax.scatter(X_[ok & s['to']] + np.random.default_rng(1).uniform(-0.15, 0.15, (ok & s['to']).sum()), y[ok & s['to']], s=5, color=CLS_COL['ant_to'], alpha=0.6, lw=0, label='ended in time-out (n=%d)' % (ok & s['to']).sum())
    ax.scatter(X_[ok & s['last']] + np.random.default_rng(2).uniform(-0.15, 0.15, (ok & s['last']).sum()), y[ok & s['last']], s=5, color=CLS_COL['win'], alpha=0.6, lw=0, label='ended in square (n=%d)' % (ok & s['last']).sum())
    xx = np.array([7.5, 13.5]); ax.plot(xx, T.maxerr * xx, 'k--', lw=0.7); ax.text(13.4, T.maxerr * 13.4 + 0.5, 'limit 0.5 x wait', fontsize=5.5, ha='right')
    ax.set_xlabel('wait drawn (s)'); ax.set_ylabel('licks the rig counted in the wait'); ax.legend(frameon=False, fontsize=5.5, loc='upper left')
    cs = T.check_segments(); ax.set_title('Time-outs follow the lick count:\n> 0.5 x wait in %.0f%% of time-out waits, never in the\nwait before a square; my count = rig count +-1 in %.1f%%' % (100 * cs['rule_to'], 100 * cs['count_within1']), loc='left', fontsize=5.8)
    # c: lick rate vs time to the end of the wait
    ax = fig.add_subplot(gs[1, 0]); label(ax, 'c')
    L = T.locate(T.lick); bins = np.arange(0, 12.01, 1.0)
    for kind, col, lab in ((True, CLS_COL['ant_to'], 'wait ends in time-out'), (False, CLS_COL['win'], 'wait ends in square')):
        cnt = np.zeros(len(bins) - 1); occ = np.zeros(len(bins) - 1)
        for j in np.where(s['to'] == kind)[0]:
            dur = s['t1'][j] - s['t0'][j]; occ += np.clip(dur - bins[:-1], 0, 1.0); m = L['seg'] == j; cnt += np.histogram(s['t1'][j] - T.lick[m], bins=bins)[0]
        ok_ = occ >= 30; ax.plot(bins[:-1][ok_] + 0.5, (cnt / np.maximum(occ, 1e-9))[ok_], color=col, marker='o', ms=2.5, label=lab)
    ax.invert_xaxis(); ax.set_xlabel('time before the end of the wait (s)'); ax.set_ylabel('lick rate (Hz)'); ax.set_ylim(0, 1.8); ax.axhline(T.maxerr, color='k', lw=0.5, ls=':')
    ax.legend(frameon=False, fontsize=5.5, loc='upper left'); ax.set_title('No ramp towards the wait end', loc='left', fontsize=6.5)
    # d: lick rate aligned to onset: hit/miss, level
    ax = fig.add_subplot(gs[1, 1]); label(ax, 'd')
    bins = np.arange(-4, 9.01, 0.5)
    def rate(sel):
        c = np.zeros(len(bins) - 1)
        for t in np.where(sel)[0]:
            c += np.histogram(T.lick[T.lick_tr == t + 1] - T.on[t], bins)[0]
        return c / max(sel.sum(), 1) / 0.5
    for sel, col, lab in ((T.hit, CLS_COL['win'], 'hit (%d)' % T.hit.sum()), (~T.hit & (T.lvl == 2), '#D98A00', 'miss, level 2 (pump +1.5 s) (%d)' % (~T.hit & (T.lvl == 2)).sum()), (~T.hit & (T.lvl == 3), '#6B7280', 'miss, level 3: no pump (%d)' % (~T.hit & (T.lvl == 3)).sum())):
        ax.plot(bins[:-1] + 0.25, rate(sel), color=col, label=lab)
    ax.axvspan(0, 3, color='#06B6D4', alpha=0.12, lw=0); ax.axvline(0, color='k', lw=0.5); ax.axvline(T.rdel, color='#D98A00', lw=0.5, ls='--')
    ax.set_xlabel('time from square onset (s)'); ax.set_ylabel('lick rate (Hz)'); ax.set_ylim(0, 3.0); ax.legend(frameon=False, fontsize=5, loc='upper right')
    ax.set_title('Licking around the square', loc='left', fontsize=6.5)
    # e: bouts per class + hit rate
    ax = fig.add_subplot(gs[1, 2]); label(ax, 'e')
    cls = ['wait_early', 'ant_to', 'ant_sq', 'win', 'post_a', 'post_b']; cnt = {c: 0 for c in cls}
    for c, v in T.bout_class_counts.items():
        k = 'win' if c in ('rew', 'win_other') else c
        if k in cnt: cnt[k] += v
    ax.barh(range(len(cls)), [cnt[c] for c in cls], color=[CLS_COL[c] for c in cls]); ax.set_yticks(range(len(cls))); ax.set_yticklabels(['early wait', 'antic. -> time-out', 'antic. -> square', 'reward window', 'post 3-6 s', 'post > 6 s'], fontsize=5.5); ax.invert_yaxis()
    for i, c in enumerate(cls): ax.text(cnt[c] + 4, i, str(cnt[c]), va='center', fontsize=5.5)
    ax.set_xlabel('lick bouts (>= 3 licks)'); ax.set_title('Lick bouts per class\n%d trials, %d hits (level 2: %.0f%%, level 3: %.0f%%)' % (T.n, T.hit.sum(), 100 * T.hit[T.lvl == 2].mean(), 100 * T.hit[T.lvl == 3].mean()), loc='left', fontsize=6)
    save(fig, 'D1_task_and_lick_classes')


# ================================================================================================= D2
def fir_stat(r, o, nm):
    b = r['out'][o]; est = b['sum'][nm]; bt = b['boot'][nm]; nu = b['null'][nm]
    return est, np.percentile(bt, [2.5, 97.5]), (est - nu.mean()) / (nu.std() + 1e-12), (1 + (np.abs(nu - nu.mean()) >= abs(est - nu.mean())).sum()) / (1 + len(nu))


def fig_d2():
    r = load('fir', PRIMARY); ex = load('extras', PRIMARY)
    fig = plt.figure(figsize=(fc.W2, 6.6)); gs = fig.add_gridspec(3, 3, hspace=0.95, wspace=0.55, left=0.07, right=0.985, top=0.94, bottom=0.07)
    lags = r['lags']
    # a: kernels of the population per class
    ax = fig.add_subplot(gs[0, 0]); label(ax, 'a')
    for c in ('wait_early', 'ant_to', 'win', 'post_a', 'post_b'):
        k = r['out']['P']['kern']['lick_' + c]; bk = r['out']['P']['boot_kern']['lick_' + c]
        ax.plot(lags, k, color=CLS_COL[c], label=CLS_LAB[c]); ax.fill_between(lags, np.percentile(bk, 2.5, 0), np.percentile(bk, 97.5, 0), color=CLS_COL[c], alpha=0.15, lw=0)
    ax.axhline(0, color='k', lw=0.4); ax.axvline(0, color='k', lw=0.4); ax.set_xlabel('delay after a lick (s)'); ax.set_ylabel('dF/F per lick'); ax.legend(frameon=False, fontsize=4.8, loc='upper center', bbox_to_anchor=(0.5, -0.28), ncol=2, columnspacing=0.8, handlelength=1.2)
    ax.set_title('Kernel per lick class\n(52 signal-like ROIs)', loc='left', fontsize=6.0)
    # b: kernel summary per class x outcome (forest)
    ax = fig.add_subplot(gs[0, 1:]); label(ax, 'b')
    cls = ['wait_early', 'ant_to', 'ant_sq', 'win', 'post_a', 'post_b']; outs = ['P', 'g1', 'g2', 'g3', 'g4', 'control']; w = 0.12
    for oi, o in enumerate(outs):
        for ci, c in enumerate(cls):
            est, ci_, z, p = fir_stat(r, o, 'lick_' + c); xpos = ci + (oi - 2.5) * w
            ax.plot([xpos, xpos], ci_, color=OUT_COL[o], lw=0.8); ax.plot(xpos, est, 'o', color=OUT_COL[o], ms=3, mfc=OUT_COL[o] if p < 0.01 else 'white', mew=0.8, label=OUT_LAB[o] if ci == 0 else None)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(len(cls))); ax.set_xticklabels([CLS_LAB[c].replace(' (', '\n(').replace(' -> ', '\n-> ') for c in cls], fontsize=5.3); ax.set_ylabel('mean kernel 0-2 s (dF/F per lick)')
    ax.legend(frameon=False, fontsize=5, ncol=3, loc='lower left'); ax.set_title('Response per lick, by lick class and ROI group (95% block-bootstrap CI; filled = shift-null p < 0.01)', loc='left', fontsize=6.3)
    # c: class differences
    ax = fig.add_subplot(gs[1, 0:2]); label(ax, 'c')
    pairs = [('ant_to', 'wait_early', 'anticip. (time-out) - early wait'), ('win', 'wait_early', 'reward window - early wait'), ('win', 'ant_to', 'reward window - anticip. (time-out)'), ('post_a', 'win', 'post 3-6 s - reward window'),
             ('post_b', 'wait_early', 'post > 6 s - early wait'), ('post_b', 'win', 'post > 6 s - reward window')]
    for oi, o in enumerate(outs):
        b = r['out'][o]
        for pi, (a, c_, lab) in enumerate(pairs):
            est = b['sum']['lick_' + a] - b['sum']['lick_' + c_]; bt = b['boot']['lick_' + a] - b['boot']['lick_' + c_]; nu = b['null']['lick_' + a] - b['null']['lick_' + c_]
            p = (1 + (np.abs(nu - nu.mean()) >= abs(est - nu.mean())).sum()) / (1 + len(nu)); xpos = pi + (oi - 2.5) * w; ci_ = np.percentile(bt, [2.5, 97.5])
            ax.plot([xpos, xpos], ci_, color=OUT_COL[o], lw=0.8); ax.plot(xpos, est, 'o', color=OUT_COL[o], ms=3, mfc=OUT_COL[o] if p < 0.01 else 'white', mew=0.8)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(len(pairs))); ax.set_xticklabels([p[2].replace(' - ', '\n- ') for p in pairs], fontsize=5.3); ax.set_ylabel('difference of kernels (dF/F per lick)')
    ax.set_title('Do the classes differ?  Class differences, same encoding', loc='left', fontsize=6.3)
    # d: ROI sign structure
    ax = fig.add_subplot(gs[1, 2]); label(ax, 'd')
    names = ['lick_wait_early', 'lick_ant_to', 'lick_win', 'lick_post_a', 'lick_post_b']
    for ki, kind in enumerate(('signal', 'noise')):
        z = ex['e1'][kind]['z']; col = fc.C_SIG if kind == 'signal' else fc.C_NOISE
        for ni, nm in enumerate(names):
            npos, nneg = (z[nm] > 3).sum(), (z[nm] < -3).sum(); ax.bar(ni + (ki - 0.5) * 0.38, -nneg / len(z[nm]) * 100, 0.36, color=col, alpha=0.5 if ki else 0.9); ax.bar(ni + (ki - 0.5) * 0.38, npos / len(z[nm]) * 100, 0.36, color=col, alpha=0.5 if ki else 0.9)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(len(names))); ax.set_xticklabels(['early\nwait', 'antic.\ntime-out', 'reward\nwindow', 'post\n3-6 s', 'post\n> 6 s'], fontsize=5.2)
    ax.set_ylabel('% ROIs z > 3 (up) / z < -3 (down)', fontsize=5.5, labelpad=1); ax.set_title('Per-ROI sign: signal-like (blue)\nvs noise-like (amber)', loc='left', fontsize=6.0)
    # e: subsession stability
    ax = fig.add_subplot(gs[2, 0]); label(ax, 'e')
    for o in outs:
        v = np.array(ex['e2'][o]); xs = np.arange(len(v)) + 1 + (outs.index(o) - 2.5) * 0.06
        ax.plot(xs, v[:, 0], '-o', color=OUT_COL[o], ms=2.5, lw=0.7); ax.vlines(xs, v[:, 1], v[:, 2], color=OUT_COL[o], lw=0.6, alpha=0.6)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks([1, 2, 3, 4]); ax.set_xlabel('subsession (30 trials, ~25 min)'); ax.set_ylabel('lick kernel 0-2 s (dF/F per lick)'); ax.set_title('Lick coupling over the session', loc='left', fontsize=6.3)
    # f: cross-session population kernels per class
    ax = fig.add_subplot(gs[2, 1:]); label(ax, 'f')
    sess = [n for n in SESS_ORDER if load('fir', n) is not None]; cls2 = ['wait_early', 'ant_to', 'ant_sq', 'win', 'post_a', 'post_b']
    for si, n in enumerate(sess):
        rr = load('fir', n)
        for ci, c in enumerate(cls2):
            if rr['n_licks'].get(c, 0) < 40: continue
            est, ci_, z, p = fir_stat(rr, 'P', 'lick_' + c); xpos = ci + (si - (len(sess) - 1) / 2) * 0.1
            ax.plot([xpos, xpos], ci_, color=CLS_COL[c], lw=0.7); ax.plot(xpos, est, fc.MARK[fc.animal(n)], color=CLS_COL[c], ms=3, mfc='white' if n in FLAGGED else CLS_COL[c], mew=0.8)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(len(cls2))); ax.set_xticklabels([CLS_LAB[c].replace(' (', '\n(').replace(' -> ', '\n-> ') for c in cls2], fontsize=5.3); ax.set_ylabel('population kernel 0-2 s (dF/F per lick)')
    ax.set_title('Population kernels, every lick session (within a class, left to right: 1506, 1423, 1520, 1744, 1752, 1317*, RB23 1704*; hollow = flagged)', loc='left', fontsize=5.6)
    save(fig, 'D2_lick_class_kernels')


# ================================================================================================= D3
CONTRASTS = [('ant_to_vs_early_to', 'anticipatory vs early wait, both in time-out waits'), ('ant_sq_vs_early_sq', 'anticipatory vs early wait, both in square waits'),
             ('ant_to_vs_wait_early', 'anticipatory (time-out) vs early wait'), ('ant_sq_vs_wait_early', 'anticipatory (square) vs early wait'), ('ant_sq_vs_ant_to', 'anticipatory: square vs time-out'),
             ('rew_vs_wait_early', 'reward bout vs early wait'), ('rew_vs_ant_to', 'reward bout vs anticipatory (time-out)'), ('post_a_vs_wait_early', 'post 3-6 s vs early wait'),
             ('post_b_vs_wait_early', 'post > 6 s vs early wait'), ('post_a_vs_rew', 'post 3-6 s vs reward bout'),
             ('pre_ant_to_vs_wait_early', 'PRE-bout (-2..0 s): anticipatory vs early wait'), ('pre_rew_vs_wait_early', 'PRE-bout: reward vs early wait'), ('pre_post_b_vs_wait_early', 'PRE-bout: post > 6 s vs early wait'),
             ('timeout_post', 'silent time-out (0..3 s after) vs random wait time'), ('timeout_pre', 'before a time-out (-3..0 s) vs random wait time'), ('free_pump_vs_nopump', 'free pump vs no pump (lick-free 0..0.8 s)')]


def get_c(m, o, key):
    if m is None: return None
    c = m['contrasts'].get(o, {}).get(key)
    if not c or not c.get('obs') or not np.isfinite(c['obs'].get('matched', np.nan)): return None
    return c


def fig_d3():
    m1 = load('matched', PRIMARY)
    fig = plt.figure(figsize=(fc.W2, 6.6)); gs = fig.add_gridspec(1, 2, wspace=1.15, left=0.20, right=0.985, top=0.95, bottom=0.07, width_ratios=[1.2, 1])
    ax = fig.add_subplot(gs[0, 0]); label(ax, 'a'); ax.yaxis.set_label_position('left')
    outs = ['P', 'g1', 'g2', 'g3', 'g4', 'control']; h = 0.11
    for ki, (key, lab) in enumerate(CONTRASTS):
        for oi, o in enumerate(outs):
            c = get_c(m1, o, key)
            if c is None: continue
            est = c['obs']['matched']; sd = c['stats']['matched']['nsd']; p = c['stats']['matched']['p']; ypos = ki + (oi - 2.5) * h
            ax.plot([est - 1.96 * sd, est + 1.96 * sd], [ypos, ypos], color=OUT_COL[o], lw=0.8); ax.plot(est, ypos, 'o', color=OUT_COL[o], ms=2.8, mfc=OUT_COL[o] if p < 0.01 else 'white', mew=0.7, label=OUT_LAB[o] if ki == 0 else None)
    ax.axvline(0, color='k', lw=0.4); ax.set_yticks(range(len(CONTRASTS))); ax.set_yticklabels([c[1] for c in CONTRASTS], fontsize=5.6); ax.invert_yaxis(); ax.set_xlabel('behaviour-matched difference of the 30 s high-passed dF/F (event - matched control)')
    ax.legend(frameon=False, fontsize=5, ncol=3, loc='lower right'); ax.set_title('RB24_261109_1506: matched contrasts (+-1.96 x shift-null SD; filled = p < 0.01)', loc='left', fontsize=6.3)
    # replication: population, per session
    ax = fig.add_subplot(gs[0, 1]); label(ax, 'b')
    REPL = ['trial_3_6', 'trial_6_10', 'ant_to_vs_early_to', 'ant_sq_vs_early_sq', 'ant_sq_vs_wait_early', 'rew_vs_wait_early', 'post_a_vs_wait_early', 'post_b_vs_wait_early', 'pre_ant_sq_vs_wait_early', 'pre_rew_vs_wait_early', 'hit_vs_miss_early']
    LABD = dict(CONTRASTS); LABD['hit_vs_miss_early'] = 'hit - miss, lick-free 0..0.8 s'; LABD['trial_3_6'] = 'trial time: 3-6 s after the square vs wait times'; LABD['trial_6_10'] = 'trial time: 6-10 s after the square vs wait times'
    keys = REPL
    sess = [n for n in SESS_ORDER if load('matched', n) is not None]
    for ki, key in enumerate(keys):
        for si, n in enumerate(sess):
            mm = load('matched_tr' if key.startswith('trial_') else 'matched', n); c = get_c(mm, 'P', key)
            if c is None or c['obs']['n_matched'] < 8 or c['obs']['n_matched'] < 0.4 * c['obs']['n_focal']: continue
            est = c['obs']['matched']; sd = c['stats']['matched']['nsd']; ypos = ki + (si - (len(sess) - 1) / 2) * 0.1; col = SESS_COL[n]
            ax.plot([est - 1.96 * sd, est + 1.96 * sd], [ypos, ypos], color=col, lw=0.6, alpha=0.7); ax.plot(est, ypos, fc.MARK[fc.animal(n)], color=col, ms=3.2, mfc='white' if (n in FLAGGED) else col, mew=0.8, label=sh(n) + (' (flagged)' if n in FLAGGED else '') if ki == 0 else None)
    ax.axvline(0, color='k', lw=0.4); ax.set_yticks(range(len(keys))); ax.set_yticklabels([LABD.get(k, {'pre_ant_sq_vs_wait_early': 'PRE-bout: anticipatory (square) vs early wait'}.get(k, k)) for k in keys], fontsize=5.4); ax.invert_yaxis(); ax.set_xlabel('matched difference, population (dF/F)')
    h_, l_ = ax.get_legend_handles_labels(); seen = {}
    for hh, ll in zip(h_, l_): seen.setdefault(ll, hh)
    ax.legend(list(seen.values()), list(seen.keys()), frameon=False, fontsize=5, loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=3)
    ax.set_title('Same contrasts in every lick session (shown when >= 40% of events matched)', loc='left', fontsize=6.0)
    save(fig, 'D3_matched_contrasts')


# ================================================================================================= D4
def fig_d4():
    hm = load('hitmiss', PRIMARY)
    fig = plt.figure(figsize=(fc.W2, 6.0)); gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.55, left=0.07, right=0.985, top=0.94, bottom=0.09)
    ax = fig.add_subplot(gs[0, 0]); label(ax, 'a'); tc = hm['tc']['P']
    ax.plot(tc['tau'], tc['hit_bs'], color=fc.C_SIG, label='hit (%d)' % hm['n_hit']); ax.plot(tc['tau'], tc['miss_bs'], color='#D98A00', label='miss (%d)' % (hm['n'] - hm['n_hit']))
    ax.axvspan(0, 3, color='#06B6D4', alpha=0.12, lw=0); ax.axvline(0, color='k', lw=0.4); ax.axhline(0, color='k', lw=0.4); ax.set_xlabel('time from square onset (s)'); ax.set_ylabel('dF/F (baseline -3..-0.5 s subtracted)')
    ax.legend(frameon=False, fontsize=5.5); ax.set_title('Population: hit and miss trials', loc='left', fontsize=6.5)
    ax = fig.add_subplot(gs[0, 1]); label(ax, 'b')
    wins = ['pre', 'early', 'mid', 'late']; wl = {'pre': '-3..0 s', 'early': '0..0.8 s\n(slow hits)', 'mid': '0.8..3 s', 'late': '3..8 s'}
    m1 = load('matched', PRIMARY)
    for oi, (o, col) in enumerate((('P', fc.C_SIG), ('control', fc.C_NOISE))):
        for wi, wn in enumerate(wins):
            d = hm['win'][o][wn]; x = wi + (oi - 0.5) * 0.36
            ax.plot([x - 0.08] * 2, d['ci'], color=col, lw=0.8); ax.plot(x - 0.08, d['diff'], 'o', color=col, ms=3, label=('raw, ' + OUT_LAB[o]) if wi == 0 else None)
            key = {'pre': 'hit_vs_miss_pre', 'early': 'hit_vs_miss_early', 'late': 'hit_vs_miss_late'}.get(wn)
            if key:
                c = get_c(m1, o, key)
                if c: sd = c['stats']['matched']['nsd']; e = c['obs']['matched']; ax.plot([x + 0.08] * 2, [e - 1.96 * sd, e + 1.96 * sd], color=col, lw=0.8); ax.plot(x + 0.08, e, 's', color=col, ms=3, mfc='white', mew=0.8, label=('matched, ' + OUT_LAB[o]) if wi == 0 else None)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(4)); ax.set_xticklabels([wl[w] for w in wins], fontsize=5.5); ax.set_ylabel('hit - miss (dF/F)'); ax.legend(frameon=False, fontsize=4.8, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    ax.set_title('Hit - miss by window: raw (filled) vs matched (open)', loc='left', fontsize=6.0)
    ax = fig.add_subplot(gs[0, 2]); label(ax, 'c')
    sess = [n for n in SESS_ORDER if load('hitmiss', n) is not None]
    for si, n in enumerate(sess):
        h = load('hitmiss', n); pr = h['predict']; x = si
        ax.plot([x - 0.12, x + 0.12], [pr['auc_m0'], pr['auc_m1']], color='#9CA3AF', lw=0.6); ax.plot(x - 0.12, pr['auc_m0'], 'o', color='#6B7280', ms=3.5, label='movement + history' if si == 0 else None); ax.plot(x + 0.12, pr['auc_m1'], 'o', color=fc.C_SIG, ms=3.5, label='+ axon pre-state' if si == 0 else None)
        ax.text(x, 0.36, 'p=%.2f' % pr['p'], fontsize=5, ha='center')
    ax.axhline(0.5, color='k', lw=0.4); ax.set_ylim(0.33, 0.8); ax.set_xticks(range(len(sess))); ax.set_xticklabels([sh(n).replace('RB24 ', '').replace('RB23 ', 'RB23 ') for n in sess], rotation=60, ha='right', fontsize=5); ax.set_ylabel('cross-validated AUC (hit)')
    ax.legend(frameon=False, fontsize=5, loc='upper right'); ax.set_title('Axon state before the square -> outcome?', loc='left', fontsize=6.3)
    ax = fig.add_subplot(gs[1, 0:2]); label(ax, 'd')
    sc = hm['scan']; kk = sc['kind']
    for wi, wn in enumerate(wins):
        for ti, tag in enumerate(('raw', 'adj')):
            p = sc['win'][wn][tag]['p']
            for ki, (kind, col) in enumerate((('signal', fc.C_SIG), ('noise', fc.C_NOISE))):
                v = 100 * (p[kk == kind] < 0.05).mean(); x = wi + (ki - 0.5) * 0.36 + (ti - 0.5) * 0.15
                ax.bar(x, v, 0.14, color=col, alpha=0.9 if tag == 'raw' else 0.45, hatch=None if tag == 'raw' else '///', edgecolor='white', lw=0.3)
    ax.axhline(5, color='k', lw=0.5, ls=':'); ax.set_xticks(range(4)); ax.set_xticklabels([wl[w] for w in wins], fontsize=5.5); ax.set_ylabel('% of ROIs with hit - miss p < 0.05'); ax.text(3.45, 6, '5% chance', fontsize=5, ha='right')
    ax.set_title('Per-ROI hit - miss: signal-like (blue, n=%d) vs noise-like (amber, n=%d)' % ((kk == 'signal').sum(), (kk == 'noise').sum()), loc='left', fontsize=6.3)
    ax.text(0.01, 0.97, 'solid = raw, hatched = after regressing the window movement', transform=ax.transAxes, fontsize=5.5, va='top')
    ax = fig.add_subplot(gs[1, 2]); label(ax, 'e')
    for si, n in enumerate(sess):
        h = load('hitmiss', n)
        for wi, wn in enumerate(wins):
            d = h['win']['P'][wn]; x = wi + (si - (len(sess) - 1) / 2) * 0.1
            if not np.isfinite(d['diff']): continue
            ax.plot([x, x], d['ci'], color='#9CA3AF', lw=0.5); ax.plot(x, d['diff'], fc.MARK[fc.animal(n)], color=fc.C_SIG if n != PRIMARY else '#E11D48', ms=2.8, mfc='white' if n in FLAGGED else (fc.C_SIG if n != PRIMARY else '#E11D48'), mew=0.7)
    ax.axhline(0, color='k', lw=0.4); ax.set_xticks(range(4)); ax.set_xticklabels([wl[w] for w in wins], fontsize=5.3); ax.set_ylabel('hit - miss, population (raw, dF/F)'); ax.set_title('All lick sessions (red = 1506)', loc='left', fontsize=6.3)
    save(fig, 'D4_hit_vs_miss')


# ================================================================================================= tables
def write_tables():
    rows = []
    for n in SESS_ORDER:
        m = load('matched', n)
        if m is None: continue
        for o, cd in m['contrasts'].items():
            for key, c in cd.items():
                ob = c.get('obs')
                if not ob or not np.isfinite(ob.get('matched', np.nan)): rows.append([n, o, key, c.get('n_events', ''), '', '', '', '', '', '', '', '']); continue
                s = c['stats']['matched']; ss = c.get('stats_sim', {}).get('matched', {})
                rows.append([n, o, key, c['n_events'], ob['n_focal'], ob['n_matched'], round(ob['naive'], 4), round(ob['matched'], 4), round(s['nsd'], 4), round(s['z'], 2) if np.isfinite(s['z']) else '', round(s['p'], 4) if np.isfinite(s['p']) else '', round(ss.get('p', np.nan), 4) if ss else ''])
    with open(os.path.join(dd.OUTDIR, 'matched_contrasts.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['session', 'outcome', 'contrast', 'n_events', 'n_valid', 'n_matched', 'naive', 'matched', 'null_sd', 'shift_z', 'shift_p', 'behaviour_only_bench_p']); w.writerows(rows)
    rows = []
    for n in SESS_ORDER:
        r = load('fir', n)
        if r is None: continue
        for o in r['outcome_names']:
            for nm in r['out'][o]['sum']:
                est, ci, z, p = fir_stat(r, o, nm); rows.append([n, o, nm, round(est, 5), round(ci[0], 5), round(ci[1], 5), round(z, 2), round(p, 4)])
    with open(os.path.join(dd.OUTDIR, 'fir_kernels.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['session', 'outcome', 'regressor', 'kernel_0_2s_dFF_per_event', 'ci_lo', 'ci_hi', 'shift_z', 'shift_p']); w.writerows(rows)
    rows = []
    for n in SESS_ORDER:
        h = load('hitmiss', n)
        if h is None: continue
        for o in h['outcome_names']:
            for wn in ('pre', 'early', 'mid', 'late'):
                d = h['win'][o][wn]; rows.append([n, o, wn, h['n_hit'], h['n'] - h['n_hit'], round(float(d['diff']), 4), round(d['ci'][0], 4), round(d['ci'][1], 4)])
    with open(os.path.join(dd.OUTDIR, 'hit_vs_miss_windows.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['session', 'outcome', 'window', 'n_hit', 'n_miss', 'hit_minus_miss_dFF', 'ci_lo', 'ci_hi']); w.writerows(rows)
    print('tables written')


if __name__ == '__main__':
    which = sys.argv[1:] or ['d1', 'd2', 'd3', 'd4', 'tables']
    if 'd1' in which: fig_d1()
    if 'd2' in which: fig_d2()
    if 'd3' in which: fig_d3()
    if 'd4' in which: fig_d4()
    if 'tables' in which: write_tables()
