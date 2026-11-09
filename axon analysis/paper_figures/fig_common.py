"""Shared style, session lists and statistics helpers of the paper-style figure set (axon analysis/paper_figures).

Conventions (stated once here, repeated in every legend):
  colour = what the thing IS  : signal-like axons blue, noise-like ROIs amber, rejected suite2p ROIs grey, behaviours their own hue (never reused for anything else)
  shape  = which animal       : RB24 circle, RB23 square, RB22 triangle (so colour stays free for the data class; shape is also readable in grey-scale)
  filled / hollow             : unflagged session (counts in the pooled numbers) / flagged session (drawn, never pooled)
  unit of analysis            : session, clustered by animal -> pooled value = mean of the per-animal means, CI = hierarchical bootstrap (sessions within animal, then animals);
                                only 3 mice exist, so every pooled number is a description of these 3 mice, not a population estimate.
"""
import os, sys, csv, glob, pickle, numpy as np

REPO = r'c:\Users\rbondarenko\projects\2P_Scripts'
CS_DIR = os.path.join(REPO, 'axon analysis', 'cross_session')
CS_FIG = os.path.join(REPO, 'axon analysis', 'plots', 'cross_session')
OUT = os.path.join(REPO, 'axon analysis', 'plots', 'summary_figures')          # (this folder was called plots/paper_figures until the user renamed it on 2026-11-09)
CACHE = os.path.join(OUT, 'cache')
for d in (OUT, CACHE): os.makedirs(d, exist_ok=True)
if CS_DIR not in sys.path: sys.path.insert(0, CS_DIR)

import matplotlib, matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

MM = 1 / 25.4
W1, W2 = 89 * MM, 183 * MM                     # single / double column width (inches)
_FONT = 'Arial' if any('Arial' in f.name for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams.update({'font.family': _FONT, 'font.size': 6.5, 'axes.labelsize': 6.5, 'axes.titlesize': 7, 'axes.titleweight': 'bold', 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'legend.fontsize': 6, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.2, 'ytick.major.size': 2.2,
                     'xtick.direction': 'out', 'ytick.direction': 'out', 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9,
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 110, 'savefig.dpi': 300, 'axes.titlepad': 3, 'axes.labelpad': 2})

# ---- colours: what the thing is
C_SIG, C_NOISE, C_REJ = '#1F5FBF', '#D98A00', '#8A8F98'
INK, MUTED, FAINT = '#222222', '#6B7280', '#D1D5DB'
C_BEH = {'face': INK, 'whisk': '#2E8B57', 'pupil': '#C2185B', 'ball': '#8C6D31', 'eyelid': '#7E57C2', 'shift': '#9A6B4F', 'dpupil': '#E57AA3'}
BEH_LAB = {'face': 'face motion', 'whisk': 'whisking', 'pupil': 'pupil area', 'ball': 'locomotion (ball)', 'loco': 'locomotion (ball)', 'eyelid': 'eyelid opening', 'shift': 'brain shift', 'dpupil': 'pupil dilation rate'}
MARK = {'RB24': 'o', 'RB23': 's', 'RB22': '^'}

# ---- sessions
QC_PATH = os.path.join(CS_FIG, 'session_qc.csv')
FLAG = {   # flagged sessions: drawn everywhere, never in the pooled numbers (reason printed in the inventory and in Figure 7)
    'AnimalRB23_261027_1306': 'face / whisking / pupil tracking inconsistent (face-whisk r 0.09), 4 signal-like ROIs',
    'AnimalRB23_261029_1704': 'face / whisking / pupil tracking inconsistent (face-whisk r 0.13), 3 signal-like ROIs',
    'AnimalRB24_261023_1527': 'interrupted, 6 min: circular null has few independent lags',
    'AnimalRB24_261104_1640': '7 min, 5-plane: circular null has few independent lags, no recorded grey',
    'AnimalRB24_261106_1529': '11-06: very little behavioural range, ball reads running 79%',
    'AnimalRB24_261106_1317': '11-06: very little behavioural range, ball reads running 82%',
    'AnimalRB24_261106_1604': '11-06: very little behavioural range, ball reads running 71%',
}


# sessions with a spout + reward (lick task) and with a RECORDED >= 120 s grey around the loom blocks (from the cross-session analyses 3 + 4)
TASK_SESSIONS = {'AnimalRB24_261029_1752', 'AnimalRB24_261030_1744', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261106_1317', 'AnimalRB23_261029_1704'}
GREY_SESSIONS = {'AnimalRB24_261023_1538', 'AnimalRB24_261027_1534', 'AnimalRB24_261106_1604', 'AnimalRB23_261023_1642', 'AnimalRB23_261027_1306', 'AnimalRB22_261027_1138'}


def animal(nm): return nm.replace('Animal', '')[:4]
def short(nm): a, d, t = nm.replace('Animal', '').split('_'); return '%s %s-%s %s' % (a, d[2:4], d[4:6], t)
def is_flagged(nm): return nm in FLAG
def session_order(names): return sorted(names, key=lambda n: (animal(n)[::-1], n))      # RB22, RB23, RB24; date inside


def _qc():
    rows = {}
    for r in csv.DictReader(open(QC_PATH)): rows[r['session']] = r
    return rows


QC = _qc()
ALL = session_order(QC)
UNFLAGGED = [n for n in ALL if not is_flagged(n)]
FLAGGED = [n for n in ALL if is_flagged(n)]


def qcf(name, key):
    try: return float(QC[name][key])
    except Exception: return np.nan


def load_state():
    res = {}
    for f in sorted(glob.glob(os.path.join(CS_FIG, 'state_part*.pkl'))): res.update(pickle.load(open(f, 'rb')))
    return res


def load_pickle(name, folder=CS_FIG): return pickle.load(open(os.path.join(folder, name), 'rb'))


# ---- drawing helpers
def sdot(ax, x, y, name, ms=4.6, color=INK, hollow_color=None, zorder=3, **kw):
    """session marker: shape = animal, filled = unflagged, hollow = flagged"""
    fl = is_flagged(name); mk = MARK[animal(name)]
    if fl: return ax.plot(x, y, mk, mfc='white', mec=hollow_color or color, mew=0.8, ms=ms, zorder=zorder, **kw)
    return ax.plot(x, y, mk, mfc=color, mec='white', mew=0.3, ms=ms, zorder=zorder, **kw)


def plabel(ax, letter, dx=-16, dy=3):
    ax.annotate(letter, xy=(0, 1), xycoords='axes fraction', xytext=(dx, dy), textcoords='offset points', fontsize=9, fontweight='bold', va='bottom', ha='left', annotation_clip=False)


def clean(ax, left=True, bottom=True):
    if not left: ax.spines['left'].set_visible(False); ax.tick_params(axis='y', left=False, labelleft=False)
    if not bottom: ax.spines['bottom'].set_visible(False); ax.tick_params(axis='x', bottom=False, labelbottom=False)


def scalebar(ax, x0, y0, dx=None, dy=None, lx=None, ly=None, color=INK, fs=5.5):
    """scale bar(s) in DATA units anchored at (x0, y0); lx / ly are the labels"""
    if dx: ax.plot([x0, x0 + dx], [y0, y0], color=color, lw=0.9, clip_on=False, solid_capstyle='butt'); ax.text(x0 + dx / 2, y0, lx, ha='center', va='top', fontsize=fs, color=color)
    if dy: ax.plot([x0, x0], [y0, y0 + dy], color=color, lw=0.9, clip_on=False, solid_capstyle='butt'); ax.text(x0, y0 + dy / 2, ly, ha='right', va='center', fontsize=fs, color=color, rotation=90)


def savefig(fig, name):
    for ext in ('pdf', 'png'): fig.savefig(os.path.join(OUT, name + '.' + ext), facecolor='white', bbox_inches=None if ext == 'png' else None)
    plt.close(fig); print('saved', name)


def fmt_p(p):
    if p is None or not np.isfinite(p): return 'p n/a'
    return 'p < 0.001' if p < 1e-3 else 'p = %.3f' % p if p < 0.1 else 'p = %.2f' % p


# ---- statistics
def zs(v): v = np.asarray(v, float); return (v - np.nanmean(v)) / np.nanstd(v)


def hboot(vals, names, n_boot=10000, seed=0):
    """pooled effect = mean over animals of the per-animal mean of the session values; CI = hierarchical bootstrap (animals with replacement, sessions within animal with replacement).
    returns dict(est, lo, hi, per_animal{animal: (mean, n)}, n_sessions, n_animals, k_pos = sessions with value > 0)"""
    vals = np.asarray(vals, float); ok = np.isfinite(vals); vals = vals[ok]; an = np.array([animal(n) for n, o in zip(names, ok) if o])
    if len(vals) == 0: return dict(est=np.nan, lo=np.nan, hi=np.nan, per_animal={}, n_sessions=0, n_animals=0, k_pos=0)
    ua = sorted(set(an)); groups = {a: vals[an == a] for a in ua}
    est = float(np.mean([g.mean() for g in groups.values()]))
    rng = np.random.default_rng(seed); bs = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.choice(ua, len(ua), replace=True)
        bs[b] = np.mean([rng.choice(groups[a], len(groups[a]), replace=True).mean() for a in pick])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(est=est, lo=float(lo), hi=float(hi), per_animal={a: (float(g.mean()), len(g)) for a, g in groups.items()}, n_sessions=len(vals), n_animals=len(ua), k_pos=int((vals > 0).sum()))


def sign_flip_p(vals, n_perm=20000, seed=0):
    """exact-ish sign-flip test over sessions (H0: symmetric around 0), two-sided; with n <= 14 sessions the p floor is 2 / 2^n"""
    v = np.asarray(vals, float); v = v[np.isfinite(v)]
    if len(v) < 3: return np.nan
    obs = abs(v.mean()); rng = np.random.default_rng(seed)
    if len(v) <= 14:
        import itertools
        signs = np.array(list(itertools.product([-1, 1], repeat=len(v)))); m = np.abs((signs * v).mean(1))
        return float((m >= obs - 1e-12).mean())
    s = rng.choice([-1, 1], (n_perm, len(v))); return float((np.abs((s * v).mean(1)) >= obs - 1e-12).mean())


# ---- results ledger (collected by every figure script, drawn / exported by fig7)
import json
LEDGER_PATH = os.path.join(CACHE, 'ledger.json')


def ledger_add(figure, key, label, est, lo, hi, n_sessions, n_animals, verdict, note='', unit=''):
    """verdict: 'supported' | 'ambiguous' | 'not seen'. Replaces an existing entry with the same (figure, key)."""
    if verdict == 'supported' and n_animals < 2: verdict = 'ambiguous'; note = (note + ' | one animal only').strip(' |')
    L = json.load(open(LEDGER_PATH)) if os.path.exists(LEDGER_PATH) else []
    L = [r for r in L if not (r['figure'] == figure and r['key'] == key)]
    f = lambda v: None if v is None or not np.isfinite(v) else float(v)
    L.append(dict(figure=figure, key=key, label=label, est=f(est), lo=f(lo), hi=f(hi), n_sessions=int(n_sessions), n_animals=int(n_animals), verdict=verdict, note=note, unit=unit))
    json.dump(L, open(LEDGER_PATH, 'w'), indent=1)


def pooled_unflagged(vals_by_session, n_boot=10000):
    """hboot over the unflagged sessions of a dict session -> value"""
    ns = [n for n in UNFLAGGED if n in vals_by_session and np.isfinite(vals_by_session[n])]
    return hboot([vals_by_session[n] for n in ns], ns, n_boot=n_boot)
