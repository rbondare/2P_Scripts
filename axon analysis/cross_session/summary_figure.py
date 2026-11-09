"""One-page summary of what the RB22-24 data support (built from the cross-session analyses 1-4).
usage: python summary_figure.py   -> plots/cross_session/SUMMARY_what_the_data_support.png / .pdf
Colour = what the thing is (same as the S / L / G / T figures): signal-like axons blue, axons with the concurrent movement removed hollow, rejected suite2p ROIs grey,
behaviour dark ink, animals RB24 blue / RB23 red / RB22 teal (session labels only).
"""
import os, csv, pickle, sys, textwrap, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import spearmanr
import xs_common as xc
from stim_response import pooled
FIG = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session'
plt.rcParams.update({'font.family': 'Arial' if any('Arial' in f.name for f in matplotlib.font_manager.fontManager.ttflist) else 'DejaVu Sans', 'font.size': 8, 'axes.linewidth': 0.75, 'axes.spines.top': False, 'axes.spines.right': False,
                     'xtick.direction': 'out', 'ytick.direction': 'out', 'figure.dpi': 110, 'savefig.dpi': 220, 'axes.titlesize': 8.6, 'axes.titleweight': 'bold'})
C_SIG, C_REJ, INK, MUTED = '#1F5FBF', '#8A8F98', '#222222', '#6B7280'
C_ANIMAL = {'RB24': '#1F5FBF', 'RB23': '#D1495B', 'RB22': '#1B998B'}
def short(nm): a, d, t = nm.replace('Animal', '').split('_'); return '%s  %s-%s  %s' % (a, d[2:4], d[4:6], t)
def animal(nm): return nm.replace('Animal', '')[:4]

# ---------------------------------------------------------------- data
qc = {r['session']: r for r in csv.DictReader(open(os.path.join(FIG, 'session_qc.csv')))}
names = sorted(qc, key=lambda n: (animal(n)[::-1], n))
def fl(r, k):
    try: return float(r[k])
    except Exception: return np.nan
FLAG = {'AnimalRB23_261027_1306', 'AnimalRB23_261029_1704', 'AnimalRB24_261106_1529', 'AnimalRB24_261106_1317', 'AnimalRB24_261106_1604', 'AnimalRB24_261023_1527', 'AnimalRB24_261104_1640'}
stim = pickle.load(open(os.path.join(FIG, 'stim_response_results.pkl'), 'rb')); task = pickle.load(open(os.path.join(FIG, 'task_results.pkl'), 'rb')); grey = pickle.load(open(os.path.join(FIG, 'grey_results.pkl'), 'rb'))
n_clear = sum(1 for n in names if fl(qc[n], 'r_face') >= 0.4 and fl(qc[n], 'p_face') < 0.01 and n not in FLAG)
cnt = {a: sum(1 for n in names if animal(n) == a) for a in ('RB24', 'RB23', 'RB22')}

fig = plt.figure(figsize=(17.5, 14.8))
gs = fig.add_gridspec(5, 13, height_ratios=[1.25, 1.25, 1.3, 1.45, 1.0], width_ratios=[1] * 7 + [0.95] + [1] * 5, hspace=0.7, wspace=0.55, left=0.03, right=0.99, top=0.925, bottom=0.015)
fig.text(0.005, 0.985, 'RB22 - RB24 LC axons in the superior colliculus: what the data support', fontsize=14, fontweight='bold', va='top')
fig.text(0.005, 0.962, '%d sessions, 3 mice (RB24 %d, RB23 %d, RB22 %d sessions with accepted axon ROIs); 1-plane and 5-plane imaging; exploratory, one pass. Signal-like = ROIs whose trace passes the V3 axon criterion; '
         '"movement removed" = a behaviour-only model (face + whisking + ball) subtracted.' % (len(names), cnt['RB24'], cnt['RB23'], cnt['RB22']), fontsize=8.4, color=MUTED, va='top')

# ---------------------------------------------------------------- A  example
S = xc.Session('AnimalRB24_261021_1630', with_rejected=True, verbose=False)
P = S.population('sig'); Psm = xc.nsmooth(P, 1.0 / xc.DT); face = S.beh['face']; pup = S.beh['pupil']
rng = np.random.default_rng(2); R = S.rej['dff']; rej = np.nanmean(R[rng.choice(len(R), int(S.sig.sum()), replace=False)], 0); rejsm = xc.nsmooth(rej, 1.0 / xc.DT)
w = int(240 / xc.DT); best, bi = -1, 0
for i in range(0, S.nb - w, int(20 / xc.DT)):
    seg = Psm[i:i + w]; f = face[i:i + w]
    if np.isfinite(seg).mean() > 0.95 and np.isfinite(f).mean() > 0.95:
        sc = np.nanstd(seg) * np.nanstd(f)
        if sc > best: best, bi = sc, i
tt = S.t[bi:bi + w] - S.t[bi]
sub = gs[0, 0:7].subgridspec(3, 1, hspace=0.12, height_ratios=[1.6, 1, 1])
a1 = fig.add_subplot(sub[0]); a1.plot(tt, rejsm[bi:bi + w], color=C_REJ, lw=0.8, label='random set of rejected ROIs'); a1.plot(tt, Psm[bi:bi + w], color=C_SIG, lw=1.3, label='signal-like axons (mean)')
a1.set_ylim(-0.15, 2.5); a1.set_ylabel('dF/F', fontsize=7.5); a1.legend(frameon=False, fontsize=6.8, loc='upper left', ncol=2); a1.set_xticklabels([]); a1.set_title('A   The shared arousal / movement signal  (RB24 10-21 1630, 4 min excerpt)', loc='left')
a2 = fig.add_subplot(sub[1], sharex=a1); a2.plot(tt, face[bi:bi + w], color=INK, lw=0.9); a2.set_ylabel('face motion\n(log)', fontsize=7.5); a2.set_xticklabels([])
a3 = fig.add_subplot(sub[2], sharex=a1); a3.plot(tt, pup[bi:bi + w], color=MUTED, lw=0.9); a3.set_ylabel('pupil\narea', fontsize=7.5); a3.set_xlabel('seconds', fontsize=7.5)
ok = np.isfinite(Psm) & np.isfinite(face); r_all = np.corrcoef(xc.bandpass(P)[ok], xc.bandpass(face)[ok])[0, 1]
a1.text(0.99, 0.97, 'whole session: r(axons, face motion) = %.2f' % r_all, transform=a1.transAxes, fontsize=7.2, va='top', ha='right', color=INK)

# ---------------------------------------------------------------- B, C  coupling across sessions and sign test
sub = gs[0:3, 8:13].subgridspec(1, 2, wspace=0.08, width_ratios=[1.25, 1])
ab = fig.add_subplot(sub[0]); ac = fig.add_subplot(sub[1], sharey=ab); y = np.arange(len(names))
for i, n in enumerate(names):
    r = qc[n]; flagged = n in FLAG; al = 0.5 if flagged else 1.0
    if flagged: ab.axhspan(i - 0.5, i + 0.5, color='#F3F4F6', lw=0, zorder=0); ac.axhspan(i - 0.5, i + 0.5, color='#F3F4F6', lw=0, zorder=0)
    for k, mk in (('r_face', 'o'), ('r_whisk', 's'), ('r_pupil', '^')):
        v = fl(r, k)
        if np.isfinite(v): ab.plot(v, i, mk, color=INK, ms=4.6, alpha=al, mfc=INK if (k != 'r_face' or fl(r, 'p_face') < 0.01) else 'white', mew=0.8)
    sg, rj = fl(r, 'roi_net_sign_signal_face'), fl(r, 'roi_net_sign_rejected_face')
    if np.isfinite(rj): ac.plot([sg, rj], [i, i], color='#D1D5DB', lw=1.2, zorder=1)
    if np.isfinite(sg): ac.plot(sg, i, 'o', color=C_SIG, ms=5, alpha=al, zorder=3)
    if np.isfinite(rj): ac.plot(rj, i, 'o', color=C_REJ, ms=5, alpha=al, zorder=3)
ab.axvline(0, color=MUTED, lw=0.5); ac.axvline(0, color=MUTED, lw=0.5); ab.set_yticks(y); ab.set_yticklabels([short(n) + (' \u2020' if n in FLAG else '') for n in names], fontsize=6.6); ab.set_ylim(len(names) - 0.4, -0.6)
for t, n in zip(ab.get_yticklabels(), names): t.set_color(C_ANIMAL[animal(n)])
plt.setp(ac.get_yticklabels(), visible=False); ac.tick_params(axis='y', length=0)
ab.set_xlabel('peak r of the axon mean\nwith each behaviour (+-10 s)', fontsize=7); ac.set_xlabel('net sign of single-ROI\ncoupling to face motion\n(+1 = all the same way)', fontsize=7); ab.set_xlim(-0.35, 1.0); ac.set_xlim(-0.5, 1.1)
ab.set_title('B   Coupling to behaviour', loc='left'); ac.set_title('C   ROI sign', loc='left')
from matplotlib.lines import Line2D
axl = fig.add_subplot(gs[3:5, 8:13]); axl.axis('off')
hl = [Line2D([], [], marker='o', ls='', color=INK, ms=5, label='B: face motion (filled = p < 0.01 vs shifted data)'), Line2D([], [], marker='s', ls='', color=INK, ms=5, label='B: whisking'), Line2D([], [], marker='^', ls='', color=INK, ms=5, label='B: pupil'),
      Line2D([], [], marker='o', ls='', color=C_SIG, ms=5, label='C: signal-like ROIs, net sign of the single-ROI coupling to face motion'), Line2D([], [], marker='o', ls='', color=C_REJ, ms=5, label='C: rejected ROIs (random mix of signs = no shared signal)')]
axl.legend(handles=hl, frameon=False, fontsize=6.8, loc='upper left', bbox_to_anchor=(0.0, 1.0))
axl.text(0.0, 0.74, textwrap.fill('%d of %d sessions are clearly coupled (r >= 0.4, p < 0.01 against circularly shifted data); the 2 short sessions (1527, 1640) go the same way. \u2020 shaded rows: low behavioural range (11-06), inconsistent face / whisking / pupil tracking (RB23 1306, 1704) or very short sessions.' % (n_clear, len(names)), 62), fontsize=6.8, color=MUTED, va='top', transform=axl.transAxes, linespacing=1.3)

# ---------------------------------------------------------------- D  stimulus responses pooled
ad = fig.add_subplot(gs[1, 0:3]); keys = [('flash', 'flash'), ('loom_expand', 'expanding\nloom'), ('loom_recede', 'receding\nloom'), ('bar', 'moving\nbar'), ('square', 'task\nsquare')]
for i, (k, lab) in enumerate(keys):
    for j, (ser, col, mfc) in enumerate((('P', C_SIG, C_SIG), ('res_mv', C_SIG, 'white'))):
        d = pooled(stim, k, ser)
        if d is None: continue
        x = i + (j - 0.5) * 0.3; ad.plot([x, x], d['null95'], color='#D1D5DB', lw=5, solid_capstyle='butt', zorder=1)
        ad.plot(x, d['mean'], 'o', color=col, mfc=mfc, ms=7, mew=1.3, zorder=3)
        if j == 0: ad.text(i, -0.055, 'n=%d' % d['n_sessions'], ha='center', fontsize=6.4, color=MUTED)
ad.axhline(0, color=MUTED, lw=0.5); ad.set_xticks(range(len(keys))); ad.set_xticklabels([k[1] for k in keys], fontsize=7); ad.set_ylim(-0.065, 0.12); ad.set_xlim(-0.6, len(keys) - 0.4)
ad.set_ylabel('response, mean dF/F (30 s high-passed)', fontsize=7); ad.set_title('D   Stimulus responses of the axon mean', loc='left')
ad.plot([], [], 'o', color=C_SIG, ms=6, label='raw'); ad.plot([], [], 'o', color=C_SIG, mfc='white', ms=6, mew=1.3, label='movement removed'); ad.plot([], [], '-', color='#D1D5DB', lw=5, label='95% of onsets shifted in time')
ad.legend(frameon=False, fontsize=6.3, loc='upper right')

# ---------------------------------------------------------------- E  loom amplitude vs evoked pupil
ae = fig.add_subplot(gs[1, 3:7]); xs_, ys_ = [], []
for n, r in stim.items():
    d = r.get('keys', {}).get('loom_expand')
    if d and 'b_pupil' in d['amp']:
        a, b = d['amp']['P'], d['amp']['b_pupil']; ok = np.isfinite(a) & np.isfinite(b); xs_ += list(b[ok]); ys_ += list(a[ok])
xs_, ys_ = np.array(xs_), np.array(ys_); rho, p = spearmanr(xs_, ys_)
ae.plot(xs_, ys_, 'o', color=C_SIG, ms=3.6, alpha=0.55, mew=0); m, c0 = np.polyfit(xs_, ys_, 1); xx = np.linspace(xs_.min(), xs_.max(), 10); ae.plot(xx, m * xx + c0, color=INK, lw=1.0)
ae.axhline(0, color=MUTED, lw=0.4); ae.axvline(0, color=MUTED, lw=0.4); ae.set_xlabel('evoked pupil change per trial (z, 0-4 s)', fontsize=7); ae.set_ylabel('axon response per trial (dF/F)', fontsize=7); ae.set_title('E   Expanding loom: response size follows the arousal change', loc='left')
ae.text(0.03, 0.96, 'Spearman rho = %+.2f, %s\n%d trials, %d sessions' % (rho, 'p < 0.0001' if p < 1e-4 else 'p = %.4f' % p, len(xs_), sum(1 for r in stim.values() if r.get('keys', {}).get('loom_expand'))), transform=ae.transAxes, fontsize=7, va='top')
RHO = rho

# ---------------------------------------------------------------- F  task
sub = gs[2, 0:3].subgridspec(1, 2, wspace=0.75, width_ratios=[1, 1])
af = fig.add_subplot(sub[0])
for nm, ls in (('AnimalRB24_261106_1317', '-'), ('AnimalRB24_261030_1744', '--')):
    ss = task[nm]['subs']; x = [s_['sub'] for s_ in ss]
    af.plot(x, [100 * s_['hit'] / ss[0]['hit'] for s_ in ss], ls, color=INK, marker='o', ms=4, lw=1.2); af.plot(x, [100 * s_['P_sd'] / ss[0]['P_sd'] for s_ in ss], ls, color=C_SIG, marker='o', ms=4, lw=1.2)
af.axhline(100, color=MUTED, lw=0.4, ls=':'); af.set_xticks([1, 2, 3, 4]); af.set_xlabel('subsession', fontsize=7); af.set_ylabel('% of subsession 1', fontsize=7); af.set_ylim(40, 160)
af.plot([], [], '-', color=INK, marker='o', ms=4, label='hit rate'); af.plot([], [], '-', color=C_SIG, marker='o', ms=4, label='axon fluctuation SD'); af.plot([], [], '-', color=MUTED, lw=1.2, label='1317'); af.plot([], [], '--', color=MUTED, lw=1.2, label='1744')
af.legend(frameon=False, fontsize=6, loc='upper right', ncol=2, columnspacing=0.8, handlelength=1.4); af.set_title('F   Lick task', loc='left')
ag = fig.add_subplot(sub[1]); tn = [n for n in task if not task[n].get('flag') and 'hit_miss' in task[n]]
for i, n in enumerate(tn):
    h = task[n]['hit_miss']; ag.plot([h['P_hp response']['diff'], h['P_res response']['diff']], [i, i], color='#D1D5DB', lw=1.2, zorder=1)
    ag.plot(h['P_hp response']['diff'], i, 'o', color=C_SIG, ms=5.5, zorder=3); ag.plot(h['P_res response']['diff'], i, 'o', color=C_SIG, mfc='white', ms=5.5, mew=1.2, zorder=3)
ag.axvline(0, color=INK, lw=0.5); ag.set_yticks(range(len(tn))); ag.set_yticklabels([short(n).replace('RB24  ', '') for n in tn], fontsize=6.0); ag.set_ylim(len(tn) - 0.4, -0.6)
ag.set_xlabel('hit - miss response (dF/F)\nfilled raw, hollow movement removed', fontsize=6.8)

# ---------------------------------------------------------------- G  grey
ah = fig.add_subplot(gs[2, 4:7]); con = grey['contrasts']
rows = [(('sig', 'mean'), 'axons: mean dF/F', C_SIG), (('sig', 'sd'), 'axons: fluctuation', C_SIG), (('sig', 'rate'), 'axons: event rate', C_SIG), (('rej', 'sd'), 'rejected: fluctuation', C_REJ),
        (('beh', 'face'), 'face motion', INK), (('beh', 'whisk'), 'whisking', INK), (('beh', 'pupil'), 'pupil area', INK), (('beh', 'eyelid'), 'eyelid opening', INK)]
for i, (k, lab, col) in enumerate(rows):
    c = con[k]['loom']
    if not c: continue
    ps = list(c['per_session'].values()); ah.plot(ps, [i] * len(ps), '.', color=col, ms=5, alpha=0.45); ah.plot(c['ci'], [i, i], color=col, lw=2.4, solid_capstyle='round'); ah.plot(c['mean'], i, 'o', color=col, ms=6.5, mec='white', mew=0.7)
    ah.text(1.01, i, 'p=%.2f' % c['p'], transform=ah.get_yaxis_transform(), fontsize=6.4, va='center', color=INK if c['p'] < 0.05 else MUTED)
ah.axvline(0, color=INK, lw=0.6); ah.set_yticks(range(len(rows))); ah.set_yticklabels([r[1] for r in rows], fontsize=7); ah.set_ylim(len(rows) - 0.4, -0.6)
for t, r in zip(ah.get_yticklabels(), rows): t.set_color(r[2])
ah.set_xlabel('after minus before a loom block (SD of the metric)\ndots = sessions, bar = 95% CI over sessions', fontsize=7); ah.set_title('G   Grey screen before vs after looming (6 sessions, 21 blocks)', loc='left')

# ---------------------------------------------------------------- H  how much weight (text strip)
ax = fig.add_subplot(gs[3, 0:8]); ax.axis('off'); ax.set_title('H   How much weight each result can carry', loc='left')
items = [('STRONG', '#1F3F8F', 'white', 'Shared arousal / movement signal in the axons.', '%d of %d sessions, 3 mice: signal-like ROIs share ONE signal (PC1 45-80%% of variance) that follows face motion, whisking and ball at lag 0 (pupil +1.2 s, slow kinetics); the ROIs the experimenter rejected do not (sign test, panel C).' % (n_clear, len(names))),
         ('MODERATE', '#4F7FD0', 'white', 'Axons respond to flashes beyond movement.', '+0.07 dF/F (6 sessions, p 0.001), +0.044 after removing movement (p 0.001). Flashes also constrict the pupil (-1.6 SD), so a residual arousal / light-reflex contribution is not excluded.'),
         ('WEAK', '#8FB3E8', INK, 'Expanding-loom response.', '+0.043 dF/F (8 sessions, p 0.005), +0.019 after removing movement (p 0.06); per trial it scales with the evoked pupil change (rho %+.2f). Receding disc, moving bar, task square: nothing.' % RHO),
         ('NOT SEEN', '#E5E7EB', INK, 'Anything specific to looming, reward, task engagement, or the grey after looming.', 'Hit - miss differences vanish once movement is removed; in 1317 the hit rate falls to 59% of its start while axon fluctuations stay flat; 21 loom blocks show no axon change (behaviour calms, axons do not follow).'),
         ('OPEN', '#FFFFFF', INK, 'The three 11-06 sessions (1529, 1317, 1604) show no coupling.', 'Behaviour barely varies there (face-motion SD 0.17-0.24 vs 0.4-0.8 earlier) and the shared signal is weaker (PC1 0.17-0.58): rig / FOV / expression change, or just a quiet, habituated animal?')]
for i, (chip, bg, fg, head, body) in enumerate(items):
    yy = 0.9 - i * 0.2
    ax.add_patch(FancyBboxPatch((0.0, yy - 0.12), 0.095, 0.14, boxstyle='round,pad=0.0,rounding_size=0.02', fc=bg, ec=INK if bg == '#FFFFFF' else bg, lw=0.8, transform=ax.transAxes))
    ax.text(0.0475, yy - 0.05, chip, transform=ax.transAxes, ha='center', va='center', fontsize=7, fontweight='bold', color=fg)
    ax.text(0.11, yy, head, transform=ax.transAxes, fontsize=8, fontweight='bold', va='center', color=INK)
    ax.text(0.11, yy - 0.04, textwrap.fill(body, 150), transform=ax.transAxes, fontsize=6.9, va='top', color=INK, linespacing=1.25)

# ---------------------------------------------------------------- I  limits
ax = fig.add_subplot(gs[4, 0:8]); ax.axis('off')
ax.set_title('I   What limits every conclusion above (fix before the next experiments)', loc='left')
cols = [('No reference channel', 'Only GCaMP is imaged: axial (z) motion and neuropil cannot be separated from signal. The 5-plane z-registration failed (planes almost identical at 3 um steps).'),
        ('Little behavioural variance', 'Habituated animals move little (11-06); the state signal needs arousal transitions to be measured. The pupil is confounded by screen luminance (flashes, breaks).'),
        ('Weak logs on the rig', 'Ball log truncated in open-loop blocks, 2nd ball sensor dead, 3.5% projector sync loss, lick line not on an aux channel, pump TTL not logged, rig-code reward rule unconfirmed.'),
        ('Design', 'No no-loom control block, 5-20 trials per class, fixed block order, grey screens sometimes declared but not shown, ROI curation not blind to the traces, axon identity / histology not established.')]
for i, (h, b) in enumerate(cols):
    x0 = i * 0.255
    ax.text(x0, 0.86, h, fontsize=7.8, fontweight='bold', va='top', color=INK, transform=ax.transAxes); ax.text(x0, 0.62, textwrap.fill(b, 38), fontsize=6.9, va='top', color=INK, transform=ax.transAxes, linespacing=1.25)
fig.savefig(os.path.join(FIG, 'SUMMARY_what_the_data_support.png'), facecolor='white', bbox_inches='tight'); fig.savefig(os.path.join(FIG, 'SUMMARY_what_the_data_support.pdf'), facecolor='white', bbox_inches='tight')
print('saved SUMMARY (n_clear = %d of %d, loom rho %.2f)' % (n_clear, len(names), RHO))
