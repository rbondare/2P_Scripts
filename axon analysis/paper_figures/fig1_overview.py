"""Figure 1 - what was recorded: data streams, FOV + ROI classes, one full session, inventory of all sessions.
usage: python fig1_overview.py   -> plots/summary_figures/Fig1_overview.pdf/.png (+ cache/inventory.pkl so that the figure rebuilds without the Z: drive for the inventory)
"""
import os, re, pickle, numpy as np
from fig_common import *
import xs_common as xc
from scipy.ndimage import binary_erosion
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

Z = r'Z:/Group Members/Rima/DATA_2P/'
EX_1P, EX_5P = 'AnimalRB24_261021_1630', 'AnimalRB24_261104_1651'


# ------------------------------------------------------------------ data
def inventory_extra():
    """sessions recorded but never aggregated: what the disk shows (TIFF count, suite2p run, ROIs with iscell > 0)"""
    f = os.path.join(CACHE, 'inventory.pkl')
    if os.path.exists(f) and not os.environ.get('REFRESH_INV'): return pickle.load(open(f, 'rb'))
    agg = {os.path.basename(p).replace('_preprocessed.mat', '') for p in __import__('glob').glob(r'Z:/Group Members/Rima/Aggregated/*_preprocessed.mat')}
    rows = []
    for s in sorted(os.listdir(Z)):
        if not any(a in s for a in ('RB22', 'RB23', 'RB24')) or s in agg: continue
        p = Z + s; tifs = [x for x in os.listdir(p) if x.endswith('.tif')]; has = os.path.isdir(p + '/suite2p/plane0'); acc = -1; nroi = -1; fs = np.nan
        if has:
            try:
                ic = np.load(p + '/suite2p/plane0/iscell.npy'); nroi = len(ic); acc = int((ic[:, 0] > 0).sum())
                fs = float(np.load(p + '/suite2p/plane0/ops.npy', allow_pickle=True).item().get('fs', np.nan))
            except Exception: pass
        rows.append(dict(session=s, tifs=len(tifs), suite2p=has, rois=nroi, accepted=acc, fs=fs))
    pickle.dump(rows, open(f, 'wb')); return rows


def parse_stim(name):
    s = QC[name]['stimuli']; d = {}
    for m in re.finditer(r'(\w+) x(\d+)', s): d[m.group(1)] = int(m.group(2))
    return d


def fov_panel(ax, name, S, label_rois=None):
    p = Z + name + '/suite2p/plane0/'
    ops = np.load(p + 'ops.npy', allow_pickle=True).item(); stat = np.load(p + 'stat.npy', allow_pickle=True); ic = np.load(p + 'iscell.npy')[:, 0] > 0
    img = ops.get('max_proj'); yr, xr = ops['yrange'], ops['xrange']
    if img is None: img = ops['meanImg']
    full = np.zeros((ops['Ly'], ops['Lx']), np.float32); full[yr[0]:yr[1], xr[0]:xr[1]] = img
    sub = full[yr[0]:yr[1], xr[0]:xr[1]]; lo, hi = np.percentile(sub, [5, 99.6]); g = np.clip((sub - lo) / (hi - lo), 0, 1) ** 0.7
    rgb = np.dstack([g, g, g]) * 0.85      # grey-scale base, dark background (axons bright)
    H, W = sub.shape
    def mask(i):
        m = np.zeros((ops['Ly'], ops['Lx']), bool); m[stat[i]['ypix'], stat[i]['xpix']] = True; return m[yr[0]:yr[1], xr[0]:xr[1]]
    acc = np.where(ic)[0]; rej = np.where(~ic)[0]
    for i in rej:
        m = mask(i); edge = m & ~binary_erosion(m); rgb[edge] = 0.55 * rgb[edge] + 0.45 * np.array([0.62, 0.64, 0.70])
    cols = {True: np.array([0.12, 0.45, 1.0]), False: np.array([1.0, 0.62, 0.0])}
    for k, i in enumerate(acc):
        m = mask(i); edge = m & ~binary_erosion(m); c = cols[bool(S.sig[k])] if not S.neuropil[k] else np.array([0.62, 0.64, 0.70])
        rgb[m] = 0.55 * rgb[m] + 0.45 * c; rgb[edge] = c
    ax.imshow(rgb, extent=(xr[0], xr[1], yr[1], yr[0]), interpolation='nearest', aspect='equal'); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    px = 1.0
    return dict(ops=ops, stat=stat, ic=ic, acc=acc, rej=rej, xr=xr, yr=yr)


def main():
    S1 = xc.Session(EX_1P, with_rejected=True, verbose=False); S5 = xc.Session(EX_5P, with_rejected=False, verbose=False)
    st = load_state(); extra = inventory_extra()
    fig = plt.figure(figsize=(W2, 262 * MM))
    gs = fig.add_gridspec(4, 1, height_ratios=[31, 42, 62, 118], hspace=0.20, left=0.10, right=0.955, top=0.985, bottom=0.015)

    # ================================================================ a  data streams
    ax = fig.add_subplot(gs[0]); ax.set_xlim(0, 100); ax.set_ylim(0, 37); ax.axis('off'); plabel(ax, 'a', dx=-35)
    def box(x, y, w, h, txt, fc='#F3F4F6', ec=FAINT, fs=6, bold=False, color=INK):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.0,rounding_size=1.2', fc=fc, ec=ec, lw=0.7)); ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=fs, fontweight='bold' if bold else 'normal', color=color, linespacing=1.25)
    def arrow(x0, y0, x1, y1): ax.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='-|>', lw=0.7, color=MUTED, shrinkA=0, shrinkB=0, mutation_scale=6))
    ax.text(0, 36.4, 'Axon signal', fontsize=7, fontweight='bold', va='top'); ax.text(37, 36.4, 'Behaviour (face camera + DLC, ball)', fontsize=7, fontweight='bold', va='top'); ax.text(70, 36.4, 'Events (projector, rig)', fontsize=7, fontweight='bold', va='top')
    box(0, 21, 31, 10, 'two-photon imaging of LC axons\nin the superior colliculus\n1 plane @ 30 Hz  or  5 planes @ 5 Hz', fc='#E8F0FC', ec=C_SIG)
    box(0, 10.5, 31, 7.5, 'suite2p ROIs; accepted by hand\nor rejected (control population)')
    box(0, 0, 31, 7.5, 'dF/F ; 3 ROI classes:\nsignal-like | noise-like | rejected')
    arrow(15.5, 21, 15.5, 18); arrow(15.5, 10.5, 15.5, 7.5)
    ax.text(37, 1.4, 'all streams are aligned on the 2P clock;\nwhich channels exist differs per session (see f)', fontsize=5.6, color=MUTED, va='bottom')
    chips = [('pupil area', C_BEH['pupil']), ('eyelid opening', C_BEH['eyelid']), ('pupil position', MUTED), ('face motion energy', C_BEH['face']), ('whisking index', C_BEH['whisk']), ('locomotion (ball)', C_BEH['ball'])]
    for k, (t, c) in enumerate(chips):
        r, cc = divmod(k, 2); x = 37 + cc * 16.2; y = 25 - r * 8.5
        ax.add_patch(FancyBboxPatch((x, y), 15.2, 6.4, boxstyle='round,pad=0.0,rounding_size=1.2', fc='white', ec=c, lw=1.0)); ax.text(x + 7.6, y + 3.2, t, ha='center', va='center', fontsize=5.8, color=INK)
    ev = ['full-field flash', 'moving bar', 'looming disc (expanding / receding)', 'grey screens (before / after loom)', 'square (passive or lick task)', 'spout: licks, pump / reward']
    for k, t in enumerate(ev):
        ax.add_patch(FancyBboxPatch((70, 27.0 - k * 5.0), 30, 4.1, boxstyle='round,pad=0.0,rounding_size=1.0', fc='white', ec=FAINT, lw=0.7)); ax.text(85, 29.05 - k * 5.0, t, ha='center', va='center', fontsize=5.8)

    # ================================================================ b, c  FOV + d ROI classes
    sub = gs[1].subgridspec(1, 3, width_ratios=[1.15, 1.15, 1.0], wspace=0.06)
    a1 = fig.add_subplot(sub[0]); a2 = fig.add_subplot(sub[1]); a3 = fig.add_subplot(sub[2])
    f1 = fov_panel(a1, EX_1P, S1); f2 = fov_panel(a2, EX_5P, S5)
    plabel(a1, 'b', dx=-14, dy=12); plabel(a2, 'c', dx=-10, dy=12)
    a1.set_title('1 plane, 30 Hz\n%s  (%d min)' % (short(EX_1P), round(S1.dur / 60)), fontsize=6.2, fontweight='normal', loc='left'); a2.set_title('5 planes summed, 5 Hz\n%s  (%d min)' % (short(EX_5P), round(S5.dur / 60)), fontsize=6.2, fontweight='normal', loc='left')
    a1.plot([f1['xr'][0] + 15, f1['xr'][0] + 15 + 100], [f1['yr'][1] - 12] * 2, color='white', lw=1.4); a1.text(f1['xr'][0] + 65, f1['yr'][1] - 17, '100 px', color='white', ha='center', va='bottom', fontsize=5.2)
    # example ROIs of the three classes on the 1-plane FOV
    rr = S1.rej; si = int(np.argmax([np.corrcoef(S1.dff[k][np.isfinite(S1.dff[k])][:5000], np.nanmean(S1.dff[S1.sig], 0)[np.isfinite(S1.dff[k])][:5000])[0, 1] for k in np.where(S1.sig)[0]]))
    k_sig = int(np.where(S1.sig)[0][si]); k_noise = int(np.where(S1.noise)[0][0]); rng = np.random.default_rng(3); k_rej = int(rng.choice(np.arange(len(rr['idx']))[rr['F0'] > np.percentile(rr['F0'], 50)]))
    cent = {k_sig: S1.X['centroid'][k_sig], k_noise: S1.X['centroid'][k_noise]}
    for lab, (y0, x0), col in (('1', cent[k_sig][:2], C_SIG), ('2', cent[k_noise][:2], C_NOISE)):
        a1.text(x0 + 10, y0 - 8, lab, color=col, fontsize=7, fontweight='bold', ha='left', va='bottom', bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.75))
    sr = f1['stat'][rr['idx'][k_rej]]; a1.text(np.median(sr['xpix']) + 10, np.median(sr['ypix']) - 8, '3', color=C_REJ, fontsize=7, fontweight='bold', ha='left', va='bottom', bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.75))
    a1.legend(handles=[Line2D([], [], color=C_SIG, lw=2, label='accepted, signal-like'), Line2D([], [], color=C_NOISE, lw=2, label='accepted, noise-like'), Line2D([], [], color=C_REJ, lw=1.2, label='rejected by suite2p / curation')], loc='upper left', bbox_to_anchor=(0.0, -0.02), frameon=False, fontsize=5.6, ncol=3, columnspacing=0.8, handlelength=1.2)
    # d  example traces
    plabel(a3, 'd', dx=-16, dy=2); a3.axis('off'); tmax = 180.0; sel = (S1.t - S1.t[0] >= 600) & (S1.t - S1.t[0] < 600 + tmax)
    tt = S1.t[sel] - S1.t[sel][0]; trs = {1: (xc.nsmooth(S1.dff[k_sig], 1 / xc.DT)[sel], C_SIG), 2: (xc.nsmooth(S1.dff[k_noise], 1 / xc.DT)[sel], C_NOISE), 3: (xc.nsmooth(rr['dff'][k_rej], 1 / xc.DT)[sel], C_REJ)}
    top = 0.0
    for k in (1, 2, 3):
        y, c = trs[k]; lo, hi = np.nanpercentile(y, [0.5, 99.8]); rngk = max(hi - lo, 0.6); base = top - rngk
        a3.plot(tt, y - lo + base, color=c, lw=0.8); a3.text(-8, base + 0.35 * rngk, str(k), color=c, fontsize=7, fontweight='bold', ha='right', va='center'); top = base - 0.35
    a3.set_xlim(-14, tmax + 22); a3.set_ylim(top - 0.5, 0.2); scalebar(a3, tmax - 30, top - 0.15, dx=30, lx='30 s'); a3.plot([tmax + 8] * 2, [-1.5, -0.5], color=INK, lw=0.9, clip_on=False); a3.text(tmax + 11, -1.0, '1 dF/F', fontsize=5.5, va='center', ha='left')
    a3.set_title('example dF/F (3 min)', fontsize=6.2, fontweight='normal', loc='left')

    # ================================================================ e  one full session
    S = S5; ge = gs[2].subgridspec(6, 1, height_ratios=[0.5, 7, 1.6, 1.6, 1.6, 1.6], hspace=0.12); tt = S.t - S.t[0]
    sig = np.where(S.sig)[0]; noi = np.where(~S.sig & ~S.neuropil)[0]; pop = np.nanmean(S.dff[sig], 0); okp = np.isfinite(pop)
    def order(ix):
        Z = S.dff[ix]; sc = []
        for r in Z: ok = np.isfinite(r) & okp; sc.append(np.corrcoef(r[ok], pop[ok])[0, 1] if ok.sum() > 50 else -1)
        return ix[np.argsort(sc)[::-1]]
    ix = np.concatenate([order(sig), order(noi)]); M = np.array([xc.nsmooth(S.dff[k], 1.0 / xc.DT) for k in ix])
    axs = fig.add_subplot(ge[0]); axs.set_xlim(tt[0], tt[-1]); axs.axis('off')
    on = np.asarray(S.X['on_t'], float) - S.t[0] if np.size(S.X['on_t']) else np.array([])
    for o in on: axs.plot([o, o], [0, 1], color=INK, lw=0.5)
    axs.text(tt[0] - 5, 0.5, 'n = %d' % len(on), fontsize=5.6, color=MUTED, ha='right', va='center')
    ar = fig.add_subplot(ge[1], sharex=axs); im = ar.imshow(M, aspect='auto', cmap='Blues', vmin=-0.05, vmax=np.nanpercentile(M[:len(sig)], 99.5), extent=(tt[0], tt[-1], len(ix), 0), interpolation='nearest')
    ar.axhline(len(sig), color=INK, lw=0.6); ar.set_ylabel('ROI'); ar.set_yticks([len(sig) / 2, len(sig) + len(noi) / 2]); ar.set_yticklabels(['signal-like\n(%d)' % len(sig), 'noise-like\n(%d)' % len(noi)], fontsize=5.6); ar.tick_params(axis='y', length=0); ar.tick_params(axis='x', labelbottom=False)
    plabel(axs, 'e', dx=-38, dy=-2); axs.text(1.0, 1.5, '%s: all accepted ROIs, dF/F (1 s smoothing); ticks above = square onsets' % short(EX_5P), transform=axs.transAxes, ha='right', va='bottom', fontsize=6.2, color=MUTED)
    cax = ar.inset_axes([1.005, 0.55, 0.008, 0.4]); cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=5, length=1.5); cb.set_label('dF/F', fontsize=5.5, labelpad=1)
    rows = [('face', 'face motion (log)', S.beh.get('face')), ('whisk', 'whisking (log)', S.beh.get('whisk')), ('pupil', 'pupil area', S.beh.get('pupil')), ('ball', 'ball speed', S.beh.get('loco'))]
    for k, (key, lab, y) in enumerate(rows):
        a = fig.add_subplot(ge[2 + k], sharex=axs)
        if y is not None: a.plot(tt, y, color=C_BEH[key], lw=0.6)
        a.set_ylabel(lab, fontsize=5.6, rotation=0, ha='right', va='center', labelpad=2); a.set_yticks([]); a.spines['left'].set_visible(False)
        if k < 3: a.tick_params(axis='x', labelbottom=False, bottom=False); a.spines['bottom'].set_visible(False)
    a.set_xlabel('time in session (s)'); a.set_xlim(tt[0], tt[-1])

    # ================================================================ f  inventory of the 20 analysed sessions
    SHORT_WHY = {'AnimalRB23_261027_1306': 'tracking channels disagree (face-whisk r 0.09)', 'AnimalRB23_261029_1704': 'tracking channels disagree (face-whisk r 0.13)', 'AnimalRB24_261023_1527': 'interrupted, 6 min',
                 'AnimalRB24_261104_1640': '7 min, short null', 'AnimalRB24_261106_1529': '11-06: little behavioural range', 'AnimalRB24_261106_1317': '11-06: little behavioural range', 'AnimalRB24_261106_1604': '11-06: little behavioural range'}
    names = ALL; n = len(names); gf = gs[3].subgridspec(2, 1, height_ratios=[n + 4.5, len(extra) * 0.8 + 4.5], hspace=0.36)
    sub = gf[0].subgridspec(1, 7, width_ratios=[1.0, 0.85, 1.0, 0.42, 0.2, 1.25, 1.15], wspace=0.08)
    y = np.arange(n); axl = fig.add_subplot(sub[0]); axd = fig.add_subplot(sub[1], sharey=axl); axr = fig.add_subplot(sub[2], sharey=axl); axv = fig.add_subplot(sub[3], sharey=axl); axb = fig.add_subplot(sub[4], sharey=axl)
    axs_ = fig.add_subplot(sub[5], sharey=axl); axt2 = fig.add_subplot(sub[6], sharey=axl)
    for a in (axl, axd, axr, axv, axb, axs_, axt2): a.set_ylim(n - 0.5, -3.3); a.patch.set_visible(False)
    axl.set_xlim(0, 1); axl.axis('off'); plabel(axl, 'f', dx=-24, dy=-6)
    for i, nm in enumerate(names):
        fl = is_flagged(nm); c = MUTED if fl else INK
        axl.plot(0.035, i, MARK[animal(nm)], mfc='white' if fl else INK, mec=MUTED if fl else INK, mew=0.8, ms=4.2, clip_on=False)
        axl.text(0.12, i, short(nm) + (' †' if fl else ''), ha='left', va='center', fontsize=6, color=c)
        for a in (axl, axd, axr, axv, axb, axs_, axt2):
            if i % 2 == 0: a.axhspan(i - 0.5, i + 0.5, color='#F5F5F6', lw=0, zorder=0)
    axl.text(0.0, -2.95, 'analysed (aggregated) sessions, n = %d' % n, fontsize=6.4, fontweight='bold', ha='left', va='center')
    axl.text(0.0, -2.05, '† flagged: hollow, never pooled', fontsize=5.6, color=MUTED, ha='left', va='center')
    mins = [float(QC[nm]['minutes']) for nm in names]; axd.barh(y, mins, color=MUTED, height=0.62); axd.set_xlim(0, 125); axd.set_xticks([0, 50, 100]); axd.tick_params(axis='y', left=False, labelleft=False)
    axd.text(0, -1.0, 'minutes', fontsize=6.2, fontweight='bold', ha='left', va='bottom')
    for i, m in enumerate(mins): axd.text(m + 2, i, '%.0f' % m, fontsize=4.8, color=MUTED, va='center')
    ns = np.array([qcf(nm, 'signal_like') for nm in names]); nn = np.array([qcf(nm, 'noise_like') for nm in names]); npl = np.array([qcf(nm, 'neuropil_limited') for nm in names])
    axr.barh(y, ns, color=C_SIG, height=0.62); axr.barh(y, nn, left=ns, color=C_NOISE, height=0.62); axr.barh(y, npl, left=ns + nn, color=C_REJ, height=0.62); axr.set_xlim(0, 300); axr.set_xticks([0, 100, 200]); axr.tick_params(axis='y', left=False, labelleft=False)
    axr.text(0, -1.0, 'accepted ROIs', fontsize=6.2, fontweight='bold', ha='left', va='bottom')
    axr.legend(handles=[Rectangle((0, 0), 1, 1, color=C_SIG, label='signal-like'), Rectangle((0, 0), 1, 1, color=C_NOISE, label='noise-like'), Rectangle((0, 0), 1, 1, color=C_REJ, label='Fc ≤ 0')], frameon=False, fontsize=5.2, loc='lower right', bbox_to_anchor=(1.0, 0.0), handlelength=0.8, borderpad=0.2, labelspacing=0.25)
    chs = ['pupil', 'face', 'whisk', 'eyelid']
    for j, ch in enumerate(chs):
        for i, nm in enumerate(names):
            has = ch in st[nm]['behaviours']; axv.plot(j, i, 'o', ms=3.2, mfc=C_BEH[ch] if has else 'none', mec=C_BEH[ch] if has else FAINT, mew=0.6)
    axv.set_xlim(-0.6, 3.6); axv.set_xticks(range(4)); axv.set_xticklabels(['pupil', 'face', 'whisk', 'lid'], rotation=90, fontsize=5.4); axv.tick_params(axis='y', left=False, labelleft=False); axv.text(-0.5, -1.0, 'video', fontsize=6.2, fontweight='bold', ha='left', va='bottom')
    for i, nm in enumerate(names):
        used = QC[nm]['ball'].startswith('used'); axb.plot(0, i, 'o', ms=3.2, mfc=C_BEH['ball'] if used else 'none', mec=C_BEH['ball'] if used else FAINT, mew=0.6)
    axb.set_xlim(-0.6, 0.6); axb.set_xticks([0]); axb.set_xticklabels(['ball ok'], rotation=90, fontsize=5.4); axb.tick_params(axis='y', left=False, labelleft=False)
    cols = [('flash', 'flash', lambda d, nm: d.get('flash', 0) + d.get('fff', 0)), ('bar', 'bar', lambda d, nm: d.get('bar', 0)), ('loom_expand', 'loom exp.', lambda d, nm: d.get('loom_expand', 0)), ('loom_recede', 'loom rec.', lambda d, nm: d.get('loom_recede', 0)),
            ('square', 'square', lambda d, nm: 0 if nm in TASK_SESSIONS else d.get('square', 0)), ('task', 'lick task', lambda d, nm: d.get('square', 0) if nm in TASK_SESSIONS else 0), ('grey', 'grey', lambda d, nm: 1 if nm in GREY_SESSIONS else 0)]
    DX = 1.7
    for j, (k, lab, fn) in enumerate(cols):
        for i, nm in enumerate(names):
            v = fn(parse_stim(nm), nm)
            if v > 0: axs_.plot(j * DX, i, 's', ms=3.4, mfc=INK if not is_flagged(nm) else 'white', mec=INK if not is_flagged(nm) else MUTED, mew=0.6)
            if v > 0 and k != 'grey': axs_.text(j * DX + 0.22, i, str(v), fontsize=4.4, color=MUTED, va='center', ha='left')
    axs_.set_xlim(-0.7, DX * (len(cols) - 1) + 0.9); axs_.set_xticks([j * DX for j in range(len(cols))]); axs_.set_xticklabels([c[1] for c in cols], rotation=90, fontsize=5.4); axs_.tick_params(axis='y', left=False, labelleft=False)
    axs_.text(-0.6, -1.0, 'stimuli (n trials)', fontsize=6.2, fontweight='bold', ha='left', va='bottom')
    axt2.set_xlim(0, 1); axt2.axis('off'); axt2.text(0, -1.0, 'why flagged', fontsize=6.2, fontweight='bold', ha='left', va='bottom')
    for i, nm in enumerate(names):
        if is_flagged(nm): axt2.text(0, i, SHORT_WHY[nm], fontsize=4.9, color=MUTED, va='center', ha='left')
    for xx, mk, lab in ((0.035, 'o', 'RB24'), (0.34, 's', 'RB23'), (0.66, '^', 'RB22')): axl.plot(xx, -1.15, mk, color=INK, ms=4, clip_on=False); axl.text(xx + 0.05, -1.15, lab, fontsize=5.6, va='center')

    # ================================================================ g  recorded but not aggregated
    ex = sorted(extra, key=lambda r: (animal(r['session'])[::-1], r['session'])); m = len(ex); axg = gf[1].subgridspec(1, 3, width_ratios=[1.0, 1.0, 1.9], wspace=0.07)
    ga = fig.add_subplot(axg[0]); gb = fig.add_subplot(axg[1], sharey=ga); gc = fig.add_subplot(axg[2], sharey=ga)
    for a in (ga, gb, gc): a.set_ylim(m - 0.5, -3.0); a.patch.set_visible(False)
    ga.axis('off'); ga.set_xlim(0, 1); plabel(ga, 'g', dx=-24, dy=-6)
    ga.text(0.0, -2.7, 'recorded, NOT aggregated (n = %d)' % m, fontsize=6.4, fontweight='bold', ha='left', va='center'); ga.text(0.0, -1.7, 'what exists on disk', fontsize=5.6, color=MUTED, ha='left', va='center')
    for i, r in enumerate(ex):
        ga.plot(0.035, i, MARK[animal(r['session'])], mfc=INK, mec='white', mew=0.3, ms=4.0, clip_on=False); ga.text(0.12, i, short(r['session']), ha='left', va='center', fontsize=6)
        for a in (ga, gb, gc):
            if i % 2 == 0: a.axhspan(i - 0.5, i + 0.5, color='#F5F5F6', lw=0, zorder=0)
    gb.barh(np.arange(m), [r['tifs'] for r in ex], color=MUTED, height=0.62); gb.set_xlim(0, 130); gb.set_xticks([0, 50, 100]); gb.tick_params(axis='y', left=False, labelleft=False); gb.text(0, -0.9, 'TIFF files', fontsize=6.2, fontweight='bold', ha='left', va='bottom')
    for i, r in enumerate(ex):
        if r['suite2p'] and r['accepted'] >= 0: gc.barh(i, r['accepted'], color=C_REJ, height=0.62); gc.text(r['accepted'] + 1.5, i, str(r['accepted']), fontsize=5.2, color=MUTED, va='center')
        else: gc.text(1.0, i, 'no suite2p result', fontsize=5.2, color=MUTED, va='center')
    gc.set_xlim(0, 70); gc.tick_params(axis='y', left=False, labelleft=False); gc.text(0, -0.9, 'suite2p ROIs with iscell > 0 (classification not curated)', fontsize=6.2, fontweight='bold', ha='left', va='bottom')
    savefig(fig, 'Fig1_overview')


if __name__ == '__main__':
    main()
