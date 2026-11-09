"""Figure 7 - what did NOT work or is ambiguous.
a  RB24 11-06 1529: 4 min of the session type that shows no coupling (compare Fig. 2a)   b  behavioural range (SD of log face motion) per recording day   c  coupling vs behavioural range across sessions
d  range-restriction control: window-wise coupling vs the behavioural SD of the window (all sessions)   e  do the behavioural channels agree with each other?
f  how much more than random rejected ROIs do behaviour models explain?   g  coupling of the axon mean to brain shift (lateral motion) vs to face motion
usage: python fig7_negative.py -> plots/summary_figures/Fig7_negative.pdf/.png ; ledger entries
"""
import datetime as dt
import numpy as np
from fig_common import *
import xs_common as xc
from scipy.stats import spearmanr


def day_of(n): d = n.split('_')[1]; return (dt.date(2000 + int(d[:2]), int(d[2:4]), int(d[4:6])) - dt.date(2026, 10, 17)).days


def windowed(st, n, win_s=180.0):
    r = st[n]; P = np.asarray(r['P'], float); f = np.asarray(r['beh_grid']['face'], float)
    if len(P) != len(f): return []
    bp, bf = xc.bandpass(P), xc.bandpass(f); w = int(win_s / xc.DT); out = []
    for a in range(0, len(P) - w + 1, w):
        x, y, raw = bp[a:a + w], bf[a:a + w], f[a:a + w]; ok = np.isfinite(x) & np.isfinite(y)
        if ok.mean() < 0.9 or np.nanstd(raw) == 0: continue
        out.append((float(np.nanstd(raw)), float(np.corrcoef(x[ok], y[ok])[0, 1])))
    return out


def main():
    st = load_state(); fig = plt.figure(figsize=(W2, 150 * MM))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.36, left=0.085, right=0.985, top=0.955, bottom=0.125)
    top = gs[0].subgridspec(1, 3, wspace=0.34, width_ratios=[1.9, 1.0, 1.0])
    # ------------------------------------------------ a  raw excerpt 11-06 1529
    nm = 'AnimalRB24_261106_1529'; r = st[nm]; t = np.asarray(r['t'], float); P = xc.nsmooth(np.asarray(r['P'], float), 1.0 / xc.DT); bg = r['beh_grid']; sl = (t >= t[0] + 300) & (t < t[0] + 540); tt = t[sl] - t[sl][0]
    ga = top[0].subgridspec(4, 1, height_ratios=[1.6, 1, 1, 1], hspace=0.14)
    for k, (y, col, yl) in enumerate(((P[sl], C_SIG, 'axon mean\ndF/F'), (np.asarray(bg['face'])[sl], C_BEH['face'], 'face motion\n(log)'), (np.asarray(bg['whisk'])[sl], C_BEH['whisk'], 'whisking\n(log)'), (np.asarray(bg['pupil'])[sl], C_BEH['pupil'], 'pupil'))):
        a = fig.add_subplot(ga[k], sharex=fig.axes[-1] if k else None); a.plot(tt, y, color=col, lw=1.0 if k == 0 else 0.6); a.set_yticks([]); a.spines['left'].set_visible(False); a.set_ylabel(yl, fontsize=5.4, rotation=0, ha='right', va='center', labelpad=3)
        if k < 3: a.tick_params(axis='x', labelbottom=False, bottom=False); a.spines['bottom'].set_visible(False)
        if k == 0: plabel(a, 'a', dx=-34, dy=10); a.set_title('%s (flagged): 4 min. Face-motion SD %.2f vs %.2f in 1630 (Fig. 2a); session r(axon, face) = %.2f' % (short(nm), qcf(nm, 'beh_face_sd_log'), qcf('AnimalRB24_261021_1630', 'beh_face_sd_log'), qcf(nm, 'r_face')), fontsize=5.2, loc='left', fontweight='normal', pad=3)
    a.set_xlabel('time (s)'); a.set_xlim(0, tt[-1])
    # ------------------------------------------------ b  behavioural range per day
    ab = fig.add_subplot(top[1]); rng = np.random.default_rng(8); xs_, ys_, rs_ = [], [], []
    for n in ALL:
        sd = qcf(n, 'beh_face_sd_log'); x = day_of(n) + rng.uniform(-0.3, 0.3); sdot(ab, x, sd, n, color=INK, hollow_color=MUTED, ms=3.8)
    ab.axvspan(19.5, 20.5, color=C_NOISE, alpha=0.12, lw=0); ab.set_xticks([4, 10, 12, 13, 18, 19, 20]); ab.set_xticklabels(['10-21', '10-27', '10-29', '10-30', '11-04', '11-05', '11-06'], fontsize=5.0, rotation=90)
    ab.set_ylabel('SD of log face motion\n(behavioural range)'); ab.set_ylim(0, 0.85); plabel(ab, 'b', dx=-24, dy=10)
    # ------------------------------------------------ c  coupling vs range
    ac = fig.add_subplot(top[2]); xs_, ys_ = [], []
    for n in ALL:
        sd = qcf(n, 'beh_face_sd_log'); pr = st[n]['xc'].get('face', {}).get('peak_r', np.nan)
        if not np.isfinite(sd) or not np.isfinite(pr): continue
        sdot(ac, sd, pr, n, color=C_SIG, hollow_color=MUTED, ms=3.8); xs_.append(sd); ys_.append(pr)
    rho, p = spearmanr(xs_, ys_); ac.set_xlabel('SD of log face motion'); ac.set_ylabel('peak r (axon mean vs face motion)'); ac.axhline(0, color=MUTED, lw=0.4); ac.set_ylim(-0.4, 1.0); plabel(ac, 'c', dx=-24, dy=10)
    ac.text(0.03, 0.97, 'all %d sessions: ρ = %+.2f, %s' % (len(xs_), rho, fmt_p(p)), transform=ac.transAxes, fontsize=5.0, va='top')
    ledger_add('Fig7', 'r_vs_range', 'across sessions: peak r(axon, face) vs behavioural range (SD of log face motion), Spearman', rho, np.nan, np.nan, len(xs_), 3, 'ambiguous', note='%s; 11-06 sessions are at the low end of both axes - range restriction or real absence cannot be separated by this panel alone' % fmt_p(p), unit='rho')

    # ------------------------------------------------ bottom row
    bot = gs[1].subgridspec(1, 4, wspace=0.50, width_ratios=[1.25, 0.9, 0.9, 0.9]); ad = fig.add_subplot(bot[0]); ae = fig.add_subplot(bot[1]); af = fig.add_subplot(bot[2]); ag = fig.add_subplot(bot[3])
    # d  windowed coupling vs window SD
    allw = []
    for n in ALL:
        w = windowed(st, n); fl106 = n.split('_')[1] == '261106'
        for sd, rr in w:
            ad.plot(sd, rr, 'o', mfc=C_NOISE if fl106 else C_SIG, mec='none', ms=2.8, alpha=0.7 if fl106 else 0.45, zorder=3 if fl106 else 2)
            if not is_flagged(n): allw.append((sd, rr))
    allw = np.array(allw); bins = np.percentile(allw[:, 0], np.linspace(0, 100, 7)); mids = []; meds = []
    for i in range(6):
        m = (allw[:, 0] >= bins[i]) & (allw[:, 0] <= bins[i + 1]); mids.append(np.median(allw[m, 0])); meds.append(np.median(allw[m, 1]))
    ad.plot(mids, meds, color=INK, lw=1.4, marker='o', ms=3, zorder=5); ad.axhline(0, color=MUTED, lw=0.4); ad.set_xlabel('SD of log face motion within the 3 min window'); ad.set_ylabel('r(axon mean, face), 3 min window'); ad.set_ylim(-0.7, 1.0); plabel(ad, 'd', dx=-24, dy=10)
    rho_w, p_w = spearmanr(allw[:, 0], allw[:, 1]); ad.text(0.03, 0.03, 'unflagged windows: ρ = %+.2f (n = %d)\nblack: median per range bin; orange: 11-06 windows' % (rho_w, len(allw)), transform=ad.transAxes, fontsize=4.8, va='bottom')
    ledger_add('Fig7', 'window_r_vs_range', 'within sessions: 3 min window coupling vs window behavioural range (Spearman, unflagged windows)', rho_w, np.nan, np.nan, len(allw), 3, 'ambiguous', note='windows are not independent; tests whether low-variance stretches of coupled sessions look like 11-06', unit='rho')
    # e  tracking consistency
    for n in ALL:
        fw = qcf(n, 'beh_face_whisk_r'); pr = st[n]['xc'].get('face', {}).get('peak_r', np.nan)
        if np.isfinite(fw) and np.isfinite(pr):
            sdot(ae, fw, pr, n, color=C_SIG, hollow_color=MUTED, ms=3.8)
            if fw < 0.3: ae.annotate(short(n).replace('RB23 ', ''), (fw, pr), xytext=(3, 3), textcoords='offset points', fontsize=4.4, color=MUTED)
    ae.axhline(0, color=MUTED, lw=0.4); ae.set_xlabel('face motion vs whisking, r\n(do two tracking channels agree?)'); ae.set_ylabel('peak r (axon mean, face motion)'); ae.set_xlim(0, 1.0); ae.set_ylim(-0.4, 1.0); plabel(ae, 'e', dx=-24, dy=10)
    # f  ridge R2 axons minus random rejected sets
    vals = {}
    for i, n in enumerate(ALL):
        rr = st[n].get('ridge')
        if not rr or not rr.get('rej'): continue
        d = rr['sig']['full'] - np.nanmean(rr['rej']); sdot(af, d, len(ALL) - i, n, color=C_SIG, hollow_color=MUTED, ms=3.6)
        if not is_flagged(n): vals[n] = d
    h = pooled_unflagged(vals, n_boot=4000); af.axvline(0, color=MUTED, lw=0.5); af.axvspan(h['lo'], h['hi'], color=C_SIG, alpha=0.15, lw=0); af.axvline(h['est'], color=INK, lw=1.2)
    af.set_yticks([]); af.spines['left'].set_visible(False); af.set_xlabel('R² axons − R² random rejected sets' + chr(10) + 'unflagged %+.2f [%.2f, %.2f]; > 0: axons better' % (h['est'], h['lo'], h['hi'])); af.set_ylabel('sessions'); plabel(af, 'f', dx=-24, dy=10)
    ledger_add('Fig7', 'ridge_R2_minus_rejected', 'held-out R2 of behaviour model: signal-like axons minus random rejected ROI sets of equal size', h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'ambiguous' if h['lo'] < 0.05 else 'supported', note='a large share of the variance is shared with ROIs that were rejected: common-mode (arousal-locked motion / neuropil)', unit='R2')
    # g  shift coupling
    for n in ALL:
        s_ = st[n]['xc'].get('shift', {}).get('peak_r', np.nan); fc = st[n]['xc'].get('face', {}).get('peak_r', np.nan)
        if np.isfinite(s_) and np.isfinite(fc): sdot(ag, fc, s_, n, color=C_BEH['shift'], hollow_color=MUTED, ms=3.8)
    ag.plot([-0.4, 1.0], [-0.4, 1.0], color=MUTED, lw=0.4, ls=':'); ag.axhline(0, color=MUTED, lw=0.4); ag.set_xlabel('peak r (axon, face motion)'); ag.set_ylabel('peak r (axon mean, brain shift)'); ag.set_xlim(-0.4, 1.0); ag.set_ylim(-0.4, 1.0); plabel(ag, 'g', dx=-24, dy=10)
    vs = {n: st[n]['xc']['shift']['peak_r'] for n in UNFLAGGED if 'shift' in st[n]['xc']}; hs = pooled_unflagged(vs, n_boot=4000)
    ag.text(0.03, 0.97, 'unflagged: shift r %.2f [%.2f, %.2f]' % (hs['est'], hs['lo'], hs['hi']), transform=ag.transAxes, fontsize=4.8, va='top')
    ledger_add('Fig7', 'shift_coupling', 'peak r of the axon mean with the rigid brain shift (registration offsets, lateral motion)', hs['est'], hs['lo'], hs['hi'], hs['n_sessions'], hs['n_animals'], 'ambiguous', note='arousal moves the brain too; without a reference channel lateral / axial motion cannot be separated from signal', unit='r')
    fig.text(0.085, 0.004, 'Orange = 11-06 windows. Windows: 3 min, behaviour and axons band-passed 1-30 s, at least 90% valid frames. Marker shape = animal (circle RB24, square RB23, triangle RB22);\n'
             'filled = unflagged, hollow = flagged (not pooled); unflagged pooled values = mean of per-animal means with 95% hierarchical bootstrap CI. c, d: is the lack of coupling on 11-06 a lack of behavioural range?\n'
             'e-g: where the evidence for a pure arousal signal is weaker (tracking agreement between channels, common-mode with rejected ROIs, brain shift).', fontsize=4.7, color=MUTED)
    savefig(fig, 'Fig7_negative')


if __name__ == '__main__':
    main()
