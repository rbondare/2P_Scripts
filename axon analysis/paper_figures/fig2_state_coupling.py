"""Figure 2 - the axons follow face motion / whisking / pupil / ball.
a, b raw 4 min excerpts (1-plane session 1630 without ball; 1744 with ball)   c cross-correlation (averages over lags), d spectral coherence
e peak r per session x behaviour, f single-ROI sign test, g held-out ridge R2 vs controls, h alone vs unique contribution of each behaviour
usage: python fig2_state_coupling.py  -> plots/summary_figures/Fig2_state_coupling.pdf/.png ; ledger entries
"""
import os, pickle, numpy as np
from fig_common import *
import xs_common as xc
from event_triggered import load_events

EX_A, EX_B = 'AnimalRB24_261021_1630', 'AnimalRB24_261030_1744'
BEHS = ['face', 'whisk', 'pupil', 'eyelid', 'loco']
LAB = {'face': 'face motion', 'whisk': 'whisking', 'pupil': 'pupil area', 'eyelid': 'eyelid', 'loco': 'ball speed'}
CK = {'face': 'face', 'whisk': 'whisk', 'pupil': 'pupil', 'eyelid': 'eyelid', 'loco': 'ball'}       # colour keys


def excerpt(name, win_s=240.0):
    f = os.path.join(CACHE, 'excerpt_%s.pkl' % name)
    if os.path.exists(f): return pickle.load(open(f, 'rb'))
    S = xc.Session(name, with_rejected=True, verbose=False); P = S.population('sig'); face = S.beh['face']; Psm = xc.nsmooth(P, 1.0 / xc.DT); w = int(win_s / xc.DT); best, bi = -1, 0
    ball = S.beh.get('loco') if S.loco_ok else None
    for i in range(0, S.nb - w, int(10 / xc.DT)):
        seg, ff = Psm[i:i + w], face[i:i + w]
        if np.isfinite(seg).mean() > 0.97 and np.isfinite(ff).mean() > 0.97 and (ball is None or np.isfinite(ball[i:i + w]).mean() > 0.97):
            sc = np.nanstd(seg) * np.nanstd(ff)
            if sc > best: best, bi = sc, i
    sl = slice(bi, bi + w); sig = np.where(S.sig)[0]; pm = Psm[sl]; rows = []
    for k in sig:
        r = xc.nsmooth(S.dff[k], 1.0 / xc.DT)[sl]; ok = np.isfinite(r) & np.isfinite(pm); rows.append(np.corrcoef(r[ok], pm[ok])[0, 1] if ok.sum() > 50 else -1)
    M = np.array([xc.nsmooth(S.dff[k], 1.0 / xc.DT)[sl] for k in sig[np.argsort(rows)[::-1]]])
    rng = np.random.default_rng(2); R = S.rej['dff']; rej = np.nanmean(R[rng.choice(len(R), len(sig), replace=False)], 0)
    d = dict(t=S.t[sl] - S.t[sl][0], M=M, P=Psm[sl], rej=xc.nsmooth(rej, 1.0 / xc.DT)[sl], face=face[sl], whisk=S.beh['whisk'][sl], pupil=S.beh['pupil'][sl], ball=None if ball is None else ball[sl], n_sig=int(len(sig)), n_rej=int(len(R)),
             r_face=float(np.corrcoef(xc.bandpass(P)[np.isfinite(P) & np.isfinite(face)], xc.bandpass(face)[np.isfinite(P) & np.isfinite(face)])[0, 1]))
    pickle.dump(d, open(f, 'wb')); return d


def raw_panel(fig, spec, name, letter, title):
    d = excerpt(name); t = d['t']; gs = spec.subgridspec(6, 1, height_ratios=[2.4, 1.9, 1.0, 1.0, 1.0, 1.0], hspace=0.12)
    a0 = fig.add_subplot(gs[0]); a0.imshow(d['M'], aspect='auto', cmap='Blues', vmin=-0.05, vmax=np.nanpercentile(d['M'], 99.5), extent=(t[0], t[-1], len(d['M']), 0), interpolation='nearest')
    a0.set_yticks([]); a0.tick_params(axis='x', labelbottom=False, bottom=False); a0.spines['bottom'].set_visible(False); a0.spines['left'].set_visible(False)
    a0.set_ylabel('%d axon\nROIs' % d['n_sig'], fontsize=5.6, rotation=0, ha='right', va='center', labelpad=3); plabel(a0, letter, dx=-30, dy=8); a0.set_title(title + ' - whole-session r(axons, face motion) = %.2f' % d['r_face'], fontsize=6.0, loc='left', fontweight='normal', pad=2)
    a1 = fig.add_subplot(gs[1], sharex=a0); a1.plot(t, d['rej'], color=C_REJ, lw=0.7); a1.plot(t, d['P'], color=C_SIG, lw=1.0)
    a1.set_ylabel('mean\ndF/F', fontsize=5.6, rotation=0, ha='right', va='center', labelpad=3); a1.tick_params(axis='x', labelbottom=False, bottom=False); a1.spines['bottom'].set_visible(False)
    a1.text(1.0, 1.02, 'rejected ROIs (same n)', color=C_REJ, fontsize=5.2, transform=a1.transAxes, ha='right', va='bottom'); a1.text(0.72, 1.02, 'signal-like axons', color=C_SIG, fontsize=5.2, transform=a1.transAxes, ha='right', va='bottom')
    for k, (key, lab) in enumerate((('face', 'face motion\n(log)'), ('whisk', 'whisking\n(log)'), ('pupil', 'pupil\narea'), ('ball', 'ball speed\n(log)'))):
        a = fig.add_subplot(gs[2 + k], sharex=a0); y = d.get(key)
        if y is None: a.text(0.5, 0.5, 'ball not usable in this session', transform=a.transAxes, ha='center', va='center', fontsize=5.2, color=MUTED); a.set_yticks([])
        else: a.plot(t, y, color=C_BEH[key], lw=0.6); a.set_yticks([])
        a.set_ylabel(lab, fontsize=5.4, rotation=0, ha='right', va='center', labelpad=3); a.spines['left'].set_visible(False)
        if k < 3: a.tick_params(axis='x', labelbottom=False, bottom=False); a.spines['bottom'].set_visible(False)
    a.set_xlabel('time (s)'); a.set_xlim(t[0], t[-1])
    return d


def main():
    st = load_state(); ev = {n: load_events(n) for n in ALL}
    fig = plt.figure(figsize=(W2, 222 * MM))
    gs = fig.add_gridspec(3, 1, height_ratios=[62, 40, 62], hspace=0.40, left=0.105, right=0.985, top=0.975, bottom=0.085)
    top = gs[0].subgridspec(1, 2, wspace=0.30)
    raw_panel(fig, top[0], EX_A, 'a', '%s: 4 min excerpt (no usable ball)' % short(EX_A)); raw_panel(fig, top[1], EX_B, 'b', '%s: 4 min excerpt (with ball)' % short(EX_B))

    # ------------------------------------------------ c  cross-correlation curves, d spectral coherence
    mid = gs[1].subgridspec(1, 5, wspace=0.42, width_ratios=[1, 1, 1, 1, 1.15]); ac = []
    for j, b in enumerate(('face', 'whisk', 'pupil', 'loco')):
        a = fig.add_subplot(mid[j], sharey=ac[0] if ac else None); ac.append(a)
        curves = []; rejs = []; noi = []
        for n in ALL:
            d = st[n]['xc'].get(b)
            if d is None: continue
            fl = is_flagged(n); lg = d['lags_s']
            a.plot(lg, d['r'], color=FAINT if fl else C_SIG, lw=0.5 if fl else 0.6, alpha=1 if fl else 0.45, zorder=2, ls='--' if fl else '-')
            if not fl:
                curves.append(d['r'])
                if d.get('r_rej_curves') is not None: rejs.append(np.asarray(d['r_rej_curves']).mean(0)); rejs_all = np.asarray(d['r_rej_curves'])
                if d.get('r_noise') is not None: noi.append(d['r_noise'])
        lg = st[ALL[7]]['xc']['face']['lags_s']
        if curves: a.plot(lg, np.mean(curves, 0), color=C_SIG, lw=1.6, zorder=4)
        if rejs: m = np.mean(rejs, 0); s = np.std(rejs, 0); a.fill_between(lg, m - 2 * s, m + 2 * s, color=C_REJ, alpha=0.35, lw=0, zorder=1); a.plot(lg, m, color=C_REJ, lw=1.0, zorder=3)
        if noi: a.plot(lg, np.mean(noi, 0), color=C_NOISE, lw=1.0, zorder=3)
        a.axvline(0, color=MUTED, lw=0.4, ls=':'); a.axhline(0, color=MUTED, lw=0.4); a.set_xlim(-10, 10); a.set_xlabel('lag (s)  [behaviour later →]')
        a.set_title('%s  (n = %d)' % (LAB[b].split(' (')[0], len(curves)), fontsize=6.2, loc='left', fontweight='bold', color=C_BEH[CK[b]])
        if j == 0: a.set_ylabel('r  (axon mean vs behaviour)'); plabel(a, 'c', dx=-22, dy=6)
        else: a.tick_params(labelleft=False)
    ac[0].set_ylim(-0.4, 0.95)
    # d  coherence
    ad = fig.add_subplot(mid[4]); plabel(ad, 'd', dx=-22, dy=6)
    for b in ('face', 'whisk', 'pupil', 'loco'):
        C = []; fc = None; Cn = []
        for n in UNFLAGGED:
            e = ev[n]; d = e['coh'].get(b) if e else None
            if d is None: continue
            fc = d['fc']; C.append(d['c']); Cn.append(d['c_null']) if d['c_null'] is not None else None
        if C:
            C = np.array(C); ad.plot(fc, np.nanmean(C, 0), color=C_BEH[CK[b]], lw=1.3, label=LAB[b].split(' (')[0]); ad.fill_between(fc, np.nanpercentile(C, 25, 0), np.nanpercentile(C, 75, 0), color=C_BEH[CK[b]], alpha=0.12, lw=0)
    allnull = np.array([x for b in ('face', 'whisk', 'pupil') for n in UNFLAGGED if ev[n] and ev[n]['coh'].get(b) and ev[n]['coh'][b]['c_null'] is not None for x in [ev[n]['coh'][b]['c_null']]])
    ad.fill_between(fc, 0, np.nanmean(allnull, 0), color=C_REJ, alpha=0.4, lw=0, label='shifted behaviour'); ad.set_xscale('log'); ad.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ad.set_xlim(0.01, 2); ad.set_ylim(0, 1); ad.set_xlabel('frequency (Hz)'); ad.set_ylabel('coherence', labelpad=1)
    ad.legend(frameon=False, fontsize=4.9, loc='upper right', handlelength=1.2, borderpad=0.1, labelspacing=0.2); ad.set_title('spectral coherence', fontsize=6.2, loc='left', fontweight='bold')

    # ------------------------------------------------ e..h statistics
    bot = gs[2].subgridspec(1, 4, wspace=0.50, width_ratios=[1.35, 1.0, 1.05, 1.2]); ae = fig.add_subplot(bot[0]); af = fig.add_subplot(bot[1]); ag = fig.add_subplot(bot[2]); ah = fig.add_subplot(bot[3])
    rng = np.random.default_rng(1); ledger = {}
    for j, b in enumerate(BEHS):
        vals = {}; sig_n = 0; tot = 0
        for n in ALL:
            d = st[n]['xc'].get(b)
            if d is None: continue
            x = j + rng.uniform(-0.22, 0.22); sdot(ae, x, d['peak_r'], n, color=C_BEH[CK[b]], hollow_color=MUTED, ms=3.8)
            if not is_flagged(n): vals[n] = d['peak_r']; tot += 1; sig_n += int(d['p'] < 0.01 and d['peak_r'] > 0)
        h = pooled_unflagged(vals, n_boot=4000)
        if h['n_sessions']:
            ae.plot([j - 0.33, j + 0.33], [h['est']] * 2, color=INK, lw=1.4, zorder=5); ae.plot([j, j], [h['lo'], h['hi']], color=INK, lw=0.9, zorder=5)
            ae.text(j, 1.0, '%d/%d' % (sig_n, tot), ha='center', va='bottom', fontsize=5, color=MUTED)
            ledger[b] = h; ledger_add('Fig2', 'peak_r_' + b, 'peak r axon mean vs %s (±10 s)' % LAB[b], h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'supported' if h['lo'] > 0.2 else 'ambiguous', note='sessions with p<0.01 vs shifted data: %d/%d' % (sig_n, tot), unit='r')
    ae.axhline(0, color=MUTED, lw=0.4); ae.set_xticks(range(len(BEHS))); ae.set_xticklabels([LAB[b].replace(' ', '\n') for b in BEHS], fontsize=5.4); ae.set_ylim(-0.4, 1.0); ae.set_xlim(-0.6, len(BEHS) - 0.4); ae.set_ylabel('peak r within ±10 s'); plabel(ae, 'e', dx=-22, dy=6)
    # f  net sign of single-ROI couplings (face, whisk, pupil)
    for j, b in enumerate(('face', 'whisk', 'pupil')):
        vs = {}; vr = {}
        for n in ALL:
            r = st[n]['roi'].get(b)
            if not r or 'sig' not in r or 'rej' not in r: continue
            s_ = 2 * np.mean(r['sig']['r'] > 0) - 1; r_ = 2 * np.mean(r['rej']['r'] > 0) - 1
            xx = j * 1.0
            af.plot([xx - 0.17, xx + 0.17], [s_, r_], color=FAINT, lw=0.5, zorder=1); sdot(af, xx - 0.17, s_, n, color=C_SIG, hollow_color=C_SIG, ms=3.3); sdot(af, xx + 0.17, r_, n, color=C_REJ, hollow_color=C_REJ, ms=3.3)
            if not is_flagged(n): vs[n] = s_; vr[n] = r_
        hs = pooled_unflagged(vs, n_boot=3000); hr = pooled_unflagged(vr, n_boot=3000)
        af.plot([j - 0.3, j - 0.04], [hs['est']] * 2, color=INK, lw=1.4, zorder=5); af.plot([j + 0.04, j + 0.3], [hr['est']] * 2, color=INK, lw=1.4, zorder=5)
        ledger_add('Fig2', 'roi_sign_' + b, 'net sign of single-ROI coupling to %s: signal-like ROIs' % LAB[b], hs['est'], hs['lo'], hs['hi'], hs['n_sessions'], hs['n_animals'], 'supported', note='rejected ROIs: %.2f' % hr['est'], unit='net sign (+1 = all positive)')
    af.axhline(0, color=MUTED, lw=0.4); af.set_xticks(range(3)); af.set_xticklabels(['face', 'whisk', 'pupil'], fontsize=5.6); af.set_xlim(-0.55, 2.55); af.set_ylim(-0.6, 1.1); af.set_ylabel('net sign of single-ROI r\n(+1 = all positive)'); plabel(af, 'f', dx=-22, dy=6)
    af.plot([], [], 'o', color=C_SIG, ms=3.3, label='signal-like'); af.plot([], [], 'o', color=C_REJ, ms=3.3, label='rejected ROIs'); af.legend(frameon=False, fontsize=4.9, loc='lower left', handletextpad=0.1, borderpad=0.1)
    # g  ridge R2 vs controls
    for n in ALL:
        r = st[n].get('ridge')
        if not r: continue
        full = r['sig']['full']; noi = r.get('noise', np.nan); rej = np.nanmean(r['rej']) if r.get('rej') else np.nan; nul = np.nanmean(r['null_full']) if r.get('null_full') else np.nan
        xs_ = [0, 1, 2, 3]; ys_ = [full, noi, rej, nul]; ok = [np.isfinite(v) for v in ys_]; xx = [x for x, o in zip(xs_, ok) if o]; yy = [y for y, o in zip(ys_, ok) if o]
        ag.plot(xx, yy, color=FAINT, lw=0.5, zorder=1)
        for x, y, c in zip(xs_, ys_, (C_SIG, C_NOISE, C_REJ, MUTED)):
            if np.isfinite(y): sdot(ag, x, y, n, color=c, hollow_color=c, ms=3.3)
    vals = {n: st[n]['ridge']['sig']['full'] for n in UNFLAGGED if st[n].get('ridge')}; h = pooled_unflagged(vals, n_boot=4000)
    valr = {n: np.nanmean(st[n]['ridge']['rej']) for n in UNFLAGGED if st[n].get('ridge') and st[n]['ridge'].get('rej')}; hr = pooled_unflagged(valr, n_boot=4000)
    ag.plot([-0.3, 0.3], [h['est']] * 2, color=INK, lw=1.4, zorder=5); ag.plot([1.7, 2.3], [hr['est']] * 2, color=INK, lw=1.4, zorder=5)
    ledger_add('Fig2', 'ridge_R2', 'held-out R2 of the axon mean from behaviour (ridge, lags -1..+3 s)', h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'ambiguous', note='common-mode: random rejected-ROI sets reach %.2f' % hr['est'], unit='R2')
    ag.axhline(0, color=MUTED, lw=0.4); ag.set_xticks(range(4)); ag.set_xticklabels(['axons', 'noise-\nlike', 'rejected\nsets', 'shifted\nnull'], fontsize=5.2); ag.set_ylim(-0.12, 0.85); ag.set_ylabel('held-out R² (behaviour → axon mean)'); plabel(ag, 'g', dx=-22, dy=6)
    # h  alone vs unique per family (axon mean)
    for j, fam in enumerate(('face', 'whisk', 'pupil', 'eyelid', 'loco')):
        for k, kind in enumerate(('alone', 'unique')):
            vals = {}
            for n in ALL:
                r = st[n].get('ridge')
                if not r or fam not in r['sig'][kind]: continue
                v = r['sig'][kind][fam]
                if not np.isfinite(v): continue
                x = j + (k - 0.5) * 0.4 + rng.uniform(-0.07, 0.07); sdot(ah, x, v, n, color=C_BEH[CK[fam]], hollow_color=C_BEH[CK[fam]], ms=3.3, alpha=1.0 if kind == 'alone' else 0.38)
                if not is_flagged(n): vals[n] = v
            hh = pooled_unflagged(vals, n_boot=2000)
            if hh['n_sessions']: ah.plot([j + (k - 0.5) * 0.4 - 0.16, j + (k - 0.5) * 0.4 + 0.16], [hh['est']] * 2, color=INK, lw=1.3, zorder=5)
    ah.axhline(0, color=MUTED, lw=0.4); ah.set_xticks(range(5)); ah.set_xticklabels(['face', 'whisk', 'pupil', 'lid', 'ball'], fontsize=5.4); ah.set_xlim(-0.6, 4.6); ah.set_ylabel('held-out R²'); plabel(ah, 'h', dx=-22, dy=6)
    ah.text(0.0, 1.02, 'alone | unique (pale)', transform=ah.transAxes, fontsize=5.2, color=MUTED, va='bottom')
    # legend line for shapes
    fig.text(0.105, 0.004, 'c: thin blue = sessions (dashed grey = flagged), thick blue = mean of unflagged sessions; grey = random sets of rejected ROIs (mean +- 2 SD over sessions); amber = noise-like ROIs; ball: 4 unflagged sessions.\n'
             'e: n/N above each column = unflagged sessions with p < 0.01 against circularly shifted data. Marker shape = animal (circle RB24, square RB23, triangle RB22); filled = unflagged session (pooled), hollow = flagged (never pooled).\n'
             'Black bars = mean of the per-animal means, +- 95% hierarchical bootstrap CI (sessions within animal, then animals); n = 13 unflagged sessions from 3 mice (ball: 4 sessions, 1 mouse).', fontsize=4.8, color=MUTED)
    savefig(fig, 'Fig2_state_coupling')


if __name__ == '__main__':
    main()
