"""Figure 9 - axon calcium vs whisking: overlaid traces and the shuffle test.
columns a-d: three example sessions (one per animal)   rows: raster + overlay of the axon mean on whisking (z) | scatter | shuffle null of the peak r | cross-correlation with its null band
bottom: e observed r vs shuffle null for every session   f single-ROI r with whisking (signal-like vs rejected vs shuffle)   g r and partial r (face / pupil removed) for whisking, face, pupil
Statistic = the one of Figure 2 / session_qc.csv: both signals band-passed 1-30 s, circular cross-correlation, peak |r| within +-10 s; null = the same statistic for the behaviour circularly shifted (>= 120 s from zero lag).
usage: python fig9_whisking.py -> plots/summary_figures/Fig9_whisking.pdf/.png ; ledger entries
"""
import os, pickle, numpy as np
from fig_common import *
import xs_common as xc
from scipy.stats import gaussian_kde

DT = xc.DT; MAXLAG = int(10 / DT); SEP = int(120 / DT)
EXS = [('AnimalRB24_261021_1630', 'strong'), ('AnimalRB22_261027_1138', 'strong'), ('AnimalRB23_261023_1642', 'weak')]
C_W = C_BEH['whisk']


def bp(x): return xc.bandpass(np.asarray(x, float))


def wrap9(s, w=205):
    import textwrap
    return textwrap.fill(s.replace('\n', ' '), w)


def typical_window(P, w, win_s=180.0):
    """window whose local r (band-passed) is the median of all 3 min windows of the session - a typical, not the best, excerpt"""
    n = int(win_s / DT); starts = np.arange(0, len(P) - n, int(15 / DT)); rs = []; sds = []
    for i in starts:
        a, b = P[i:i + n], w[i:i + n]; ok = np.isfinite(a) & np.isfinite(b)
        if ok.mean() < 0.85: rs.append(np.nan); sds.append(np.nan); continue
        rs.append(np.corrcoef(a[ok], b[ok])[0, 1]); sds.append(np.nanstd(b[ok]))
    rs = np.array(rs); sds = np.array(sds); ok = np.isfinite(rs); cand = np.where(ok)[0]
    k = cand[np.argmin(np.abs(rs[cand] - np.nanmedian(rs[cand])))]
    return starts[k], n, rs[k], float(np.nanmedian(rs[cand]))


def example(name):
    f = os.path.join(CACHE, 'whisk_ex_%s.pkl' % name)
    if os.path.exists(f): return pickle.load(open(f, 'rb'))
    st = load_state()[name]; P = np.asarray(st['P'], float); W = np.asarray(st['beh_grid']['whisk'], float); fP, fW = bp(P), bp(W)
    i0, n, rloc, rmed = typical_window(fP, fW)
    S = xc.Session(name, with_rejected=False, verbose=False); sig = np.where(S.sig)[0]
    M = np.array([bp(S.dff[k]) for k in sig]); ok = np.isfinite(fP) & np.isfinite(fW)
    rr = np.array([np.corrcoef(m[ok], fW[ok])[0, 1] for m in M]); order = np.argsort(rr)[::-1]
    sl = slice(i0, i0 + n)
    _, r_lag, full = xc.xcorr_fft(fP, fW, MAXLAG); (lg, rpk), nul = xc.peak_stats(full, MAXLAG, SEP)
    nn = len(full); null_r0 = full[SEP:nn - SEP]
    d = dict(t=S.t[sl] - S.t[sl][0], M=M[order][:, sl], P=fP[sl], W=fW[sl], P_dff=P[sl], n_sig=int(len(sig)), Wnan=~np.isfinite(W[sl]), r_loc=float(rloc), r_med=float(rmed), r0=float(full[0]), peak=float(rpk), peak_lag=float(lg * DT), nul=nul, null_r0=null_r0.astype(np.float32),
             lags_s=np.arange(-MAXLAG, MAXLAG + 1) * DT, r_lag=r_lag, x=fW[ok][::10], y=fP[ok][::10], sdP=float(np.nanstd(fP)), sd_dff=float(np.nanstd(P)))
    pickle.dump(d, open(f, 'wb')); return d


def partial_stats(name, st):
    """r and partial r (the other two behaviours regressed out of both signals) at lag 0 for whisking, face, pupil (pupil taken 1.2 s later), with a circular-shift null"""
    f = os.path.join(CACHE, 'whisk_partial_%s.pkl' % name)
    if os.path.exists(f): return pickle.load(open(f, 'rb'))
    r = st[name]; bg = r['beh_grid']; out = {}
    if not r.get('usable') or any(b not in bg for b in ('whisk', 'face', 'pupil')): pickle.dump(out, open(f, 'wb')); return out
    fP = bp(r['P']); B = {'whisk': bp(bg['whisk']), 'face': bp(bg['face']), 'pupil': np.roll(bp(bg['pupil']), -int(round(1.2 / DT)))}
    ok = np.isfinite(fP) & np.all([np.isfinite(v) for v in B.values()], 0); n = int(ok.sum())
    if n < 600: pickle.dump(out, open(f, 'wb')); return out
    y = fP[ok]; Z = {k: v[ok] for k, v in B.items()}
    def resid(v, others):
        X = np.stack([o - o.mean() for o in others] + [np.ones_like(v)], 1); return v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
    rng = np.random.default_rng(5); shifts = rng.integers(SEP, n - SEP, 600)
    for b in ('whisk', 'face', 'pupil'):
        oth = [Z[k] for k in Z if k != b]; ry, rb = resid(y, oth), resid(Z[b], oth)
        raw = float(np.corrcoef(y, Z[b])[0, 1]); par = float(np.corrcoef(ry, rb)[0, 1])
        nr = np.array([np.corrcoef(y, np.roll(Z[b], s))[0, 1] for s in shifts]); np_ = np.array([np.corrcoef(ry, np.roll(rb, s))[0, 1] for s in shifts])
        out[b] = dict(raw=raw, partial=par, null_raw=float(np.percentile(np.abs(nr), 95)), null_par=float(np.percentile(np.abs(np_), 95)), p_par=float((1 + (np.abs(np_) >= abs(par)).sum()) / (1 + len(np_))))
    pickle.dump(out, open(f, 'wb')); return out


def main():
    st = load_state(); fig = plt.figure(figsize=(W2, 244 * MM))
    gs = fig.add_gridspec(2, 1, height_ratios=[168, 58], hspace=0.24, left=0.085, right=0.985, top=0.968, bottom=0.115)
    top = gs[0].subgridspec(5, 3, height_ratios=[1.0, 1.5, 2.1, 1.5, 1.7], hspace=0.66, wspace=0.28)
    for j, (nm, kind) in enumerate(EXS):
        d = example(nm); t = d['t']; dd = st[nm]['xc']['whisk']; rej = [p[1] for p in dd['rej_peaks']]
        # row 0  raster of the signal-like ROIs (band-passed, sorted by their r with whisking)
        a0 = fig.add_subplot(top[0, j]); M = d['M']; a0.imshow(M, aspect='auto', cmap='Blues', vmin=0, vmax=np.nanpercentile(M, 99), extent=(t[0], t[-1], len(M), 0), interpolation='nearest')
        a0.set_yticks([]); a0.tick_params(labelbottom=False, bottom=False); a0.spines['bottom'].set_visible(False); a0.spines['left'].set_visible(False)
        a0.set_title('%s  (%s, %s)' % (short(nm), 'RB24' if 'RB24' in nm else 'RB23' if 'RB23' in nm else 'RB22', 'strong coupling' if kind == 'strong' else 'weak coupling'), fontsize=6.0, loc='left', pad=3)
        a0.set_ylabel('%d ROIs' % d['n_sig'], fontsize=5.4, rotation=0, ha='right', va='center', labelpad=2)
        if j == 0: plabel(a0, 'a', dx=-30, dy=6)
        # row 1  overlay (z)
        a1 = fig.add_subplot(top[1, j], sharex=a0); zP = d['P'] / d['sdP']; zW = d['W'] / np.nanstd(d['W'])
        a1.plot(t, zW, color=C_W, lw=0.9, label='whisking'); a1.plot(t, zP, color=C_SIG, lw=1.0, label='axon mean')
        nm_ = d['Wnan']; dd_ = np.diff(np.r_[0, nm_.astype(int), 0]); s_, e_ = np.where(dd_ == 1)[0], np.where(dd_ == -1)[0]          # camera not recording (between stimulus blocks): no whisking data, the axon trace continues
        for a_, b_ in zip(s_, e_):
            if (b_ - a_) * DT >= 1.5: a1.axvspan(t[a_], t[min(b_, len(t) - 1)], color=FAINT, alpha=0.55, lw=0, zorder=0); a1.text(0.5 * (t[a_] + t[min(b_, len(t) - 1)]), 0.03 if t[a_] < t[-1] - 40 else 0.97, 'no video', transform=a1.get_xaxis_transform(), ha='center', va='bottom' if t[a_] < t[-1] - 40 else 'top', fontsize=4.6, color=MUTED)
        a1.axhline(0, color=MUTED, lw=0.3); a1.set_ylabel('z (band-passed)', fontsize=5.4); a1.set_xlabel('time (s)', fontsize=5.6); a1.set_xlim(t[0], t[-1])
        y0, y1 = a1.get_ylim(); a1.set_ylim(y0 - 0.9, y1); a1.plot([t[-1] - 22, t[-1] - 2], [y0 - 0.55] * 2, color=INK, lw=0.9, solid_capstyle='butt'); a1.text(t[-1] - 12, y0 - 0.45, '20 s', ha='center', va='bottom', fontsize=5.2)
        a1.text(0.0, 1.03, 'whisking', color=C_W, fontsize=5.2, transform=a1.transAxes, va='bottom', fontweight='bold'); a1.text(0.24, 1.03, 'axon mean', color=C_SIG, fontsize=5.2, transform=a1.transAxes, va='bottom', fontweight='bold')
        a1.text(1.0, 1.03, 'typical window: r = %.2f' % d['r_loc'], transform=a1.transAxes, ha='right', va='bottom', fontsize=4.7, color=MUTED)
        a1.tick_params(labelbottom=False, bottom=False); a1.spines['bottom'].set_visible(False); a1.set_xlabel('')
        if j == 0: plabel(a1, 'b', dx=-30, dy=6)
        # row 2  scatter, whole session
        a2 = fig.add_subplot(top[2, j]); x = d['x'] / np.nanstd(d['x']); y = d['y'] / d['sdP']
        a2.plot(x, y, 'o', ms=1.3, mfc=C_SIG, mec='none', alpha=0.22, rasterized=True); m_, c_ = np.polyfit(x, y, 1); xx = np.array([np.percentile(x, 1), np.percentile(x, 99)]); a2.plot(xx, m_ * xx + c_, color=INK, lw=1.0)
        a2.set_xlabel('whisking (z, band-passed)'); a2.set_ylabel('axon mean (z)'); a2.set_xlim(np.percentile(x, 0.5) - 0.3, np.percentile(x, 99.5) + 0.3); a2.set_ylim(np.percentile(y, 0.5) - 0.3, np.percentile(y, 99.5) + 0.3)
        a2.text(0.04, 0.96, 'r (lag 0) = %.2f\nwhole session, every 2 s' % d['r0'], transform=a2.transAxes, va='top', fontsize=5.4)
        if j == 0: plabel(a2, 'c', dx=-30, dy=6)
        # row 3  shuffle null of the peak statistic
        a3 = fig.add_subplot(top[3, j]); nul = d['nul']; hi = max(0.85, abs(d['peak']) + 0.08); bins = np.linspace(0, hi, 40)
        a3.hist(nul, bins=bins, color=C_REJ, alpha=0.75, lw=0, density=True, label='time-shifted whisking')
        a3.axvline(abs(d['peak']), color=C_SIG, lw=1.6, label='observed'); ymax = a3.get_ylim()[1]; a3.plot([], [], 'v', ms=2.8, color=C_REJ, mec=INK, mew=0.3, label='random rejected-ROI sets')
        for rr_ in rej: a3.plot(abs(rr_), ymax * 0.06, 'v', ms=2.6, color=C_REJ, mec=INK, mew=0.3, clip_on=False)
        a3.set_xlabel('peak |r| within ±10 s'); a3.set_ylabel('density'); a3.set_xlim(0, hi); a3.set_yticks([])
        a3.set_title('observed %.2f; shuffle 95%% %.2f; %s' % (abs(d['peak']), np.percentile(nul, 95), fmt_p(dd['p'])), fontsize=5.0, loc='left', fontweight='normal', pad=2)
        if j == 0: plabel(a3, 'd', dx=-30, dy=6); a3.legend(frameon=False, fontsize=4.6, loc='center', handlelength=1.0, borderpad=0.1, bbox_to_anchor=(0.55, 0.6))
        # row 4  cross-correlation with the shuffle band
        a4 = fig.add_subplot(top[4, j]); lg = d['lags_s']; lo_, hi_ = np.percentile(d['null_r0'], [2.5, 97.5])
        a4.axhspan(lo_, hi_, color=C_REJ, alpha=0.35, lw=0, label='shuffle 95%'); a4.plot(lg, d['r_lag'], color=C_SIG, lw=1.3, label='axons')
        rc = dd.get('r_rej_curves')
        if rc is not None: a4.plot(lg, np.mean(rc, 0), color=C_REJ, lw=1.0, label='rejected ROIs')
        a4.axvline(0, color=MUTED, lw=0.4, ls=':'); a4.axhline(0, color=MUTED, lw=0.3); a4.set_xlim(-10, 10); a4.set_xlabel('lag (s)  [whisking later →]'); a4.set_ylabel('r')
        if j == 0: plabel(a4, 'e', dx=-30, dy=6); a4.legend(frameon=False, fontsize=4.6, loc='upper left', handlelength=1.0, borderpad=0.1, labelspacing=0.2)

    # ------------------------------------------------ f, g, h
    bot = gs[1].subgridspec(1, 3, wspace=0.40, width_ratios=[1.35, 1.0, 1.35]); af = fig.add_subplot(bot[0]); ag = fig.add_subplot(bot[1]); ah = fig.add_subplot(bot[2]); rng = np.random.default_rng(3)
    # f  observed vs null per session
    names = [n for n in ALL if st[n].get('usable') and 'whisk' in st[n]['xc']]; vals = {}; ksig = 0
    order = sorted(names, key=lambda n: st[n]['xc']['whisk']['peak_r'])
    for i, n in enumerate(order):
        d = st[n]['xc']['whisk']; fl = is_flagged(n)
        af.plot([i, i], [-d['null95'], d['null95']], color=C_REJ, lw=3.2, alpha=0.55, solid_capstyle='butt', zorder=1); sdot(af, i, d['peak_r'], n, color=C_W, hollow_color=MUTED, ms=4.2)
        if not fl: vals[n] = d['peak_r']; ksig += int(d['p'] < 0.01 and d['peak_r'] > 0)
    h = pooled_unflagged(vals, n_boot=6000); af.axhline(0, color=MUTED, lw=0.4)
    af.set_xticks(range(len(order))); af.set_xticklabels([short(n) + ('*' if is_flagged(n) else '') for n in order], rotation=90, fontsize=4.2); af.set_ylim(-0.5, 0.9); af.set_xlim(-0.7, len(order) - 0.3); af.set_ylabel('peak r with whisking (±10 s)'); plabel(af, 'f', dx=-24, dy=6)
    af.text(0.02, 0.97, 'dots: observed; grey bar: 95%% of shuffled |r| of that session; * = flagged\nunflagged: pooled r = %.2f [%.2f, %.2f]; %d/%d sessions p < 0.01 and r > 0' % (h['est'], h['lo'], h['hi'], ksig, h['n_sessions']), transform=af.transAxes, fontsize=4.7, va='top')
    ledger_add('Fig9', 'whisk_peak_r', 'peak r axon mean vs whisking, shuffle-tested per session', h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'supported', note='%d/%d unflagged sessions p<0.01 vs time-shifted whisking; RB23 1642 and RB22 1427 are the weak sessions' % (ksig, h['n_sessions']), unit='r')
    # g  single-ROI r
    edges = np.linspace(-0.6, 1.0, 33); ctr = 0.5 * (edges[1:] + edges[:-1]); pools = {}; frac = {'sig': {}, 'rej': {}}
    for grp in ('sig', 'rej'):
        allr = []
        for n in UNFLAGGED:
            r = st[n].get('roi', {}).get('whisk', {}).get(grp)
            if not r or len(r['r']) < 3: continue
            allr += list(r['r']); frac[grp][n] = float(np.mean((np.asarray(r['r']) > 0) & (np.asarray(r['p']) < 0.01)))
        pools[grp] = np.array(allr)
    for grp, col in (('rej', C_REJ), ('sig', C_SIG)):
        cnt, _ = np.histogram(pools[grp], bins=edges); ag.fill_between(ctr, 0, 100 * cnt / cnt.sum(), step='mid', color=col, alpha=0.55 if grp == 'sig' else 0.5, lw=0.0, label=('signal-like ROIs' if grp == 'sig' else 'rejected ROIs') + ' (n = %d)' % len(pools[grp]))
        ag.step(ctr, 100 * cnt / cnt.sum(), where='mid', color=col, lw=0.9)
    ag.axvline(0, color=MUTED, lw=0.4, ls=':'); ag.set_xlabel('single-ROI r with whisking (band-passed)'); ag.set_ylabel('% of ROIs'); ag.set_ylim(0, 15.5); ag.set_xlim(-0.65, 1.05); ag.legend(frameon=False, fontsize=4.8, loc='upper left', handlelength=1.0, borderpad=0.1); plabel(ag, 'g', dx=-24, dy=6)
    hs = pooled_unflagged(frac['sig'], n_boot=4000); hr = pooled_unflagged(frac['rej'], n_boot=4000)
    ag.text(0.02, 0.80, 'ROIs with r > 0 and p < 0.01 against their own\ntime-shift null (chance 0.5%%), unflagged sessions:\nsignal-like %.0f%%, rejected %.0f%%' % (100 * hs['est'], 100 * hr['est']), transform=ag.transAxes, ha='left', va='top', fontsize=4.5)
    ledger_add('Fig9', 'roi_frac_whisk', 'fraction of signal-like ROIs with r > 0 and p < 0.01 against time-shifted whisking', hs['est'], hs['lo'], hs['hi'], hs['n_sessions'], hs['n_animals'], 'ambiguous', note='coupling is real but not specific to the accepted axons: %.0f%% of the rejected ROIs of the same sessions pass the same test (common-mode)' % (100 * hr['est']), unit='fraction of ROIs')
    # h  r and partial r
    cols = [('whisk', C_W), ('face', C_BEH['face']), ('pupil', C_BEH['pupil'])]; labs = []; xi = 0; PS = {n: partial_stats(n, st) for n in ALL}
    for b, col in cols:
        for kind in ('raw', 'partial'):
            v = {}; nl = []
            for n in ALL:
                d = PS[n].get(b)
                if not d: continue
                x = xi + rng.uniform(-0.17, 0.17); y = d[kind]; sdot(ah, x, y, n, color=col, hollow_color=MUTED, ms=3.6, alpha=1.0 if kind == 'raw' else 0.5)
                if not is_flagged(n): v[n] = y; nl.append(d['null_raw' if kind == 'raw' else 'null_par'])
            hh = pooled_unflagged(v, n_boot=3000)
            if hh['n_sessions']:
                ah.fill_between([xi - 0.34, xi + 0.34], -np.mean(nl), np.mean(nl), color=C_REJ, alpha=0.35, lw=0, zorder=0); ah.plot([xi - 0.3, xi + 0.3], [hh['est']] * 2, color=INK, lw=1.4, zorder=5); ah.plot([xi, xi], [hh['lo'], hh['hi']], color=INK, lw=0.9, zorder=5)
                if kind == 'partial': ledger_add('Fig9', 'partial_r_' + b, 'partial r of axon mean with %s (the other two behaviours regressed out, lag 0)' % b, hh['est'], hh['lo'], hh['hi'], hh['n_sessions'], hh['n_animals'], 'supported' if hh['lo'] > 0.1 else 'ambiguous', note='raw r %.2f; band = 95%% of time-shifted |r|' % pooled_unflagged({n: PS[n][b]['raw'] for n in UNFLAGGED if PS[n].get(b)}, n_boot=2000)['est'], unit='r')
            labs.append('%s\n%s' % ({'whisk': 'whisk', 'face': 'face', 'pupil': 'pupil'}[b], 'r' if kind == 'raw' else 'partial')); xi += 1
        xi += 0.5
    pos = [0, 1, 2.5, 3.5, 5, 6]; ah.set_xticks([p for p in pos]); ah.set_xticklabels(labs, fontsize=4.9); ah.axhline(0, color=MUTED, lw=0.4); ah.set_ylim(-0.3, 0.95); ah.set_xlim(-0.6, 6.6); ah.set_ylabel('r at lag 0 (pupil +1.2 s)'); plabel(ah, 'h', dx=-24, dy=6)
    ah.text(0.0, 1.03, 'partial: the other two behaviours regressed out of both signals;\ngrey = 95% of time-shifted |r|', transform=ah.transAxes, fontsize=4.6, color=MUTED, va='bottom')
    fig.text(0.085, 0.004, wrap9('Axon mean = mean best dF/F of the signal-like ROIs; both signals band-passed 1-30 s (slow drift and fast noise removed) for all statistics. a, b: ROI raster (band-passed dF/F, white = 0, sorted by r\n'
             'with whisking) and overlaid z-scored traces of a typical 3 min window (median windowed r among windows with >= 85% valid whisking samples, not the best); grey band = camera not recording (10-16 s between stimulus blocks), the axon trace is continuous. c-e: whole session.\n'
             'd: null = whisking circularly shifted by >= 120 s (348 offsets; same statistic: peak |r| within +-10 s); triangles = the same statistic for 20 random sets of rejected ROIs of equal size\n'
             '(they share part of the signal: common-mode, see Figs 2g, 7f). f-h: marker shape = animal (circle RB24, square RB23, triangle RB22), filled = unflagged (pooled), hollow = flagged (*).\n'
             'Black = mean of per-animal means with 95% hierarchical bootstrap CI. Examples: one session per animal; RB23 1642 is a weak session shown on purpose.'), fontsize=4.7, color=MUTED)
    savefig(fig, 'Fig9_whisking')


if __name__ == '__main__':
    main()
