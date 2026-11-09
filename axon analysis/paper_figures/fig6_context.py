"""Figure 6 - context: the lick task and the grey screens before / after looming.
a  6 min of RB24 10-30 1744 (lick task): licks, square onsets (hit / miss), axon mean, face motion, pupil   b  engagement vs axon fluctuation through two sessions (1744, and 1317 where engagement collapses)
c  square-aligned hit vs miss averages (raw / movement removed)   d  120 s of grey before and after loom blocks (unflagged grey sessions)
e  hit - miss response per session   f  change of hit rate vs change of axon fluctuation across subsessions   g  grey contrasts: across a loom block, and loom minus plain grey
usage: python fig6_context.py -> plots/summary_figures/Fig6_context.pdf/.png ; ledger entries
"""
import numpy as np, scipy.io as sio
from fig_common import *
import xs_common as xc
from grey_compare import contrasts

EX = 'AnimalRB24_261030_1744'; EX2 = 'AnimalRB24_261106_1317'
TASK_UNFL = [n for n in UNFLAGGED if n in TASK_SESSIONS]


def load_task_mat(n): return {k: v for k, v in sio.loadmat(r'C:\Users\rbondarenko\2P_cache\xs\%s_task.mat' % n, squeeze_me=True).items() if not k.startswith('__')}


def rolling_hit(on, hit, w=10):
    k = np.ones(w) / w; h = np.convolve(hit, k, mode='same'); h[:w // 2] = np.nan; h[-(w // 2):] = np.nan; return h


def main():
    st = load_state(); T = load_pickle('task_results.pkl'); G = load_pickle('grey_results.pkl'); res = G['results']
    fig = plt.figure(figsize=(W2, 222 * MM))
    gs = fig.add_gridspec(3, 1, height_ratios=[50, 52, 62], hspace=0.46, left=0.085, right=0.985, top=0.965, bottom=0.07)

    # ------------------------------------------------ a  raw excerpt
    top = gs[0].subgridspec(1, 2, wspace=0.18, width_ratios=[1.6, 1.0]); ga = top[0].subgridspec(4, 1, height_ratios=[1.1, 1.6, 1.0, 1.0], hspace=0.12)
    r = st[EX]; t = np.asarray(r['t'], float); P = xc.nsmooth(np.asarray(r['P'], float), 1.0 / xc.DT); bg = r['beh_grid']; tk = load_task_mat(EX); t0, t1 = 120.0, 480.0; sl = (t >= t0) & (t < t1); tt = t[sl] - t0
    a0 = fig.add_subplot(ga[0]); lk = tk['lick_t']; lk = lk[(lk >= t0) & (lk < t1)] - t0
    a0.vlines(lk, 0, 1, color=INK, lw=0.3); on = np.asarray(tk['on_t'], float); se = np.asarray(tk['st_end'], float); hit = np.asarray(tk['hit'], float)
    for o, e, h in zip(on, se, hit):
        if t0 <= o < t1: a0.axvspan(o - t0, e - t0, color=C_REJ, alpha=0.35, lw=0); a0.plot(o - t0, 1.18, 'v', color=INK, mfc=INK if h else 'white', ms=3.2, mew=0.7, clip_on=False)
    a0.set_ylim(0, 1.3); a0.set_yticks([]); a0.spines['left'].set_visible(False); a0.set_ylabel('licks', fontsize=5.4, rotation=0, ha='right', va='center', labelpad=3); a0.tick_params(axis='x', labelbottom=False, bottom=False); a0.spines['bottom'].set_visible(False)
    plabel(a0, 'a', dx=-28, dy=14); a0.set_title('%s: lick task, 6 min  (\u25bc = square onset: filled hit, open miss; grey = square on screen)' % short(EX), fontsize=5.8, loc='left', fontweight='normal', pad=6)
    for k, (y, col, yl, lw) in enumerate(((P[sl], C_SIG, 'axon mean\ndF/F', 1.0), (np.asarray(bg['face'])[sl], C_BEH['face'], 'face\nmotion', 0.6), (np.asarray(bg['pupil'])[sl], C_BEH['pupil'], 'pupil', 0.6))):
        a = fig.add_subplot(ga[1 + k], sharex=a0); a.plot(tt, y, color=col, lw=lw); a.set_yticks([]); a.spines['left'].set_visible(False); a.set_ylabel(yl, fontsize=5.4, rotation=0, ha='right', va='center', labelpad=3)
        if k < 2: a.tick_params(axis='x', labelbottom=False, bottom=False); a.spines['bottom'].set_visible(False)
    a.set_xlabel('time (s)'); a.set_xlim(0, t1 - t0)

    # ------------------------------------------------ b  engagement timelines (% of the first 10 min / 20 trials)
    gb = top[1].subgridspec(2, 1, hspace=0.50)
    for q, nm in enumerate((EX, EX2)):
        a = fig.add_subplot(gb[q]); d = T[nm]; tl = d['timeline']; tr = d['trials']; tmin = np.asarray(tl['t'], float) / 60; on = np.asarray(tr['on'], float) / 60; h = np.asarray(tr['hit'], float)
        rh = rolling_hit(on, h, 12); base = np.nanmean(h[:20]); lr = np.asarray(tl['lick_rate'], float); sd = np.asarray(tl['P_sd'], float)
        a.plot(on, 100 * rh / base, color=INK, lw=1.2, label='hit rate'); a.plot(tmin, 100 * lr / np.nanmean(lr[:10]), color=C_REJ, lw=1.0, label='lick rate'); a.plot(tmin, 100 * sd / np.nanmean(sd[:10]), color=C_SIG, lw=1.2, label='axon fluctuation SD')
        for s0 in np.asarray(tr['sub_t0'], float)[1:]: a.axvline(s0 / 60, color=MUTED, lw=0.4, ls=':')
        a.axhline(100, color=MUTED, lw=0.3); a.set_ylim(0, 230); a.set_ylabel('% of start', fontsize=5.6); a.set_title('%s%s' % (short(nm), '  \u2020 flagged' if is_flagged(nm) else ''), fontsize=5.6, loc='left', fontweight='normal', pad=2)
        if q == 0: a.legend(frameon=False, fontsize=4.8, loc='upper right', ncol=3, handlelength=1.0, columnspacing=0.8, borderpad=0.1, bbox_to_anchor=(1.0, 1.12)); plabel(a, 'b', dx=-24, dy=14); a.tick_params(labelbottom=False)
        else: a.set_xlabel('time in session (min)')

    # ------------------------------------------------ c, d  averages
    mid = gs[1].subgridspec(2, 2, wspace=0.34, hspace=0.20, width_ratios=[1.0, 1.55]); rel = np.asarray(T[EX]['trace_rel'], float)
    ac1 = fig.add_subplot(mid[0, 0]); ac2 = fig.add_subplot(mid[1, 0], sharex=ac1)
    for a, kh, km, lab in ((ac1, 'trace_hit', 'trace_miss', 'axon mean, raw'), (ac2, 'trace_res_hit', 'trace_res_miss', 'axon mean, movement removed')):
        H = []; M = []
        for n in TASK_SESSIONS:
            d = T[n]; fl = is_flagged(n)
            h_ = np.nanmean(d[kh], 0); m_ = np.nanmean(d[km], 0) if len(d[km]) else np.full_like(h_, np.nan)
            a.plot(rel, h_, color=FAINT if fl else C_SIG, lw=0.5, alpha=1 if fl else 0.45, ls='--' if fl else '-'); a.plot(rel, m_, color=FAINT if fl else C_REJ, lw=0.5, alpha=1 if fl else 0.7, ls='--' if fl else '-')
            if not fl: H.append(h_); M.append(m_)
        a.plot(rel, np.nanmean(H, 0), color=C_SIG, lw=1.6, label='hit'); a.plot(rel, np.nanmean(M, 0), color=C_REJ, lw=1.6, label='miss'); a.axvline(0, color=INK, lw=0.5, ls=':'); a.axhline(0, color=MUTED, lw=0.4); a.set_ylabel(lab, fontsize=5.4)
    ac1.tick_params(labelbottom=False); ac2.set_xlabel('time from square onset (s)'); ac1.legend(frameon=False, fontsize=4.8, loc='upper left', handlelength=1.0, borderpad=0.1); plabel(ac1, 'c', dx=-24, dy=8)
    ac1.set_title('%d unflagged task sessions (thin: sessions, dashed grey: 1317 / RB23 1704)' % len(TASK_UNFL), fontsize=5.2, loc='left', fontweight='normal')
    ad1 = fig.add_subplot(mid[0, 1]); ad2 = fig.add_subplot(mid[1, 1], sharex=ad1); GU = [n for n in UNFLAGGED if n in res and 'peri' in res[n]]
    nb = 12; edges = np.linspace(0, 120, nb + 1); ctr = 0.5 * (edges[1:] + edges[:-1]); allm = {}
    for key, a, col, lab in (('sig', ad1, C_SIG, 'signal-like axons'), ('rej', ad1, C_REJ, 'rejected ROIs'), ('b_face', ad2, C_BEH['face'], 'face motion'), ('b_pupil', ad2, C_BEH['pupil'], 'pupil'), ('b_whisk', ad2, C_BEH['whisk'], 'whisking')):
        ms = []
        for n in GU:
            pk = res[n]['peri']
            if key not in pk['pre']: continue
            rl = pk['rel']
            def bins(M): return np.array([[np.nanmean(M[j][(rl >= edges[i]) & (rl < edges[i + 1])]) if np.isfinite(M[j][(rl >= edges[i]) & (rl < edges[i + 1])]).any() else np.nan for i in range(nb)] for j in range(len(M))])
            bp_, bq_ = bins(pk['pre'][key]), bins(pk['post'][key]); allv = np.concatenate([bp_.ravel(), bq_.ravel()]); mu, sd = np.nanmean(allv), np.nanstd(allv)
            if not np.isfinite(sd) or sd == 0: continue
            ms.append(np.concatenate([(np.nanmean(bp_, 0) - mu) / sd, (np.nanmean(bq_, 0) - mu) / sd]))
        if ms:
            m = np.nanmean(ms, 0)
            for sl_, off in ((slice(0, nb), -120), (slice(nb, 2 * nb), 14)): a.plot(ctr + off, m[sl_], color=col, lw=1.5, label=lab if off == -120 else None)
            if key == 'sig':
                for q_ in ms:
                    for sl_, off in ((slice(0, nb), -120), (slice(nb, 2 * nb), 14)): a.plot(ctr + off, q_[sl_], color=C_SIG, lw=0.4, alpha=0.35)
    for a in (ad1, ad2): a.axvspan(0, 14, color='#E5E7EB', lw=0); a.axhline(0, color=MUTED, lw=0.4, ls=':')
    ad1.text(7, 0.97, 'loom block', transform=ad1.get_xaxis_transform(), ha='center', va='top', fontsize=4.8, color=MUTED); ad1.tick_params(labelbottom=False); ad1.set_ylabel('z (session grey)', fontsize=5.6); ad2.set_ylabel('z (session grey)', fontsize=5.6)
    ad2.set_xticks([-120, -60, 0, 14, 74, 134]); ad2.set_xticklabels(['\u2212120', '\u221260', '0', '', '+60', '+120']); ad2.set_xlabel('s before the loom block  |  s after it', fontsize=5.8)
    ad1.legend(frameon=False, fontsize=4.8, loc='upper left', ncol=2, handlelength=1.0, borderpad=0.1); ad2.legend(frameon=False, fontsize=4.8, loc='upper right', ncol=3, handlelength=1.0, borderpad=0.1, columnspacing=0.8); plabel(ad1, 'd', dx=-24, dy=8)
    ad1.set_title('recorded grey screens, %d unflagged sessions (%d loom blocks)' % (len(GU), sum(len(res[n]['peri_blocks']) for n in GU)), fontsize=5.2, loc='left', fontweight='normal')

    # ------------------------------------------------ e..g statistics
    bot = gs[2].subgridspec(1, 4, wspace=0.50, width_ratios=[0.8, 1.0, 0.95, 0.75]); ae = fig.add_subplot(bot[0]); af = fig.add_subplot(bot[1]); ag = fig.add_subplot(bot[2]); ag2 = fig.add_subplot(bot[3], sharey=ag); rng = np.random.default_rng(6)
    for q, (key, col, lab) in enumerate((('P_hp response', C_SIG, 'raw'), ('P_res response', C_MV if False else '#8FB3E8', 'movement\nremoved'))):
        vals = {}
        for n in TASK_SESSIONS:
            d = T[n]['hit_miss'][key]['diff']; sdot(ae, q + rng.uniform(-0.18, 0.18), d, n, color=col, hollow_color=MUTED, ms=3.8)
            if not is_flagged(n): vals[n] = d
        h = pooled_unflagged(vals, n_boot=4000)
        if h['n_sessions']:
            ae.plot([q - 0.3, q + 0.3], [h['est']] * 2, color=INK, lw=1.4, zorder=5); ae.plot([q, q], [h['lo'], h['hi']], color=INK, lw=0.9, zorder=5)
            ledger_add('Fig6', 'hit_minus_miss_' + ('raw' if q == 0 else 'mv'), 'axon response to the square, hit minus miss trials (%s)' % ('raw' if q == 0 else 'movement removed'), h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'supported' if (q == 0 and h['lo'] > 0) else 'not seen' if h['lo'] < 0 < h['hi'] else 'ambiguous', note='lick-task sessions; one animal (RB24) only; hit trials contain licking, miss trials do not', unit='dF/F (high-passed)')
    ae.axhline(0, color=MUTED, lw=0.4); ae.set_xticks([0, 1]); ae.set_xticklabels(['raw', 'movement\nremoved'], fontsize=5.4); ae.set_xlim(-0.6, 1.6); ae.set_ylabel('hit \u2212 miss response (dF/F)'); plabel(ae, 'e', dx=-24, dy=8)
    # f  engagement vs axon fluctuation
    xs_, ys_ = [], []
    for n in TASK_SESSIONS:
        subs = T[n]['subs']
        if len(subs) < 2: continue
        hx = np.array([100 * s_['hit'] / subs[0]['hit'] for s_ in subs]); sy = np.array([100 * s_['P_sd'] / subs[0]['P_sd'] for s_ in subs]); fl = is_flagged(n)
        af.plot(hx, sy, '-', color=FAINT if fl else MUTED, lw=0.7, zorder=1)
        for j, (x, y) in enumerate(zip(hx, sy)):
            if j == 0: continue
            sdot(af, x, y, n, color=C_SIG, hollow_color=MUTED, ms=4.0)
            if not fl: xs_.append(x); ys_.append(y)
    af.axhline(100, color=MUTED, lw=0.3, ls=':'); af.axvline(100, color=MUTED, lw=0.3, ls=':'); af.set_xlabel('hit rate, % of subsession 1'); af.set_ylabel('axon fluctuation SD, % of subsession 1'); plabel(af, 'f', dx=-24, dy=8)
    af.text(0.03, 0.97, 'later subsessions\n(lines start at 100 / 100);\n1317 \u2020: hit rate falls\nto 59%, axon SD stays', transform=af.transAxes, fontsize=4.6, color=MUTED, va='top', ha='left')
    ledger_add('Fig6', 'engagement_vs_axonSD', 'axon fluctuation SD in later subsessions of the lick task (% of subsession 1)', float(np.mean(ys_)), float(np.percentile(ys_, 2.5)), float(np.percentile(ys_, 97.5)), len(TASK_UNFL), 1, 'not seen', note='mean hit rate in the same subsessions %.0f%% of subsession 1: it fell only in 1744 (77 to 50%%) and in flagged 1317 (to 59%%), axon SD did not follow; range = spread of subsession values, not a CI' % np.mean(xs_), unit='% of subsession 1')
    # g  grey contrasts (two axes sharing the rows)
    GUN = {n: res[n] for n in GU}; c_un = contrasts(GUN); c_all = contrasts({n: res[n] for n in res if 'windows' in res[n]})
    rows = [('sig', 'mean', 'axons: mean dF/F', C_SIG), ('sig', 'sd', 'axons: fluctuation SD', C_SIG), ('sig', 'rate', 'axons: event rate', C_SIG), ('rej', 'sd', 'rejected ROIs: SD', C_REJ), ('beh', 'face', 'face motion', C_BEH['face']), ('beh', 'whisk', 'whisking', C_BEH['whisk']), ('beh', 'pupil', 'pupil area', C_BEH['pupil']), ('beh', 'eyelid', 'eyelid opening', C_BEH['eyelid'])]
    cols = (('loom', 'after − before\na loom block', ag), ('diff', 'minus the same time\nover plain grey', ag2))
    for q, (ck, ttl, a) in enumerate(cols):
        for i, (pop, m, lab, col) in enumerate(rows):
            call = c_all.get((pop, m), {}).get(ck); cun = c_un.get((pop, m), {}).get(ck)
            if call:
                for n, v in call['per_session'].items():
                    if np.isfinite(v): a.plot(v, i + rng.uniform(-0.17, 0.17), MARK[animal(n)], mfc='white' if is_flagged(n) else col, mec=MUTED if is_flagged(n) else 'white', mew=0.7 if is_flagged(n) else 0.3, ms=3.2, zorder=3)
            if cun:
                a.plot([cun['ci'][0], cun['ci'][1]], [i, i], color=INK, lw=1.2, zorder=4); a.plot(cun['mean'], i, '|', color=INK, ms=7, mew=1.6, zorder=5)
                a.text(2.95, i, 'p = %.2f' % cun['p'], fontsize=4.8, va='center', color=INK if cun['p'] < 0.05 else MUTED)
                if (pop, m) in (('sig', 'mean'), ('sig', 'sd'), ('sig', 'rate')): ledger_add('Fig6', 'grey_%s_%s_%s' % (ck, pop, m), 'grey screen change %s: axons %s (SD of metric)' % ('across a loom block' if ck == 'loom' else '(loom minus plain grey)', m), cun['mean'], cun['ci'][0], cun['ci'][1], cun['n'], len({animal(n) for n in cun['per_session']}), 'not seen' if cun['p'] > 0.05 else 'ambiguous', note='sign-flip p = %.2f over %d sessions; no no-loom control block' % (cun['p'], cun['n']), unit='SD of metric')
        a.axvline(0, color=MUTED, lw=0.5); a.set_xlim(-2.4, 3.6); a.set_xticks([-2, -1, 0, 1, 2]); a.set_title(ttl, fontsize=5.4, loc='left', fontweight='bold'); a.set_xlabel('SD of the metric', fontsize=5.6)
    ag.set_yticks(range(len(rows))); ag.set_yticklabels([r[2] for r in rows], fontsize=5.4); ag.set_ylim(len(rows) - 0.4, -0.7); ag2.tick_params(axis='y', left=False, labelleft=False)
    for t_, r_ in zip(ag.get_yticklabels(), rows): t_.set_color(r_[3])
    plabel(ag, 'g', dx=-62, dy=14)
    fig.text(0.085, 0.004, 'Hit trials contain licking, miss trials do not: any raw hit - miss difference is also a movement difference. Pooled task sessions: 1752, 1744, 1423 (11-05), 1520 (one animal); 1317 and RB23 1704 flagged.\n'
             'Movement removed = face + whisking + ball ridge model. Grey: 100 s windows (10 s guard) before / after recorded 120 s grey screens; unflagged sessions RB24 1538, 1534, RB23 1642, RB22 1138;\n'
             'g: dots = sessions (hollow = flagged 1604, 1306), black bar = bootstrap 95% CI over the unflagged sessions, p = exact sign-flip over sessions; there is no no-loom control block.\n'
             'Marker shape = animal (circle RB24, square RB23, triangle RB22), filled = unflagged, hollow = flagged (not pooled).', fontsize=4.7, color=MUTED)
    savefig(fig, 'Fig6_context')


if __name__ == '__main__':
    main()
