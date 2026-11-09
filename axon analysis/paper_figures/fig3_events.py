"""Figure 3 - arousal events and state tuning.
a  three example epochs (face-motion bout, pupil dilation, locomotion onset)  b  event-triggered averages (axon mean z, behaviour z) with the circular-shift null
c  tuning curves (axon mean vs decile of behaviour)  d  response size per session vs null  e  active - quiet, running - still (dF/F) vs random rejected ROIs
usage: python fig3_events.py -> plots/summary_figures/Fig3_events.pdf/.png ; ledger entries
"""
import numpy as np
from fig_common import *
import xs_common as xc
from event_triggered import load_events, LAGS

EVT = [('face_bout', 'face-motion\nbout onset', 'face'), ('whisk_bout', 'whisking\nbout onset', 'whisk'), ('pupil_dilation', 'pupil dilation\nonset', 'pupil'), ('loco_onset', 'locomotion\nonset (ball)', 'ball'), ('saccade', 'saccade-like eye\nmovement', 'dpupil')]
BEH_KEY = {'face': 'face', 'whisk': 'whisk', 'pupil': 'pupil', 'loco': 'loco'}


def example_event(st, ev, name, key, win=(-15, 40)):
    e = ev[name]['events'][key]; r = st[name]; t = np.asarray(r['t'], float); P = np.asarray(r['P'], float); idx = e['onset_idx']
    k = int(np.argmin(np.abs(e['R_events'] - np.nanmean(e['R_events']))))           # the event whose response is closest to the session mean: representative, not the best one
    i0 = int(idx[k]); a = i0 + int(win[0] / xc.DT); b = i0 + int(win[1] / xc.DT); sl = slice(a, b + 1); tt = (np.arange(a, b + 1) - i0) * xc.DT; bg = r['beh_grid']
    out = dict(t=tt, P=xc.nsmooth(P, 1.0 / xc.DT)[sl], face=np.asarray(bg['face'])[sl], pupil=np.asarray(bg['pupil'])[sl], ball=np.asarray(bg['loco'])[sl] if 'loco' in bg else None, n=len(idx), k=k + 1, R=float(e['R_events'][k]))
    return out


def main():
    st = load_state(); ev = {n: load_events(n) for n in ALL}
    fig = plt.figure(figsize=(W2, 232 * MM))
    gs = fig.add_gridspec(4, 1, height_ratios=[38, 70, 36, 50], hspace=0.42, left=0.085, right=0.985, top=0.955, bottom=0.075)

    # ------------------------------------------------ a  examples
    ga = gs[0].subgridspec(1, 3, wspace=0.28)
    ex = [('AnimalRB24_261021_1630', 'face_bout', 'face-motion bout'), ('AnimalRB24_261029_1752', 'pupil_dilation', 'pupil dilation'), ('AnimalRB24_261105_1423', 'loco_onset', 'locomotion onset')]
    for j, (nm, key, lab) in enumerate(ex):
        d = example_event(st, ev, nm, key); sub = ga[j].subgridspec(4, 1, height_ratios=[1.8, 1, 1, 1], hspace=0.12); t = d['t']
        axs = []
        for k, (y, col, yl) in enumerate(((d['P'], C_SIG, 'mean\ndF/F'), (d['face'], C_BEH['face'], 'face\nmotion'), (d['pupil'], C_BEH['pupil'], 'pupil'), (d['ball'], C_BEH['ball'], 'ball'))):
            a = fig.add_subplot(sub[k], sharex=axs[0] if axs else None); axs.append(a)
            if y is None: a.text(0.5, 0.5, 'no usable ball', transform=a.transAxes, ha='center', va='center', fontsize=5, color=MUTED)
            else: a.plot(t, y, color=col, lw=0.9 if k == 0 else 0.6)
            a.axvline(0, color=INK, lw=0.5, ls='--'); a.set_yticks([]); a.spines['left'].set_visible(False); a.set_ylabel(yl, fontsize=5.2, rotation=0, ha='right', va='center', labelpad=2)
            if k < 3: a.tick_params(axis='x', labelbottom=False, bottom=False); a.spines['bottom'].set_visible(False)
        a.set_xlabel('time from onset (s)', fontsize=5.8); a.set_xlim(t[0], t[-1])
        axs[0].set_title('%s: %s\n(event %d of %d, closest to the session mean)' % (short(nm), lab, d['k'], d['n']), fontsize=5.6, loc='left', fontweight='normal', pad=2)
        if j == 0: plabel(axs[0], 'a', dx=-30, dy=12)

    # ------------------------------------------------ b  event-triggered averages
    gb = gs[1].subgridspec(2, len(EVT), hspace=0.12, wspace=0.34, height_ratios=[1.35, 1])
    lags = LAGS; axb = []
    for j, (key, lab, bk) in enumerate(EVT):
        a1 = fig.add_subplot(gb[0, j], sharey=axb[0] if axb else None); a2 = fig.add_subplot(gb[1, j], sharex=a1); axb.append(a1)
        nt = []; means = []; ns = 0
        for n in ALL:
            e = ev[n]['events'].get(key) if ev[n] else None
            if not e or e['n'] < 3: continue
            fl = is_flagged(n); m = e['Pz'][0]
            a1.plot(lags, m, color=FAINT if fl else C_SIG, lw=0.5, alpha=1 if fl else 0.45, ls='--' if fl else '-', zorder=2)
            if not fl: means.append(m); nt.append(e['null_traces_z'].astype(float)); ns += 1
        if means:
            L = min(len(x) for x in nt); pooled = np.mean([x[:L] for x in nt], 0)          # mean over sessions of independent per-session shift draws = null of the pooled mean
            a1.fill_between(lags, np.percentile(pooled, 2.5, 0), np.percentile(pooled, 97.5, 0), color=C_REJ, alpha=0.5, lw=0, zorder=1); a1.plot(lags, np.mean(means, 0), color=C_SIG, lw=1.7, zorder=4)
        a1.axvline(0, color=INK, lw=0.5, ls=':'); a1.axhline(0, color=MUTED, lw=0.4); a1.set_title('%s\n(n = %d sessions)' % (lab, ns), fontsize=6, loc='left', fontweight='bold'); a1.tick_params(labelbottom=False)
        if j == 0: a1.set_ylabel('axon mean (z)'); plabel(a1, 'b', dx=-24, dy=16)
        else: a1.tick_params(labelleft=False)
        for bkey, col in (('face', C_BEH['face']), ('whisk', C_BEH['whisk']), ('pupil', C_BEH['pupil']), ('loco', C_BEH['ball'])):
            Y = [ev[n]['events'][key][bkey][0] for n in UNFLAGGED if ev[n] and key in ev[n]['events'] and ev[n]['events'][key].get('n', 0) >= 3 and ev[n]['events'][key].get(bkey) is not None]
            if Y: a2.plot(lags, np.nanmean(Y, 0), color=col, lw=1.0, label={'face': 'face', 'whisk': 'whisk', 'pupil': 'pupil', 'loco': 'ball'}[bkey])
        a2.axvline(0, color=INK, lw=0.5, ls=':'); a2.axhline(0, color=MUTED, lw=0.4); a2.set_xlabel('time from onset (s)', fontsize=5.8); a2.set_xlim(-10, 20)
        if j == 0: a2.set_ylabel('behaviour (z)'); a2.legend(frameon=False, fontsize=4.8, loc='upper left', handlelength=1.0, borderpad=0.1, labelspacing=0.15)
        else: a2.tick_params(labelleft=False)
    axb[0].set_ylim(-0.7, 2.6)

    # ------------------------------------------------ c  tuning curves
    gc = gs[2].subgridspec(1, 4, wspace=0.30); axc = []
    for j, (bk, lab) in enumerate((('face', 'face motion'), ('whisk', 'whisking'), ('pupil', 'pupil area (+1.2 s)'), ('loco', 'ball speed'))):
        a = fig.add_subplot(gc[j], sharey=axc[0] if axc else None); axc.append(a); Y = []
        for n in ALL:
            t_ = ev[n]['tuning'].get(bk) if ev[n] else None
            if not t_: continue
            fl = is_flagged(n); a.plot(t_['x'], t_['yz'], color=FAINT if fl else C_BEH[{'loco': 'ball'}.get(bk, bk)], lw=0.5, alpha=1 if fl else 0.5, ls='--' if fl else '-')
            if not fl: Y.append(t_['yz'])
        if Y: a.plot((np.arange(10) + 0.5) * 10, np.nanmean(Y, 0), color=C_BEH[{'loco': 'ball'}.get(bk, bk)], lw=1.8)
        a.axhline(0, color=MUTED, lw=0.4); a.set_xlabel('%s, percentile' % lab, fontsize=5.8); a.set_title('n = %d' % len(Y), fontsize=5.6, loc='left', fontweight='normal')
        if j == 0: a.set_ylabel('axon mean (z)'); plabel(a, 'c', dx=-24, dy=6)
        else: a.tick_params(labelleft=False)

    # ------------------------------------------------ d  response size per session, e  states
    gd = gs[3].subgridspec(1, 2, wspace=0.28, width_ratios=[1.55, 1.0]); ad = fig.add_subplot(gd[0]); ae = fig.add_subplot(gd[1]); rng = np.random.default_rng(3)
    for j, (key, lab, bk) in enumerate(EVT):
        vals = {}; k_sig = 0; tot = 0; nullw = []; nr = []
        for n in ALL:
            e = ev[n]['events'].get(key) if ev[n] else None
            if not e or e['n'] < 3: continue
            sd = e['sd']; x = j + rng.uniform(-0.2, 0.2); sdot(ad, x, e['Rz'], n, color=C_BEH[bk], hollow_color=MUTED, ms=3.6)
            if not is_flagged(n): vals[n] = e['R']; tot += 1; k_sig += int(e['p'] < 0.01 and e['R'] > 0); nullw.append(np.percentile(e['null_R'] / sd, [2.5, 97.5])); nr.append(e['null_R'] / sd)
        if nullw: nw = np.mean(nullw, 0); ad.fill_between([j - 0.38, j + 0.30], nw[0], nw[1], color=C_REJ, alpha=0.4, lw=0, zorder=0)
        h = pooled_unflagged(vals, n_boot=4000)
        if h['n_sessions']:
            zs_ = {n: ev[n]['events'][key]['Rz'] for n in vals}; hz = pooled_unflagged(zs_, n_boot=4000)
            L = min(len(x) for x in nr); pn = np.mean([x[:L] for x in nr], 0); ppool = float((1 + (np.abs(pn) >= abs(np.mean(list(zs_.values())))).sum()) / (1 + L)); lo_n, hi_n = np.percentile(pn, [2.5, 97.5])
            ad.fill_between([j + 0.34, j + 0.46], lo_n, hi_n, color=INK, alpha=0.35, lw=0, zorder=0)
            ad.plot([j - 0.33, j + 0.27], [hz['est']] * 2, color=INK, lw=1.4, zorder=5); ad.plot([j - 0.03, j - 0.03], [hz['lo'], hz['hi']], color=INK, lw=0.9, zorder=5); ad.text(j, 2.55, '%d/%d' % (k_sig, tot), ha='center', va='bottom', fontsize=5, color=MUTED)
            verdict = 'supported' if hz['lo'] > 0 and k_sig >= 0.6 * tot and tot >= 6 else 'ambiguous' if hz['lo'] > 0 or k_sig >= 0.4 * tot else 'not seen'
            win = ev[[n for n in vals][0]]['events'][key]['win']
            ledger_add('Fig3', 'event_' + key, 'axon response (%g..%g s minus %g..%g s, session z) to %s' % (win[0][0], win[0][1], win[1][0], win[1][1], lab.replace('\n', ' ')), hz['est'], hz['lo'], hz['hi'], hz['n_sessions'], hz['n_animals'], verdict,
                       note='dF/F %.2f [%.2f, %.2f]; sessions with p<0.01: %d/%d; pooled circular-shift p = %.3f' % (h['est'], h['lo'], h['hi'], k_sig, tot, ppool), unit='z')
    ad.axhline(0, color=MUTED, lw=0.4); ad.set_xticks(range(len(EVT))); ad.set_xticklabels([e[1] for e in EVT], fontsize=5.2); ad.set_xlim(-0.6, len(EVT) - 0.4); ad.set_ylim(-0.8, 2.5); ad.set_ylabel('response (session z)'); plabel(ad, 'd', dx=-24, dy=8)
    ad.text(0.0, 1.07, 'n/N = unflagged sessions with p < 0.01 (circular shift); light grey = per-session 95% null range, dark grey = null of the pooled mean', transform=ad.transAxes, fontsize=4.8, color=MUTED, va='bottom')
    # e  states
    for j, (key, lab, colr) in enumerate((('face_active_minus_quiet', 'active − quiet\n(face-motion terciles)', C_BEH['face']), ('run_minus_still', 'running − still\n(ball)', C_BEH['ball']))):
        vals = {}
        for n in ALL:
            s = st[n].get('states', {})
            if key not in s: continue
            v = s[key][0]; x = j + rng.uniform(-0.2, 0.2); sdot(ae, x - 0.12, v, n, color=colr, hollow_color=MUTED, ms=3.6)
            if not is_flagged(n): vals[n] = v
            if key == 'face_active_minus_quiet' and s.get(key + '_rej'):
                rv = np.nanmean(s[key + '_rej']); sdot(ae, x + 0.28, rv, n, color=C_REJ, hollow_color=C_REJ, ms=3.2)
        h = pooled_unflagged(vals, n_boot=4000)
        if h['n_sessions']:
            ae.plot([j - 0.4, j - 0.0], [h['est']] * 2, color=INK, lw=1.4, zorder=5); ae.plot([j - 0.2, j - 0.2], [h['lo'], h['hi']], color=INK, lw=0.9, zorder=5)
            ledger_add('Fig3', 'state_' + key, '%s (dF/F)' % lab.replace('\n', ' '), h['est'], h['lo'], h['hi'], h['n_sessions'], h['n_animals'], 'supported' if h['lo'] > 0 else 'ambiguous', unit='dF/F')
    ae.axhline(0, color=MUTED, lw=0.4); ae.set_xticks([0.08, 1.08]); ae.set_xticklabels(['active − quiet\n(face-motion terciles)', 'running − still\n(ball)'], fontsize=5.2); ae.set_xlim(-0.6, 1.6); ae.set_ylabel('difference in mean dF/F'); plabel(ae, 'e', dx=-24, dy=8)
    ae.plot([], [], 'o', color=C_REJ, ms=3.2, label='random rejected ROIs'); ae.legend(frameon=False, fontsize=4.9, loc='upper right', handletextpad=0.1, borderpad=0.1)
    fig.text(0.085, 0.004, 'Onsets: hysteresis crossings of the smoothed behaviour after ≥ 5 s (pupil 8 s) of quiescence, excluding 2 s before to 6 s after any stimulus / trial onset; the triggering behaviour rises at 0 by construction.\n'
             'Response = mean(0..5 s) − mean(−6..−1 s) in session z; saccade-like events (eye displacement > 4° in 100 ms from 50 Hz DLC eye position, 3.5-28 per min, detection noise-limited): mean(0..2 s) − mean(−3..−1 s).\n'
             'Null: the same onsets circularly shifted by a random common offset (>= 60 s), 500 draws; null of the pooled mean = mean over sessions of independent draws.\n'
             'Marker shape = animal, filled = unflagged, hollow = flagged (not pooled); black bars = mean of per-animal means +- 95% hierarchical bootstrap CI.', fontsize=4.8, color=MUTED)
    savefig(fig, 'Fig3_events')


if __name__ == '__main__':
    main()
