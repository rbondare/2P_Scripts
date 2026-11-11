"""Grey BREAKS between the subsessions of the lick-task sessions: what the animal does and what the LC axons do (task_deep_dive, 2026-11-11).

usage: python dd_break.py <part> [session ...]      parts: table | aligned | state | lta | matched | all      (pickles + CSVs in OUTDIR; REFRESH=1 recomputes)

BACKGROUND.  Between two subsessions (3-4 per session) the acquisition is restarted: the last trial ends, then 55-100 s of imaging with NO trial (the 'break'), then an imaging gap of 0.3-7 s and the
next subsession.  The Arduino lick log of the last trial keeps recording through the break (dd_common.Task._breaks gives the window and the log-capacity rule), so licks in the break are known on the
2P clock.  After the last trial of 1317 / 1423 / 1744 / 1752 there is the same kind of period ('tail', 96-135 s).  Nothing is ever rewarded in a break.  The PUPIL is 40-140% larger in every break
and falls back at its end: the visual condition of the break is not the grey of the trials (luminance), whatever the screen shows.

PARTS.
 table    per break: duration, licks / bouts / rate, rate by time in the break, pupil / face / whisking / jaw level relative to the task waits, lick rhythm vs the waits        -> break_table.csv
 aligned  break onset / offset aligned averages of the axon population (slow-baseline dF/F: Fc / running median over +-300 s, so that the 10th-percentile baseline of the usual dF/F cannot absorb a
          90 s level change), noise-like control, pupil, face, whisking, jaw, lick rate; window means vs the pre-break window; null = the same statistic at random alignment points of the task  -> break_aligned.csv
 state    inside the break vs equally long task stretches: SD of the 30 s high-passed population, its coupling to face / whisking / jaw / pupil dilation rate / lick rate                 -> break_state.csv
 lta      lick-bout-triggered average (population, jaw) for break bouts vs wait / anticipatory / reward / post bouts
 matched  behaviour-matched bout contrasts: break bouts vs the other classes (+ the same on the 'pre' window), break windows of 10 s vs random task windows with equal movement (engine of
          cross_session/matched_behaviour.py via dd_matched.py)                                                                                                       -> cache/matchedbrk_<session>.pkl
"""
import os, sys, csv, pickle, time, zlib, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dd_common as dd
import xs_common as xc
from xs_common import DT
from scipy.ndimage import percentile_filter

OUTDIR = dd.OUTDIR; CACHE = os.path.join(OUTDIR, 'cache'); os.makedirs(CACHE, exist_ok=True)
SESS = [s for s in dd.TASK_SESSIONS if s != 'AnimalRB24_261105_1520']           # 1520 has no break (3 s between its two subsessions)
UNFLAGGED = [s for s in SESS if s not in dd.FLAGGED]
hp30 = lambda x: x - xc.nsmooth(x, 30.0 / DT)
zs = lambda v: (v - np.nanmean(v)) / np.nanstd(v)
N_RAND = 400


def load(name):
    S, T = dd.load_session(name); return S, T


def slow_dff(S, idx):
    """population dF/F with a baseline that is slow compared with a break: Fc / running median over +-300 s - 1, per ROI (computed on 1 s bins, upsampled), mean over the ROIs of `idx`"""
    Fc = S.Fc[idx].astype(float); k = int(round(1.0 / DT)); n = (Fc.shape[1] // k) * k
    out = []
    for r in Fc:
        b = np.nanmean(r[:n].reshape(-1, k), 1)
        if (~np.isfinite(b)).all(): continue
        bad = ~np.isfinite(b)
        if bad.any(): b = np.interp(np.arange(len(b)), np.where(~bad)[0], b[~bad])
        med = percentile_filter(b, 50, size=601, mode='nearest')
        if np.nanmedian(med) < 5: continue                              # no light: dF/F meaningless
        d = np.repeat(b / med - 1, k); d = np.concatenate([d, np.full(Fc.shape[1] - n, d[-1])])
        d[~np.isfinite(r)] = np.nan; out.append(d)
    return np.nanmean(np.stack(out), 0)


def wmean(v, e0, a, b):
    i1 = int(np.floor((a - e0) / DT)); i2 = int(np.floor((b - e0) / DT))
    if i1 < 0 or i2 > len(v) or i2 <= i1: return np.nan
    x = v[i1:i2]; return float(np.nanmean(x)) if np.isfinite(x).sum() >= 0.5 * (i2 - i1) else np.nan


def sess_series(S, T):
    """series on the 5 Hz grid used by the break analyses"""
    sd = {'P': slow_dff(S, np.where(S.sig)[0])}
    if S.noise.sum() >= 3: sd['ctl'] = slow_dff(S, np.where(S.noise)[0])
    for b in ('pupil', 'face', 'whisk', 'jaw'):
        if b in S.beh: sd[b] = zs(S.beh[b])
    sd['lick'] = xc.nsmooth(S.beh['lick'], 2.0 / DT) if 'lick' in S.beh else None
    return sd


# =============================================================================================================  table
def part_table():
    rows = []; hdr = ['session', 'kind', 'break', 'dur_s', 'licks', 'lick_rate_hz', 'bouts', 'bout_rate_per_min', 'rate_0_20s', 'rate_20_40s', 'rate_40_60s', 'rate_60_end', 'wait_rate_hz', 'ili_median_ms_break', 'ili_median_ms_wait',
                      'licks_per_bout_break', 'licks_per_bout_wait', 'pupil_ratio', 'face_ratio', 'whisk_ratio', 'jaw_ratio', 'pupil_pre', 'pupil_in', 'pupil_post', 'log_capped']
    for s in SESS:
        S, T = load(s)
        if not T.brk: continue
        sg = T.seg; wait_rate = sg['n_mine'].sum() / (sg['t1'] - sg['t0']).sum()
        wm = np.zeros(S.nb, bool)
        for t in range(T.n): wm |= (S.t >= T.tr0[t] + 1) & (S.t < T.on[t])
        lv = {b: np.nanmean(S.beh[b][wm]) for b in ('pupil', 'face', 'whisk', 'jaw') if b in S.beh}
        # lick rhythm: inter-lick intervals < 0.5 s (inside bouts), waits vs breaks
        inw = np.zeros(len(T.lick), bool)
        for a, b in zip(sg['t0'], sg['t1']): inw |= (T.lick >= a) & (T.lick < b)
        def ili(m):
            x = np.diff(T.lick[m]); x = x[(x < 0.5) & (x > 0)]; return 1e3 * np.median(x) if len(x) else np.nan
        inb = T.in_break(T.lick) >= 0
        bw = [(T.bout_on[(T.bout_on >= a) & (T.bout_on < b)], T.bout_n[(T.bout_on >= a) & (T.bout_on < b)]) for a, b in zip(sg['t0'], sg['t1'])]
        nbw = np.concatenate([x[1] for x in bw]) if bw else np.array([])
        for j, x in enumerate(T.brk):
            L = x['t1'] - x['t0']; m = (T.lick >= x['t0']) & (T.lick < x['t1']); bm = (T.bout_on >= x['t0']) & (T.bout_on < x['t1'])
            r = lambda a, b: ((T.lick >= x['t0'] + a) & (T.lick < min(x['t0'] + b, x['t1']))).sum() / max(min(b, L) - a, 1e-9) if L > a else np.nan
            e0 = S.edges[0]; lvl = lambda b: wmean(S.beh[b], e0, x['t0'] + 10, x['t1'] - 5) / lv[b] if b in lv else np.nan
            pp = [wmean(S.beh['pupil'], e0, x['t0'] - 12, x['t0'] - 2), wmean(S.beh['pupil'], e0, x['t0'] + 10, x['t1'] - 5), wmean(S.beh['pupil'], e0, x['t1'] + 8, x['t1'] + 25)]
            rows.append([s, x['kind'], x['b'], round(L, 1), int(m.sum()), round(m.sum() / L, 3), int(bm.sum()), round(60 * bm.sum() / L, 2), round(r(0, 20), 2), round(r(20, 40), 2), round(r(40, 60), 2), round(r(60, 1e9), 2), round(wait_rate, 3),
                         round(ili(m), 0), round(ili(inw), 0), round(float(np.mean(T.bout_n[bm])), 2) if bm.any() else np.nan, round(float(nbw.mean()), 2) if len(nbw) else np.nan,
                         round(lvl('pupil'), 3), round(lvl('face'), 3), round(lvl('whisk'), 3), round(lvl('jaw'), 3), round(pp[0], 0), round(pp[1], 0), round(pp[2], 0) if np.isfinite(pp[2]) else np.nan, int(x['capped'])])
    with open(os.path.join(OUTDIR, 'break_table.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(hdr); w.writerows(rows)
    print('%-12s %-4s %2s %5s %5s %6s %5s | %5s %5s %5s %5s | wait %5s | pupil/face/whisk/jaw ratio' % ('session', 'kind', 'b', 'dur', 'licks', 'Hz', 'bouts', '0-20', '20-40', '40-60', '60+', 'Hz'))
    for r in rows: print('%-12s %-4s %2d %5.0f %5d %6.2f %5d | %5.2f %5.2f %5.2f %5.2f | wait %5.2f | %.2f %.2f %.2f %.2f' % (r[0][-11:], r[1], r[2], r[3], r[4], r[5], r[6], r[8], r[9], r[10], r[11], r[12], r[17], r[18], r[19], r[20]))
    return rows


# =============================================================================================================  aligned
WINS = {'onset': {'pre': (-30.0, -5.0), 'early': (0.0, 10.0), 'mid': (10.0, 40.0), 'late': (40.0, 'L-5')}, 'offset': {'pre': (-30.0, -5.0), 'post': (0.0, 20.0)}}


def _win(w, L):
    a, b = w; return a, (min(L - 5.0, 80.0) if b == 'L-5' else b)


def _anchor_stats(S, v, t_anchor, kind, L, e0):
    out = {}
    for nm, w in WINS[kind].items():
        a, b = _win(w, L); out[nm] = wmean(v, e0, t_anchor + a, t_anchor + b)
    return out


def part_aligned():
    res = {}; rng = np.random.default_rng(11)
    for s in SESS:
        S, T = load(s)
        if not T.brk: continue
        sd = sess_series(S, T); e0 = S.edges[0]; ser = {k: v for k, v in sd.items() if v is not None}
        brks = [x for x in T.brk if x['kind'] == 'brk']
        # task stretches usable as random anchors: the whole analysis window [-30, +100] s must not touch a break / tail
        blk = [(x['t0'] - 40.0, x['t1'] + 40.0) for x in T.brk]
        def free(t0, t1): return all(t1 < a or t0 > b for a, b in blk) and t0 > S.t[0] + 5 and t1 < S.t[-1] - 5
        out = {'onset': {}, 'offset': {}, 'n_brk': len(brks), 'L': [x['t1'] - x['t0'] for x in brks]}
        for kind in ('onset', 'offset'):
            for nm, v in ser.items():
                ev = [_anchor_stats(S, v, x['t0'] if kind == 'onset' else x['t1'], kind, x['t1'] - x['t0'], e0) for x in brks]
                nul = []
                for _ in range(N_RAND):
                    L = float(rng.choice(out['L']))
                    for _try in range(50):
                        ta = rng.uniform(S.t[0] + 40, S.t[-1] - 120)
                        if free(ta - 40, ta + L + 40): break
                    else: continue
                    nul.append(_anchor_stats(S, v, ta, kind, L, e0))
                out[kind][nm] = dict(event=ev, null=nul)
        res[s] = out
        print(s, 'done', flush=True)
    pickle.dump(res, open(os.path.join(CACHE, 'break_aligned.pkl'), 'wb'))
    return res


def summarise_aligned(res=None, pool=UNFLAGGED, n_boot=5000, seed=3):
    """pooled change vs pre-window (per break, mean over breaks of the pooled sessions, sessions weighted by their number of breaks), the same on random task alignments"""
    res = res or pickle.load(open(os.path.join(CACHE, 'break_aligned.pkl'), 'rb')); rng = np.random.default_rng(seed); rows = []
    for kind, wins in (('onset', ('early', 'mid', 'late')), ('offset', ('post',))):
        for nm in ('P', 'ctl', 'pupil', 'face', 'whisk', 'jaw', 'lick'):
            for w in wins:
                for grp, names in (('pooled', pool), ('flagged', [s for s in res if s not in pool])):
                    ev, nulls, ns = [], [], []
                    for s in names:
                        if s not in res or nm not in res[s][kind]: continue
                        d = res[s][kind][nm]; e = [x[w] - x['pre'] for x in d['event']]; e = [x for x in e if np.isfinite(x)]
                        nl = np.array([x[w] - x['pre'] for x in d['null']]); nl = nl[np.isfinite(nl)]
                        if len(e) == 0 or len(nl) < 30: continue
                        ev += e; nulls.append(nl); ns.append(len(e))
                    if not ev: continue
                    obs = float(np.mean(ev))
                    # null: draw as many random anchors per session as it has breaks, pool the same way
                    nd = np.array([np.mean(np.concatenate([rng.choice(nl, n) for nl, n in zip(nulls, ns)])) for _ in range(n_boot)])
                    # CI of the observed pooled mean: bootstrap over breaks inside sessions
                    bs = np.array([np.mean(rng.choice(ev, len(ev))) for _ in range(n_boot)])
                    p = float((1 + (np.abs(nd - nd.mean()) >= abs(obs - nd.mean())).sum()) / (1 + n_boot))
                    rows.append([kind, nm, w, grp, len(ev), round(obs, 4), round(np.percentile(bs, 2.5), 4), round(np.percentile(bs, 97.5), 4), round(float(nd.mean()), 4), round(float(nd.std()), 4), round(p, 4),
                                 int(np.sum(np.array(ev) > 0))])
    with open(os.path.join(OUTDIR, 'break_aligned.csv'), 'w', newline='') as fh:
        w_ = csv.writer(fh); w_.writerow(['alignment', 'series', 'window', 'group', 'n_breaks', 'change_vs_pre', 'ci_lo', 'ci_hi', 'null_mean', 'null_sd', 'p_vs_random_anchors', 'n_positive']); w_.writerows(rows)
    print('%-7s %-6s %-6s %-8s %3s %9s %20s %8s %7s %6s' % ('align', 'series', 'window', 'group', 'n', 'change', '95% CI (breaks)', 'null sd', 'p', 'n>0'))
    for r in rows: print('%-7s %-6s %-6s %-8s %3d %+9.3f [%+.3f, %+.3f] %8.3f %7.3f %4d' % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[9], r[10], r[11]))
    return rows


# =============================================================================================================  state (SD and coupling inside the break)
def part_state():
    res = {}; rng = np.random.default_rng(21)
    for s in SESS:
        S, T = load(s)
        if not T.brk: continue
        e0 = S.edges[0]; P = hp30(S.population('sig')); C = hp30(S.population('noise')) if S.noise.sum() >= 3 else None
        beh = {b: hp30(zs(S.beh[b])) for b in ('face', 'whisk', 'jaw') if b in S.beh}
        if 'dpupil' in S.beh: beh['dpupil'] = hp30(S.beh['dpupil'] / np.nanstd(S.beh['dpupil']))
        beh['lick'] = hp30(xc.nsmooth(S.beh['lick'], 1.0 / DT))
        blk = [(x['t0'] - 20.0, x['t1'] + 20.0) for x in T.brk]
        def stat(a, b):
            i1, i2 = int(np.floor((a - e0) / DT)), int(np.floor((b - e0) / DT))
            if i1 < 0 or i2 > S.nb: return None
            o = {'sd_P': np.nanstd(P[i1:i2])}
            if C is not None: o['sd_ctl'] = np.nanstd(C[i1:i2])
            for k, v in beh.items():
                ok = np.isfinite(P[i1:i2]) & np.isfinite(v[i1:i2])
                o['r_' + k] = np.corrcoef(P[i1:i2][ok], v[i1:i2][ok])[0, 1] if ok.sum() > 50 else np.nan
                o['sd_' + k] = np.nanstd(v[i1:i2])
            return o
        brks = [x for x in T.brk if x['kind'] == 'brk']
        ev = [stat(x['t0'] + 10, x['t1'] - 5) for x in brks]; ev = [e for e in ev if e]
        Ls = [x['t1'] - x['t0'] - 15 for x in brks]; nul = []
        for _ in range(N_RAND):
            L = float(rng.choice(Ls))
            for _try in range(50):
                ta = rng.uniform(S.t[0] + 40, S.t[-1] - L - 40)
                if all(ta + L < a or ta > b for a, b in blk): break
            else: continue
            o = stat(ta, ta + L)
            if o: nul.append(o)
        res[s] = dict(event=ev, null=nul); print(s, len(ev), 'breaks,', len(nul), 'null windows', flush=True)
    pickle.dump(res, open(os.path.join(CACHE, 'break_state.pkl'), 'wb')); return res


def summarise_state(res=None, pool=UNFLAGGED, n_boot=5000, seed=4):
    res = res or pickle.load(open(os.path.join(CACHE, 'break_state.pkl'), 'rb')); rng = np.random.default_rng(seed); rows = []
    keys = [k for k in res[pool[0]]['event'][0].keys() if k.startswith('r_') or k.startswith('sd_')]
    for k in keys:
        for grp, names in (('pooled', pool), ('flagged', [s for s in res if s not in pool])):
            ev, nulls, ns = [], [], []
            for s in names:
                if s not in res or k not in res[s]['event'][0]: continue
                e = np.array([x[k] for x in res[s]['event']], float); nl = np.array([x[k] for x in res[s]['null']], float); nl = nl[np.isfinite(nl)]
                if k.startswith('sd_'):
                    base = np.nanmean(nl); e = e / base; nl = nl / base                      # ratio to the session's own task level
                e = e[np.isfinite(e)]
                if len(e) == 0 or len(nl) < 30: continue
                ev += list(e); nulls.append(nl); ns.append(len(e))
            if not ev: continue
            obs = float(np.mean(ev)); nd = np.array([np.mean(np.concatenate([rng.choice(nl, n) for nl, n in zip(nulls, ns)])) for _ in range(n_boot)]); bs = np.array([np.mean(rng.choice(ev, len(ev))) for _ in range(n_boot)])
            p = float((1 + (np.abs(nd - nd.mean()) >= abs(obs - nd.mean())).sum()) / (1 + n_boot))
            rows.append([k, grp, len(ev), round(obs, 4), round(np.percentile(bs, 2.5), 4), round(np.percentile(bs, 97.5), 4), round(float(nd.mean()), 4), round(float(nd.std()), 4), round(p, 4)])
    with open(os.path.join(OUTDIR, 'break_state.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['statistic', 'group', 'n_breaks', 'in_break', 'ci_lo', 'ci_hi', 'task_mean', 'task_sd_of_pooled_mean', 'p_vs_random_task_windows']); w.writerows(rows)
    print('%-10s %-8s %3s %8s %18s %8s %7s %7s   (sd_* = ratio to the task level, r_* = correlation with the population)' % ('statistic', 'group', 'n', 'break', '95% CI', 'task', 'sd', 'p'))
    for r in rows: print('%-10s %-8s %3d %+8.3f [%+.3f, %+.3f] %+8.3f %7.3f %7.3f' % tuple(r[:2] + [r[2], r[3], r[4], r[5], r[6], r[7], r[8]]))
    return rows


# =============================================================================================================  lick-bout-triggered averages
TAU = np.arange(-5.0, 10.0 + DT / 2, DT)


def part_lta():
    res = {}
    for s in SESS:
        S, T = load(s)
        if not T.brk: continue
        P = hp30(S.population('sig')); C = hp30(S.population('noise')) if S.noise.sum() >= 3 else None; J = hp30(zs(S.beh['jaw'])) if 'jaw' in S.beh else None; e0 = S.edges[0]
        out = {}
        for c in ('brk', 'tail', 'wait_early', 'ant_to', 'ant_sq', 'rew', 'post_a', 'post_b'):
            b = T.bouts(c)
            if len(b) < 5: continue
            i0 = np.floor((b - e0) / DT).astype(int)
            def tr(v):
                ii = i0[:, None] + np.round(TAU / DT).astype(int)[None, :]; ok = (ii >= 0) & (ii < S.nb); o = np.full(ii.shape, np.nan); o[ok] = v[ii[ok]]; return o
            d = {'n': len(b), 'P': tr(P)}
            if C is not None: d['ctl'] = tr(C)
            if J is not None: d['jaw'] = tr(J)
            out[c] = d
        res[s] = out; print(s, {c: v['n'] for c, v in out.items()}, flush=True)
    pickle.dump(res, open(os.path.join(CACHE, 'break_lta.pkl'), 'wb')); return res


# =============================================================================================================  matched contrasts
def break_contrasts(R, T):
    import dd_matched as dm
    mb = dm.mb; tb = dm.tb; out = []; p = R.par; nb = R.nb
    g0, g1 = [int(round(g / DT)) for g in p['guard']]
    cl = {c: T.bouts(c) for c in ('brk', 'tail', 'wait_early', 'ant_to', 'ant_sq', 'rew', 'post_a', 'post_b')}
    for foc in ('brk', 'tail'):
        for ref in ('wait_early', 'ant_to', 'ant_sq', 'rew', 'post_a', 'post_b'):
            if len(cl[foc]) < dm.MIN_EV or len(cl[ref]) < dm.MIN_EV: continue
            for tag, spec in (('', (0.0, 3.0)), ('pre_', (-2.0, 0.0))):
                out.append(dict(key='%s%s_vs_%s' % (tag, foc, ref), kind='vs_class', spec=spec, focal=tb(R, cl[foc]), ref=tb(R, cl[ref])))
    # time in the break: the first 20 s may still be post-trial licking (lick rate is highest there), later bouts are 'pure' break licking
    bt0 = np.array([T.brk[j]['t0'] for j in T.in_break(T.bouts('brk'))]) if len(cl['brk']) else np.zeros(0)
    since = cl['brk'] - bt0
    for key, sel in (('brkLate', since >= 20.0), ('brkEarly', since < 20.0)):
        if sel.sum() < dm.MIN_EV or len(cl['wait_early']) < dm.MIN_EV: continue
        for tag, spec in (('', (0.0, 3.0)), ('pre_', (-2.0, 0.0))):
            out.append(dict(key='%s%s_vs_wait_early' % (tag, key), kind='vs_class', spec=spec, focal=tb(R, cl['brk'][sel]), ref=tb(R, cl['wait_early'])))
    # whole break windows (10 s) vs random task windows with equal movement: the pool excludes every window that touches a break / tail or a square onset guard
    brk_ev = np.concatenate([np.arange(int(np.floor((x['t0'] - R.e0) / DT)), int(np.floor((x['t1'] - R.e0) / DT)), int(round(1.0 / DT))) for x in T.brk]) if T.brk else np.zeros(0, int)
    for kind in ('brk', 'tail'):
        foc = []
        for x in T.brk:
            if x['kind'] != kind: continue
            t = x['t0'] + 10.0
            while t + 10.0 <= x['t1'] - 5.0: foc.append(t); t += 10.0
        if len(foc) >= dm.MIN_EV:
            out.append(dict(key='%s_state_10s' % kind, kind='vs_pool', spec=(0.0, 10.0), focal=tb(R, np.array(foc)), blockers=[(R.allon, g0, g1), (brk_ev, 0, int(round(5.0 / DT)))]))
    return out


def part_matched(names):
    import dd_matched as dm, matched_behaviour as mb
    from dd_fir import groups_from_log
    for nm in names:
        f = os.path.join(CACHE, 'matchedbrk_%s.pkl' % nm)
        if os.path.exists(f) and not os.environ.get('REFRESH'): print('cached', nm); continue
        t0 = time.time(); R, T = dm.make_rec(nm)
        if not T.brk: print('no break', nm); continue
        out = dict(name=nm, n_sig=R.n_sig, names=R.names, par=R.par, contrasts={}, counts=T.bout_class_counts)
        cons = break_contrasts(R, T); hp = lambda x: x - xc.nsmooth(x, R.par['hp_s'] / DT); outs = {'P': R.P}
        if R.S.noise.sum() >= 3: outs['control'] = hp(R.S.population('noise'))
        print('%s: behaviours %s, %d contrasts' % (nm, R.names, len(cons)), flush=True)
        for on, y in outs.items():
            R.set_outcome(y); out['contrasts'][on] = {}; R.par['n_sur'] = dm.N_SUR if on == 'P' else 0
            for c in cons:
                rng = np.random.default_rng(zlib.crc32((nm + on + c['key']).encode()))
                try: rr = mb.run_contrast(R, c, dm.N_NULL, rng, with_extras=(on == 'P'))
                except Exception as e: rr = dict(key=c['key'], note='failed: %s' % e, n_events=len(c['focal']))
                out['contrasts'][on][c['key']] = rr
                o = rr.get('obs')
                if o and np.isfinite(o.get('matched', np.nan)):
                    s = rr['stats']['matched']; ss = rr.get('stats_sim', {}).get('matched', {})
                    print('  %-7s %-26s n %3d valid %3d matched %3d | naive %+.3f matched %+.3f (shift z %+.1f p %.3f | bench z %+.1f p %.3f)' % (on, c['key'], rr['n_events'], o['n_focal'], o['n_matched'], o['naive'], o['matched'], s['z'], s['p'], ss.get('z', np.nan), ss.get('p', np.nan)), flush=True)
                else: print('  %-7s %-26s n %3d: %s' % (on, c['key'], rr.get('n_events', 0), rr.get('note', '')), flush=True)
        out['sec'] = time.time() - t0; pickle.dump(out, open(f, 'wb')); print('  %.0f s' % out['sec'], flush=True)


if __name__ == '__main__':
    part = sys.argv[1] if len(sys.argv) > 1 else 'all'; names = sys.argv[2:] or SESS
    if part in ('table', 'all'): part_table()
    if part in ('aligned', 'all'):
        if os.environ.get('REFRESH') or not os.path.exists(os.path.join(CACHE, 'break_aligned.pkl')): part_aligned()
        summarise_aligned()
    if part in ('state', 'all'):
        if os.environ.get('REFRESH') or not os.path.exists(os.path.join(CACHE, 'break_state.pkl')): part_state()
        summarise_state()
    if part in ('lta', 'all'): part_lta()
    if part in ('matched', 'all'): part_matched(names)
