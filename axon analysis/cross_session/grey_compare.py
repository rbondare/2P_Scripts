"""Cross-session analysis 3: the RECORDED GREY screens before and after the looming blocks (spontaneous activity in a controlled visual state).

usage: python grey_compare.py [session ...]      (default GREY_SESSIONS; output grey_results.pkl + grey_windows.csv next to the figures)

Every looming block of these sessions is preceded and followed by 120 s of static grey (specParams.pre_grey_t / post_grey_t; not stimulus frames, so b_t0 / b_t1 of the
aggregated file exclude them). A grey counts as recorded if the block is >= 121 s away from its neighbouring blocks / the imaging edge (stimulus gap = grey + wait).
Windows (100 s, 10 s guard against the last / first trial): PRE = [t0-110, t0-10], POST = [t1+10, t1+110].
Per window and population (signal-like mean = 'sig', noise-like mean = 'noise', random sets of rejected suite2p ROIs of the size of the signal-like set = 'rej'):
  mean best dF/F | SD of the band-passed (1-30 s) trace | event rate (peaks > 2 SD of the pooled grey windows, >= 3 s apart) | fraction of time above 1 SD
  behaviour (z over the session): face motion, whisking, pupil, eyelid, locomotion (ball only if plausible) ; coupling r / slope of sig vs face motion in the window
  'res' = sig after a behaviour-only ridge model (face + whisking + locomotion, lags -1..+3 s) trained on the grey windows of the OTHER blocks of the session (leave-one-block-out)
Contrasts (standardised by the SD of the metric over all windows of the session, then mean over sessions; exact sign-flip test over sessions + bootstrap over sessions):
  across a loom block  : POST(k) - PRE(k)                       (one loom block and ~190 s elapsed)
  across plain grey    : PRE(k+1) - POST(k)  for consecutive loom blocks (no loom, ~170 s elapsed: the time-in-session / drift reference)
  difference           : across loom - across grey              (what the loom block adds beyond elapsed time; the only contrast that is not confounded with time)
Epoch comparison: coupling of the signal-like mean with face motion in spontaneous block / grey windows / loom blocks / other stimulus blocks of the same session.
"""
import sys, os, pickle, itertools, time, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xs_common as xc
from xs_common import DT
from stim_response import fit_ridge, predict
from scipy.signal import find_peaks

FIG = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session'
GREY_SESSIONS = ['AnimalRB24_261023_1538', 'AnimalRB24_261027_1534', 'AnimalRB24_261106_1604', 'AnimalRB23_261023_1642', 'AnimalRB23_261027_1306', 'AnimalRB22_261027_1138']
PRE_W, POST_W, GREY_LEN, GUARD = (-110.0, -10.0), (10.0, 110.0), 120.0, 1.0
BEH = ['face', 'whisk', 'pupil', 'eyelid', 'loco']
POPS = ['sig', 'noise', 'rej']
METRICS = ['mean', 'sd', 'rate', 'active']
rng = np.random.default_rng(7)


def zs(v):
    v = np.asarray(v, float); s = np.nanstd(v); return (v - np.nanmean(v)) / (s if s > 0 else 1.0)


def blocks_of(S):
    X = S.X; keys = np.atleast_1d(np.asarray(X['blk_key'])).astype(str); sp = np.atleast_2d(np.asarray(X['blk_span'], float)); it = np.atleast_1d(np.asarray(X['blk_interrupted'], float)).astype(bool)
    return [dict(k=i, key=keys[i], t0=sp[i, 0], t1=sp[i, 1], interrupted=bool(it[i])) for i in range(len(keys)) if np.isfinite(sp[i]).all()]


def grey_blocks(S, blks):
    """loom blocks with a recorded pre / post grey (gap to the neighbouring blocks / imaging edge >= GREY_LEN + 1 s)"""
    out = []; t_first, t_last = S.tca[0], S.tca[-1]
    for b in blks:
        if not b['key'].startswith('loom') or b['interrupted']: continue
        prev = [c['t1'] for c in blks if c['t1'] < b['t0'] - 1e-6]; nxt = [c['t0'] for c in blks if c['t0'] > b['t1'] + 1e-6]
        gap_pre = b['t0'] - (max(prev) if prev else t_first); gap_post = (min(nxt) if nxt else t_last) - b['t1']
        b = dict(b, gap_pre=gap_pre, gap_post=gap_post, pre_ok=gap_pre >= GREY_LEN + 1, post_ok=gap_post >= GREY_LEN + 1)
        out.append(b)
    return out


class Traces:
    def __init__(self, S, grey_mask):
        self.S = S; nb = S.nb; self.e0 = S.edges[0]
        self.raw = {}; self.bp = {}
        self.raw['sig'] = S.population('sig'); self.raw['noise'] = S.population('noise') if S.noise.sum() >= 3 else None
        self.raw['rej'] = None; self.rej_sets = []
        if S.rej is not None and len(S.rej['dff']) >= 10:
            R = S.rej['dff']; m = min(int(S.sig.sum()), len(R))
            self.rej_sets = [np.nanmean(R[rng.choice(len(R), m, replace=False)], 0) for _ in range(20)]
        for k, v in self.raw.items():
            if v is not None: self.bp[k] = xc.bandpass(v)
        self.rej_bp = [xc.bandpass(v) for v in self.rej_sets]
        self.beh = {b: zs(S.beh[b]) for b in BEH if b in S.beh and (b != 'loco' or S.loco_ok) and np.isfinite(S.beh[b]).mean() > 0.3}
        self.beh_bp = {b: xc.bandpass(S.beh[b]) for b in self.beh}
        # event thresholds from the pooled grey windows
        self.thr = {}; self.ev = {}
        for k in list(self.bp) + [('rej', i) for i in range(len(self.rej_bp))]:
            x = self.bp[k] if not isinstance(k, tuple) else self.rej_bp[k[1]]
            sd = np.nanstd(x[grey_mask & np.isfinite(x)]) if (grey_mask & np.isfinite(x)).sum() > 100 else np.nanstd(x)
            xs_ = xc.nsmooth(x, 1.0 / DT); pk, _ = find_peaks(np.nan_to_num(xs_), height=2 * sd, distance=int(3 / DT))
            ev = np.zeros(len(x), bool); ev[pk] = True; self.ev[k] = ev; self.thr[k] = sd

    def idx(self, a, b):
        i1 = int(round((a - self.e0) / DT)); i2 = int(round((b - self.e0) / DT)); return max(i1, 0), min(i2, self.S.nb)


def win_metrics(T, a, b):
    i1, i2 = T.idx(a, b); n = i2 - i1
    out = dict(valid=False)
    if n < 0.9 * (b - a) / DT: return out
    ps = T.raw['sig'][i1:i2]
    if np.isfinite(ps).mean() < 0.8: return out
    out['valid'] = True; sl = slice(i1, i2)
    def pm(raw, bp, ev, thr):
        x = raw[sl]; xb = bp[sl]
        return dict(mean=float(np.nanmean(x)), sd=float(np.nanstd(xb)), rate=float(ev[sl].sum() / ((b - a) / 60.0)), active=float(np.nanmean(xb > thr)))
    out['sig'] = pm(T.raw['sig'], T.bp['sig'], T.ev['sig'], T.thr['sig'])
    out['noise'] = pm(T.raw['noise'], T.bp['noise'], T.ev['noise'], T.thr['noise']) if T.raw['noise'] is not None else None
    if T.rej_sets:
        ms = [pm(T.rej_sets[i], T.rej_bp[i], T.ev[('rej', i)], T.thr[('rej', i)]) for i in range(len(T.rej_sets))]
        out['rej'] = {m: float(np.mean([d[m] for d in ms])) for m in METRICS}
    else: out['rej'] = None
    out['beh'] = {b: float(np.nanmean(T.beh[b][sl])) for b in T.beh}
    cp = {}
    xb = T.bp['sig'][sl]
    for b in ('face', 'whisk'):
        if b in T.beh_bp:
            yb = T.beh_bp[b][sl]; ok = np.isfinite(xb) & np.isfinite(yb)
            if ok.sum() > 100 and np.nanstd(yb[ok]) > 0: cp['r_' + b] = float(np.corrcoef(xb[ok], yb[ok])[0, 1]); cp['slope_' + b] = float(np.cov(xb[ok], yb[ok])[0, 1] / np.var(yb[ok]))
    out['coupling'] = cp
    return out


def design(T, fam, lags):
    return np.concatenate([np.stack([np.roll(T.beh_bp[b], -l) for l in lags], 1) for b in fam], 1)


def stat_group(vals_by_session, nboot=4000):
    """vals: {session: standardised diff (mean over its blocks)} -> mean over sessions, exact sign-flip p, bootstrap CI over sessions"""
    v = np.array([x for x in vals_by_session.values() if np.isfinite(x)])
    if len(v) < 3: return None
    m = v.mean(); signs = np.array(list(itertools.product([-1, 1], repeat=len(v)))); nul = (signs * v).mean(1)
    p = float((np.abs(nul) >= abs(m) - 1e-12).mean())
    bs = np.array([v[rng.integers(0, len(v), len(v))].mean() for _ in range(nboot)])
    return dict(n=len(v), mean=float(m), p=p, ci=(float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))), per_session=dict(vals_by_session))


def analyse(name):
    S = xc.Session(name, with_rejected=True, verbose=True); blks = blocks_of(S); gb = grey_blocks(S, blks)
    res = dict(name=name, dur_min=S.dur / 60, n_sig=int(S.sig.sum()), n_noise=int(S.noise.sum()), loco_ok=S.loco_ok, beh=[b for b in BEH if b in S.beh and (b != 'loco' or S.loco_ok)])
    res['blocks'] = [dict(b) for b in gb]
    ok_b = [b for b in gb if b['pre_ok'] and b['post_ok']]
    print('  loom blocks %d, with both greys recorded: %d (%s)' % (len(gb), len(ok_b), ', '.join('blk%d gaps %.0f/%.0f s' % (b['k'] + 1, b['gap_pre'], b['gap_post']) for b in gb)))
    if len(ok_b) < 2: res['why'] = 'fewer than 2 loom blocks with both greys'; return res
    # grey mask for event thresholds
    gm = np.zeros(S.nb, bool)
    for b in ok_b:
        for (a, c) in ((b['t0'] + PRE_W[0], b['t0'] + PRE_W[1]), (b['t1'] + POST_W[0], b['t1'] + POST_W[1])):
            i1 = int(round((a - S.edges[0]) / DT)); i2 = int(round((c - S.edges[0]) / DT)); gm[max(i1, 0):max(i2, 0)] = True
    T = Traces(S, gm)
    # ---- windows
    wins = []
    for j, b in enumerate(ok_b):
        for kind, (a, c) in (('pre', (b['t0'] + PRE_W[0], b['t0'] + PRE_W[1])), ('post', (b['t1'] + POST_W[0], b['t1'] + POST_W[1]))):
            w = win_metrics(T, a, c); w.update(block=j, blk_k=b['k'], kind=kind, t0=a, t1=c, loom_dir=b['key']); wins.append(w)
    # ---- behaviour-adjusted sig (leave-one-block-out ridge on the grey windows)
    fam = [b for b in ('face', 'whisk', 'loco') if b in T.beh_bp]; lags = (np.arange(-1.0, 3.01, 0.5) / DT).round().astype(int); res['adj_family'] = fam
    Pbp = T.bp['sig']
    if fam:
        Xd = design(T, fam, lags)
        for j in range(len(ok_b)):
            tr_mask = np.zeros(S.nb, bool)
            for w in wins:
                if w['block'] != j and w['valid']:
                    i1, i2 = T.idx(w['t0'], w['t1']); tr_mask[i1:i2] = True
            model = fit_ridge(np.where(tr_mask[:, None], Xd, np.nan), np.where(tr_mask, Pbp, np.nan))
            if model is None: continue
            Phat = predict(model, Xd); resid = Pbp - Phat
            for w in wins:
                if w['block'] == j and w['valid']:
                    i1, i2 = T.idx(w['t0'], w['t1']); w['res'] = dict(mean=float(np.nanmean(resid[i1:i2])), sd=float(np.nanstd(resid[i1:i2])), pred_mean=float(np.nanmean(Phat[i1:i2])))
    res['windows'] = wins
    # ---- peri-block traces (for the time-course figure): -120..0 before t0 and 0..120 after t1 on the 0.2 s grid, session units
    n120 = int(GREY_LEN / DT); tr = dict(rel=np.arange(n120) * DT, pre={}, post={})
    series = {'sig': T.raw['sig'], 'noise': T.raw['noise'], 'rej': (np.mean(T.rej_sets, 0) if T.rej_sets else None)}
    series.update({'b_' + b: T.beh[b] for b in T.beh})
    for nm, v in series.items():
        if v is None: continue
        tr['pre'][nm] = []; tr['post'][nm] = []
        for b in ok_b:
            i1, i2 = T.idx(b['t0'] - GREY_LEN, b['t0']); seg = np.full(n120, np.nan); sl = v[i1:i2]; seg[:len(sl)] = sl[:n120]; tr['pre'][nm].append(seg)
            i1, i2 = T.idx(b['t1'], b['t1'] + GREY_LEN); seg = np.full(n120, np.nan); sl = v[i1:i2]; seg[:len(sl)] = sl[:n120]; tr['post'][nm].append(seg)
        tr['pre'][nm] = np.array(tr['pre'][nm]); tr['post'][nm] = np.array(tr['post'][nm])
    res['peri'] = tr; res['peri_blocks'] = [dict(k=b['k'], key=b['key'], dur=b['t1'] - b['t0']) for b in ok_b]
    # ---- epoch comparison: coupling and fluctuation size in spontaneous block / grey / loom blocks / other stimulus blocks
    ep = {}
    def mask_of(spans, trim=0.0):
        m = np.zeros(S.nb, bool)
        for (a, c) in spans:
            i1, i2 = T.idx(a + trim, c - trim); m[i1:i2] = True
        return m
    other = [(b['t0'], b['t1']) for b in blks if (not b['key'].startswith('loom')) and b['key'] not in ('spont',) and not b['interrupted']]
    spont = [(b['t0'], b['t1']) for b in blks if b['key'] == 'spont' and not b['interrupted']]
    loom = [(b['t0'], b['t1']) for b in gb if not b['interrupted']]
    greys = []
    for b in ok_b: greys += [(b['t0'] + PRE_W[0], b['t0'] + PRE_W[1]), (b['t1'] + POST_W[0], b['t1'] + POST_W[1])]
    for nm, spans in (('spont block', spont), ('grey (before / after loom)', greys), ('loom blocks', loom), ('other stimulus blocks', other)):
        if not spans: continue
        m = mask_of(spans); xb = T.bp['sig']; d = dict(minutes=float(m.sum() * DT / 60), sd=float(np.nanstd(xb[m])), mean=float(np.nanmean(T.raw['sig'][m])))
        for b in ('face', 'whisk'):
            if b in T.beh_bp:
                yb = T.beh_bp[b]; ok = m & np.isfinite(xb) & np.isfinite(yb)
                if ok.sum() > 300 and np.nanstd(yb[ok]) > 0: d['r_' + b] = float(np.corrcoef(xb[ok], yb[ok])[0, 1]); d['slope_' + b] = float(np.cov(xb[ok], yb[ok])[0, 1] / np.var(yb[ok])); d['sd_' + b] = float(np.nanstd(yb[ok]))
        for b in ('face', 'whisk', 'pupil', 'loco'):
            if b in T.beh: d['mean_' + b] = float(np.nanmean(T.beh[b][m]))
        ep[nm] = d
    res['epochs'] = ep
    return res


def contrasts(results):
    """standardised contrasts over sessions"""
    rows = {}
    keys = []
    for pop in POPS:
        for m in METRICS: keys.append((pop, m))
    keys += [('res', 'mean'), ('res', 'sd')] + [('beh', b) for b in BEH] + [('coupling', 'r_face'), ('coupling', 'slope_face'), ('coupling', 'r_whisk')]
    def val(w, k):
        pop, m = k
        if not w['valid']: return np.nan
        if pop in POPS: d = w.get(pop); return np.nan if d is None else d[m]
        if pop == 'res': d = w.get('res'); return np.nan if d is None else d[m]
        if pop == 'beh': return w['beh'].get(m, np.nan)
        if pop == 'coupling': return w['coupling'].get(m, np.nan)
    for k in keys:
        loom_d, grey_d, dd, pre_first, post_last = {}, {}, {}, {}, {}
        for nm, r in results.items():
            if 'windows' not in r: continue
            W = r['windows']; allv = np.array([val(w, k) for w in W]); sc = np.nanstd(allv)
            if not np.isfinite(sc) or sc == 0: continue
            nb = max(w['block'] for w in W) + 1
            P = {(w['block'], w['kind']): val(w, k) for w in W}
            ld = [(P.get((j, 'post'), np.nan) - P.get((j, 'pre'), np.nan)) / sc for j in range(nb)]
            gd = {j: (P.get((j + 1, 'pre'), np.nan) - P.get((j, 'post'), np.nan)) / sc for j in range(nb - 1) if r['peri_blocks'][j + 1]['k'] == r['peri_blocks'][j]['k'] + 1}
            ddj = [ld[j] - gd[j] for j in gd if np.isfinite(ld[j]) and np.isfinite(gd[j])]
            if np.isfinite(ld).any(): loom_d[nm] = float(np.nanmean(ld))
            if gd and np.isfinite(list(gd.values())).any(): grey_d[nm] = float(np.nanmean(list(gd.values())))
            if ddj: dd[nm] = float(np.mean(ddj))
            if np.isfinite(P.get((0, 'pre'), np.nan)) and np.isfinite(P.get((nb - 1, 'post'), np.nan)): pre_first[nm] = float((P[(nb - 1, 'post')] - P[(0, 'pre')]) / sc)
        rows[k] = dict(loom=stat_group(loom_d), grey=stat_group(grey_d), diff=stat_group(dd), first_last=stat_group(pre_first))
    return rows


def write_csv(results, path):
    rows = []
    for nm, r in results.items():
        for w in r.get('windows', []):
            if not w['valid']: continue
            row = dict(session=nm, block=w['blk_k'] + 1, kind=w['kind'], direction=w['loom_dir'], t_start=round(w['t0'], 1), t_end=round(w['t1'], 1))
            for pop in POPS:
                d = w.get(pop)
                if d:
                    for m in METRICS: row['%s_%s' % (pop, m)] = round(d[m], 4)
            if 'res' in w: row['res_mean'] = round(w['res']['mean'], 4); row['res_sd'] = round(w['res']['sd'], 4)
            for b, v in w['beh'].items(): row['z_' + b] = round(v, 3)
            for k, v in w['coupling'].items(): row[k] = round(v, 3)
            rows.append(row)
    cols = []
    for r in rows:
        for k in r:
            if k not in cols: cols.append(k)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); [w.writerow(r) for r in rows]


if __name__ == '__main__':
    names = sys.argv[1:] or GREY_SESSIONS; results = {}
    for nm in names:
        t0 = time.time()
        try: results[nm] = analyse(nm)
        except Exception as e:
            import traceback; traceback.print_exc(); results[nm] = dict(name=nm, why='error: %s' % e)
        print('  -> %s (%.0f s)' % (results[nm].get('why', 'ok'), time.time() - t0), flush=True)
    con = contrasts(results)
    pickle.dump(dict(results=results, contrasts=con), open(os.path.join(FIG, 'grey_results.pkl'), 'wb'))
    write_csv(results, os.path.join(FIG, 'grey_windows.csv'))
    print('\nstandardised contrasts (mean over sessions; sign-flip p; bootstrap CI):')
    for k, v in con.items():
        s = lambda d: 'n/a' if d is None else '%+.2f [%+.2f, %+.2f] p=%.3f n=%d' % (d['mean'], d['ci'][0], d['ci'][1], d['p'], d['n'])
        print('  %-22s loom %s | grey %s | diff %s' % ('%s/%s' % k, s(v['loom']), s(v['grey']), s(v['diff'])))
    print('GREY_DONE')
