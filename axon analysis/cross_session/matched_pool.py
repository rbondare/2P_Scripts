"""Pooling of the per-session results of matched_behaviour.py (session = unit, animal = cluster, unflagged sessions only; flagged sessions are listed but never pooled).

load(cache)                  -> {session: result}
pooled(res, key, sessions)   -> pooled matched effect (mean of animal means, hierarchical bootstrap CI), the same minus the behaviour-only benchmark ("excess"),
                                p from the pooled shift null and from the pooled benchmark (independent per-session draws, animal-weighted mean of the draws)
group_diff(res, key, A, B)   -> difference between two groups of sessions (e.g. task - passive) with an exact label-permutation p over sessions
table(res)                   -> flat rows for the csv
"""
import os, sys, glob, pickle, csv, itertools, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'paper_figures'))
import matched_behaviour as mb
from fig_common import hboot, animal, is_flagged, FLAG

TASK = ['AnimalRB24_261029_1752', 'AnimalRB24_261030_1744', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261106_1317', 'AnimalRB23_261029_1704', 'AnimalRB24_261109_1506']
PASSIVE = ['AnimalRB24_261027_1456', 'AnimalRB24_261104_1651', 'AnimalRB23_261027_1423', 'AnimalRB24_261106_1529']
VISUAL = ['AnimalRB24_261021_1630', 'AnimalRB24_261023_1538', 'AnimalRB24_261027_1534', 'AnimalRB23_261023_1642', 'AnimalRB22_261027_1138', 'AnimalRB22_261023_1427',
          'AnimalRB24_261023_1527', 'AnimalRB24_261104_1640', 'AnimalRB24_261106_1604', 'AnimalRB23_261027_1306']


def load(cache='cache'):
    res = {}
    for f in sorted(glob.glob(os.path.join(mb.OUT, cache, '*.pkl'))):
        r = pickle.load(open(f, 'rb')); res[r['name']] = r
    return res


def get(res, s, key):
    c = res.get(s, {}).get('contrasts', {}).get(key)
    if c and c.get('obs') and c.get('stats') and np.isfinite(c['obs'].get('matched', np.nan)): return c
    return None


def amean(vals, names):
    """mean of the per-animal means (the pooled estimate of fig_common.hboot)"""
    vals = np.asarray(vals, float); an = np.array([animal(n) for n in names]); return float(np.mean([vals[an == a].mean() for a in sorted(set(an))]))


def _scale(res, s, norm): return res[s]['sd_P'] if norm else 1.0


MIN_COV, MIN_MATCHED = 0.5, 10        # a session x contrast is 'testable' only if >= 50% of its event windows found a behaviour match and >= 10 pairs exist


def testable(c):
    return c is not None and c['obs']['n_matched'] >= MIN_MATCHED and c['obs']['n_matched'] / max(c['obs']['n_focal'], 1) >= MIN_COV


def pooled(res, key, sessions, field='matched', norm=False, include_flagged=False, n_boot=5000, only_testable=True, min_sessions=2):
    ss = [s for s in sessions if get(res, s, key) is not None and (include_flagged or not is_flagged(s)) and (testable(get(res, s, key)) or not only_testable)
          and np.isfinite(get(res, s, key)['stats'][field]['p'])]          # a degenerate shift null (< 8 distinct draws) gives p = nan: that session is not pooled
    if len(ss) < min_sessions: return None
    C = [get(res, s, key) for s in ss]; sc = np.array([_scale(res, s, norm) for s in ss])
    vals = np.array([c['obs'][field] for c in C]) / sc; out = dict(sessions=ss, vals=vals, key=key, field=field, norm=norm)
    h = hboot(vals, ss, n_boot=n_boot); out.update(est=h['est'], lo=h['lo'], hi=h['hi'], k_pos=h['k_pos'], n_animals=h['n_animals'], per_animal=h['per_animal'])
    out['naive'] = amean([c['obs']['naive'] for c in C] / sc, ss); out['reg'] = amean([c['obs']['reg'] for c in C] / sc, ss)
    # pooled shift null: independent draws per session
    K = min(len(c['null'][field]) for c in C)
    if K >= 50:
        nm = np.array([c['stats'][field]['nmean'] for c in C]) / sc; N = np.stack([c['null'][field][:K] / s_ for c, s_ in zip(C, sc)], 1) - nm[None, :]
        nul = np.array([amean(N[j], ss) for j in range(K)]); oc = amean(vals - nm, ss)
        out.update(p_shift=float((1 + (np.abs(nul) >= abs(oc)).sum()) / (1 + K)), eff_shift=oc, null_sd=float(nul.std()), null_lo=float(np.percentile(nul, 2.5)), null_hi=float(np.percentile(nul, 97.5)))
    # pooled behaviour-only benchmark
    if all('bench' in c and field in c.get('stats_sim', {}) for c in C):
        Kb = min(len(c['bench'][field]) for c in C)
        if Kb >= 20:
            bm = np.array([c['stats_sim'][field]['nmean'] for c in C]) / sc; B = np.stack([c['bench'][field][:Kb] / s_ for c, s_ in zip(C, sc)], 1) - bm[None, :]
            nulb = np.array([amean(B[j], ss) for j in range(Kb)]); exv = vals - bm; ob = amean(exv, ss); hb = hboot(exv, ss, n_boot=n_boot)
            out.update(bench=amean(bm, ss), excess=ob, excess_lo=hb['lo'], excess_hi=hb['hi'], excess_vals=exv, p_sim=float((1 + (np.abs(nulb) >= abs(ob)).sum()) / (1 + Kb)),
                       bench_sd=float(nulb.std()), excess_kpos=int((exv > 0).sum()))
    out['coverage'] = float(np.mean([c['n_matched'] / max(c['n_focal'], 1) for c in C]))
    return out


def group_diff(res, key, A, B, field='excess', norm=False):
    """mean(A) - mean(B) of the per-session value (matched effect, or matched minus behaviour-only benchmark) with an exact label-permutation p over sessions (unflagged only)"""
    def vals(L):
        o = []
        for s in L:
            c = get(res, s, key)
            if c is None or is_flagged(s) or not testable(c): continue
            v = c['obs']['matched'] - (c['stats_sim']['matched']['nmean'] if field == 'excess' and 'matched' in c.get('stats_sim', {}) else 0.0); o.append((s, v / _scale(res, s, norm)))
        return o
    a, b = vals(A), vals(B)
    if len(a) < 2 or len(b) < 2: return None
    allv = np.array([v for _, v in a + b]); na = len(a); obs = allv[:na].mean() - allv[na:].mean(); cnt = 0; tot = 0
    for comb in itertools.combinations(range(len(allv)), na):
        m = np.zeros(len(allv), bool); m[list(comb)] = True; d = allv[m].mean() - allv[~m].mean(); tot += 1; cnt += abs(d) >= abs(obs) - 1e-12
    return dict(diff=float(obs), p=float(cnt / tot), n_a=len(a), n_b=len(b), a=a, b=b, min_p=1.0 / tot)


def table(res):
    rows = []
    for s in sorted(res):
        r = res[s]
        for k, c in r.get('contrasts', {}).items():
            o = c.get('obs');
            if not o: continue
            st = c.get('stats', {}); ss = c.get('stats_sim', {})
            g = lambda d, kk, f: d.get(kk, {}).get(f, np.nan)
            rows.append(dict(session=s, flagged=int(is_flagged(s)), task=int(s in mb.TASK_SESSIONS), contrast=k, label=mb.LABEL.get(k, k), n_events=c['n_events'], n_valid=o['n_focal'], n_matched=o['n_matched'],
                             coverage=o['n_matched'] / max(o['n_focal'], 1), n_pool=o['n_ref'], sd_P=r.get('sd_P', np.nan), naive=o['naive'], reg=o['reg'], matched=o['matched'], adj=o['adj'],
                             shift_nmean=g(st, 'matched', 'nmean'), shift_nsd=g(st, 'matched', 'nsd'), shift_z=g(st, 'matched', 'z'), shift_p=g(st, 'matched', 'p'),
                             bench_mean=g(ss, 'matched', 'nmean'), bench_sd=g(ss, 'matched', 'nsd'), sim_z=g(ss, 'matched', 'z'), sim_p=g(ss, 'matched', 'p'),
                             excess=o['matched'] - g(ss, 'matched', 'nmean'), smd_before_max=float(np.nanmax(np.abs(o['smd_before']))) if np.size(o.get('smd_before')) else np.nan,
                             smd_after_max=float(np.nanmax(np.abs(o['smd_after']))) if np.size(o.get('smd_after')) else np.nan, behaviours=','.join(r.get('behaviours', []))))
    return rows


def write_csv(res, path):
    rows = table(res)
    if not rows: return
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
        for r in rows: w.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
