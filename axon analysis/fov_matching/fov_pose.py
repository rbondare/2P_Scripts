"""Step 3b: global placement of the sessions in one frame from the pairwise candidates (cache/pairwise.pkl) -> plots/fov_matching/fov_pose.csv + cache/poses.pkl

The global maximum of a pair's correlation map is sometimes a wrong peak (sparse axons, partial overlap). So a session is only placed where the candidate positions
predicted from SEVERAL already-placed sessions agree (consensus within TOL_UM), or where a single very strong, well-overlapping, above-null peak exists.
Growth starts from the best-connected session; sessions that never reach consensus stay UNPLACED (different FOV, or too weak) - they are not forced into the frame.
After growth all placed pairs vote again in a robust least-squares refinement (edges = the candidate peak closest to the current prediction; Huber-style reweighting).
Pose P_s = position (um) of session s's pixel (0, 0) centre in the reference frame (reference = the session placed first, at (0, 0)).
"""
import os, sys, pickle, csv, numpy as np
from fov_common import *

TOL_UM = 2.5
R_STRONG, EXC_STRONG, OV_STRONG = 0.75, 0.25, 0.5      # single-neighbour acceptance: peak r, excess over the phase-randomised null, overlap fraction
R_MULTI, EXC_MULTI = 0.45, 0.05                        # per-neighbour minimum in a >= 2-neighbour consensus


def load():
    d = pickle.load(open(os.path.join(CACHE, 'pairwise.pkl'), 'rb'))
    return d['res']


def peaks_between(res, a, b):
    """list of (t, peak) with t = position of b's origin in a's frame (for (a, b) in either stored order)"""
    if (a, b) in res: return [(np.array(p['t']), p) for p in res[(a, b)]['peaks']], res[(a, b)]
    r = res[(b, a)]; return [(-np.array(p['t']), p) for p in r['peaks']], r


def votes(res, poses, s):
    """candidate positions of s predicted by every placed member and peak: list of (pos(2), member, weight, r, excess, overlap)"""
    out = []
    for m, Pm in poses.items():
        pk, info = peaks_between(res, m, s)
        for t, p in pk:
            exc = p['r'] - info['null_max']
            out.append((Pm + t, m, max(exc, 0.0) * max(p['overlap'], 0.2), p['r'], exc, p['overlap']))
    return out


def best_cluster(v):
    """greedy cluster of the votes within TOL_UM; support = sum over members of the best vote of that member in the cluster"""
    best = None
    for c in sorted(v, key=lambda x: -x[2])[:60]:
        cen = c[0]
        mem = {}
        for x in v:
            if np.hypot(*(x[0] - cen)) <= TOL_UM and (x[1] not in mem or x[2] > mem[x[1]][2]): mem[x[1]] = x
        # refine the centre as the weighted mean
        w = np.array([x[2] + 1e-6 for x in mem.values()]); pos = np.array([x[0] for x in mem.values()])
        cen2 = (pos * w[:, None]).sum(0) / w.sum()
        ok = [x for x in mem.values() if x[3] >= R_MULTI and x[4] >= EXC_MULTI]
        support = float(sum(x[2] for x in ok))
        strong = [x for x in mem.values() if x[3] >= R_STRONG and x[4] >= EXC_STRONG and x[5] >= OV_STRONG]
        cand = dict(pos=cen2, support=support, n=len(ok), members=list(mem.keys()), strong=len(strong), rbest=max(x[3] for x in mem.values()))
        if best is None or (cand['strong'] > 0, cand['n'], cand['support']) > (best['strong'] > 0, best['n'], best['support']): best = cand
    return best


def anchor_degree(res):
    deg = {s: 0.0 for s in SESSIONS}
    for (a, b), info in res.items():
        p = info['peaks'][0]
        if p['r'] >= R_STRONG and p['r'] - info['null_max'] >= EXC_STRONG and p['overlap'] >= OV_STRONG:
            deg[a] += p['r']; deg[b] += p['r']
    return deg


def grow(res):
    deg = anchor_degree(res); ref = max(deg, key=deg.get)
    poses = {ref: np.zeros(2)}; log = [(ref, 'reference', 0, 0.0)]
    while True:
        cands = []
        for s in SESSIONS:
            if s in poses: continue
            v = votes(res, poses, s)
            if not v: continue
            c = best_cluster(v)
            if c is None: continue
            if c['strong'] > 0 or c['n'] >= 2: cands.append((c['strong'] > 0, c['n'], c['support'], s, c))
        if not cands: break
        cands.sort(key=lambda x: (x[1], x[2]), reverse=True)
        _, n, sup, s, c = cands[0]
        poses[s] = c['pos']; log.append((s, 'strong' if c['strong'] > 0 and c['n'] < 2 else 'consensus', c['n'], c['rbest']))
    return poses, log


def refine(res, poses, iters=4):
    names = list(poses); idx = {s: i for i, s in enumerate(names)}
    P = {s: poses[s].copy() for s in names}
    edges = []
    for a in names:
        for b in names:
            if a >= b: continue
            pk, info = peaks_between(res, a, b)
            pred = P[b] - P[a]
            best = min(pk, key=lambda x: np.hypot(*(x[0] - pred)))
            if np.hypot(*(best[0] - pred)) <= TOL_UM and best[1]['r'] >= R_MULTI:
                exc = best[1]['r'] - info['null_max']
                edges.append((a, b, best[0], max(exc, 0.02) * max(best[1]['overlap'], 0.2), best[1]['r']))
    for it in range(iters):
        n = len(names); A = np.zeros((2 * len(edges) + 2, 2 * n)); y = np.zeros(2 * len(edges) + 2); W = np.zeros(2 * len(edges) + 2)
        # re-pick the peak closest to the current prediction each iteration and apply robust weights
        rows = []
        for (a, b, t, w, r) in edges:
            pk, info = peaks_between(res, a, b); pred = P[b] - P[a]
            tt = min(pk, key=lambda x: np.hypot(*(x[0] - pred)))[0]
            resid = np.hypot(*(tt - pred)); hub = 1.0 if resid <= 1.0 else 1.0 / resid
            rows.append((a, b, tt, w * hub))
        for k, (a, b, tt, w) in enumerate(rows):
            for d in range(2):
                A[2 * k + d, 2 * idx[b] + d] = 1; A[2 * k + d, 2 * idx[a] + d] = -1; y[2 * k + d] = tt[d]; W[2 * k + d] = w
        # gauge: the first session stays at its place
        r0 = names[0]
        for d in range(2): A[2 * len(edges) + d, 2 * idx[r0] + d] = 1; y[2 * len(edges) + d] = P[r0][d]; W[2 * len(edges) + d] = 1e3
        sw = np.sqrt(W)[:, None]
        sol, *_ = np.linalg.lstsq(A * sw, y * sw[:, 0], rcond=None)
        P = {s: sol[2 * idx[s]:2 * idx[s] + 2] for s in names}
    # residuals
    res_edges = []
    for (a, b, t, w, r) in edges:
        pk, info = peaks_between(res, a, b); pred = P[b] - P[a]
        tt = min(pk, key=lambda x: np.hypot(*(x[0] - pred)))[0]
        res_edges.append((a, b, float(np.hypot(*(tt - pred))), r))
    return P, res_edges


def main():
    res = load()
    poses, log = grow(res)
    print('growth order:')
    for s, how, n, r in log: print('  %s %-9s n=%d rbest=%.2f' % (short(s), how, n, r))
    P, resid = refine(res, poses)
    ref = log[0][0]; off = P[ref].copy()
    P = {s: p - off for s, p in P.items()}
    rows = []
    for s in SESSIONS:
        if s in P:
            e = [x for x in resid if s in (x[0], x[1])]
            rr = [x[2] for x in e]; rv = [x[3] for x in e]
            conf = 'high' if sum(v >= 0.6 for v in rv) >= 3 else ('medium' if sum(v >= 0.45 for v in rv) >= 2 else 'low')
            rows.append(dict(session=s, short=short(s), placed='yes', y_um=round(P[s][0], 2), x_um=round(P[s][1], 2), n_edges=len(e), rms_resid_um=round(float(np.sqrt(np.mean(np.square(rr)))) if rr else np.nan, 2),
                             max_resid_um=round(max(rr), 2) if rr else np.nan, best_r=round(max(rv), 2) if rv else np.nan, median_r=round(float(np.median(rv)), 2) if rv else np.nan, confidence=conf))
        else:
            rows.append(dict(session=s, short=short(s), placed='no', y_um='', x_um='', n_edges=0, rms_resid_um='', max_resid_um='', best_r='', median_r='', confidence='unplaced'))
    p = os.path.join(OUT, 'fov_pose.csv')
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    pickle.dump(dict(P=P, resid=resid, log=log, ref=ref), open(os.path.join(CACHE, 'poses.pkl'), 'wb'))
    print('\nplacement (reference %s):' % short(ref))
    for r in rows: print('  %s %-8s y=%7s x=%7s edges=%2d rms=%s max=%s best_r=%s med_r=%s %s' % (r['short'], r['placed'], r['y_um'], r['x_um'], r['n_edges'], r['rms_resid_um'], r['max_resid_um'], r['best_r'], r['median_r'], r['confidence']))
    print('edges used: %d, residual > 2 um: %d' % (len(resid), sum(x[2] > 2 for x in resid)))


if __name__ == '__main__':
    main()
