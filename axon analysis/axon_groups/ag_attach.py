"""Axon grouping layer, step 4: attach REJECTED ROIs to the units of the accepted ROIs, with a graded confidence.

Attachment is an ANNOTATION: it never changes a unit trace (design principle: ROI level stays primary). For every rejected ROI j and every unit u within D_BUNDLE um:
  trace of u built WITHOUT any member whose footprint overlaps j (no shared pixels, so shared photons cannot create the correlation),
  residual r (same removal as the grouping), z / p against the pooled shift null, Benjamini-Hochberg over all (j, u) tests, r in each half of the session, footprint gap.
Confidence: high   q < 1e-3, r >= 0.4, r >= 0.25 in BOTH halves, gap <= D_FRAG
            medium q < 1e-2, r >= 0.3, r >= 0.15 in both halves, gap <= D_BUNDLE
            low    q < 5e-2, r >= 0.2
            none   otherwise
Hidden merge constituents (src == 2) are scored the same way as a sanity check (are they recovered into the unit of their parent?). This is circular for the merge used correlation after PC1
removal too, so it demonstrates that the machinery works, not that attachments are correct.
`orphan`: rejected, not attached at any level, activity-like on its own (V3 transient test or skew > 0.5) and not touching an accepted ROI: candidate for review (review_queue.py).
"""
import numpy as np
import ag_pairs as ap, ag_group as ag, ag_noise


def _half_r(P, x, j, halves):
    out = []
    for h in halves:
        a = P.Zr[j, h]; b = x[h]; sa, sb = a.std(), b.std(); out.append(float(((a - a.mean()) @ (b - b.mean())) / len(h) / (sa * sb)) if sa > 1e-9 and sb > 1e-9 else np.nan)
    return out


def attach(P, G, d_um=ag.D_BUNDLE, verbose=True):
    U = P.U; gap = U.gap_um(d_um); halves = G.halves; v = P.v
    rows = np.where(U.src > 0)[0]; units = G.units; tests = []
    for j in rows:
        for k, mem in enumerate(units):
            g = gap[j, mem]
            if not np.isfinite(g).any() or g.min() > d_um: continue
            keep = [m for m, gg in zip(mem, g) if gg > 0.01]                                   # members overlapping j are left out of the unit trace
            if not keep: continue
            x = ag._ztrace(P, keep); r, z, p = P.vs_trace(x, rows=[j]); tests.append(dict(roi=int(j), unit=k, gap=float(g.min()), r=float(r[0]), z=float(z[0]), p=float(p[0]), n_used=len(keep), n_unit=len(mem), half=_half_r(P, x, j, halves)))
    if not tests: return dict(tests=[], best={}, conf={int(j): 'none' for j in rows}, orphan={int(j): False for j in rows}, constituent_recovery=[])
    q = ap.bh_fdr([t['p'] for t in tests])
    for t, qq in zip(tests, q): t['q'] = float(qq)
    best = {}
    for t in tests:
        if t['roi'] not in best or t['z'] > best[t['roi']]['z']: best[t['roi']] = t
    pool = {}
    for j, t in best.items():
        h = t['half']; s1, s2 = (h[0], h[1]) if np.all(np.isfinite(h)) else (np.nan, np.nan)
        if t['q'] < 1e-3 and t['r'] >= 0.4 and s1 >= 0.25 and s2 >= 0.25 and t['gap'] <= ag.D_FRAG: c = 'high'
        elif t['q'] < 1e-2 and t['r'] >= 0.3 and s1 >= 0.15 and s2 >= 0.15 and t['gap'] <= ag.D_BUNDLE: c = 'medium'
        elif t['q'] < 5e-2 and t['r'] >= 0.2: c = 'low'
        else: c = 'none'
        pool[j] = c
    geom = ag_noise.footprint_geometry(U); passes = U.v3()['passes']; sk = U.get_skew()
    conf = {int(j): pool.get(int(j), 'none') for j in rows}
    orphan = {int(j): (U.src[j] == 1 and conf[int(j)] == 'none' and (passes[j] or sk[j] > 0.5) and geom['adjacent'][j] == 0) for j in rows}
    out = dict(tests=tests, best=best, conf=conf, orphan=orphan)
    if verbose:
        cnt = lambda s: {c: sum(1 for j in rows if U.src[j] == s and conf[int(j)] == c) for c in ('high', 'medium', 'low', 'none')}
        print('  attach: rejected ROIs %s | hidden constituents (overlap their parent, so the parent is left out of the unit trace; descriptive only) %s | orphans %d' % (cnt(1), cnt(2), sum(orphan.values())))
    return out
