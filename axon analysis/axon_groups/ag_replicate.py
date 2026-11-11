"""Axon grouping layer, replication check across two sessions of the SAME field of view (RB24_261105_1423 / 1520).   python ag_replicate.py [sessionA sessionB]

The only independent test of the grouping and of the attachments: ROIs are matched across sessions through the fov_matching master axons (Recordings_DB rois.csv, `master_axon`, accepted AND rejected ROIs).
 (a) units: pairs of master axons that are an accepted ROI in both sessions and within D_BUNDLE in both: is "same unit in A" predictive of "same unit in B"?
 (b) attachments: a rejected ROI of A attached to a unit; its master axon in B: is it a member of / attached to the unit that holds the matching partners in B? By confidence class.
"""
import os, sys, csv, pickle, numpy as np
import ag_data as ad, ag_group as ag

RDB = r'Z:\Group Members\Rima\Recordings_DB\rois.csv'


def master_of(name):
    m = {}
    with open(RDB, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['session'] == name and r['master_axon'] != '': m[int(r['roi'])] = int(float(r['master_axon']))
    return m


def load(name):
    with open(os.path.join(ad.CACHE, '%s_groups.pkl' % name), 'rb') as fh: R = pickle.load(fh)
    return R, ad.load_universe(name)


def unit_row_of(R, U):
    """row of the Universe -> unit id for accepted ROIs; for rejected ROIs: (unit id of the best attachment, confidence)"""
    acc = {int(U.idx[i]): int(R['unit_of'][i]) for i in range(U.N) if U.src[i] == 0}
    att = {int(U.idx[j]): (b['unit'], R['attach']['conf'][j]) for j, b in R['attach']['best'].items()}
    return acc, att


def replicate(A, B, verbose=True):
    RA, UA = load(A); RB, UB = load(B); mA, mB = master_of(A), master_of(B); accA, attA = unit_row_of(RA, UA); accB, attB = unit_row_of(RB, UB)
    idxA = {int(UA.idx[i]): i for i in range(UA.N)}; idxB = {int(UB.idx[i]): i for i in range(UB.N)}
    revB = {}; [revB.setdefault(m, []).append(r) for r, m in mB.items()]
    # (a) pairs of accepted ROIs matched to accepted ROIs of the other session
    pairs = []
    ids = [r for r in accA if r in mA and any(x in accB for x in revB.get(mA[r], []))]
    gapA, gapB = UA.gap_um(60.0), UB.gap_um(60.0); to_b = {r: [x for x in revB[mA[r]] if x in accB][0] for r in ids}
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            i, j = ids[a], ids[b]; ia, ja = idxA[i], idxA[j]; ib, jb = idxB[to_b[i]], idxB[to_b[j]]
            if ib == jb or gapA[ia, ja] > ag.D_BUNDLE or gapB[ib, jb] > ag.D_BUNDLE: continue
            pairs.append((accA[i] == accA[j], accB[to_b[i]] == accB[to_b[j]], RA['cluster_of'][ia] == RA['cluster_of'][ja], RB['cluster_of'][ib] == RB['cluster_of'][jb]))
    P = np.array(pairs, bool) if pairs else np.zeros((0, 4), bool); out = dict(n_pairs=len(P))
    if len(P):
        for lab, a_, b_ in (('unit', 0, 1), ('cluster', 2, 3)):
            sa, sb = P[:, a_], P[:, b_]; out[lab] = dict(n_same_A=int(sa.sum()), p_same_B_given_same_A=float(sb[sa].mean()) if sa.any() else np.nan, n_diff_A=int((~sa).sum()), p_same_B_given_diff_A=float(sb[~sa].mean()) if (~sa).any() else np.nan)
    # (b) attachments of rejected ROIs of A, judged in B
    res = {c: dict(n=0, member_or_attached_same_unit=0, other_unit=0, unattached=0, no_match=0) for c in ('high', 'medium', 'low', 'none')}
    for r, (u, conf) in attA.items():
        if r not in mA: continue
        partners = [x for x in RA['units'][u]]; pm = [mA[int(UA.idx[p])] for p in partners if int(UA.idx[p]) in mA]
        pb = [x for m in pm for x in revB.get(m, []) if x in accB]; unitsB = {accB[x] for x in pb}; res[conf]['n'] += 1
        cand = [x for x in revB.get(mA[r], []) if x != r or True]
        if not cand or not unitsB: res[conf]['no_match'] += 1; continue
        x = cand[0]
        if x in accB: ub = accB[x]; res[conf]['member_or_attached_same_unit' if ub in unitsB else 'other_unit'] += 1
        elif x in attB and attB[x][1] != 'none': res[conf]['member_or_attached_same_unit' if attB[x][0] in unitsB else 'other_unit'] += 1
        else: res[conf]['unattached'] += 1
    out['attach'] = res
    if verbose:
        print('%s -> %s: %d ROI pairs usable (both accepted in both sessions, within %d um)' % (A[-9:], B[-9:], len(P), ag.D_BUNDLE))
        for lab in ('unit', 'cluster'):
            if lab in out: o = out[lab]; print('  same %s in A: %d pairs -> same in B %.2f | different in A: %d pairs -> same in B %.2f' % (lab, o['n_same_A'], o['p_same_B_given_same_A'], o['n_diff_A'], o['p_same_B_given_diff_A']))
        for c, d in res.items():
            n = d['n'] - d['no_match']
            if n: print('  attachments %-6s (n %3d with a partner in B): same unit in B %.2f | other unit %.2f | unattached %.2f' % (c, n, d['member_or_attached_same_unit'] / n, d['other_unit'] / n, d['unattached'] / n))
    return out


if __name__ == '__main__':
    a, b = (sys.argv[1:3] if len(sys.argv) > 2 else ('AnimalRB24_261105_1423', 'AnimalRB24_261105_1520'))
    o1 = replicate(a, b); o2 = replicate(b, a)
    with open(os.path.join(ad.CACHE, 'replication_%s_%s.pkl' % (a[-4:], b[-4:])), 'wb') as fh: pickle.dump(dict(ab=o1, ba=o2), fh)
