"""Cross-session summary of the grouping pilot: python ag_summary.py [session ...]  (reads <session>_groups.pkl and <session>_validation_pc1.pkl from the cache)"""
import os, sys, pickle, numpy as np
import ag_data as ad

DEFAULT = ['AnimalRB24_261021_1630', 'AnimalRB22_261027_1138', 'AnimalRB23_261023_1642', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520']


def load(name, kind):
    f = os.path.join(ad.CACHE, '%s_%s.pkl' % (name, kind))
    if not os.path.exists(f): return None
    with open(f, 'rb') as fh: return pickle.load(fh)


def main(names):
    print('%-26s %4s %5s %5s | %-22s | %-30s' % ('session', 'acc', 'units', 'multi', 'tiers (ROIs)', 'attach rejected: high/med/low/none, orphans'))
    for nm in names:
        R = load(nm, 'groups'); V = load(nm, 'validation_pc1')
        if R is None: print(nm, 'no groups'); continue
        tiers = list(R['tier'].values()); multi = [g for g in R['units'] if len(g) >= 2]; conf = R['attach']['conf']; src = R['src']
        cnt = {c: sum(1 for j, v in conf.items() if src[j] == 1 and v == c) for c in ('high', 'medium', 'low', 'none')}
        print('%-26s %4d %5d %5d | %-22s | %d/%d/%d/%d, orphans %d' % (nm, (src == 0).sum(), len(R['units']), len(multi), ' '.join('%s%d' % (t, tiers.count(t)) for t in 'ABCD'), cnt['high'], cnt['medium'], cnt['low'], cnt['none'], sum(R['attach']['orphan'].values())))
        print('    sizes of multi-member units %s | tier A (averaged) units: %d | why-not-A: %s' % (sorted([len(g) for g in multi], reverse=True),
              sum(1 for d in R['unit_info'] if d.get('avg_ok')), {w: sum(1 for d in R['unit_info'] if w in d.get('why', [])) for w in ('untestable_n2', 'untestable_n3', 'heterogeneous', 'unstable', 'v3_fail')}))
        print('    noise reference: pool %d, standardised mean difference to the members (F0, npix, Fneu/F) %s' % (R['pool_size'], np.round(R['noise_smd'], 2)))
        if V is None: continue
        n0 = V['null']; print('    null runs: %s | real data %d edges' % (' | '.join('%s: %d of %d edges, %d units' % (k, v['edges'], v['candidate_pairs'], v['multi_units']) for k, v in n0.items() if k != 'real data'), n0['real data']['edges']))
        H = V['hetero']; fm = {}
        for h in H: fm.setdefault(h['mix'], []).append((h['rho'], h['false_merge'], h['flagged']))
        for mx in ('2+2', '3+3'):
            if mx in fm: print('    sim %s fragments: false merge by rho %s | flagged heterogeneous %s' % (mx, ' '.join('%.2f:%.2f' % (r, f) for r, f, _ in fm[mx]), ' '.join('%.2f:%.2f' % (r, g) for r, _, g in fm[mx] if np.isfinite(g))))
        ho = V['held_out']; rows = [r for h in ho for r in h['units']]
        if rows:
            res = np.array([r['res_r_heldout'] for r in rows]); ctl = np.array([r['res_r_ctrl_p95'] for r in rows]); gain = np.array([r['gain'] for r in rows]); cg = np.array([r['ctrl_gain_p95'] for r in rows]); cm = np.array([r['ctrl_gain_median'] for r in rows])
            print('    held-out residual r among members: median %.2f; above the near-ROI control p95 in %d of %d units | SNR gain median %.2f vs ideal %.2f, above control median in %d of %d, above control p95 in %d' % (
                np.median(res), (res > ctl).sum(), len(res), np.median(gain), np.median([r['ideal'] for r in rows]), (gain > cm).sum(), len(gain), (gain > cg).sum()))
        for h in ho:
            I = h['inflation']
            print('    fold %d  PC1 share: ROIs %.2f | units in-sample %.2f, held-out %.2f (ROI level in the held-out half %.2f) | random merges held-out %.2f ;  mean pair r ROIs %.2f | units in %.2f out %.2f | random %.2f' % (
                h['fold'], I['pc1_roi_A'], I['pc1_units_in_sample'], I['pc1_units_held_out'], I['pc1_roi_B'], I['pc1_random_merge_held_out'], I['r_roi_B'], I['r_units_in_sample'], I['r_units_held_out'], I['r_random_merge_held_out']))
        cs = {}
        for h in ho:
            for b, v in h['coupling_spread'].items(): cs.setdefault(b, []).append((v['rho'], v['rho_with_snr']))
        print('    members differ reproducibly in behaviour coupling (rho half1 vs half2 | rho with member SNR): %s' % ' ; '.join('%s %.2f | %.2f' % (b, np.mean([x[0] for x in v]), np.mean([x[1] for x in v])) for b, v in cs.items()))
        nc = V.get('noise_compare', [])
        if nc:
            for m, lab in (('abs_r_beh', '|r| behaviour'), ('snr', 'SNR'), ('ev_per_min', 'events/min')):
                above = sum(1 for r in nc if r[m]['frac_noise_ge_unit'] < 0.05); print('    units vs mean of n matched noise blobs, %s: unit above 95%% of noise draws in %d of %d units' % (lab, above, len(nc)))


if __name__ == '__main__':
    main(sys.argv[1:] or DEFAULT)
