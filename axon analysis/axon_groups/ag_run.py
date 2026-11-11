"""Axon grouping layer, runner.   python ag_run.py <session> [--level pc1] [--validate] [--refresh]

Per session (analysis layer only; nothing on Z: is written, iscell.npy / Fall.mat / aggregated files untouched):
  1 universe (ag_data)  2 pair statistics (ag_pairs)  3 homogeneity threshold from simulated noisy copies with the session's real noise (ag_validate, cached)
  4 grouping of the accepted ROIs (ag_group)  5 attachment of rejected ROIs (ag_attach)  6 noise reference (ag_noise)
  outputs  C:\\Users\\rbondarenko\\2P_cache\\axon_groups\\<session>_groups.pkl,  <session>_roi_groups.csv,  <session>_units.npz (unit dF/F on the 5 Hz grid, MATLAB-readable)
  --validate adds the report card (ag_validate.null_runs / sim_hetero, ag_report.held_out) -> <session>_validation.pkl
"""
import os, sys, csv, pickle, time, numpy as np
import ag_data as ad, ag_pairs as ap, ag_group as ag, ag_attach as aa, ag_noise as an, ag_validate as av, ag_report as ar

CACHE = ad.CACHE


def homog_threshold(U, P, refresh=False, n_rep=100):
    f = os.path.join(CACHE, '%s_homog_%s_%d.pkl' % (U.name, P.level, n_rep))
    if os.path.exists(f) and not refresh:
        with open(f, 'rb') as fh: d = pickle.load(fh)
    else:
        thr, table, rtab, succ = av.calibrate_homogeneity(U, P, n_rep=n_rep); d = dict(table=table, rtab=rtab, success=succ)
        with open(f, 'wb') as fh: pickle.dump(d, fh)
    ns = sorted(d['table']); vals = {n: float(np.quantile(d['table'][n], 0.95)) for n in ns}
    def thr(n): return vals[min(max(n, ns[0]), ns[-1])]
    thr.table = vals; return thr, d


def run(name, level='pc1', validate=False, refresh=False, verbose=True, n_rep=100):
    t0 = time.time(); U = ad.load_universe(name, refresh=refresh); P = ap.PairStats(U, level, verbose=verbose)
    thr, hd = homog_threshold(U, P, refresh=refresh, n_rep=n_rep)
    if verbose: print('  rank-1 thresholds (95th pct of simulated homogeneous units): %s' % {n: round(v, 2) for n, v in thr.table.items()})
    G = ag.Grouping(P, homog_thr=thr, verbose=verbose)
    A = aa.attach(P, G, verbose=verbose)
    attached = np.array([A['conf'].get(int(j), 'none') != 'none' for j in range(U.N)]) if A.get('conf') else np.zeros(U.N, bool)
    pool = an.noise_pool(U, exclude=attached)
    members = [i for g in G.units if len(g) >= 2 for i in g] or list(np.where(U.src == 0)[0])
    NR = an.match_noise(U, members, pool, k=20)
    if verbose: print('  noise reference: pool %d blobs, top-20 matched (standardised mean difference to the group members: %s)' % (pool.sum(), dict(zip(NR['features'], np.round(NR['smd'], 2)))))
    out = dict(name=name, level=level, params={k: v for k, v in G.params.items() if k != 'homog_thr'}, units=G.units, unit_info=G.unit_info, tier=G.tier, clusters=G.clusters_of_units, unit_of=G.unit_of, cluster_of=G.cluster_of,
               attach=A, noise_top=NR['top'], noise_draws=NR['draws'], noise_pool_rows=np.where(pool)[0], noise_smd=NR['smd'], pool_size=int(pool.sum()), idx=U.idx, src=U.src, cls=U.cls, thr=thr.table)
    os.makedirs(CACHE, exist_ok=True)
    with open(os.path.join(CACHE, '%s_groups.pkl' % name), 'wb') as fh: pickle.dump(out, fh)
    # csv: one row per ROI
    geom = an.footprint_geometry(U); passes = U.v3()['passes']
    with open(os.path.join(CACHE, '%s_roi_groups.csv' % name), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['suite2p_idx', 'source', 'v3_class', 'unit_id', 'unit_size', 'cluster_id', 'tier', 'tier_why', 'attach_unit', 'attach_conf', 'attach_r', 'attach_q', 'attach_gap_um', 'orphan', 'noise_ref_top20', 'merged_parent_idx'])
        for i in range(U.N):
            src = ('accepted', 'rejected', 'constituent')[U.src[i]]; u = int(G.unit_of[i]) if U.src[i] == 0 else ''
            d = G.unit_info[G.unit_of[i]] if U.src[i] == 0 else None; b = A['best'].get(i) if A.get('best') else None
            w.writerow([int(U.idx[i]), src, U.cls[i] if U.src[i] == 0 else ('passes' if passes[i] else 'fails'), u, d['n'] if d else '', int(G.cluster_of[i]) if U.src[i] == 0 else '', G.tier.get(i, '') if U.src[i] == 0 else '',
                        ','.join(d['why']) if d else '', b['unit'] if b else '', A['conf'].get(i, '') if U.src[i] > 0 else '', round(b['r'], 3) if b else '', '%.2g' % b['q'] if b else '', round(b['gap'], 2) if b else '',
                        int(A['orphan'].get(i, False)) if U.src[i] > 0 else '', int(i in set(NR['top'])), int(U.inmerge[i]) if U.src[i] == 2 else ''])
    unit_dff = np.array([np.nanmean(U.dff[g], 0) for g in G.units]); np.savez(os.path.join(CACHE, '%s_units.npz' % name), t=U.t, unit_dff=unit_dff, unit_members_suite2p=np.array([[int(U.idx[i]) for i in g] for g in G.units], dtype=object),
                                                                               unit_tier=np.array([G.tier[g[0]] for g in G.units]), unit_avg_ok=np.array([bool(d.get('avg_ok', False)) for d in G.unit_info]))
    if validate:
        val = dict(null=av.null_runs(U, P), homog_table={n: v.tolist() for n, v in hd['table'].items()}, success=hd['success'])
        val['noise_compare'] = ar.compare_noise(U, P, out); val['sweep'] = av.param_sweep(U)
        val['hetero'] = av.sim_hetero(U, P, thr=thr, n_rep=30); val['held_out'] = [ar.held_out(U, level, fold=k, homog_thr=thr, verbose=verbose) for k in (0, 1)]
        with open(os.path.join(CACHE, '%s_validation_%s.pkl' % (name, level)), 'wb') as fh: pickle.dump(val, fh)
    if verbose: print('  done %s in %.0f s' % (name, time.time() - t0))
    return out


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]; level = 'pc1'
    if '--level' in sys.argv: level = sys.argv[sys.argv.index('--level') + 1]; args = [a for a in args if a != level]
    nrep = int(sys.argv[sys.argv.index('--nrep') + 1]) if '--nrep' in sys.argv else 100
    args = [a for a in args if a != str(nrep) or '--nrep' not in sys.argv]
    for nm in args: run(nm, level=level, validate='--validate' in sys.argv, refresh='--refresh' in sys.argv, n_rep=nrep)
