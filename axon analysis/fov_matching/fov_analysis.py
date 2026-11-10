"""Step 5b: what differs between sessions / days (tables for the figures; all numbers also printed).

Inputs: cache/master.pkl (master axons x sessions, class + contrast), roi_match_pairs.csv, fov_inventory.csv, plots/cross_session/session_qc.csv, cache/poses_fine.pkl
 1. Leave-one-out recall: confirmed axons for session s = master axons accepted in >= 2 OTHER sessions and fully in view of s. recall = fraction accepted in s (unbiased
    sensitivity of s w.r.t. the consensus of the others). Missing ones are split by class (rejected ROI present / structure, no ROI / nothing visible).
 2. Contrast (robust z of the axon-scale DoG of the mean image, Vcorr, max_proj-mean) at the confirmed footprints: where accepted vs where missing vs background (95th pct).
 3. Persistence in the 11-04..11-06 block: in how many of the sessions each axon (fully in view of all of them) was accepted; compared with the distribution expected if each
    session had accepted a RANDOM subset of the same N axons (per-session rates p_s = n_s / N, N solved so that the expected number of axons seen at least once equals the
    observed number); heterogeneity = more axons seen always or only once than that model gives.
 4. Directional carry-over matrix P(accepted in s | accepted in t, in view of s) from the pair matching.
 5. Per-session descriptors vs counts (day, hFastZ, behavioural range, ...), exploratory.
Outputs: fov_session_summary.csv, fov_persistence.csv, fov_carryover.csv, cache/analysis.pkl
"""
import os, csv, pickle, collections, numpy as np
from fov_common import *

NOV = ['AnimalRB24_261104_1540', 'AnimalRB24_261104_1651', 'AnimalRB24_261105_1423', 'AnimalRB24_261105_1520', 'AnimalRB24_261105_1549',
       'AnimalRB24_261106_1317', 'AnimalRB24_261106_1529', 'AnimalRB24_261106_1555', 'AnimalRB24_261106_1604']      # 261104_1640 excluded from all ROI analyses (junk band ROIs)
ACCEPT_ALL = {'AnimalRB24_261106_1604'}          # iscell all-1 (user-curated 93 + 67 new): no rejected class available


def main():
    d = pickle.load(open(os.path.join(CACHE, 'master.pkl'), 'rb'))
    cl, rows, names = d['clusters'], d['rows'], d['names']; thr = d['thr']
    inv = read_inventory()
    pres = collections.defaultdict(dict)
    for r in rows: pres[r['master']][r['session']] = r
    out = []
    # ---- 1 + 2 per-session summary
    for s in names:
        conf = [i for i, c in enumerate(cl) if len([x for x in c if x != s]) >= 2 and pres[i][s]['cls'] != 'out of FOV']
        cnt = collections.Counter(pres[i][s]['cls'] for i in conf); t = len(conf)
        def zs(cls_set, key):
            v = [pres[i][s].get(key) for i in conf if pres[i][s]['cls'] in cls_set]
            v = np.array([x for x in v if x is not None and x != ''], float); return (float(np.nanmean(v)) if len(v) else np.nan, len(v))
        n_all_in = sum(pres[i][s]['cls'] != 'out of FOV' for i in range(len(cl)))
        nfr = np.load(os.path.join(s2p_dir(s), 'F.npy'), mmap_mode='r').shape[1]
        vol = 5.0017 if int(inv[s]['n_slices']) > 1 else float(inv[s]['frame_rate'])       # 5-plane sessions: suite2p frames = volumes
        n_single = sum(1 for i, c in enumerate(cl) if list(c) == [s])
        row = dict(session=s, short=short(s), day=inv[s]['tif_first_mtime'][5:10], n_accepted=int(inv[s]['n_roi_accepted']), n_roi_total=inv[s]['n_roi_total'], accept_all='yes' if s in ACCEPT_ALL else 'no',
                   n_singletons=n_single, n_confirmed_in_view=t, recall=round(cnt['accepted'] / t, 3) if t else np.nan,
                   frac_rejected_present=round(cnt['rejected ROI present'] / t, 3) if t else np.nan, frac_structure_no_roi=round(cnt['structure, no ROI'] / t, 3) if t else np.nan,
                   frac_nothing=round(cnt['nothing visible'] / t, 3) if t else np.nan,
                   z_dog_accepted=round(zs({'accepted'}, 'z_dog')[0], 2), z_dog_missing=round(zs({'rejected ROI present', 'structure, no ROI', 'nothing visible'}, 'z_dog')[0], 2),
                   z_dog_nothing=round(zs({'nothing visible'}, 'z_dog')[0], 2), z_vcorr_accepted=round(zs({'accepted'}, 'z_vcorr')[0], 2), z_vcorr_missing=round(zs({'rejected ROI present', 'structure, no ROI', 'nothing visible'}, 'z_vcorr')[0], 2),
                   z_vcorr_nothing=round(zs({'nothing visible'}, 'z_vcorr')[0], 2), bg95_dog=round(thr[s]['dog'], 2), bg95_vcorr=round(thr[s]['vcorr'], 2),
                   zoom=inv[s]['zoom'], um_per_px=inv[s]['um_per_px_x'], fastz=inv[s]['fastz_pos'], axes_z=inv[s]['axes_z'], n_frames=int(nfr), minutes=round(nfr / vol / 60, 1), face_sd=inv[s]['beh_face_sd_log'])
        out.append(row)
    with open(os.path.join(OUT, 'fov_session_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    print('per-session leave-one-out summary (confirmed = accepted in >= 2 other sessions):')
    print('  session     acc  conf_in_view recall  rejPres struct  nothing | z_dog acc/missing/nothing | z_vcorr acc/missing/nothing | bg95 dog/vcorr')
    for r in out:
        print('  %s %4d %5d      %.2f   %.2f    %.2f    %.2f   | %5s %5s %5s | %5s %5s %5s | %s %s' % (r['short'], r['n_accepted'], r['n_confirmed_in_view'], r['recall'], r['frac_rejected_present'], r['frac_structure_no_roi'], r['frac_nothing'],
                                                                                                           r['z_dog_accepted'], r['z_dog_missing'], r['z_dog_nothing'], r['z_vcorr_accepted'], r['z_vcorr_missing'], r['z_vcorr_nothing'], r['bg95_dog'], r['bg95_vcorr']))
    # ---- 3 persistence in November block
    nov = [s for s in NOV if s in names]
    common = [i for i in range(len(cl)) if all(pres[i][s]['cls'] != 'out of FOV' for s in nov)]
    X = np.array([[pres[i][s]['cls'] == 'accepted' for s in nov] for i in common], bool)
    seen = X.sum(1); keep = seen >= 1; X = X[keep]; seen = seen[keep]; S_obs = len(X)
    n_s = X.sum(0)
    print('\nPERSISTENCE in the 11-04..06 block: %d sessions, %d axons fully in view of all and accepted at least once' % (len(nov), S_obs))
    f_obs = np.bincount(seen, minlength=len(nov) + 1)[1:]
    # model: random subset per session of N axons
    def expected_seen(N):
        p = np.clip(n_s / N, 0, 1); q0 = np.prod(1 - p); return N * (1 - q0)
    lo, hi = S_obs, S_obs * 50
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if expected_seen(mid) < S_obs: lo = mid
        else: hi = mid
    N_hat = 0.5 * (lo + hi); p = np.clip(n_s / N_hat, 0, 1)
    rng = np.random.default_rng(0); sims = []
    for _ in range(400):
        Y = rng.random((int(round(N_hat)), len(nov))) < p[None, :]; k = Y.sum(1); k = k[k >= 1]; sims.append(np.bincount(k, minlength=len(nov) + 1)[1:] * (S_obs / max(len(k), 1)))
    sims = np.array(sims); f_mod = sims.mean(0); f_lo, f_hi = np.percentile(sims, [2.5, 97.5], axis=0)
    print('  sessions seen in:       ' + ' '.join('%4d' % k for k in range(1, len(nov) + 1)))
    print('  observed (axons):       ' + ' '.join('%4d' % v for v in f_obs))
    print('  random-subset model:    ' + ' '.join('%4.0f' % v for v in f_mod) + '   (N_total solved = %.0f)' % N_hat)
    print('  model 95%% band low:     ' + ' '.join('%4.0f' % v for v in f_lo))
    print('  model 95%% band high:    ' + ' '.join('%4.0f' % v for v in f_hi))
    # sensitivity: without the two accept-all sessions
    nov2 = [s for s in nov if s not in ACCEPT_ALL]; ix = [nov.index(s) for s in nov2]
    X2 = X[:, ix]; seen2 = X2.sum(1); X2 = X2[seen2 >= 1]; seen2 = seen2[seen2 >= 1]
    print('  without the accept-all sessions (%d sessions): observed %s' % (len(nov2), ' '.join('%d' % v for v in np.bincount(seen2, minlength=len(nov2) + 1)[1:])))
    with open(os.path.join(OUT, 'fov_persistence.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['n_sessions_seen', 'observed', 'model_mean', 'model_lo95', 'model_hi95'])
        for k in range(len(nov)): w.writerow([k + 1, int(f_obs[k]), round(float(f_mod[k]), 1), round(float(f_lo[k]), 1), round(float(f_hi[k]), 1)])
    # ---- 4 carry-over matrix
    pr = list(csv.DictReader(open(os.path.join(OUT, 'roi_match_pairs.csv'))))
    idx = {s: i for i, s in enumerate(names)}; C = np.full((len(names), len(names)), np.nan); Nn = np.zeros_like(C)
    for r in pr:
        a, b = r['a'], r['b']; m = int(r['n_matched'])
        na, nb = int(r['n_a_inview']), int(r['n_b_inview'])
        if na >= 5: C[idx[a], idx[b]] = m / na; Nn[idx[a], idx[b]] = na          # P(accepted in b | accepted in a, a's ROI in view of b)
        if nb >= 5: C[idx[b], idx[a]] = m / nb; Nn[idx[b], idx[a]] = nb
    with open(os.path.join(OUT, 'fov_carryover.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['from \\ to'] + [short(s) for s in names])
        for s in names: w.writerow([short(s)] + ['' if not np.isfinite(C[idx[s], idx[t]]) else round(C[idx[s], idx[t]], 3) for t in names])
    sd = []; ad = []; far = []
    day = {s: np.datetime64('2026-' + inv[s]['tif_first_mtime'][5:10]) for s in names}
    for s in nov:
        for t in nov:
            if s == t or not np.isfinite(C[idx[s], idx[t]]): continue
            gap = abs((day[s] - day[t]).astype(int))
            (sd if gap == 0 else ad if gap == 1 else far).append(C[idx[s], idx[t]])
    print('\nCARRY-OVER within the November block, P(accepted in t | accepted in s): same day %.2f (n=%d), 1 day apart %.2f (n=%d), 2 days apart %.2f (n=%d)' % (
        np.mean(sd), len(sd), np.mean(ad), len(ad), np.mean(far) if far else np.nan, len(far)))
    # same comparison restricted to pairs without accept-all sessions
    sd2 = []; ad2 = []
    for s in nov2:
        for t in nov2:
            if s == t or not np.isfinite(C[idx[s], idx[t]]): continue
            gap = abs((day[s] - day[t]).astype(int)); (sd2 if gap == 0 else ad2 if gap == 1 else []).append(C[idx[s], idx[t]])
    print('   without the accept-all sessions: same day %.2f (n=%d), 1 day apart %.2f (n=%d)' % (np.mean(sd2), len(sd2), np.mean(ad2), len(ad2)))
    pickle.dump(dict(out=out, nov=nov, X=X, seen=seen, f_obs=f_obs, f_mod=f_mod, f_lo=f_lo, f_hi=f_hi, N_hat=N_hat, C=C, Nn=Nn, names=names, common=[common[i] for i in np.where(keep)[0]]), open(os.path.join(CACHE, 'analysis.pkl'), 'wb'))


if __name__ == '__main__':
    main()
