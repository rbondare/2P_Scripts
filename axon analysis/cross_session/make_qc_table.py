"""Session QC table of the cross-session analyses: plots/cross_session/session_qc.csv  (one row per session, from state_part*.pkl and stim_response_results.pkl)"""
import os, sys, glob, pickle, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import xs_common as xc
FIG = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\cross_session'
res = {}
for f in sorted(glob.glob(os.path.join(FIG, 'state_part*.pkl'))): res.update(pickle.load(open(f, 'rb')))
stim = pickle.load(open(os.path.join(FIG, 'stim_response_results.pkl'), 'rb'))
rows = []
for nm in sorted(res, key=lambda n: (n.replace('Animal', '')[:4][::-1], n)):
    r = res[nm]; row = dict(session=nm, rois_accepted=r.get('n_roi'), signal_like=r.get('n_sig'), noise_like=r.get('n_noise'), neuropil_limited=r.get('n_neuropil'), rejected_rois_used=r.get('n_rej'),
                            minutes=round(r.get('dur_min', np.nan), 1), ball=('used' if r.get('loco_ok') else 'not used') + ' (' + str(r.get('loco_note')) + ')')
    if r.get('usable'):
        bg = r['beh_grid']; bp = {k: xc.bandpass(v.astype(float)) for k, v in bg.items() if k in ('face', 'whisk', 'pupil')}
        def rr(a, b):
            if a not in bp or b not in bp: return ''
            ok = np.isfinite(bp[a]) & np.isfinite(bp[b]); return round(float(np.corrcoef(bp[a][ok], bp[b][ok])[0, 1]), 2) if ok.sum() > 100 else ''
        row.update(beh_face_whisk_r=rr('face', 'whisk'), beh_face_pupil_r=rr('face', 'pupil'), beh_whisk_pupil_r=rr('whisk', 'pupil'), beh_face_sd_log=round(float(np.nanstd(bg['face'])), 2) if 'face' in bg else '')
        c = r['coh']['sig']; rj = [x for x in r['coh']['rej'] if x]
        row.update(pc1_signal=round(c['pc1'], 2), pc1_rejected_sets=round(np.mean([x['pc1'] for x in rj]), 2) if rj else '', same_sign_signal=round(c['same_sign'], 2), same_sign_rejected=round(np.mean([x['same_sign'] for x in rj]), 2) if rj else '')
        for b in ('pupil', 'face', 'whisk', 'loco', 'shift'):
            d = r['xc'].get(b)
            if d: row['r_%s' % b] = round(d['peak_r'], 2); row['lag_%s_s' % b] = round(d['peak_lag_s'], 1); row['p_%s' % b] = round(d['p'], 3)
            sg = r['roi'].get(b, {}).get('sig'); rg = r['roi'].get(b, {}).get('rej')
            if sg: row['roi_net_sign_signal_%s' % b] = round((sg['pos'] - sg['neg']) / sg['n'], 2)
            if rg: row['roi_net_sign_rejected_%s' % b] = round((rg['pos'] - rg['neg']) / rg['n'], 2)
        if 'ridge' in r:
            row['ridge_R2_signal'] = round(r['ridge']['sig']['full'], 2); row['ridge_R2_noise_like'] = round(r['ridge'].get('noise', np.nan), 2) if 'noise' in r['ridge'] else ''
            row['ridge_R2_rejected_sets_mean'] = round(float(np.nanmean(r['ridge']['rej'])), 2) if len(r['ridge'].get('rej', [])) else ''
    else: row['excluded'] = r.get('why')
    s = stim.get(nm, {}); row['stimuli'] = ', '.join('%s x%d' % (k, v['n']) for k, v in s.get('keys', {}).items())
    rows.append(row)
cols = []
for r in rows:
    for k in r:
        if k not in cols: cols.append(k)
with open(os.path.join(FIG, 'session_qc.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); [w.writerow(r) for r in rows]
print('saved session_qc.csv', len(rows), 'sessions')
