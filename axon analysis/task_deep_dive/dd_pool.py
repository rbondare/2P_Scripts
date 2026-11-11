"""Pool the per-session matched contrasts (population outcome P) over the lick sessions -> pooled_contrasts.csv + printed table.

Unit = session; ALL sessions are one animal (RB24; the only RB23 lick session is flagged), so the pooled numbers describe repeatability across days in ONE mouse, not across animals.
Included per contrast: unflagged sessions with >= 8 matched events and >= 40% of the events matched (the flagged 1317 / RB23 1704 are listed separately).
Pooled = mean of the session estimates; CI = bootstrap over sessions (percentile, n is 3-5: wide by construction); n_pos = sessions with a positive estimate;
Stouffer z = sum of the sessions' shift-null z / sqrt(n) (each session's z comes from its own circular-shift null, so autocorrelation is inside it).
"""
import os, sys, csv, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dd_figures as f
import dd_common as dd

KEYS = ['trial_3_6', 'trial_6_10', 'ant_to_vs_early_to', 'ant_sq_vs_early_sq', 'ant_to_vs_wait_early', 'ant_sq_vs_wait_early', 'ant_sq_vs_ant_to', 'rew_vs_wait_early', 'rew_vs_ant_to', 'post_a_vs_wait_early', 'post_b_vs_wait_early', 'post_a_vs_rew',
        'pre_ant_to_vs_early_to', 'pre_ant_sq_vs_early_sq', 'pre_ant_to_vs_wait_early', 'pre_ant_sq_vs_wait_early', 'pre_rew_vs_wait_early', 'pre_post_b_vs_wait_early',
        'timeout_post', 'timeout_pre', 'hit_vs_miss_pre', 'hit_vs_miss_early', 'hit_vs_miss_late', 'free_pump_vs_nopump']


def main():
    rng = np.random.default_rng(1); rows = []
    print('%-26s %-5s %8s %20s %6s %8s %s' % ('contrast', 'n', 'pooled', '95% CI (bootstrap)', 'n_pos', 'Stouffer', 'sessions (estimate, matched/events)'))
    for k in KEYS:
        est, zs_, names, det = [], [], [], []
        for n in f.SESS_ORDER:
            m = f.load('matched_tr' if k.startswith('trial_') else 'matched', n); c = f.get_c(m, 'P', k) if m else None
            if c is None: continue
            o = c['obs']
            if n in f.FLAGGED or o['n_matched'] < 8 or o['n_matched'] < 0.4 * o['n_focal']: continue
            est.append(o['matched']); zs_.append(c['stats']['matched']['z']); names.append(n); det.append('%s %+.3f %d/%d' % (f.sh(n)[5:], o['matched'], o['n_matched'], o['n_focal']))
        if len(est) < 2: rows.append([k, len(est), '', '', '', '', '', '; '.join(det)]); print('%-26s %-5d (too few sessions)' % (k, len(est))); continue
        est = np.array(est); zs_ = np.array(zs_); ok = np.isfinite(zs_)
        bs = [rng.choice(est, len(est)).mean() for _ in range(5000)]; ci = np.percentile(bs, [2.5, 97.5]); st = zs_[ok].sum() / np.sqrt(ok.sum())
        rows.append([k, len(est), round(est.mean(), 4), round(ci[0], 4), round(ci[1], 4), int((est > 0).sum()), round(st, 2), '; '.join(det)])
        print('%-26s %-5d %+8.3f [%+.3f, %+.3f] %4d/%d %+8.1f  %s' % (k, len(est), est.mean(), ci[0], ci[1], (est > 0).sum(), len(est), st, ' | '.join(det)))
    with open(os.path.join(dd.OUTDIR, 'pooled_contrasts.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['contrast', 'n_sessions', 'pooled_matched_dFF', 'ci_lo', 'ci_hi', 'n_positive', 'stouffer_z', 'sessions']); w.writerows(rows)


if __name__ == '__main__':
    main()
