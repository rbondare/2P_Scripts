"""FIR lick-class model with the licks of the GREY BREAKS as their own classes ('brk' = between subsessions, 'tail' = after the last trial).

usage: python dd_break_fir.py [session ...]       (one pickle per session: OUTDIR/cache/firbrk_<session>.pkl; REFRESH=1 recomputes)

Identical to dd_fir.py (same design, same covariates, same ridge / bootstrap / shift-null) except that LICK_CLASSES gets two more entries, so the kernel of a break lick is estimated next to
the kernels of wait / anticipatory / reward / post licks, with the movement covariates (jaw, face, whisking, pupil, shift, ball) taken out linearly.  The covariates are 30 s high-passed
(dd_fir.build_design), i.e. the sustained pupil dilation of the break itself is NOT a regressor: only its fast part is.
"""
import os, sys, pickle, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dd_common as dd
import dd_fir

dd_fir.LICK_CLASSES = ['wait_early', 'ant_to', 'ant_sq', 'win', 'post_a', 'post_b', 'brk', 'tail']

# BOX=1: the break itself is modelled too, so that the kernel of a break lick is the LICK-LOCKED part only.  Without it a break-wide shift of the axon level (or a transient at the break boundaries)
# would be absorbed by the break-lick kernel, because break licks are spread over the whole break.  Added: a boxcar (1 inside the break window) per kind (columns 'cov_box_brk' / 'cov_box_tail',
# not part of the kernel summaries; its weight is the sustained level of the break) and FIR kernels (-1 .. +4 s) at the break start ('brk_on', 'tail_on') and end ('brk_off', 'tail_off').
BOX = bool(os.environ.get('BOX'))
if BOX:
    _orig_build = dd_fir.build_design

    def _build_box(S, T, covs=('jaw', 'face', 'whisk', 'pupil', 'shift', 'loco')):
        X, cmap, used, cls = _orig_build(S, T, covs); cols = [X]; k = X.shape[1]; cmap = dict(cmap)
        def add(name, M):
            nonlocal k; cols.append(M); cmap[name] = slice(k, k + M.shape[1]); k += M.shape[1]
        for kind in ('brk', 'tail'):
            w = [x for x in T.brk if x['kind'] == kind]
            box = np.zeros(S.nb)
            for x in w: box[(S.t >= x['t0']) & (S.t < x['t1'])] = 1.0
            add('cov_box_' + kind, (box - box.mean())[:, None])
            if len(w):
                add(kind + '_on', dd_fir.lagged(dd_fir.event_series(np.array([x['t0'] for x in w]), S.edges), dd_fir.LAG_Y))
                add(kind + '_off', dd_fir.lagged(dd_fir.event_series(np.array([x['t1'] for x in w]), S.edges), dd_fir.LAG_Y))
        return np.concatenate(cols, 1), cmap, used, cls
    dd_fir.build_design = _build_box
PREFIX = 'firbrk2_' if BOX else 'firbrk_'

if __name__ == '__main__':
    os.makedirs(os.path.join(dd.OUTDIR, 'cache'), exist_ok=True)
    for nm in (sys.argv[1:] or dd.TASK_SESSIONS):
        f = os.path.join(dd.OUTDIR, 'cache', PREFIX + nm + '.pkl')
        if os.path.exists(f) and not os.environ.get('REFRESH'): print('cached', nm); continue
        T = dd.Task(nm)
        if not T.brk: print('no break in', nm); continue
        try:
            r = dd_fir.run_session(nm); pickle.dump(r, open(f, 'wb')); print('  %.0f s' % r['sec'], flush=True)
        except Exception:
            import traceback; traceback.print_exc()
    print('FIRBRK_DONE')
