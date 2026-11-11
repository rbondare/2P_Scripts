"""Behaviour-matched contrasts of the lick-task events (engine: cross_session/matched_behaviour.py, same matching / nulls / behaviour-only benchmark).

usage: python dd_matched.py [session ...]       (one pickle per session in OUTDIR/cache/matched_<session>.pkl; env DD_OUTCOMES='P,g1,...' DD_NNULL, DD_NSUR)

Every event of a class (a lick bout of a class, a trial of an outcome, a boundary of a wait) is paired with 1-k events of the CONTROL class / control windows of the same session whose
movement in the response window is the same (face motion, whisking, JAW motion, lick rate in the window, ball where plausible; pupil STATE just before the window), the axon response
(30 s high-passed best dF/F, mean over the window) is compared inside the pairs.  A difference = activity the measured movement does not explain.  Windows are relative to the event
(bout onset, square onset, wait end).  Contrast definitions:
  bout classes  (window 0-3 s after the bout onset, matched on the movement of that window):   ant_to vs wait_early (closeness to the wait end), rew vs wait_early, rew vs ant_to,
                post_a vs wait_early, post_b vs wait_early, post_a vs rew, ant_sq vs ant_to, ant_sq vs wait_early (when there are enough bouts)
  pre-bout      (window -2..0 s before the bout onset: preparation / anticipation signal before the movement starts) the same pairs
  hit / miss    pre (-3..0 s), early (0..0.8 s after the onset, hits with a first lick later than 0.9 s only: no movement of licking in the window in either class), late (3..8 s)
  wait end      silent time-out: window 0..3 s after a time-out boundary and -3..0 s before it, vs random wait windows (all licking states matched)
  free pump     lick-free windows after the pump: level-2 trials without a lick (pump at onset + RewardDelay) vs level-3 trials without a lick (no pump), window 0..0.8 s after onset + RewardDelay
"""
import os, sys, pickle, time, zlib, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import dd_common as dd
import xs_common as xc
from xs_common import DT
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'cross_session'))
import matched_behaviour as mb
from matched_behaviour import _cum, _blocked

OUTCOMES = os.environ.get('DD_OUTCOMES', 'P').split(',')
N_NULL = int(os.environ.get('DD_NNULL', 200)); N_SUR = int(os.environ.get('DD_NSUR', 100))
MIN_EV = 8


class DRec(mb.Rec):
    BEH_ALL = ('pupil', 'face', 'whisk', 'loco', 'jaw', 'lick')

    def set_outcome(self, P):
        self.P = P; self.cumP = _cum(P); self._sur = None; self._cache = {}; self.cumD = None


def make_rec(name, beh=('face', 'whisk', 'jaw', 'lick', 'loco'), par=None):
    S, T = dd.load_session(name)
    p = dict(beh=beh, n_sur=N_SUR, n_null=N_NULL); p.update(par or {})
    R = DRec(name, p, S=S)
    for b in ('jaw', 'lick'):
        if b in R.blag: R.blag[b] = 0           # the lag fit of a movement that does not move the axon signal would only pick a noise peak
    sq = np.array([b for b in (np.floor((T.on - R.e0) / DT).astype(int)) if 0 <= b < R.nb]); R.allon = sq
    return R, T


def tb(R, t):
    return mb.to_bin(R, t)


def contrasts(R, T):
    """list of contrast dicts for this session"""
    out = []; p = R.par; nb = R.nb
    g0, g1 = [int(round(g / DT)) for g in p['guard']]
    st_edges = np.floor((np.sort(np.unique(T.d['sub_t0']))[1:] - R.e0) / DT).astype(int) if np.size(T.d['sub_t0']) > 1 else None
    cl = {c: T.bouts(c) for c in ('wait_early', 'ant_to', 'ant_sq', 'rew', 'post_a', 'post_b')}
    pairs = [('ant_to', 'wait_early'), ('ant_sq', 'wait_early'), ('ant_sq', 'ant_to'), ('rew', 'wait_early'), ('rew', 'ant_to'), ('post_a', 'wait_early'), ('post_b', 'wait_early'), ('post_a', 'rew')]
    for a, b in pairs:
        if len(cl[a]) < MIN_EV or len(cl[b]) < MIN_EV: continue
        for tag, spec in (('', (0.0, 3.0)), ('pre_', (-2.0, 0.0))):
            out.append(dict(key='%s%s_vs_%s' % (tag, a, b), kind='vs_class', spec=spec, focal=tb(R, cl[a]), ref=tb(R, cl[b])))
    # anticipation with the CONTEXT held constant: the wait licks of waits that end in a square are sparse (< 0.5 Hz), those that end in a time-out dense; early_sq / early_to = wait bouts more than
    # ANT_S s before the end of a wait of the same kind as the anticipatory bouts they are compared with
    bl = T.bloc; wait_far = (bl['phase'] == 'wait') & (bl['tte'] > dd.ANT_S)
    early_sq = T.bout_on[wait_far & (bl['endk'] == 'sq')]; early_to = T.bout_on[wait_far & (bl['endk'] == 'to')]
    for key, a, b in (('ant_sq_vs_early_sq', cl['ant_sq'], early_sq), ('ant_to_vs_early_to', cl['ant_to'], early_to)):
        if len(a) < MIN_EV or len(b) < MIN_EV: continue
        for tag, spec in (('', (0.0, 3.0)), ('pre_', (-2.0, 0.0))):
            out.append(dict(key=tag + key, kind='vs_class', spec=spec, focal=tb(R, a), ref=tb(R, b)))
    # hit / miss
    hit = T.hit; on = T.on
    if hit.sum() >= MIN_EV and (~hit).sum() >= MIN_EV:
        out.append(dict(key='hit_vs_miss_pre', kind='vs_class', spec=(-3.0, 0.0), focal=tb(R, on[hit]), ref=tb(R, on[~hit]), strata_edges=st_edges))
        slow = hit & (T.lat1 > 0.9)
        if slow.sum() >= MIN_EV: out.append(dict(key='hit_vs_miss_early', kind='vs_class', spec=(0.0, 0.8), focal=tb(R, on[slow]), ref=tb(R, on[~hit]), strata_edges=st_edges))
        out.append(dict(key='hit_vs_miss_late', kind='vs_class', spec=(3.0, 8.0), focal=tb(R, on[hit]), ref=tb(R, on[~hit]), strata_edges=st_edges))
    # the trial's own post-square period vs random times of the session (licking states matched): is the axon level different 3-6 / 6-10 s after the square at all?
    for key, off, spec in (('trial_3_6', 3.0, (0.0, 3.0)), ('trial_6_10', 6.0, (0.0, 4.0))):
        out.append(dict(key=key, kind='vs_pool', spec=spec, focal=tb(R, on + off), blockers=[(R.allon, g0, g1)]))
    # silent time-outs
    if T.has_timeouts:
        te = T.seg['t1'][T.seg['to']]
        lickb = tb(R, T.bout_on)
        for key, spec in (('timeout_post', (0.0, 3.0)), ('timeout_pre', (-3.0, 0.0))):
            out.append(dict(key=key, kind='vs_pool', spec=spec, focal=tb(R, te), blockers=[(R.allon, g0, g1)]))
    # free pump (level 2, no lick) vs no pump (level 3, no lick), lick-free window after onset + RewardDelay
    fl = np.full(T.n, np.nan)
    for t in range(T.n):
        rel = T.lick - T.on[t]; m = (T.lick_tr == t + 1) & (rel >= 0)
        if m.any(): fl[t] = rel[m].min()
    nolick = ~(fl < T.rdel + 0.9)
    A = (T.lvl == 2) & nolick & (T.rew_kind == 2); B = (T.lvl == 3) & nolick & (T.rew_kind == 0)
    if A.sum() >= MIN_EV and B.sum() >= MIN_EV:
        out.append(dict(key='free_pump_vs_nopump', kind='vs_class', spec=(0.0, 0.8), focal=tb(R, on[A] + T.rdel), ref=tb(R, on[B] + T.rdel), strata_edges=st_edges))
    return out


def run_session(name, verbose=True):
    t0 = time.time(); R, T = make_rec(name); out = dict(name=name, n_sig=R.n_sig, names=R.names, par=R.par, contrasts={}, counts=T.bout_class_counts)
    from dd_fir import groups_from_log
    groups = groups_from_log(name, R.S); hp = lambda x: x - xc.nsmooth(x, R.par['hp_s'] / DT)
    outs = {'P': R.P}
    for gi, g in enumerate(groups): outs['g%d' % (gi + 1)] = hp(np.nanmean(R.S.dff[g], 0))
    if R.S.noise.sum() >= 3: outs['control'] = hp(R.S.population('noise'))
    cons = contrasts(R, T)
    if verbose: print('%s: behaviours %s, %d contrasts, outcomes %s' % (name, R.names, len(cons), OUTCOMES), flush=True)
    for on in OUTCOMES:
        if on not in outs: continue
        R.set_outcome(outs[on]); out['contrasts'][on] = {}; R.par['n_sur'] = N_SUR if on == 'P' else min(N_SUR, 40)
        for c in cons:
            rng = np.random.default_rng(zlib.crc32((name + on + c['key']).encode()))
            try: rr = mb.run_contrast(R, c, N_NULL, rng, with_extras=(on == 'P'))
            except Exception as e:
                rr = dict(key=c['key'], note='failed: %s' % e, n_events=len(c['focal']))
            out['contrasts'][on][c['key']] = rr
            if verbose:
                o = rr.get('obs')
                if o and np.isfinite(o.get('matched', np.nan)):
                    s = rr['stats']['matched']; ss = rr.get('stats_sim', {}).get('matched', {})
                    print('  %-4s %-26s n %3d valid %3d matched %3d | naive %+.3f matched %+.3f (shift z %+.1f p %.3f | bench z %+.1f p %.3f) | max|SMD| %.2f->%.2f' % (
                        on, c['key'], rr['n_events'], o['n_focal'], o['n_matched'], o['naive'], o['matched'], s['z'], s['p'], ss.get('z', np.nan), ss.get('p', np.nan),
                        np.nanmax(np.abs(o['smd_before'])), np.nanmax(np.abs(o['smd_after']))), flush=True)
                else: print('  %-4s %-26s n %3d: %s' % (on, c['key'], rr.get('n_events', 0), rr.get('note', '')), flush=True)
    out['sec'] = time.time() - t0
    return out


if __name__ == '__main__':
    os.makedirs(os.path.join(dd.OUTDIR, 'cache'), exist_ok=True)
    tag = os.environ.get('DD_TAG', '')
    for nm in (sys.argv[1:] or dd.TASK_SESSIONS):
        f = os.path.join(dd.OUTDIR, 'cache', 'matched%s_%s.pkl' % (tag, nm))
        if os.path.exists(f) and not os.environ.get('REFRESH'): print('cached', nm); continue
        try:
            r = run_session(nm); pickle.dump(r, open(f, 'wb')); print('  %.0f s' % r['sec'], flush=True)
        except Exception:
            import traceback; traceback.print_exc()
    print('MATCHED_DONE')
