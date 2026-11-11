"""Shared code of the lick-task deep dive (task_deep_dive/*.py): trial / wait-segment / lick-class table of one session from the export of export_task_events.py.

TASK LOGIC (read from MouseTrainingScripts/timed_Training.m, level 2-4, and confirmed on RB24_261109_1506 by the rig's own counts, check_segments()):
  every trial starts with a wait X (randi 8..13 s).  At the END of the wait the rig counts the licks of that wait (PreLicks); if PreLicks > maxErrorLicks * X (0.5 licks/s)
  it is a TIME-OUT: a new wait X' is drawn and the count restarts, the screen stays grey, TimeOutSound is 0 (silent).  Otherwise the square (3 s) follows: a lick 0.1-3 s after
  its onset = hit; after the square a 8.5 s grey.  waitTimes(trial, k) = the k-th wait of the trial, errorLicks(trial, k) = the rig's lick count of that wait.
  A time-out is therefore decided at the end of a wait from the licks of the WHOLE wait (not punished lick by lick), and nothing marks it for the animal.

LICK CLASSES (bout level, a bout = >= 3 licks with gaps < 0.5 s, first lick = bout onset; per-lick level for the FIR model):
  wait_early  bout onset in a wait segment, more than ANT_S s before the end of that wait
  ant_to      bout onset in the last ANT_S s of a wait that ENDS IN A TIME-OUT  (anticipation, licking that is too much -> time-out)
  ant_sq      bout onset in the last ANT_S s of the wait that ends in the SQUARE  (anticipation, licking below the limit)
  rew         the bout that holds the first lick 0.1-3 s after the square onset (reward trigger / collection; task script cons_on)
  win_other   other bouts that start 0.1-3 s after the onset
  post_a      bout onset 3-6 s after the onset (after the square, reward already delivered / collected)
  post_b      bout onset 6 s after the onset to the end of the trial
  brk         bout onset in a GREY BREAK between two subsessions (added 2026-11-11): the stimulus is grey, no trial runs and nothing is ever rewarded; imaging continues up to the
              acquisition gap.  The lick log of the last trial of a subsession keeps recording through the break (see Task._breaks for the window and the log-capacity caveat)
  tail        the same after the LAST trial of the session (96-135 s of imaging after it in 1317 / 1423 / 1744 / 1752; the lick log of 1506's last trial is missing)
"""
import os, sys, numpy as np, scipy.io as sio
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'cross_session'))
import xs_common as xc
XS_V2 = os.environ.get('TASK_OUT', r'C:\Users\rbondarenko\2P_cache\xs_v2')
xc.XS_DIR = XS_V2
DT = xc.DT
ANT_S = 3.0
BRK_GUARD_S = 0.5            # a break starts this long after the end of the last trial's post-square grey
BRK_MIN_S = 20.0             # shorter stretches (RB24_261105_1520: 3 s) are not breaks
LOG_CAP = 242                # rows per trial of lickData (V1 script): a trial that reaches it lost its later licks
OUTDIR = r'c:\Users\rbondarenko\projects\2P_Scripts\axon analysis\plots\task_deep_dive'
os.makedirs(OUTDIR, exist_ok=True)

# the lick-task sessions (all RB24 except the flagged RB23 one); level / RewardDelay as in the stimulus logs
TASK_SESSIONS = ['AnimalRB24_261109_1506', 'AnimalRB24_261106_1317', 'AnimalRB24_261105_1520', 'AnimalRB24_261105_1423', 'AnimalRB24_261030_1744', 'AnimalRB24_261029_1752', 'AnimalRB23_261029_1704']
FLAGGED = {'AnimalRB23_261029_1704': 'behaviour channels inconsistent, lick clock from face motion (+-0.2 s), ball truncated: shown, not pooled',
           'AnimalRB24_261106_1317': 'no behaviour coupling (11-06 sessions) and 135 ROIs incl. noise-like by design: shown, not pooled'}


def _item(v):
    return v.item() if isinstance(v, np.ndarray) and v.ndim == 0 else v


class Task:
    def __init__(self, name):
        self.name = name
        d = sio.loadmat(os.path.join(XS_V2, name + '_task.mat'), squeeze_me=True)
        T = {k: _item(v) for k, v in d.items() if not k.startswith('__')}
        self.d = T
        f = lambda k: np.atleast_1d(np.asarray(T[k], float))
        self.n = int(T['nTr']); n = self.n
        self.tr0, self.on, self.st_end, self.tr1 = f('tr_t0'), f('on_t'), f('st_end'), f('tr_t1')
        self.sub = f('sub_of').astype(int); self.lvl = f('lvl_t').astype(int); self.hit = f('hit') > 0; self.n_to = f('n_to').astype(int)
        self.lat1 = f('lat1'); self.rew_t = f('rew_t'); self.rew_kind = f('rew_kind').astype(int); self.trig = f('trig_t'); self.cons_on = f('cons_on')
        self.lw = f('lw'); self.rdel = float(T['rdel']); self.Tend = float(T['Tend'])
        self.lick = np.atleast_1d(np.asarray(T['lick_t'], float)); self.lick_tr = np.atleast_1d(np.asarray(T['lick_tr'], float)).astype(int)
        self.bout_on = f('bout_on'); self.bout_off = f('bout_off') if 'bout_off' in T else None; self.bout_n = f('bout_n')
        def mat(v):                                                      # squeeze_me turns a single column into a vector: (trials,) -> (trials, 1)
            v = np.asarray(v, float)
            if v.ndim == 0: return v.reshape(1, 1)
            return v.reshape(-1, 1) if v.ndim == 1 else v
        W = mat(T['waitTimes']); self.W = W
        self.E = mat(T['errorLicks']) if 'errorLicks' in T else None
        self.maxerr = float(T.get('sp_maxErrorLicks', np.nan)); self.minwait = float(T.get('sp_minWait', np.nan)); self.maxwait = float(T.get('sp_maxWait', np.nan))
        self.has_timeouts = bool(np.any(self.n_to > 0))
        self.t_first, self.t_last = float(T['t_first']), float(T['t_last'])
        self.sub_t0 = f('sub_t0'); self.sub_t1 = f('sub_t1')
        self._breaks(); self._segments(); self._classify_bouts()

    # ------------------------------------------------------------------ grey breaks between the subsessions (+ the tail after the last trial)
    def _breaks(self):
        """windows of the grey breaks: from the end of the last trial of a subsession (+ BRK_GUARD_S) to the start of the next imaging gap (the acquisition is restarted there; for the
        tail: the end of the imaging).  The lick log of the last trial of a subsession keeps recording through the break; if it reached its capacity (LOG_CAP rows) the window ends at its
        last lick (later licks were lost: unobservable, NOT 'no licking').  kind = 'brk' (between subsessions) | 'tail' (after the last trial)"""
        gaps = np.atleast_1d(np.asarray(self.d.get('imaging_gap_t', []), float)); n_sub = len(self.sub_t0)
        self.brk = []
        for b in range(n_sub):
            tr_in = np.where(self.sub == b + 1)[0]
            if len(tr_in) == 0: continue
            t_end = float(self.sub_t1[b]); last = tr_in[-1] + 1
            if b < n_sub - 1: g = gaps[gaps > t_end - 1.0]; t_stop = float(g[0]) if len(g) else float(self.sub_t0[b + 1]); kind = 'brk'
            else: t_stop = self.t_last; kind = 'tail'
            nl = int((self.lick_tr == last).sum()); capped = nl >= LOG_CAP
            lk = self.lick[self.lick_tr == last]
            if capped: t_stop = min(t_stop, float(lk.max()))
            if kind == 'tail' and self.n < int(self.d['n_rec']): continue      # a trial after the last logged one ran without a lick log (1506): the tail is not observed
            t0 = t_end + BRK_GUARD_S
            logged = bool(nl > 0) and (last <= self.n)
            if t_stop - t0 < BRK_MIN_S: continue
            self.brk.append(dict(b=b + 1, kind=kind, t0=t0, t1=t_stop, last_trial=int(last), n_logged_last=nl, capped=capped, log_ok=logged))
        self.brk_t0 = np.array([x['t0'] for x in self.brk]); self.brk_t1 = np.array([x['t1'] for x in self.brk]); self.brk_kind = np.array([x['kind'] for x in self.brk], dtype=object)

    def in_break(self, t_ev, kind=None):
        """index of the break window holding each event time (-1 = none); kind restricts to 'brk' or 'tail'"""
        t_ev = np.asarray(t_ev, float); out = np.full(len(t_ev), -1)
        for j, x in enumerate(self.brk):
            if kind is not None and x['kind'] != kind: continue
            out[(t_ev >= x['t0']) & (t_ev < x['t1'])] = j
        return out

    # ------------------------------------------------------------------ wait segments
    def _segments(self):
        rows = []
        for t in range(self.n):
            w = self.W[t][self.W[t] > 0]
            if len(w) == 0: continue
            b = self.tr0[t] + np.concatenate([[0], np.cumsum(w)])
            b[-1] = self.on[t]                                       # last wait ends at the square onset (frame time), earlier ends from the logged waits (<= 0.15 s off)
            for k in range(len(w)):
                rows.append((t, k, b[k], b[k + 1], w[k], k == len(w) - 1, self.E[t, k] if (self.E is not None and t < self.E.shape[0] and k < self.E.shape[1]) else np.nan))
        r = np.array(rows, float)
        self.seg = dict(trial=r[:, 0].astype(int), k=r[:, 1].astype(int), t0=r[:, 2], t1=r[:, 3], X=r[:, 4], last=r[:, 5] > 0, rig=r[:, 6])
        self.seg['to'] = ~self.seg['last']                           # a segment that is not the last one of its trial ended in a time-out
        s = self.seg
        s['n_mine'] = np.array([((self.lick >= a) & (self.lick < b) & (self.lick_tr == t + 1)).sum() for t, a, b in zip(s['trial'], s['t0'], s['t1'])])

    def check_segments(self):
        s = self.seg; ok = np.isfinite(s['rig'])
        out = dict(n_seg=len(s['t0']), n_to=int(s['to'].sum()), count_equal=float((s['n_mine'][ok] == s['rig'][ok]).mean()), count_within1=float((np.abs(s['n_mine'][ok] - s['rig'][ok]) <= 1).mean()))
        if self.maxerr == self.maxerr:
            rule = s['rig'] > self.maxerr * s['X']
            out['rule_to'] = float(rule[s['to']].mean()); out['rule_last'] = float(rule[s['last']].mean())
            rule2 = s['n_mine'] > self.maxerr * s['X']; out['rule_mine_to'] = float(rule2[s['to']].mean()); out['rule_mine_last'] = float(rule2[s['last']].mean())
        return out

    # ------------------------------------------------------------------ lick -> trial position
    def locate(self, t_ev, trial=None):
        """position of event times in the trial structure.  returns dict arrays: trial, phase ('wait' | 'win' | 'post' | 'other'), time to the end of its wait (wait only), end kind of that wait, tau = time from square onset"""
        t_ev = np.asarray(t_ev, float); n = len(t_ev)
        tri = np.full(n, -1); phase = np.full(n, 'other', dtype=object); tte = np.full(n, np.nan); endk = np.full(n, '', dtype=object); tau = np.full(n, np.nan); seg_i = np.full(n, -1)
        # trial = the one whose [tr0, next tr0) contains the event
        idx = np.searchsorted(self.tr0, t_ev, side='right') - 1
        jb = self.in_break(t_ev)
        for i in range(n):
            if jb[i] >= 0: phase[i] = 'break' if self.brk[jb[i]]['kind'] == 'brk' else 'tail'; continue
            t = idx[i]
            if t < 0 or t >= self.n: continue
            if t_ev[i] > self.tr1[t] + 0.5: continue                 # in the break after the last trial of a subsession
            tri[i] = t; tau[i] = t_ev[i] - self.on[t]
            if t_ev[i] < self.on[t]:
                phase[i] = 'wait'; ks = np.where((self.seg['trial'] == t) & (self.seg['t0'] <= t_ev[i] + 1e-9) & (self.seg['t1'] > t_ev[i]))[0]
                if len(ks): j = ks[0]; seg_i[i] = j; tte[i] = self.seg['t1'][j] - t_ev[i]; endk[i] = 'to' if self.seg['to'][j] else 'sq'
            elif t_ev[i] < self.on[t] + self.lw[1]: phase[i] = 'win'
            else: phase[i] = 'post'
        return dict(trial=tri, phase=phase, tte=tte, endk=endk, tau=tau, seg=seg_i)

    def classify(self, t_ev):
        """lick class of event times (bout onsets or single licks): see the module docstring"""
        L = self.locate(t_ev); cls = np.full(len(t_ev), 'other', dtype=object)
        w = L['phase'] == 'wait'
        cls[w & (L['tte'] > ANT_S)] = 'wait_early'
        cls[w & (L['tte'] <= ANT_S) & (L['endk'] == 'to')] = 'ant_to'
        cls[w & (L['tte'] <= ANT_S) & (L['endk'] == 'sq')] = 'ant_sq'
        tau = L['tau']; ph = L['phase']
        cls[(ph == 'win') & (tau >= self.lw[0])] = 'win_other'
        cls[(ph == 'win') & (tau < self.lw[0])] = 'win_early'            # < 0.1 s after the onset: before the lick window opens (not a hit lick)
        cls[(ph == 'post') & (tau < self.lw[1] + 3.0)] = 'post_a'
        cls[(ph == 'post') & (tau >= self.lw[1] + 3.0)] = 'post_b'
        cls[ph == 'break'] = 'brk'; cls[ph == 'tail'] = 'tail'
        return cls, L

    def _classify_bouts(self):
        b = self.bout_on
        self.bcls, self.bloc = self.classify(b)
        ro = np.isin(b, self.cons_on[np.isfinite(self.cons_on)])
        self.bcls[ro] = 'rew'
        self.bout_class_counts = {c: int((self.bcls == c).sum()) for c in np.unique(self.bcls)}

    def bouts(self, cls):
        return self.bout_on[self.bcls == cls]

    # ------------------------------------------------------------------ pump (reward) bookkeeping
    def pump_info(self):
        """pump time per trial (script rule: first lick of the square + RewardDelay, level 2 without a lick: onset + RewardDelay) and whether the mouse licked before it"""
        return self.rew_t, self.rew_kind


def lick_rate_grid(task, edges):
    """licks per second on the 5 Hz grid"""
    c, _ = np.histogram(task.lick, bins=edges)
    return c / DT


def load_session(name, with_lick=True):
    """xs session (xs_v2, with JawME) + task table; adds the lick rate (licks / s, 5 Hz grid) to S.beh['lick']"""
    S = xc.Session(name, with_rejected=False, verbose=False)
    T = Task(name)
    if with_lick: S.beh['lick'] = lick_rate_grid(T, S.edges)
    return S, T
