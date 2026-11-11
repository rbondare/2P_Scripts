"""Export the trial / lick / wait-segment table of lick-task sessions for the task deep dive (task_deep_dive/*.py).

usage:  python export_task_events.py <session> [<session> ...]   (builds a MATLAB wrapper from axon_analysis_task_analysis_V1.m, runs it headless, writes <OUT>\\<session>_task.mat)

The wrapper is the task script itself (settings + checks + trial structure + lick clock + bouts + pump times: lines 1 .. the line before 'kind_n =') followed by an export tail and
the script's local functions, so every number (lick clock fit, pump rule, bouts) is the one of the task script.  Added to what the 2026-11-08 export held:
  waitTimes        trials x segments  the wait that was drawn at the start of every trial and after every time-out  (cumsum = the time of every wait end, time-outs included)
  errorLicks       trials x segments  licks the rig counted in each wait segment (the rig's own count: validates the segment boundaries and the lick clock)
  timeoutHistory, maxErrorLicks, minWait, maxWait, TimeOutSound, RewardDelay, RewardMode ... from specParams
The xs export (traces + behaviour, now with JawME) is written by cross_session/export_xs.m into the same folder.
"""
import os, re, sys, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
V1 = os.path.join(REPO, 'axon analysis', 'axon_analysis_task_analysis_V1.m')
OUT = os.environ.get('TASK_OUT', r'C:\Users\rbondarenko\2P_cache\xs_v2')
MATLAB = r'C:\Program Files\MATLAB\R2024b\bin\matlab.exe'

TAIL = r"""
%% ---- export tail (task_deep_dive/export_task_events.py)
out_dir = '%OUT%';  if ~isfolder(out_dir); mkdir(out_dir); end
X = struct();
X.session = session;  X.nTr = nTr;  X.n_rec = n_rec;  X.sub_of = sub_of;  X.tr_t0 = tr_t0;  X.on_t = on_t;  X.st_end = st_end;  X.tr_t1 = tr_t1;  X.wait_s = wait_s;  X.n_to = n_to;
X.lvl_t = lvl_t;  X.lw = lw;  X.rdel = rdel;  X.rmode = rmode;  X.Tend = S.Tend;  X.P_log = S.P;  X.has_lick = double(has_lick);  X.subsessions = S.subsessions;  X.trials_per_sub = S.trials;
X.sub_t0 = sub_t0;  X.sub_t1 = sub_t1;  X.rew_t = rew_t;  X.trig_t = trig_t;  X.rew_kind = rew_kind;  X.hit = double(hit);  X.n_win = n_win;  X.lat1 = lat1;  X.cons_on = cons_on;  X.cons_n = cons_n;
X.n_pre = n_pre;  X.n_post = n_post;  X.late_on = late_on;  X.lick_t = lick_t;  X.lick_tr = lick_tr;  X.bout_on = bout_on;  X.bout_n = bout_n;
if exist('bout_off', 'var'); X.bout_off = bout_off; end
if exist('bout_iso', 'var'); X.bout_iso = double(bout_iso); end
X.sp_on = sp_on;  X.sp_tr = sp_tr;  X.sp_n = sp_n;  X.imaging_gap_t = imaging_gap_t;  X.t_first = tCa(1);  X.t_last = tCa(end);
X.waitTimes = double(S.waitTimes(1:nTr, :));
if isfield(S, 'errorLicks'); X.errorLicks = double(S.errorLicks(1:min(nTr, size(S.errorLicks, 1)), :)); end
if isfield(S, 'timeoutHistory'); X.timeoutHistory = double(S.timeoutHistory(:)); end
if isfield(S, 'StateOfAgent'); X.state_of_agent = double(S.StateOfAgent(:)); end
if isfield(S, 'levelHistory'); X.levelHistory = double(S.levelHistory(:)); end
sp = S.specParams;  fn = fieldnames(sp);
for i = 1:numel(fn)
    v = sp.(fn{i});
    if (isnumeric(v) || islogical(v)) && numel(v) <= 4; X.(['sp_' fn{i}]) = double(v); elseif (ischar(v) || isstring(v)) && numel(v) < 60; X.(['sp_' fn{i}]) = char(v); end
end
if has_beh && isfield(data.Behav, 'JawME'); X.has_jaw = 1; else; X.has_jaw = 0; end
if isfield(data, 'VideoParams') && isfield(data.VideoParams, 'Whisker') && isfield(data.VideoParams.Whisker, 'whisker_box'); X.whisker_version = 'front_top'; else; X.whisker_version = 'old'; end
save(fullfile(out_dir, [session '_task.mat']), '-struct', 'X', '-v7');
fprintf('EXPORT_TASK_OK %s: %d trials, %d licks, %d bouts, time-outs %d, waits %d segments\n', session, nTr, numel(lick_t), numel(bout_on), sum(n_to), size(X.waitTimes, 2));
try; export_xs(session, out_dir); catch ME; fprintf('EXPORT_XS_FAILED %s\n', ME.message); end
"""


def build(session):
    L = open(V1, encoding='utf8', errors='replace').read().split('\n')
    i_cut = next(i for i, l in enumerate(L) if l.startswith('kind_n = arrayfun'))
    i_fun = next(i for i, l in enumerate(L) if l.startswith('%% local functions'))
    head = L[:i_cut]
    head = [re.sub(r"^session\s*=\s*'[^']*';", "session   = '%s';" % session, l) if l.startswith('session') else l for l in head]
    assert any(session in l for l in head), 'session line not replaced'
    repo = REPO.replace('\\', '\\\\') if False else REPO
    pre = ["addpath('%s', '%s');" % (repo, os.path.join(repo, 'axon analysis', 'cross_session'))]
    return '\n'.join(pre + head) + '\n' + TAIL.replace('%OUT%', OUT) + '\n' + '\n'.join(L[i_fun:]) + '\n'


def run(session, wdir):
    os.makedirs(wdir, exist_ok=True)
    f = os.path.join(wdir, 'wrap_%s.m' % session.replace('Animal', ''))
    # MATLAB script names must be valid identifiers
    name = 'wrap_' + re.sub(r'\W', '_', session)
    f = os.path.join(wdir, name + '.m')
    open(f, 'w', encoding='utf8').write(build(session))
    log = os.path.join(wdir, name + '.log')
    t0 = time.time()
    cmd = [MATLAB, '-batch', "cd('%s'); %s" % (wdir, name)]
    with open(log, 'w') as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    txt = open(log, errors='replace').read()
    ok = 'EXPORT_TASK_OK' in txt
    print('%s: %s in %.0f s (rc %d)' % (session, 'OK' if ok else 'FAILED', time.time() - t0, r.returncode), flush=True)
    if not ok: print(txt[-3000:])
    return ok


if __name__ == '__main__':
    wdir = os.environ.get('TASK_WRAP_DIR', r'C:\Users\rbondarenko\2P_cache\task_wrap')
    for s in sys.argv[1:]:
        run(s, wdir)
