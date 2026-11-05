%% axon_analysis_alley_V1
% Axon analysis for the closed-loop corridor ('alleyCL') with visuomotor-mismatch events, aggregated with newAggregation/save_preprocessed.m.
% Same layout as axon_analysis_V3 / _task_analysis_V1: checks first (PASS / WARN / FAIL), then the results. Nothing is hard-coded to a ROI count,
% frame rate or number of events; alley_block picks the alleyCL block of a multi-block session.
%
% Stimulus (present_stimulus_addnew.m, case 'alleyCL'; the copy of the script that ran is stored next to the stim file):
%   - Two sine gratings (25 deg period, 100% contrast), left of x = 350 px with orientation 0 deg and right of it with 180 deg, a black vertical bar
%     (40 px) at x = 350 as the vanishing point: a corridor. Both halves drift with the mouse's FORWARD running speed on the ball (closed loop):
%     shift per frame = asind(ifi * v / 20 cm) / px2deg (+ rotation term, rotation_gain 1). No ball movement -> static gratings.
%     The stimulus-PC ball read is non-blocking and averaged with the previous frame (ForwardSpeed = (old + new) / 2).
%   - Mismatch events (options.specParams.mismatch): for 1 s the grating flow is multiplied by gain 0 (the flow FREEZES) while the mouse keeps running.
%     Events start only after >= 1 s of uninterrupted forward speed > 3 cm/s, 5-12 s after the end of the previous event (first 10 s of a trial and the
%     last 2 s are free), 20% are SHAM (scheduled, gated and logged in exactly the same way but the flow is left alone). Everything is logged per frame
%     in Frameinfo.alleyCL: speed, rotation, flow_scale, event id, ball-read time; and per event (events struct: onset / offset frame, gain, is_sham,
%     speed at onset, wait for running).
%   - 3 trials x 600 s were planned; recordings were stopped by hand (interrupted = 1).
%
% What the script does:
%   Part A  checks   0 settings | 1 load | 2 design as stored (events, gating, flow) | 3 timing and sync (imaging, projector pulses vs PC clock, event
%                    times on the 2P clock, ball vs stimulus-PC speed, closed-loop latency) | 4 behaviour and locomotion | 5 ROI quality and motion |
%                    6 ROI selection and groups | 7 session overview
%   Part B  results  8 event-triggered dF/F and behaviour: mismatch vs sham vs pseudo-events from a simulation of the scheduler | 9 per-event amplitudes, behaviour
%                    regression, change over the session | 10 natural stops / run onsets as a reference for 'visual flow stops' | 11 stimulus start,
%                    trial-2 start | 12 summary + findings table
%
% Reading the output: the question is whether the axons carry a signal that is specific to the flow freeze (a prediction error). Three controls are used:
%   sham events (same gating, nothing changes on the screen), pseudo-events from a SIMULATION OF THE SCHEDULER on the logged running speed (the gate
%   ties event times to the running: an event that had to wait fires 1 s after a run onset, while the speed is still rising; this alone makes speed,
%   whisking and every ROI that follows them 'respond' at the event, mismatch and sham alike), and a regression on the behaviour that follows the
%   freeze (the mouse may slow down, whisk, dilate: those drive LC axons without any prediction error).
% Outputs left in the workspace: chk (all checks), ev_t / ev_off_t (event times, 2P clock), ev_sham, roiq, signal_like, sel_roi, grp_rois, ETA / stats (tables),
% findings.

%% 0. Settings (edit only this block)
session     = 'AnimalRB22_261023_1427';
alley_block = 'auto';                                      % index into Stimuli, or 'auto' = first alleyCL block
agg_dir     = 'Z:\Group Members\Rima\Aggregated';
s2p_dir     = fullfile('Z:\Group Members\Rima\DATA_2P', session, 'suite2p', 'plane0');   % only used for the rigid-shift check
repo_root   = 'C:\Users\rbondarenko\projects\2P_Scripts';
save_figs   = false;                                       % true: write every figure as PNG (300 dpi) + vector PDF and the tables as CSV to fig_dir
fig_dir     = fullfile(repo_root, 'axon analysis', 'plots', 'alley_V1', session);

roi_sel     = 'signal'; % ROIs of part B: 'all' | 'signal' (signal-like) | 'good' (signal-like and not motion-coupled) | a list, e.g. [3 5 8]
sync_corr   = true;     % true: move the event times by the local (projector pulse - stimulus-PC clock) offset; false: PC-clock times as aggregated

% event-triggered analysis
ev_pre_s    = 4;        % aligned windows: -ev_pre_s ... +ev_post_s around the flow freeze
ev_post_s   = 8;
base_win    = [-2 -0.5];% baseline (s from onset) subtracted from every aligned trace
resp_win    = [0 2];    % main response window: the 1 s freeze + 1 s after it
halt_win    = [0 1];    % during the freeze
post_win    = [1 2.5];  % after the flow resumes
run_gate_cms = 3;       % the stimulus' running gate (specParams.mismatch.min_speed_cms), re-applied offline to the ball
excl_s      = 4;        % natural stops / run onsets keep at least this distance (s) from every logged event
stop_hi_mms = 100;      % natural stop: speed falls from > stop_hi_mms (mean of -2...-0.5 s) to < stop_lo_mms and stays there stop_hold_s
stop_lo_mms = 30;
stop_hold_s = 1;

% signal / ROI parameters (as in axon_analysis_V3 / task_analysis_V1)
neuropil_factor = 0.5;  % must equal dFFsettings.neuropil_factor of save_preprocessed
dff_window_s    = 15;
smooth_s        = 0.5;
max_interp_s    = 0.5;  % behaviour gaps up to this long are interpolated onto the imaging grid
slow_ratio_thr  = 4;    % ROI counts as signal-like if var(smoothed dFF) >= slow_ratio_thr x its white-noise expectation
mot_r_thr       = 0.15;
grp_r_thr       = 0.8;  % ROIs of the selection whose smoothed dF/F correlates above this are one group
n_perm          = 2000; % label permutations (mismatch vs sham)
n_boot          = 2000;
n_null          = 1000; % simulated realisations of the event scheduler (null distribution)
rng(1);

col = struct('dff', [0 0 0], 'mm', [0.85 0.1 0.1], 'sham', [0.5 0.5 0.5], 'null', [0.1 0.35 0.8], 'pupil', [0.85 0.15 0.15], 'whisk', [0.1 0.55 0.25], ...
    'face', [0.9 0.55 0.05], 'eye', [0.5 0.3 0.7], 'loc', [0.1 0.35 0.8], 'rigid', [0.5 0.5 0.5]);
beh_lab = struct('LocSpeed', 'running speed (mm/s)', 'PupilArea', 'pupil area (px^2)', 'WhiskIndex', 'whisking index', 'FacePadME', 'face motion', 'EyeLidOpening', 'eyelid opening', 'RigidShift', 'rigid shift (px)');
beh_col = struct('LocSpeed', col.loc, 'PupilArea', col.pupil, 'WhiskIndex', col.whisk, 'FacePadME', col.face, 'EyeLidOpening', col.eye, 'RigidShift', col.rigid);
newfig = @(nm, w, h) figure('Name', nm, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 w h], 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 8, ...
    'DefaultAxesBox', 'off', 'DefaultAxesTickDir', 'out', 'DefaultAxesLineWidth', 0.75, 'DefaultLineLineWidth', 0.75, 'DefaultTextFontName', 'Arial', 'DefaultLegendBox', 'off');

%% 1. Load file and basic integrity
addpath(repo_root, fullfile(repo_root, 'analysis_functions'), fullfile(repo_root, 'axon analysis'));
if save_figs && ~isfolder(fig_dir); mkdir(fig_dir); end
savef = @(fh, nm) save_fig(fh, fullfile(fig_dir, nm));

sev      = {'PASS', 'WARN', 'FAIL'};
grade    = @(ok, lvl) sev{1 + (~all(ok)) * lvl};                  % lvl 1: WARN when not ok, lvl 2: FAIL when not ok
addchk   = @(c, name, st, msg) [c; {name, st, msg}];
lst      = @(v) ['[' strtrim(sprintf('%g ', v)) ']'];
printchk = @(c, n0) cellfun(@(a, b, d) fprintf('  [%s] %-36s %s\n', b, a, d), c(n0+1:end, 1), c(n0+1:end, 2), c(n0+1:end, 3));
chk = cell(0, 3);  n0 = 0;

data_file = fullfile(agg_dir, [session '_preprocessed.mat']);
assert(isfile(data_file), 'Aggregated file not found: %s', data_file);
fmeta  = whos('-file', data_file);
inFile = {fmeta.name};
want   = {'CaData','TimeCa','Stimuli','Triggers','Behav','LocomotionCal','LocomotionNonCal','VideoParams','AnimalInfo','SI_Info'};
fprintf('Loading %s ...\n', data_file);
data = load(data_file, want{ismember(want, inFile)});

core = {'CaData','TimeCa','Stimuli','Triggers'};
miss = core(~ismember(core, inFile));
chk  = addchk(chk, 'file: core variables', grade(isempty(miss), 2), sprintf('%d/%d present %s', numel(core)-numel(miss), numel(core), strjoin(miss, ',')));
assert(isempty(miss), 'Missing core variables: %s', strjoin(miss, ', '));

types = {data.Stimuli.type};
alley_idx = find(strcmp(types, 'alleyCL'));
assert(~isempty(alley_idx), 'No alleyCL block in %s (blocks: %s)', session, strjoin(types, ', '));
if ischar(alley_block) || isstring(alley_block); kb = alley_idx(1); else; kb = alley_block; end
assert(strcmp(types{kb}, 'alleyCL'), 'Stimuli(%d) is %s, not alleyCL', kb, types{kb});
S  = data.Stimuli(kb);
Fi = S.Frameinfo;
has_cl = isfield(Fi, 'alleyCL') && ~isempty(Fi.alleyCL);
chk = addchk(chk, 'stimulus: closed-loop log', grade(has_cl, 2), sprintf('alleyCL block %d of %d blocks (%d alleyCL), Frameinfo.alleyCL present: %d', kb, numel(types), numel(alley_idx), has_cl));
assert(has_cl, 'This block has no Frameinfo.alleyCL (stimulus script without the closed-loop log)');
C   = Fi.alleyCL;
spp = S.specParams;
ev  = C.events;
assert(~isempty(ev), 'No mismatch events were logged in this block');
nEv = numel(ev);

nz = find(~cellfun(@isempty, {data.CaData.Ca_dFF}));
assert(~isempty(nz), 'CaData has no calcium traces');
ca   = data.CaData(nz(1));
dFF  = double(ca.Ca_dFF);  F = ca.Ca_F;  Fneu = ca.Ca_F_neuropil;  F0 = ca.Ca_baseline;
[nROI, nT] = size(dFF);
tCa  = double(data.TimeCa(1, :));
cent = ca.Ca_centroid_voxel;
dtc  = diff(tCa);  fs = 1 / median(dtc);
chk = addchk(chk, 'Ca arrays', grade(isequal(size(F), size(dFF)) && numel(tCa) == nT && size(cent, 1) == nROI, 2), sprintf('%d ROIs x %d frames, %d time stamps, %d centroids', nROI, nT, numel(tCa), size(cent, 1)));
chk = addchk(chk, 'dFF: finite', grade(all(isfinite(dFF(:))), 2), sprintf('%.4f%% NaN/Inf', 100 * mean(~isfinite(dFF(:)))));
chk = addchk(chk, 'imaging: monotonic, no gaps', grade(all(dtc > 0) && max(dtc) < 1.5 * median(dtc), 1), sprintf('%.3f Hz, max dt / median dt = %.2f, %d gap(s) > 1.5 dt', fs, max(dtc) / median(dtc), sum(dtc > 1.5 * median(dtc))));
has_beh = isfield(data, 'Behav') && ~isempty(data.Triggers.TimeCamera);
has_loc = isfield(data, 'LocomotionCal') && ~isempty(data.LocomotionCal) && isfield(data.LocomotionCal, 'Triggers_sec');
has_si  = isfield(data, 'SI_Info');
if isfield(data, 'AnimalInfo')
    ai = data.AnimalInfo;  im_date = '';
    if isfield(ai, 'Dates') && ~isempty(ai.Dates.Imaging); im_date = char(ai.Dates.Imaging(end)); end
    fprintf('\n%s | %s | imaged %s | stimulus %s (block %d)\n', ai.AnimalName, ai.Genotype, im_date, S.addnote, kb);
end
fprintf('%d ROIs, %d frames, %.0f s\n\n', nROI, nT, tCa(end) - tCa(1));
printchk(chk, n0);  n0 = size(chk, 1);

%% 2. Design as stored in the stimulus file
cols = C.perframe_cols;  cx = @(nm) find(strcmp(cols, nm), 1);
nfr  = size(C.perframe, 1);
ifi  = S.ifi;
v_cl     = C.perframe(:, cx('forward_speed_cms'));             % online forward speed on the stimulus PC (cm/s)
rot_cl   = C.perframe(:, cx('rotation_deg_per_frame'));
scale_cl = C.perframe(:, cx('flow_scale'));
id_cl    = C.perframe(:, cx('mismatch_event_id'));
tball_ms = C.perframe(:, cx('ball_read_ms'));
tcode_ms = C.perframe(:, cx('corridor_code_ms'));
flow_cl  = C.flow_px_per_frame;                                 % delivered shift per frame and region (px): [left right patch]
mmp = C.mismatch;
ev_on_f  = double([ev.onset_frame])';  ev_off_f = double([ev.offset_frame])';
ev_off_f(isnan(ev_off_f)) = nfr + 1;
ev_sham  = logical([ev.is_sham])';  ev_M = ~ev_sham;
ev_gain  = double([ev.gain])';  ev_trial = double([ev.trial])';  ev_speed_on = double([ev.speed_at_onset_cms])';
ev_wait  = double([ev.wait_s])';  ev_dur_pl = double([ev.dur_planned_s])';  ev_trunc = logical([ev.truncated])';
ev_off_f = min(ev_off_f, nfr);

tsf = double(S.TimeStimulusFrame(:));  tsf = tsf(1:nfr);            % 2P clock (s) of every presented frame, as aggregated
assert(all(isfinite(tsf)), 'TimeStimulusFrame has NaN inside the presented frames');
sfc = double(Fi.subtrial_frame_count(1:nfr));
tr_start = [1; find(diff(sfc) < 0) + 1];  tr_end = [tr_start(2:end) - 1; nfr];  nTrP = numel(tr_start);
tr_of_frame = cumsum(ismember((1:nfr)', tr_start));
b_t0 = tsf(tr_start);  b_t1 = tsf(tr_end);

n_real = sum(ev_M);  n_sham = sum(ev_sham);
chk = addchk(chk, 'design: stimulus parameters', grade(true, 1), sprintf('%s | trial %g s x %d planned | orientations %s (rows 1-3 drift direction of left / right / patch, 4-6 forward-flow weight) | alley origin %g px | grating period %g deg | wall distance %g cm | rotation gain %g', ...
    S.addnote, S.stimulus_trial_t, S.trials, mat2str(spp.grating_orientations'), spp.alley_origin, spp.grating_spatial_period, spp.grat_rad_distance_cm, spp.rotation_gain));
chk = addchk(chk, 'design: run length', grade(~S.interrupted && nTrP == S.trials, 1), sprintf('%d of %d planned trials started, %.0f s presented (planned %g s), interrupted = %d', nTrP, S.trials, nfr * ifi, S.trials * S.stimulus_trial_t, S.interrupted));
chk = addchk(chk, 'design: mismatch settings', grade(true, 1), sprintf('duration %s s, flow gain %s, interval %g-%g s, first %g s free, last %g s free, running gate %d (> %g cm/s for %g s), sham fraction %.2f', ...
    mat2str(mmp.dur), mat2str(mmp.gain), mmp.interval(1), mmp.interval(2), mmp.first_delay, mmp.end_guard, mmp.require_running, mmp.min_speed_cms, mmp.min_run_s, mmp.sham_fraction));
chk = addchk(chk, 'events: counts', grade(n_real >= 20 && n_sham >= 8, 1), sprintf('%d events (%d mismatch with gain %s, %d sham), %d per trial %s, %d truncated; rate %.1f / min. Few shams limit the sham comparison (%d shams; 20%% of %d = %.0f expected)', ...
    nEv, n_real, mat2str(unique(ev_gain(ev_M))'), n_sham, round(mean(accumarray(ev_trial, 1))), lst(accumarray(ev_trial, 1)'), sum(ev_trunc), 60 * nEv / (nfr * ifi), n_sham, nEv, mmp.sham_fraction * nEv));
% logged per-frame columns must agree with the event table
in_ev = false(nfr, 1);  ok_id = true;
for k = 1:nEv
    fr = ev_on_f(k):(ev_off_f(k) - 1);  in_ev(fr) = true;
    ok_id = ok_id && all(id_cl(fr) == ev(k).id) && (ev_on_f(k) == 1 || id_cl(ev_on_f(k) - 1) == 0);
end
in_mm = false(nfr, 1);  in_sh = in_mm;
for k = 1:nEv
    if ev_sham(k); in_sh(ev_on_f(k):ev_off_f(k) - 1) = true; else; in_mm(ev_on_f(k):ev_off_f(k) - 1) = true; end
end
ok_sc = all(scale_cl(in_mm) == 0) && all(scale_cl(~in_mm) == 1);
chk = addchk(chk, 'events: per-frame log consistent', grade(ok_id && ok_sc && all(~isnan(id_cl)), 2), sprintf('event id column marks exactly the frames onset..offset-1 of all %d events (%d frames, %.2f%% of the run); flow_scale is 0 on the frames of the %d mismatch events and 1 elsewhere (sham events keep 1)', nEv, sum(in_ev), 100 * mean(in_ev), n_real));
ev_dur_ms = 1e3 * (tsf(ev_off_f) - tsf(ev_on_f));
chk = addchk(chk, 'events: duration', grade(all(abs(ev_dur_ms / 1e3 - ev_dur_pl) < 0.05), 1), sprintf('%d frames each (%s), PC-clock duration %.0f-%.0f ms (planned %g s; the event ends at the first frame after the time elapsed)', median(ev_off_f - ev_on_f), lst(unique(ev_off_f - ev_on_f)'), min(ev_dur_ms), max(ev_dur_ms), ev_dur_pl(1)));
gap = tsf(ev_on_f(2:end)) - tsf(ev_off_f(1:end-1));  gap_ok = gap(ev_trial(2:end) == ev_trial(1:end-1));
chk = addchk(chk, 'events: spacing', grade(min(gap_ok) >= mmp.interval(1) - 0.05, 1), sprintf('end of event to next onset: %.1f-%.1f s, median %.1f s (planned wait %g-%g s + wait for running); first onset %.1f s after the run start (free %g s)', min(gap_ok), max(gap_ok), median(gap_ok), mmp.interval(1), mmp.interval(2), tsf(ev_on_f(1)) - tsf(1), mmp.first_delay));
% the stimulus' gate, re-applied to the logged online speed
w_run = round(mmp.min_run_s / ifi);
gate_on = arrayfun(@(f) all(v_cl(max(1, f - w_run + 1):f) > mmp.min_speed_cms), ev_on_f);
chk = addchk(chk, 'events: gate satisfied at onset', grade(mean(gate_on) > 0.95, 1), sprintf('%d/%d events have the logged forward speed > %g cm/s in the %g s before the first frame of the event (the speed at the first frame itself is %.1f-%.1f cm/s, median %.1f); wait for running at arming: %d events with 0 s, max %.1f s', ...
    sum(gate_on), nEv, mmp.min_speed_cms, mmp.min_run_s, min(ev_speed_on), max(ev_speed_on), median(ev_speed_on), sum(ev_wait == 0), max(ev_wait)));
% delivered flow and the closed-loop rule
px2deg = C.px2deg_factor;
w_fwd  = spp.grating_orientations(end-2:end);  w_fwd = w_fwd(:)';
fwd_pred = asind(ifi * v_cl * w_fwd / spp.grat_rad_distance_cm) / px2deg;     % forward part of the shift, [nfr x 3]
fwd_del  = (flow_cl(:, 1) + flow_cl(:, 2)) / 2;  rot_del = (flow_cl(:, 1) - flow_cl(:, 2)) / 2;
free = ~in_ev & isfinite(fwd_del) & isfinite(fwd_pred(:, 1));
rr = corr(fwd_del(free), fwd_pred(free, 1));
chk = addchk(chk, 'flow: delivered vs closed-loop rule', grade(rr > 0.99, 1), sprintf('delivered shift (mean of left / right) vs asind(ifi*v/%g cm)/px2deg of the logged speed: r = %.3f; during mismatch events mean |shift| = %.4f px/frame (outside events %.2f), during sham events %.2f px/frame (logged speed %.1f cm/s)', ...
    spp.grat_rad_distance_cm, rr, mean(abs(fwd_del(in_mm)), 'omitnan'), mean(abs(fwd_del(free))), mean(abs(fwd_del(in_sh)), 'omitnan'), mean(v_cl(in_sh), 'omitnan')));
rot_frac = std(rot_del(free)) / std(fwd_del(free) + rot_del(free));
r_rf = corr(rot_del(free), fwd_del(free));
bl = polyfit(v_cl(free), flow_cl(free, 1), 1);  br = polyfit(v_cl(free), flow_cl(free, 2), 1);
chk = addchk(chk, 'flow: rotation term', grade(rot_frac < 0.25 && abs(r_rf) < 0.3, 1), sprintf(['flow = forward term + rotation term (left and right halves drift in opposite directions): SD %.2f vs %.2f px/frame (rotation = %.0f%% of the SD of the left shift); ' ...
    'rotation is correlated with the forward term, r = %.2f (a rotation estimate that follows running comes from the dead second ball sensor); slope of the delivered shift on speed: left %.3f, right %.3f px/frame per cm/s (gain ratio right / left %.2f); mean |shift| left %.2f, right %.2f px/frame; the two halves move against each other (opposite sign) in %.0f%% of the frames'], ...
    std(rot_del(free)), std(fwd_del(free)), 100 * rot_frac, r_rf, bl(1), br(1), br(1) / bl(1), mean(abs(flow_cl(free, 1))), mean(abs(flow_cl(free, 2))), 100 * mean(sign(flow_cl(free, 1)) .* sign(flow_cl(free, 2)) < 0)));
long_fr = diff(tsf) > 1.5 * ifi;
chk = addchk(chk, 'stimulus: frame timing', grade(C.mean_frame_interval_s < 1.05 * ifi && mean(long_fr) < 0.01, 1), sprintf('mean interval %.3f ms (ifi %.3f), %d frames > 1.5 ifi (%.3f%%), longest %.0f ms, corridor code per frame %.2f ms (max %.1f ms), ball read %.2f ms (max %.1f ms)', ...
    1e3 * C.mean_frame_interval_s, 1e3 * ifi, sum(long_fr), 100 * mean(long_fr), 1e3 * max(diff(tsf)), median(tcode_ms), max(tcode_ms), median(tball_ms), max(tball_ms)));
printchk(chk, n0);  n0 = size(chk, 1);

%% 3. Timing and synchronisation
% 3a  projector pulses against the stimulus-PC clock. save_preprocessed fell back to the PC clock (TimeStimulusFrame = frame time - Clockstart) when the
%     red-frame pulses deviated from a straight line; here every red frame is paired with its nearest pulse (unique while |offset| < half a red period),
%     which also gives the local offset at every event.
nth   = S.red_sync_nth_frame;
red_f = (1:nth:nfr)';
P     = double(data.Triggers.TimeProjector(:));
P     = P(P > tsf(1) - 1 & P < tsf(nfr) + 1);
jn    = interp1(P, (1:numel(P))', tsf(red_f), 'nearest', 'extrap');
d_red = P(jn) - tsf(red_f);                                       % pulse minus PC-clock time of the red frame (s)
dup   = [false; diff(jn) == 0];
d_ok  = abs(d_red) < 0.45 * nth * ifi & ~dup;
tr_red = tr_of_frame(red_f);
fit_p = nan(nTrP, 2);  fit_res = nan(size(d_red));
for k = 1:nTrP
    ii = d_ok & tr_red == k;
    if nnz(ii) > 10; fit_p(k, :) = polyfit(tsf(red_f(ii)), d_red(ii), 1); fit_res(ii) = d_red(ii) - polyval(fit_p(k, :), tsf(red_f(ii))); end
end
ri = S.RedSyncInfo;
chk = addchk(chk, 'sync: as aggregated', grade(S.RedFrameSynchronized, 1), sprintf('method %s, %d pulses for %d red frames (%d lost): %s. The aggregated TimeStimulusFrame is the PC clock', ri.method, ri.n_pulses, ri.n_expected, ri.n_lost, ri.reason));
chk = addchk(chk, 'sync: pulses paired with red frames', grade(mean(d_ok) > 0.97, 1), sprintf('%d/%d red frames have a pulse within +-%.0f ms (%d without; lost pulses and the first frames after a stall)', sum(d_ok), numel(d_ok), 1e3 * 0.45 * nth * ifi, sum(~d_ok)));
chk = addchk(chk, 'sync: PC clock vs pulses', grade(max(abs(d_red(d_ok))) < 0.033, 1), sprintf('pulse - PC time %.1f..%.1f ms (median %.1f); within a trial a straight line with %s ms/min (residual SD %.2f ms, max %.1f ms)', ...
    1e3 * min(d_red(d_ok)), 1e3 * max(d_red(d_ok)), 1e3 * median(d_red(d_ok)), lst(round(60e3 * fit_p(:, 1)', 2)), 1e3 * std(fit_res, 'omitnan'), 1e3 * max(abs(fit_res))));
if nTrP > 1
    t_chg = b_t0(2);
    step_ms = 1e3 * (polyval(fit_p(2, :), t_chg) - polyval(fit_p(1, :), t_chg));
    chk = addchk(chk, 'sync: step at the trial change', grade(abs(step_ms) < 33, 1), sprintf('the offset steps by %.1f ms between the fits of trial 1 and trial 2 at t = %.1f s (one extra / missing stimulus frame is %.1f ms)', step_ms, t_chg, 1e3 * ifi));
end
loc_off = nan(nEv, 1);  loc_off_off = loc_off;
for k = 1:nEv
    ii = abs(red_f - ev_on_f(k)) <= 8 * nth & d_ok;   if nnz(ii) >= 5; loc_off(k) = median(d_red(ii)); end
    ii = abs(red_f - ev_off_f(k)) <= 8 * nth & d_ok;  if nnz(ii) >= 5; loc_off_off(k) = median(d_red(ii)); end
end
chk = addchk(chk, 'sync: local offset at the events', grade(~any(isnan(loc_off)), 1), sprintf('%d/%d events have >= 5 pulses within +-8 red frames; offset %.1f..%.1f ms. Using it: %d, so the event times are accurate to ~1 ms instead of up to %.0f ms', ...
    sum(isfinite(loc_off)), nEv, 1e3 * min(loc_off), 1e3 * max(loc_off), sync_corr, 1e3 * max(abs(loc_off))));
ev_t     = tsf(ev_on_f)  + sync_corr * fillmissing_zero(loc_off);
ev_off_t = tsf(ev_off_f) + sync_corr * fillmissing_zero(loc_off_off);
chk = addchk(chk, 'events: inside imaging', grade(all(ev_t - ev_pre_s >= tCa(1) & ev_t + ev_post_s <= tCa(end)), 1), sprintf('event windows [-%g, +%g s] inside the imaging record %.1f-%.1f s: %d/%d', ev_pre_s, ev_post_s, tCa(1), tCa(end), sum(ev_t - ev_pre_s >= tCa(1) & ev_t + ev_post_s <= tCa(end)), nEv));
acq = data.Triggers.Stimulus_acq_timestamps;
chk = addchk(chk, 'imaging: acquisitions', grade(true, 1), sprintf('%d acquisition(s) [start end] s: %s; stimulus frames %.1f-%.1f s', size(acq, 1), mat2str(round(acq, 1)), tsf(1), tsf(nfr)));

% 3b  ball: complete log? second sensor?
loc_ok = false;  tt = [];  fw = [];  dtl = NaN;  spd_g = nan(1, nT);  spd_fast = [];
if has_loc
    L  = data.LocomotionCal(min(kb, numel(data.LocomotionCal)));
    nV = size(L.Values, 1);  nTr_ball = numel(L.Triggers_sec);
    loc_ok = abs(nV - nTr_ball) <= 2 && nV > 1000;
    chk = addchk(chk, 'ball: readings vs triggers', grade(loc_ok, 1), sprintf('%d readings, %d ball triggers (%.1f Hz); Triggers_sec %.1f-%.1f s. A reading count far below the trigger count = truncated log (unusable)', nV, nTr_ball, nTr_ball / (L.Triggers_sec(end) - L.Triggers_sec(1)), L.Triggers_sec(1), L.Triggers_sec(end)));
    if isfield(data, 'LocomotionNonCal') && isfield(data.LocomotionNonCal, 'Values')
        V = double(data.LocomotionNonCal(min(kb, numel(data.LocomotionNonCal))).Values(:, 1:4));  nzf = mean(V ~= 0, 1);
    else
        V = double(L.Values(:, 1:4));  nzf = mean(V ~= 0, 1);
    end
    chk = addchk(chk, 'ball: sensors alive', grade(all(nzf > 0.2), 1), sprintf('fraction of non-zero readings per column [x1 y1 x2 y2]: %s. A dead second sensor leaves forward speed from the first sensor only; sideways / rotation estimates are then unreliable', mat2str(round(nzf, 3))));
end
if loc_ok
    n   = min([numel(L.Triggers_sec), numel(L.Forward_mm)]);
    tt  = double(L.Triggers_sec(1:n))';  fw = double(L.Forward_mm(1:n))';  dtl = median(diff(tt));
    spd_fast = movmean(fw, max(1, round(0.1 / dtl))) / dtl;                    % mm/s, 0.1 s
    spd_s    = movmean(fw, max(1, round(0.25 / dtl))) / dtl;                   % mm/s, 0.25 s (statistics)
    ii = tCa >= tt(1) & tCa <= tt(end);
    spd_g(ii) = interp1(tt, spd_s, tCa(ii), 'linear');
    % 3c  ball vs the speed the stimulus PC used: correlation, scale and the closed-loop latency (ball motion -> speed seen by the stimulus code)
    lags = -0.1:0.01:0.4;  r_lag = nan(size(lags));
    for i = 1:numel(lags)
        a = interp1(tt, spd_fast, tsf - lags(i), 'linear', NaN);
        okk = isfinite(a) & ~in_ev;  r_lag(i) = corr(v_cl(okk), a(okk));
    end
    [r_best, ib] = max(r_lag);  lag_best = lags(ib);
    a = interp1(tt, spd_fast, tsf - lag_best, 'linear', NaN);  okk = isfinite(a) & ~in_ev;
    slope = (a(okk) \ (10 * v_cl(okk)));
    chk = addchk(chk, 'ball: stimulus-PC speed vs ball', grade(r_best > 0.8, 1), sprintf('logged online forward speed vs the ball speed on the 2P clock (0.1 s): r = %.2f at a lag of %.0f ms (r = %.2f at 0 ms); stimulus speed = %.2f x ball speed; closed-loop latency ball -> speed used for the flow ~ %.0f ms plus display (1-2 frames, ~%.0f ms)', ...
        r_best, 1e3 * lag_best, r_lag(lags == 0), slope, 1e3 * lag_best, 2e3 * ifi));
    in_run = tCa >= tsf(1) & tCa <= tsf(nfr);
    chk = addchk(chk, 'locomotion: speed', grade(median(spd_g(in_run), 'omitnan') < 250, 1), sprintf('forward speed during the stimulus: median %.0f mm/s, mean %.0f, p90 %.0f, p99 %.0f; > %g cm/s %.0f%% of the time; backwards (< -30 mm/s) %.1f%%. WARN above a median of 250 mm/s: fast for a mouse on a ball, check the ball calibration (the closed-loop gain depends on it)', ...
        median(spd_g(in_run), 'omitnan'), mean(spd_g(in_run), 'omitnan'), prctile(spd_g(in_run), 90), prctile(spd_g(in_run), 99), run_gate_cms, 100 * mean(spd_g(in_run) > 10 * run_gate_cms, 'omitnan'), 100 * mean(spd_g(in_run) < -30, 'omitnan')));
else
    chk = addchk(chk, 'locomotion: present', grade(false, 1), 'ball log unusable or missing: running-speed based analyses (pseudo-events, regression, natural stops) are skipped');
end
printchk(chk, n0);  n0 = size(chk, 1);

%% 4. Behaviour on the imaging grid (arousal proxies)
beh_ca = struct();  beh_names = {};
if loc_ok; beh_ca.LocSpeed = spd_g;  beh_names{end+1} = 'LocSpeed'; end
tcam = [];
if has_beh
    tcam = double(data.Triggers.TimeCamera(:))';
    Bh = data.Behav;
    nCam = numel(tcam);
    bf = fieldnames(Bh);  brow = cellfun(@(f) size(Bh.(f), 1), bf);
    chk = addchk(chk, 'Behav: rows == camera triggers', grade(all(brow == nCam), 2), sprintf('%d camera triggers (%.1f Hz), %d Behav fields, mismatching: %s', nCam, 1 / median(diff(tcam)), numel(bf), strjoin(bf(brow ~= nCam)', ',')));
    in_stim = tcam >= tsf(1) & tcam <= tsf(nfr);
    for nm = {'PupilArea', 'WhiskIndex', 'FacePadME', 'EyeLidOpening'}
        nm = nm{1};
        if ~isfield(Bh, nm); continue; end
        x = double(Bh.(nm)(:))';
        if strcmp(nm, 'PupilArea') && isfield(Bh, 'PupilBad'); x(logical(Bh.PupilBad(:))') = NaN; end
        nf = mean(~isfinite(x(in_stim)));
        chk = addchk(chk, ['behaviour: ' nm], grade(nf < 0.10, 1), sprintf('%.1f%% missing during the stimulus, median %.3g, p1-p99 %.3g-%.3g', 100 * nf, median(x(in_stim), 'omitnan'), prctile(x(in_stim), 1), prctile(x(in_stim), 99)));
        ok = isfinite(x);  tv = tcam(ok);
        xi = interp1(tv, x(ok), tCa, 'linear', NaN);
        ix = interp1(tv, 1:numel(tv), tCa, 'nearest', 'extrap');
        xi(abs(tv(ix) - tCa) > max_interp_s) = NaN;
        beh_ca.(nm) = xi;  beh_names{end+1} = nm; %#ok<SAGROW>
    end
    if isfield(Bh, 'WhiskOn')
        chk = addchk(chk, 'whisking: bouts resolvable', grade(mean(Bh.WhiskOn(in_stim)) < 0.9, 1), sprintf('whisking in %.0f%% of the stimulus frames: %s', 100 * mean(Bh.WhiskOn(in_stim)), 'WhiskIndex is used as a continuous measure'));
    end
end
printchk(chk, n0);  n0 = size(chk, 1);

%% 5. ROI signal quality and motion (checks)
Fc    = F - neuropil_factor * Fneu;
w_sm  = max(1, round(smooth_s * fs));
dsm   = movmean(dFF, w_sm, 2);
med_s = median(dsm, 2);  dev = dsm - med_s;  neg = dev < 0;
sig_neg  = max(sqrt(sum(dev .^ 2 .* neg, 2) ./ max(sum(neg, 2), 1)), eps('single'));
dd       = diff(dFF, 1, 2);
hf_noise = 1.4826 * median(abs(dd - median(dd, 2)), 2) / sqrt(2);
slow_ratio = var(dsm, 0, 2) ./ (hf_noise .^ 2 / w_sm);
up = dev > 3 * sig_neg;  dn = dev < -3 * sig_neg;
ev_pos = sum(diff([false(nROI, 1) up], 1, 2) == 1, 2);
ev_neg = sum(diff([false(nROI, 1) dn], 1, 2) == 1, 2);
frac_Fc_le0 = mean(Fc <= 0, 2);
rows = 1:min(nROI, 30);
Fcr  = double(Fc(rows, :)) - double(F0(rows, :));
den  = movmedian(Fcr, 2 * ceil(dff_window_s / median(dtc)) + 1, 2);
re_err = median(abs(Fcr ./ den - 1 - dFF(rows, :)), 2);
chk = addchk(chk, 'dFF: reproduces from F/Fneu/baseline', grade(max(re_err) < 1e-3, 1), sprintf('median |recomputed - stored| <= %.1e (neuropil %.2f, +-%g s)', max(re_err), neuropil_factor, dff_window_s));
den_bad = false(nROI, 1);  den_bad(rows) = any(den <= 0, 2);
bad_F       = frac_Fc_le0 > 0.01 | den_bad;
noise_like  = slow_ratio < slow_ratio_thr | ev_pos < 2 * max(ev_neg, 1);
signal_like = ~bad_F & ~noise_like;
roiq = table((1:nROI)', cent(:, 1), cent(:, 2), arrayfun(@(c) c.npix, ca.Ca_cellstats(:)), median(F, 2), 100 * frac_Fc_le0, std(dFF, 0, 2), hf_noise, slow_ratio, ev_pos, ev_neg, signal_like, ...
    'VariableNames', {'roi', 'y', 'x', 'npix', 'F_med', 'pct_Fc_le0', 'std_dff', 'hf_noise', 'slow_ratio', 'ev_pos', 'ev_neg', 'signal_like'});
disp(roiq);
chk = addchk(chk, 'ROIs: signal-like', grade(sum(signal_like) >= max(3, 0.25 * nROI), 1), sprintf('%d/%d ROIs pass (slow_ratio >= %g, ev_pos >= 2*ev_neg, Fc > 0). Noise-like: %s. Pure noise gives std(dFF) ~0.71 under this dF/F definition: median std %.2f, so single-ROI signal-to-noise is low and group means are used', ...
    sum(signal_like), nROI, slow_ratio_thr, lst(find(noise_like & ~bad_F)'), median(std(dFF, 0, 2))));
% motion
rigid = nan(1, nT);  have_motion = false;  mot_flag = false(nROI, 1);  motion_r = nan(nROI, 1);
fall_file = fullfile(s2p_dir, 'Fall.mat');
px_um = NaN;  if has_si; px_um = abs(diff(data.SI_Info.hRoiManager.imagingFovUm([1 2]))) / data.SI_Info.hRoiManager.pixelsPerLine; end
if isfile(fall_file)
    o = load(fall_file, 'ops');  xo = double(o.ops.xoff(:))';  yo = double(o.ops.yoff(:))';
    chk = addchk(chk, 'motion: ops length', grade(numel(xo) == nT, 1), sprintf('%d motion samples vs %d frames', numel(xo), nT));
    if numel(xo) == nT
        have_motion = true;
        rigid = hypot(xo - median(xo), yo - median(yo));
        mot_lim_px = 1 / px_um;  if ~isfinite(mot_lim_px); mot_lim_px = 6; end
        chk = addchk(chk, 'motion: rigid shift', grade(prctile(rigid, 99) < mot_lim_px, 1), sprintf('p99 %.1f px (%.2f um), max %.1f px, x %d..%d, y %d..%d px; WARN above %.1f px (1 um)', prctile(rigid, 99), prctile(rigid, 99) * px_um, max(rigid), min(xo), max(xo), min(yo), max(yo), mot_lim_px));
        mv = movmean(rigid, w_sm);  mz = (mv - mean(mv)) / std(mv);
        dz = (dsm - mean(dsm, 2)) ./ std(dsm, 0, 2);
        motion_r = dz * mz' / (nT - 1);
        null_r = zeros(nROI, 200);
        for p = 1:200; null_r(:, p) = dz * circshift(mz, randi([round(30 * fs), nT - round(30 * fs)]))' / (nT - 1); end
        mot_flag = abs(motion_r) > max(prctile(abs(null_r), 99, 2), mot_r_thr);
        chk = addchk(chk, 'motion: ROIs coupled to the rigid shift', grade(~any(mot_flag & signal_like), 1), sprintf('%d ROIs follow the rigid shift beyond the shuffle (|r| > p99 and > %.2f); signal-like among them: %s. Median |r| %.2f. The mouse runs most of the time, so coupling to running and to the shift cannot be separated by selection; both are regressors in part B', sum(mot_flag), mot_r_thr, lst(find(mot_flag & signal_like)'), median(abs(motion_r))));
        beh_ca.RigidShift = rigid;  beh_names{end+1} = 'RigidShift';
    end
else
    chk = addchk(chk, 'motion: suite2p ops', grade(false, 1), sprintf('not checked: %s not found', fall_file));
end
printchk(chk, n0);  n0 = size(chk, 1);

%% 6. ROI selection and groups
good_roi = find(signal_like & ~mot_flag)';
if ischar(roi_sel) || isstring(roi_sel)
    switch lower(char(roi_sel))
        case 'all';    sel_roi = 1:nROI;
        case 'signal'; sel_roi = find(signal_like)';
        case 'good';   sel_roi = good_roi;  if isempty(sel_roi); sel_roi = find(signal_like)'; end
        otherwise;     error('roi_sel must be ''all'', ''signal'', ''good'' or a list of ROI numbers');
    end
else
    sel_roi = unique(double(roi_sel(:)'));
end
if isempty(sel_roi); sel_roi = 1:nROI; end
nSel = numel(sel_roi);
% ROIs on one axon share their events: group by correlation of the smoothed dF/F (group 1 = largest)
if nSel >= 2; gid = cluster(linkage(pdist(dsm(sel_roi, :), 'correlation'), 'average'), 'cutoff', 1 - grp_r_thr, 'criterion', 'distance'); else; gid = ones(nSel, 1); end
gid = gid(:)';  gsz = accumarray(gid(:), 1)';  gmn = accumarray(gid(:), sel_roi(:), [], @min)';
[~, gord] = sort(-1000 * gsz + gmn);  [~, gid] = ismember(gid, gord);
ng = max(gid);
grp_rois = arrayfun(@(g) sel_roi(gid == g), 1:ng, 'UniformOutput', false);
grp_lab  = cellfun(@(r) ['ROI ' regexprep(num2str(r), '\s+', ',')], grp_rois, 'UniformOutput', false);
% control: the noise-like ROIs outside the selection (what they show at an event is common-mode: motion, neuropil, z-shift)
ctl_roi = find(noise_like & ~bad_F & ~ismember((1:nROI)', sel_roi))';
grp_mean = cell2mat(cellfun(@(r) mean(dFF(r, :), 1), grp_rois(:), 'UniformOutput', false));
grp_all  = [mean(dFF(sel_roi, :), 1); grp_mean];   glab_all = [{sprintf('mean of %d selected ROIs', nSel)}, grp_lab];
if ~isempty(ctl_roi); grp_all = [grp_all; mean(dFF(ctl_roi, :), 1)]; glab_all = [glab_all, {['control: noise-like ROI ' regexprep(num2str(ctl_roi), '\s+', ',')]}]; end
ngr = size(grp_all, 1);
chk = addchk(chk, 'ROI selection', grade(true, 1), sprintf('%d/%d ROIs selected (roi_sel = %s): %s; %d groups (r > %g): %s; control group %s', nSel, nROI, strtrim(evalc('disp(roi_sel)')), lst(sel_roi), ng, grp_r_thr, strjoin(grp_lab, ' | '), lst(ctl_roi)));
printchk(chk, n0);  n0 = size(chk, 1);

%% 7. Session overview
nB = numel(beh_names);
fh = newfig('Session overview', 18, 15);
tl = tiledlayout(fh, 2 + nB, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, [strrep(session, '_', '\_') ' | alleyCL'], 'FontSize', 9, 'FontWeight', 'normal');
ax = gobjects(0);
ax(1) = nexttile(tl);  hold on;
for k = 1:nTrP; patch([b_t0(k) b_t1(k) b_t1(k) b_t0(k)], [0 0 1 1], [0.93 0.93 0.93], 'EdgeColor', 'none'); end
plot(repmat(ev_t(ev_M)', 2, 1), repmat([0; 1], 1, n_real), 'Color', col.mm, 'LineWidth', 0.6);
plot(repmat(ev_t(ev_sham)', 2, 1), repmat([0; 1], 1, n_sham), 'Color', [0 0 0], 'LineWidth', 0.6);
ylim([0 1]);  yticks([]);  ylabel('events');
text(0.005, 0.5, 'red = mismatch (flow freezes 1 s), black = sham', 'Units', 'normalized', 'FontSize', 6);
ax(2) = nexttile(tl);  plot(tCa, mean(dFF(sel_roi, :), 1), 'Color', [0.7 0.7 0.7]);  hold on;  plot(tCa, movmean(mean(dFF(sel_roi, :), 1), round(fs)), 'Color', col.dff);  ylabel('mean dF/F');
for i = 1:numel(beh_names)
    ax(end+1) = nexttile(tl); %#ok<SAGROW>
    plot(tCa, beh_ca.(beh_names{i}), 'Color', beh_col.(beh_names{i}));  ylabel(beh_lab.(beh_names{i}), 'FontSize', 6);
end
xlabel('session time (s)');
linkaxes(ax, 'x');  xlim(ax(1), [tCa(1) tCa(end)]);  set(ax(1:end-1), 'XTickLabel', []);
if save_figs; savef(fh, '01_overview'); end

%% 8. Event-triggered dF/F and behaviour: mismatch vs sham vs simulated-scheduler pseudo-events
% rows analysed: every ROI, the selection mean / groups / control, and the behaviour traces
Z = [dFF; grp_all];  row_name = [arrayfun(@(r) sprintf('ROI %d', r), 1:nROI, 'UniformOutput', false), glab_all];
row_roi = 1:nROI;  row_grp = nROI + (1:ngr);
for i = 1:numel(beh_names); Z = [Z; beh_ca.(beh_names{i})]; row_name{end+1} = beh_names{i}; end %#ok<SAGROW>
row_beh = nROI + ngr + (1:numel(beh_names));
nY = size(Z, 1);
tg = (-ev_pre_s * fs:ev_post_s * fs) / fs;
bi = tg >= base_win(1) & tg < base_win(2);
ri_ = tg >= resp_win(1) & tg <= resp_win(2);  hi_ = tg >= halt_win(1) & tg < halt_win(2);  pi_ = tg >= post_win(1) & tg < post_win(2);
A = align_rows(Z, tCa, ev_t, tg);                              % [nTg x nEv x nY]
A = A - mean(A(bi, :, :), 1, 'omitnan');
amp_of = @(AA, wi) reshape(mean(AA(wi, :, :), 1, 'omitnan'), size(AA, 2), size(AA, 3));      % [nEv x nY]
AMP  = amp_of(A, ri_);  AMPh = amp_of(A, hi_);  AMPp = amp_of(A, pi_);

% Pseudo-events by SIMULATING THE SCHEDULER on the logged speed. The gate ties the event times to the running: when the mouse is not running steadily at the planned
% time, the event waits and fires exactly min_run_s after the next run onset, i.e. while the speed is rising. Pseudo-events at random gate-satisfying times, or matched
% on the speed history (tried first), do not reproduce that and show a 'response' in speed, whisking and every ROI that follows them, in the real and in the sham events
% alike. Here the stimulus rule (first_delay, random interval after the end of the previous event, speed > gate for >= min_run_s, last allowed start) is run on the same
% logged online speed with new random intervals: n_null realisations, n_real events of each realisation are the pseudo-mismatch events. Under 'the freeze does nothing'
% this is an exact randomisation test (as long as the freeze does not change the running the scheduler reads).
have_null = true;
w_run = round(mmp.min_run_s / ifi) + 1;                                           % frames of uninterrupted running before an event
okrun = movmin(v_cl, [w_run - 1 0]) > mmp.min_speed_cms;  okrun(isnan(okrun)) = false;
nxt = inf(nfr, 1);  nxt(okrun) = find(okrun);
nxt = flipud(cummin(flipud(nxt)));                                                % first frame >= this one that satisfies the gate (Inf if none)
dur_fr = round(median(ev_off_f - ev_on_f));
last_f = tr_start + round((S.stimulus_trial_t - mmp.end_guard - max(mmp.dur)) / ifi);
sim_ev = cell(n_null, 1);  sim_wait = sim_ev;
for d = 1:n_null; [sim_ev{d}, sim_wait{d}] = sim_scheduler(nxt, tr_start, tr_end, last_f, ifi, mmp, dur_fr); end
n_sim = cellfun(@numel, sim_ev);
% 2P time of every stimulus frame with the local projector offset (the PC-clock time of a frame + offset interpolated between the red frames)
if sync_corr && nnz(d_ok) > 10
    off_all = interp1(tsf(red_f(d_ok)), movmedian(d_red(d_ok), 9), tsf, 'linear', 'extrap');
else
    off_all = zeros(nfr, 1);
end
tf_all = tsf + off_all;
fr_idx = min(max(round(interp1(tCa, 1:nT, tf_all, 'linear', 'extrap')), 1), nT);
% amplitude (response window minus baseline) of every row for every stimulus frame as event time, from cumulative sums
Zn = Z;  Zn(~isfinite(Zn)) = 0;
cs = [zeros(nY, 1), cumsum(Zn, 2)];  cn = [zeros(nY, 1), cumsum(double(isfinite(Z)), 2)];
rs = round(resp_win * fs);  bs_ = round(base_win * fs);
amp_f = win_mean(cs, cn, fr_idx, rs(1), rs(2)) - win_mean(cs, cn, fr_idx, bs_(1), bs_(2) - 1);      % [nfr x nY]
obs_f = mean(amp_f(ev_on_f(ev_M), :), 1, 'omitnan');
nM = n_real;
null_amp = nan(n_null, nY);
for d = 1:n_null
    e = sim_ev{d};  if numel(e) < nM; continue; end
    null_amp(d, :) = mean(amp_f(e(randperm(numel(e), nM)), :), 1, 'omitnan');
end
n_tr = min(300, n_null);  null_tr = nan(n_tr, numel(tg), nY);
for d = 1:n_tr
    e = sim_ev{d};  if numel(e) < nM; continue; end
    Ad = align_rows(Z, tCa, tf_all(e(randperm(numel(e), nM))), tg);  Ad = Ad - mean(Ad(bi, :, :), 1, 'omitnan');
    null_tr(d, :, :) = reshape(mean(Ad, 2, 'omitnan'), 1, numel(tg), nY);
end
% does the simulation reproduce the selection of the real events?
vm = @(ef) mean(v_cl(min(max(ef(:) + (-w_run:w_run), 1), nfr)), 1, 'omitnan');
vr = vm(ev_on_f);  vs_ = vm(vertcat(sim_ev{1:min(200, n_null)}));
pre_i = 1:w_run;  post_i = w_run + 1 + (0:w_run);
sw_ = vertcat(sim_wait{1:min(200, n_null)});
chk = addchk(chk, 'pseudo-events: scheduler simulation', grade(abs(mean(ev_wait == 0) - mean(cellfun(@(w) mean(w == 0), sim_wait))) < 0.15 && abs(mean(n_sim) - nEv) < 0.25 * nEv, 1), ...
    sprintf('%d realisations of the scheduler on the logged speed: %.0f events per run (real %d), %.0f%% fire at once (real %.0f%%), the others wait %.1f s for running (real %.1f s, median); speed in the 1 s before / after the first frame: simulated %.1f / %.1f cm/s, real %.1f / %.1f cm/s. The gate is satisfied in %.0f%% of the frames while the speed is > %g cm/s in %.0f%%', ...
    n_null, mean(n_sim), nEv, 100 * mean(cellfun(@(w) mean(w == 0), sim_wait)), 100 * mean(ev_wait == 0), median(sw_(sw_ > 0)), median(ev_wait(ev_wait > 0)), ...
    mean(vs_(pre_i)), mean(vs_(post_i)), mean(vr(pre_i)), mean(vr(post_i)), 100 * mean(okrun), mmp.min_speed_cms, 100 * mean(v_cl > mmp.min_speed_cms, 'omitnan')));
chk = addchk(chk, 'pseudo-events: fast vs exact amplitude', grade(max(abs(obs_f - mean(AMP(ev_M, :), 1, 'omitnan')) ./ max(std(AMP(ev_M, :), 0, 1, 'omitnan') / sqrt(n_real), eps)) < 0.5, 1), sprintf('amplitudes from cumulative sums (used for the null) vs the interpolated ones: max difference %.2f standard errors', max(abs(obs_f - mean(AMP(ev_M, :), 1, 'omitnan')) ./ max(std(AMP(ev_M, :), 0, 1, 'omitnan') / sqrt(n_real), eps))));
printchk(chk, n0);  n0 = size(chk, 1);
% statistics per row: mismatch mean (bootstrap CI), sham mean, mismatch - sham (label permutation), pseudo-event null (z, p)
st = repmat(struct('mm', NaN, 'mm_lo', NaN, 'mm_hi', NaN, 'sham', NaN, 'diff', NaN, 'p_perm', NaN, 'null_mu', NaN, 'null_sd', NaN, 'z', NaN, 'p_null', NaN), nY, 1);
mm_i = find(ev_M);  sh_i = find(ev_sham);
for y = 1:nY
    x = AMP(mm_i, y);
    st(y).mm = mean(x, 'omitnan');
    ci = boot_ci_mean(x, n_boot);  st(y).mm_lo = ci(1);  st(y).mm_hi = ci(2);
    st(y).sham = mean(AMP(sh_i, y), 'omitnan');
    st(y).diff = st(y).mm - st(y).sham;
    pd = zeros(n_perm, 1);
    for p = 1:n_perm; r = randperm(nEv); pd(p) = mean(AMP(r(1:n_real), y), 'omitnan') - mean(AMP(r(n_real+1:end), y), 'omitnan'); end
    st(y).p_perm = (1 + sum(abs(pd) >= abs(st(y).diff))) / (n_perm + 1);
    if have_null
        nv = null_amp(isfinite(null_amp(:, y)), y);
        st(y).null_mu = mean(nv);  st(y).null_sd = std(nv);
        st(y).z = (obs_f(y) - st(y).null_mu) / st(y).null_sd;
        st(y).p_null = (1 + sum(abs(nv - st(y).null_mu) >= abs(obs_f(y) - st(y).null_mu))) / (numel(nv) + 1);
    end
end
stats = struct2table(st);
stats = addvars(stats, row_name(:), 'Before', 1, 'NewVariableNames', 'name');
fprintf('\nEvent-triggered amplitude (mean over %g-%g s minus baseline %g...%g s), mismatch n = %d, sham n = %d, simulated-scheduler pseudo-event draws %d (null_mu, null_sd, z, p_null):\n', resp_win(1), resp_win(2), base_win(1), base_win(2), n_real, n_sham, sum(isfinite(null_amp(:, 1))));
disp(stats(setdiff(1:nY, row_roi), {'name', 'mm', 'mm_lo', 'mm_hi', 'sham', 'diff', 'p_perm', 'null_mu', 'null_sd', 'z', 'p_null'}));
roi_stats = stats(row_roi, :);
[~, ord] = sort(roi_stats.p_null);  q_fdr = nan(nROI, 1);  q_fdr(ord) = bh_fdr(roi_stats.p_null(ord));  roi_stats.q_fdr = q_fdr;
fprintf('Single ROIs (selected = %s):\n', lst(sel_roi));
disp(roi_stats(:, {'name', 'mm', 'mm_lo', 'mm_hi', 'sham', 'diff', 'p_perm', 'z', 'p_null', 'q_fdr'}));

% windows: during the freeze and after it
for wnm = {'halt', 'post'}
    AM = eval(['AMP' iff(strcmp(wnm{1}, 'halt'), 'h', 'p')]);
    fprintf('\n%s window (%s s): mismatch vs sham', wnm{1}, mat2str(eval([wnm{1} '_win'])));
    fprintf('\n  %-44s %9s %9s %9s %8s\n', 'row', 'mismatch', 'sham', 'diff', 'p_perm');
    for y = [row_grp row_beh]
        pd = zeros(n_perm, 1);
        for p = 1:n_perm; r = randperm(nEv); pd(p) = mean(AM(r(1:n_real), y), 'omitnan') - mean(AM(r(n_real+1:end), y), 'omitnan'); end
        dobs = mean(AM(mm_i, y), 'omitnan') - mean(AM(sh_i, y), 'omitnan');
        fprintf('  %-44s %9.3f %9.3f %9.3f %8.3f\n', row_name{y}(1:min(end, 44)), mean(AM(mm_i, y), 'omitnan'), mean(AM(sh_i, y), 'omitnan'), dobs, (1 + sum(abs(pd) >= abs(dobs))) / (n_perm + 1));
    end
end

% figure 03: dF/F of the groups
fh = newfig('Event-triggered dF/F', 18, 4 + 3.2 * ceil(ngr / 2));
tl = tiledlayout(fh, ceil(ngr / 2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s | flow freeze (red, n=%d), sham (grey, n=%d), simulated-scheduler pseudo-events (blue band = 95%% of the null mean)', strrep(session, '_', '\_'), n_real, n_sham), 'FontSize', 8, 'FontWeight', 'normal');
for g = 1:ngr
    a = nexttile(tl);  hold(a, 'on');  y = row_grp(g);
    patch(a, [0 1 1 0], [-1 -1 1 1] * 10, [1 0.9 0.9], 'EdgeColor', 'none');
    if have_null; qq = prctile(null_tr(:, :, y), [2.5 97.5], 1); fill(a, [tg fliplr(tg)], [qq(1, :) fliplr(qq(2, :))], col.null, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
    h1 = sem_trace(a, tg, squeeze(A(:, sh_i, y)), col.sham, true);
    h2 = sem_trace(a, tg, squeeze(A(:, mm_i, y)), col.mm, true);
    mtr = mean(A(:, mm_i, y), 2, 'omitnan');
    ylim(a, 1.3 * [min(min(mtr), -0.15) max(max(mtr), 0.15)]);  xlim(a, [-ev_pre_s ev_post_s]);  yline(a, 0, ':');
    title(a, sprintf('%s | resp %.3f [%.3f, %.3f], null z %.1f', glab_all{g}(1:min(end, 34)), st(y).mm, st(y).mm_lo, st(y).mm_hi, st(y).z), 'FontSize', 6.5, 'FontWeight', 'normal');
    if g > ngr - 2; xlabel(a, 'time from freeze onset (s)'); end
    if mod(g, 2) == 1; ylabel(a, 'dF/F'); end
    if g == 1; legend(a, [h2 h1], {'mismatch', 'sham'}, 'FontSize', 6, 'Location', 'northeast'); end
end
if save_figs; savef(fh, '03_eta_dff'); end

% figure 04: behaviour around the freeze
if ~isempty(row_beh)
    fh = newfig('Event-triggered behaviour', 18, 4 + 3.2 * ceil(numel(row_beh) / 2));
    tl = tiledlayout(fh, ceil(numel(row_beh) / 2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, [strrep(session, '_', '\_') ' | behaviour around the freeze: mismatch (red), sham (grey), simulated-scheduler pseudo-events (blue)'], 'FontSize', 8, 'FontWeight', 'normal');
    for b = 1:numel(row_beh)
        a = nexttile(tl);  hold(a, 'on');  y = row_beh(b);
        patch(a, [0 1 1 0], [-1 -1 1 1] * 1e5, [1 0.9 0.9], 'EdgeColor', 'none');
        if have_null; qq = prctile(null_tr(:, :, y), [2.5 97.5], 1); fill(a, [tg fliplr(tg)], [qq(1, :) fliplr(qq(2, :))], col.null, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
        sem_trace(a, tg, squeeze(A(:, sh_i, y)), col.sham, true);
        sem_trace(a, tg, squeeze(A(:, mm_i, y)), col.mm, true);
        lo = min(mean(A(:, mm_i, y), 2)); hi = max(mean(A(:, mm_i, y), 2));  pad = 0.4 * max(hi - lo, eps);
        if all(isfinite([lo hi])); ylim(a, [lo - pad hi + pad]); end
        xlim(a, [-ev_pre_s ev_post_s]);  yline(a, 0, ':');
        title(a, sprintf('%s | resp %.3g [%.3g, %.3g], null z %.1f', beh_lab.(row_name{y}), st(y).mm, st(y).mm_lo, st(y).mm_hi, st(y).z), 'FontSize', 6.5, 'FontWeight', 'normal');
        if b > numel(row_beh) - 2; xlabel(a, 'time from freeze onset (s)'); end
    end
    if save_figs; savef(fh, '04_eta_behaviour'); end
end

% figure 05: single ROIs
fh = newfig('Event-triggered dF/F per ROI', 18, 9);
tl = tiledlayout(fh, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
a = nexttile(tl);
M_roi = squeeze(mean(A(:, mm_i, row_roi), 2, 'omitnan'))';                         % [nROI x nTg]
sd_roi = squeeze(std(A(bi, :, row_roi), 0, [1 2], 'omitnan'));                  % single-event baseline SD per ROI
imagesc(a, tg, 1:nROI, M_roi ./ (max(sd_roi(:), eps) / sqrt(n_real)));  axis(a, 'xy');  clim(a, [-4 4]);  colormap(a, div_map());  cb = colorbar(a);  ylabel(cb, 'mean / (baseline SD / sqrt(n))');
xline(a, [0 1], 'k:');  xlabel(a, 'time from freeze onset (s)');  ylabel(a, 'ROI');  yticks(a, 1:nROI);
title(a, 'mismatch mean per ROI in standard errors of the baseline noise', 'FontSize', 6.5, 'FontWeight', 'normal');
a = nexttile(tl);  hold(a, 'on');
errorbar(a, 1:nROI, roi_stats.mm, roi_stats.mm - roi_stats.mm_lo, roi_stats.mm_hi - roi_stats.mm, 'o', 'Color', col.mm, 'MarkerFaceColor', col.mm, 'MarkerSize', 3);
if have_null; errorbar(a, (1:nROI) + 0.25, roi_stats.null_mu, 1.96 * roi_stats.null_sd, 's', 'Color', col.null, 'MarkerSize', 3); end
plot(a, 1:nROI, roi_stats.sham, 'x', 'Color', col.sham);
plot(a, find(signal_like), zeros(1, nnz(signal_like)) - 0.25, 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
yline(a, 0, ':');  xlabel(a, 'ROI');  ylabel(a, 'dF/F, response 0-2 s');  xlim(a, [0 nROI + 1]);
title(a, {'mismatch (red, 95% CI), simulated null (blue, 95%),', 'sham (x); triangles = signal-like ROIs'}, 'FontSize', 6.5, 'FontWeight', 'normal');
if save_figs; savef(fh, '05_eta_per_roi'); end

%% 9. Per-event amplitudes: behaviour regression, dependence on speed and on the order in the session
y_main = row_grp(1);                                                                % mean of the selected ROIs
bnames = intersect({'LocSpeed', 'PupilArea', 'WhiskIndex', 'FacePadME'}, row_name(row_beh), 'stable');
Xb = nan(nEv, numel(bnames));
for b = 1:numel(bnames); xb = AMP(:, strcmp(row_name, bnames{b})); Xb(:, b) = (xb - mean(xb, 'omitnan')) / std(xb, 'omitnan'); end
reg_names = [{'mismatch (vs sham)'}, strcat({'d '}, bnames)];
reg_out = table();
ex = cell(0, 5);  exadd = @(c, a, b, v, ci, n) [c; {a, b, v, ci, n}];
ci_str = @(x) sprintf('[%.3f, %.3f]', x(1), x(2));
if ~isempty(bnames)
    fprintf('\nPer-event regression of the %s response (0-%g s) on mismatch (vs sham) and the behaviour change in the same window:\n', glab_all{1}, resp_win(2));
    for y = row_grp
        [bb, ci_, pp] = ols_boot(AMP(:, y), [ones(nEv, 1) double(ev_M) Xb], n_boot);
        rn = row_name{y}(1:min(end, 30));
        reg_out = [reg_out; table(repmat({rn}, numel(reg_names), 1), reg_names(:), bb(2:end), ci_(2:end, 1), ci_(2:end, 2), pp(2:end), 'VariableNames', {'row', 'regressor', 'coef', 'lo', 'hi', 'p'})]; %#ok<AGROW>
    end
    disp(reg_out);
end
% response vs speed at onset, and vs order in the session
r_spd = corr(ev_speed_on(ev_M), AMP(ev_M, y_main), 'Type', 'Spearman', 'Rows', 'complete');
r_ord = corr(find(ev_M), AMP(ev_M, y_main), 'Type', 'Spearman', 'Rows', 'complete');
half = ev_t(ev_M) < median(ev_t(ev_M));
fprintf('\n%s: Spearman vs speed at onset %.2f, vs event number %.2f; first half %.3f, second half %.3f (mismatch events); trial 1 %.3f (n=%d), trial 2 %.3f (n=%d)\n', glab_all{1}, r_spd, r_ord, ...
    mean(AMP(mm_i(half), y_main)), mean(AMP(mm_i(~half), y_main)), mean(AMP(mm_i(ev_trial(mm_i) == 1), y_main)), sum(ev_trial(mm_i) == 1), mean(AMP(mm_i(ev_trial(mm_i) == 2), y_main)), sum(ev_trial(mm_i) == 2));

% selection check: events that fired at once vs events that had to wait for running (those fire 1 s into a run onset, with the speed still rising)
w0 = find(ev_M & ev_wait == 0);  w1 = find(ev_M & ev_wait > 0.3);
yl_ = [y_main, row_beh(ismember(row_name(row_beh), {'LocSpeed', 'WhiskIndex'}))];
fprintf('Mismatch events that fired at once (wait 0 s, n=%d) vs after waiting > 0.3 s for running (n=%d); response 0-%g s:\n', numel(w0), numel(w1), resp_win(2));
for y = yl_
    fprintf('  %-28s %+9.3g vs %+9.3g\n', row_name{y}(1:min(end, 28)), mean(AMP(w0, y), 'omitnan'), mean(AMP(w1, y), 'omitnan'));
    ex = exadd(ex, '9 selection', [row_name{y} ': waited vs at once'], mean(AMP(w1, y), 'omitnan') - mean(AMP(w0, y), 'omitnan'), '', numel(w1)); %#ok<SAGROW>
end

fh = newfig('Per-event amplitudes', 18, 8);
tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
a = nexttile(tl);  hold(a, 'on');
swarm(a, 1, AMP(mm_i, y_main), col.mm);  swarm(a, 2, AMP(sh_i, y_main), col.sham);
if have_null; swarm(a, 3, null_amp(:, y_main), col.null); end
plot(a, [0.6 1.4], mean(AMP(mm_i, y_main)) * [1 1], 'k', 'LineWidth', 1.5);  plot(a, [1.6 2.4], mean(AMP(sh_i, y_main)) * [1 1], 'k', 'LineWidth', 1.5);
xticks(a, 1:3);  xticklabels(a, {sprintf('mismatch (%d)', n_real), sprintf('sham (%d)', n_sham), 'null: simulated scheduler (means)'});  xtickangle(a, 20);  ylabel(a, 'dF/F response 0-2 s');  yline(a, 0, ':');
title(a, glab_all{1}, 'FontSize', 7, 'FontWeight', 'normal');
a = nexttile(tl);  hold(a, 'on');
plot(a, ev_speed_on(ev_M), AMP(ev_M, y_main), 'o', 'Color', col.mm, 'MarkerSize', 3);  plot(a, ev_speed_on(ev_sham), AMP(ev_sham, y_main), 'x', 'Color', col.sham);
xlabel(a, 'speed at onset (cm/s)');  ylabel(a, 'dF/F response');  title(a, sprintf('\\rho = %.2f', r_spd), 'FontSize', 7, 'FontWeight', 'normal');
a = nexttile(tl);  hold(a, 'on');
plot(a, find(ev_M), AMP(ev_M, y_main), 'o', 'Color', col.mm, 'MarkerSize', 3);  plot(a, find(ev_sham), AMP(ev_sham, y_main), 'x', 'Color', col.sham);
xlabel(a, 'event number');  title(a, sprintf('\\rho = %.2f', r_ord), 'FontSize', 7, 'FontWeight', 'normal');
if save_figs; savef(fh, '06_per_event'); end

%% 9b. Arousal coupling over the whole stimulus period (positive control: can this ROI selection follow anything?)
% Correlation of the smoothed group traces with the smoothed behaviour traces at lags of -3...+3 s (positive lag = behaviour AFTER the axon trace); the null is the maximum
% |r| over the lags after circularly shifting the behaviour trace by >= 30 s.
idx_s = find(tCa >= tsf(1) & tCa <= tsf(nfr));  n_s = numel(idx_s);
lag_n = round((-3:0.5:3) * fs);
Gz = grp_all(:, idx_s);  Gz = movmean(Gz, w_sm, 2);  Gz = (Gz - mean(Gz, 2)) ./ std(Gz, 0, 2);
coup = table();
for b = 1:numel(beh_names)
    nm = beh_names{b};
    yb_ = beh_ca.(nm)(idx_s);  yb_ = fillmissing(movmean(yb_, w_sm, 'omitnan'), 'linear', 'EndValues', 'nearest');  yz = (yb_ - mean(yb_)) / std(yb_);
    r_lag = zeros(ngr, numel(lag_n));
    for li = 1:numel(lag_n); r_lag(:, li) = Gz * circshift(yz, -lag_n(li))' / (n_s - 1); end
    nmax = zeros(ngr, 100);
    for p = 1:100
        ys = circshift(yz, randi([round(30 * fs), n_s - round(30 * fs)]));
        rr_ = zeros(ngr, numel(lag_n));
        for li = 1:numel(lag_n); rr_(:, li) = Gz * circshift(ys, -lag_n(li))' / (n_s - 1); end
        nmax(:, p) = max(abs(rr_), [], 2);
    end
    [~, ib_] = max(abs(r_lag), [], 2);
    coup = [coup; table(repmat({nm}, ngr, 1), glab_all(:), r_lag(:, lag_n == 0), r_lag(sub2ind(size(r_lag), (1:ngr)', ib_)), lag_n(ib_)' / fs, prctile(nmax, 99, 2), ...
        'VariableNames', {'behaviour', 'group', 'r_lag0', 'r_best', 'lag_s', 'null_p99_absr'})]; %#ok<AGROW>
end
coup.group = cellfun(@(s) s(1:min(end, 30)), coup.group, 'UniformOutput', false);
fprintf('\nCoupling of the group traces to behaviour over the whole stimulus period (smoothed %g s):\n', smooth_s);
disp(coup);
sig_c = abs(coup.r_best) > coup.null_p99_absr;
chk = addchk(chk, 'arousal coupling (positive control)', grade(any(sig_c(strcmp(coup.group, coup.group{1}))), 1), sprintf('the %s follows behaviour beyond the shuffle (|r| at the best lag > null p99) for: %s. If none, this selection is blind to arousal and a missing mismatch response says little', glab_all{1}, strjoin(coup.behaviour(strcmp(coup.group, coup.group{1}) & sig_c)', ', ')));
printchk(chk, n0);  n0 = size(chk, 1);

%% 10. Reference: natural flow stops and run onsets (the flow stops / starts with the mouse, no prediction error)
nat = struct('name', {'natural stop', 'run onset'}, 't', {[], []});
if loc_ok
    fsp = spd_g;  st_ok = false(1, nT);  ro_ok = false(1, nT);
    below = fsp < stop_lo_mms;
    idx = find(diff([false below]) == 1);
    for i = idx
        j2 = min(nT, i + round(stop_hold_s * fs));  jb = max(1, i - round(2 * fs)):max(1, i - round(0.5 * fs));
        if all(below(i:j2)) && mean(fsp(jb), 'omitnan') > stop_hi_mms && i > round(2 * fs); st_ok(i) = true; end
    end
    above = fsp > stop_hi_mms;
    idx = find(diff([false above]) == 1);
    for i = idx
        jb = max(1, i - round(1.5 * fs)):max(1, i - 1);  j2 = min(nT, i + round(1 * fs));
        if i > round(1.5 * fs) && all(below(jb)) && mean(fsp(i:j2), 'omitnan') > stop_hi_mms; ro_ok(i) = true; end
    end
    for q = 1:2
        tq = tCa(find(iff(q == 1, st_ok, ro_ok)));
        okq = tq > tsf(1) + ev_pre_s & tq < tsf(nfr) - ev_post_s & min(abs(tq(:) - [ev_t; ev_off_t]'), [], 2)' > excl_s;
        nat(q).t = tq(okq);
    end
    fprintf('\nReference events away from the freezes (> %g s): %d natural stops (speed %g -> < %g mm/s, held %g s), %d run onsets\n', excl_s, numel(nat(1).t), stop_hi_mms, stop_lo_mms, stop_hold_s, numel(nat(2).t));
    chk = addchk(chk, 'reference: natural stops / run onsets', grade(numel(nat(1).t) >= 10, 1), sprintf('%d natural stops and %d run onsets away from the logged events; the mouse runs %.0f%% of the time, so stops are rare. Thresholds: %g -> %g mm/s for the stop', numel(nat(1).t), numel(nat(2).t), 100 * mean(spd_g(tCa >= tsf(1) & tCa <= tsf(nfr)) > 10 * run_gate_cms, 'omitnan'), stop_hi_mms, stop_lo_mms));
    An = cell(1, 2);
    for q = 1:2
        if numel(nat(q).t) < 3; continue; end
        An{q} = align_rows(Z, tCa, nat(q).t, tg);  An{q} = An{q} - mean(An{q}(bi, :, :), 1, 'omitnan');
        ampq = amp_of(An{q}, ri_);
        c2 = boot_ci_mean(ampq(:, row_grp(1)), n_boot);
        fprintf('  %-13s (n=%d): %s response 0-%g s %+.3f [%+.3f, %+.3f]; mismatch %+.3f [%+.3f, %+.3f]\n', nat(q).name, numel(nat(q).t), glab_all{1}, resp_win(2), mean(ampq(:, row_grp(1)), 'omitnan'), c2(1), c2(2), st(row_grp(1)).mm, st(row_grp(1)).mm_lo, st(row_grp(1)).mm_hi);
        ex = exadd(ex, '10 reference', [nat(q).name ': ' glab_all{1}], mean(ampq(:, row_grp(1)), 'omitnan'), ci_str(c2), numel(nat(q).t)); %#ok<SAGROW>
    end
    fh = newfig('Natural stops and run onsets', 18, 4 + 3.2 * ceil((ngr + 1) / 2));
    tl = tiledlayout(fh, ceil((ngr + 1) / 2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, [strrep(session, '_', '\_') ' | flow freeze while running (red) vs natural stop (orange) vs run onset (green)'], 'FontSize', 8, 'FontWeight', 'normal');
    for g = [1:ngr, 0]
        a = nexttile(tl);  hold(a, 'on');
        if g > 0; y = row_grp(g); ttl = glab_all{g}(1:min(end, 34)); else; y = find(strcmp(row_name, 'LocSpeed')); ttl = 'running speed (mm/s)'; end
        patch(a, [0 1 1 0], [-1 -1 1 1] * 1e5, [1 0.9 0.9], 'EdgeColor', 'none');
        sem_trace(a, tg, squeeze(A(:, mm_i, y)), col.mm, true);
        if ~isempty(An{1}); sem_trace(a, tg, squeeze(An{1}(:, :, y)), [0.9 0.5 0.05], true); end
        if ~isempty(An{2}); sem_trace(a, tg, squeeze(An{2}(:, :, y)), [0.1 0.6 0.2], true); end
        xlim(a, [-ev_pre_s ev_post_s]);  yline(a, 0, ':');  title(a, ttl, 'FontSize', 6.5, 'FontWeight', 'normal');
    end
    if save_figs; savef(fh, '07_reference_events'); end
end

%% 11. Stimulus start and trial change
% t0 = first presented frame (gratings appear after wait_init s of grey), later trials: the position of the gratings is reset (new_trial), end: last frame
ep_t = [tsf(1); b_t0(2:end); tsf(nfr) + ifi];  ep_name = [{'stimulus start'}, arrayfun(@(k) sprintf('trial %d start', k), 2:nTrP, 'UniformOutput', false), {'stimulus end'}];
ep_ok = ep_t - 2 >= tCa(1) & ep_t + 5 <= tCa(end);
tg_e = (-2 * fs:5 * fs) / fs;  bi_e = tg_e >= -2 & tg_e < -0.5;  ri_e = tg_e >= 0 & tg_e <= 3;
m_all = grp_all(1, :);
rand_amp = nan(2000, 1);
for p = 1:numel(rand_amp)
    i0 = randi([round(2 * fs) + 1, nT - round(5 * fs) - 1]);  rand_amp(p) = mean(m_all(i0 + (0:round(3 * fs)))) - mean(m_all(i0 - round(2 * fs):i0 - round(0.5 * fs)));
end
fprintf('\nEpoch responses of %s (mean 0-3 s minus baseline -2...-0.5 s; null = %d random times in the session: SD %.3f):\n', glab_all{1}, numel(rand_amp), std(rand_amp));
fh = newfig('Epochs', 18, 5);  tl = tiledlayout(fh, 1, numel(ep_t), 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(ep_t)
    a = nexttile(tl);  hold(a, 'on');
    if ~ep_ok(k); title(a, [ep_name{k} ': outside imaging'], 'FontSize', 6.5, 'FontWeight', 'normal'); continue; end
    xk = interp1(tCa, m_all, ep_t(k) + tg_e);  xk = xk - mean(xk(bi_e));
    plot(a, tg_e, xk, 'Color', col.dff);  yyaxis(a, 'right');
    if loc_ok; plot(a, tg_e, interp1(tCa, spd_g, ep_t(k) + tg_e), 'Color', col.loc); end
    yyaxis(a, 'left');  xline(a, 0, ':');
    amp_k = mean(xk(ri_e));  fprintf('  %-34s %+.3f dF/F, z = %.1f\n', ep_name{k}, amp_k, amp_k / std(rand_amp));
    title(a, sprintf('%s: %+.3f (z %.1f)', ep_name{k}, amp_k, amp_k / std(rand_amp)), 'FontSize', 6.5, 'FontWeight', 'normal');
    xlabel(a, 'time (s)');
end
if save_figs; savef(fh, '08_epochs'); end

%% 12. Summary and findings table
yb = row_grp(1);  has_ctl = ~isempty(ctl_roi);  if has_ctl; yc = row_grp(end); end
fprintf('\n===== %s | alleyCL block %d | %d mismatch, %d sham events =====\n', session, kb, n_real, n_sham);
fprintf('Primary test (%s, dF/F 0-%g s): mismatch %+.3f [%+.3f, %+.3f]; sham %+.3f; mismatch - sham %+.3f (p_perm = %.3f)', glab_all{1}, resp_win(2), st(yb).mm, st(yb).mm_lo, st(yb).mm_hi, st(yb).sham, st(yb).diff, st(yb).p_perm);
if have_null; fprintf('; simulated-scheduler pseudo-events %+.3f +- %.3f (z = %.1f, p = %.3f)', st(yb).null_mu, st(yb).null_sd, st(yb).z, st(yb).p_null); end
if has_ctl; fprintf('\nControl (noise-like ROIs): %+.3f [%+.3f, %+.3f]\n', st(yc).mm, st(yc).mm_lo, st(yc).mm_hi); else; fprintf('\n'); end
for b = 1:numel(row_beh)
    y = row_beh(b);
    fprintf('Behaviour %-14s 0-%g s: mismatch %+.3g [%+.3g, %+.3g], sham %+.3g, matched pseudo %+.3g +- %.2g (z = %.1f)\n', row_name{y}, resp_win(2), st(y).mm, st(y).mm_lo, st(y).mm_hi, st(y).sham, st(y).null_mu, st(y).null_sd, st(y).z);
end
ex = exadd(ex, '8 mismatch', [glab_all{1} ': response 0-2 s'], st(yb).mm, ci_str([st(yb).mm_lo st(yb).mm_hi]), n_real);
ex = exadd(ex, '8 sham', [glab_all{1} ': response 0-2 s'], st(yb).sham, '', n_sham);
ex = exadd(ex, '8 mismatch - sham', glab_all{1}, st(yb).diff, sprintf('p_perm = %.3f', st(yb).p_perm), n_real);
if have_null; ex = exadd(ex, '8 simulated-scheduler pseudo-events', glab_all{1}, st(yb).null_mu, sprintf('SD %.3f, z = %.1f, p = %.3f', st(yb).null_sd, st(yb).z, st(yb).p_null), nM); end
for g = 2:ngr
    y = row_grp(g);  ex = exadd(ex, '8 groups', glab_all{g}, st(y).mm, sprintf('%s, sham %.3f, p_perm %.3f, null z %.1f', ci_str([st(y).mm_lo st(y).mm_hi]), st(y).sham, st(y).p_perm, st(y).z), n_real);
end
for b = 1:numel(row_beh)
    y = row_beh(b);  ex = exadd(ex, '8 behaviour', row_name{y}, st(y).mm, sprintf('%s, sham %.3g, null z %.1f', ci_str([st(y).mm_lo st(y).mm_hi]), st(y).sham, st(y).z), n_real);
end
ex = exadd(ex, '9 order', 'Spearman amplitude vs speed at onset', r_spd, '', n_real);
ex = exadd(ex, '9 order', 'Spearman amplitude vs event number', r_ord, '', n_real);
findings = cell2table(ex, 'VariableNames', {'analysis', 'statistic', 'value', 'ci_or_p', 'n'});
if save_figs
    writetable(findings, fullfile(fig_dir, 'findings.csv'));  writetable(stats, fullfile(fig_dir, 'stats_all_rows.csv'));  writetable(roi_stats, fullfile(fig_dir, 'stats_per_roi.csv'));
    writetable(roiq, fullfile(fig_dir, 'roi_quality.csv'));
    qc = cell2table(chk, 'VariableNames', {'check', 'status', 'message'});  writetable(qc, fullfile(fig_dir, 'checks.csv'));
end
fprintf('\nChecks: %d PASS, %d WARN, %d FAIL\n', sum(strcmp(chk(:, 2), 'PASS')), sum(strcmp(chk(:, 2), 'WARN')), sum(strcmp(chk(:, 2), 'FAIL')));

%% local functions
function [ev_f, wait_s] = sim_scheduler(nxt, tr_start, tr_end, last_f, ifi, mmp, dur_fr)
% one realisation of the mismatch scheduler (alleyCL_mismatch_step) on the logged online speed: nxt(f) = first frame >= f at which the gate (speed > min_speed
% for min_run_s) is satisfied (Inf if none). Per trial: the first event may be armed first_delay s after the trial start; an armed event fires at the first gate-
% satisfying frame (if before the last allowed start); the next one is armed an interval U(interval) after the end of the event.
ev_f = zeros(0, 1);  wait_s = zeros(0, 1);
for k = 1:numel(tr_start)
    p = tr_start(k) + round(mmp.first_delay / ifi);
    while p <= tr_end(k)
        f = nxt(p);
        if ~isfinite(f) || f > tr_end(k) || f > last_f(k); break; end
        ev_f(end+1, 1) = f;  wait_s(end+1, 1) = (f - p) * ifi; %#ok<AGROW>
        p = f + dur_fr + round((mmp.interval(1) + diff(mmp.interval) * rand) / ifi);
    end
end
end

function m = win_mean(cs, cn, idx, lo, hi)
% mean of every row over the samples idx+lo ... idx+hi from the cumulative sums cs / counts cn (leading zero column); [numel(idx) x nRows], NaN outside the record
a = idx(:) + lo;  b = idx(:) + hi;
ok = a >= 1 & b <= size(cs, 2) - 1;
m = nan(numel(idx), size(cs, 1));
m(ok, :) = ((cs(:, b(ok) + 1) - cs(:, a(ok))) ./ (cn(:, b(ok) + 1) - cn(:, a(ok))))';
end

function v = iff(c, a, b)
% inline conditional
if c; v = a; else; v = b; end
end

function save_fig(fh, base)
% PNG (300 dpi) and vector PDF of a figure
exportgraphics(fh, [base '.png'], 'Resolution', 300);
exportgraphics(fh, [base '.pdf'], 'ContentType', 'vector');
end

function n = nM_(mask)
% number of true elements (struct field initialiser)
n = sum(mask);
end

function idx = repmat_frames(on, off)
% all frame indices of the intervals on(k) ... off(k)-1
idx = [];
for k = 1:numel(on); idx = [idx; (on(k):off(k)-1)']; end %#ok<AGROW>
end

function x = fillmissing_zero(x)
% NaN -> 0 (no offset known: the PC-clock time is kept)
x(isnan(x)) = 0;
end

function A = align_rows(X, t, te, tg)
% X [nY x nT] on times t; returns [numel(tg) x numel(te) x nY], NaN where the window leaves the record
Tq = te(:)' + tg(:);
A = nan(numel(tg), numel(te), size(X, 1));
for y = 1:size(X, 1); A(:, :, y) = interp1(t, X(y, :), Tq, 'linear', NaN); end
end

function ci = boot_ci_mean(x, nb)
% 95% bootstrap interval of the mean of x (NaN dropped)
x = x(isfinite(x));
if numel(x) < 3; ci = [NaN NaN]; return; end
ci = prctile(mean(x(randi(numel(x), numel(x), nb)), 1), [2.5 97.5]);
end

function q = bh_fdr(p)
% Benjamini-Hochberg q values of ascending sorted p
n = numel(p);  q = p(:) .* n ./ (1:n)';
for i = n-1:-1:1; q(i) = min(q(i), q(i+1)); end
q = min(q, 1);
end

function h = sem_trace(ax, t, M, c, band)
% mean over the columns of M [time x events] (line, returned handle), +- SEM as a cloud if band
t = t(:);  mu = mean(M, 2, 'omitnan');  se = std(M, 0, 2, 'omitnan') / sqrt(max(1, sum(any(isfinite(M), 1))));
hold(ax, 'on');
ok = isfinite(mu) & isfinite(se);
if band; fill(ax, [t(ok); flipud(t(ok))], [mu(ok) + se(ok); flipud(mu(ok) - se(ok))], c, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); end
h = plot(ax, t, mu, 'Color', c, 'LineWidth', 1.2);
end

function swarm(ax, x0, y, c)
% simple jittered dot plot
y = y(:);  jit = (rand(numel(y), 1) - 0.5) * 0.5;
plot(ax, x0 + jit, y, '.', 'Color', c, 'MarkerSize', 6);
end

function cm = div_map()
% blue - white - red
cm = [linspace(0.1, 1, 128)' linspace(0.3, 1, 128)' linspace(0.8, 1, 128)'; ones(128, 1) linspace(1, 0.15, 128)' linspace(1, 0.15, 128)'];
end

function [b, ci, p] = ols_boot(y, X, n_boot)
% least squares y = X b with a case-resampling bootstrap: b [p x 1], ci [p x 2] (2.5 / 97.5 %), p [p x 1] two-sided bootstrap p (floor 1/n_boot)
ok = all(isfinite([y(:), X]), 2);  y = y(ok);  X = X(ok, :);  n = numel(y);
ws = warning('off', 'all');
b = X \ y(:);
bs = zeros(n_boot, size(X, 2));
for i = 1:n_boot; r = randi(n, n, 1); bs(i, :) = (X(r, :) \ y(r))'; end
warning(ws);
ci = prctile(bs, [2.5 97.5])';
p = max(2 * min(mean(bs > 0, 1), mean(bs < 0, 1)), 1 / n_boot)';
end
