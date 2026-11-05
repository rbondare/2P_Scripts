%% axon_analysis_spont_V1
% Spontaneous activity of axons (SC) in a session aggregated with newAggregation/save_preprocessed.m: what drives the shared signal of the axons, and
% is it neural at all? Works on the 'spontaneous' block (grey screen, no stimulus) of any such file, also of a session that was stopped early
% (AnimalRB24_261023_1527: the 300 s spontaneous block is the only real block; use axon_analysis_V3 with drop_dead_blocks for its checks).
%
% Sections
%   1 load, epoch, ROI quality (same definitions as axon_analysis_V3: signal-like, motion-coupled -> 'good')
%   2 behaviour on a common 10 Hz grid: pupil, whisking, face motion, locomotion, eye speed (saccades), brain motion (suite2p rigid shift)
%   3 population structure: ROI correlation, axon groups (bundle = one axon), PCA of the signal-like ROIs, positions
%   4 correlation with each behaviour: lag 0, best lag within +-corr_max_lag_s, PARTIAL r (the other behaviours regressed out), circular-shift null
%   5 encoding model: held-out R^2 of the dF/F from lagged behaviours (blocked CV ridge); single, unique and shared variance of each behaviour
%   6 event-triggered averages: locomotion onset, isolated whisking bouts (with / without running), saccades; and behaviour around dF/F events
%   6b states (quiet / whisking without running / running) and the large events of the shared signal: behaviour at each, kinetics
%   7 controls: brain motion (speed, displacement, non-rigid) regressed out, neuropil and raw fluorescence at running onset / PC1 peaks, per-ROI coupling to motion; activity left over after the behaviour
%   8 figures, 9 tables
%
% Statistics: one session = one animal, 300 s. The units are the ROI GROUPS (ROIs with r > grp_r_thr share one axon), behaviours are slow and
% strongly correlated, so every p value comes from circular shifts of the behaviour inside the epoch (shifts >= min_shift_s, ALL admissible shifts),
% which keeps the autocorrelation of both signals. Exploratory.
%
% Locomotion: the ball log of open-loop blocks is truncated (see memory ballcam-truncation-rb22-24), so loc_source = 'auto' uses the ball log only if it
% is complete for the epoch and otherwise a speed PROXY from the body video (locomotion_proxy/<session>_locproxy.mat: t, speed_mm_s; log-speed r ~0.9
% against the real ball speed in same-session reference data, over-estimates slow speeds).
%
% Outputs left in the workspace: spont_tbl (correlation, partial r, best lag per unit x behaviour), enc_tbl (held-out R^2), eta_tbl (event-triggered
% amplitudes and latencies), roiq, sel_roi, signal_like, groups.

%% 0. Settings (edit only this block)
session   = 'AnimalRB24_261023_1527';
agg_dir   = 'Z:\Group Members\Rima\Aggregated';
s2p_dir   = fullfile('Z:\Group Members\Rima\DATA_2P', session, 'suite2p', 'plane0');          % Fall.mat: rigid / non-rigid motion (section 7)
repo_root = 'C:\Users\rbondarenko\projects\2P_Scripts';
proxy_dir = fullfile(repo_root, 'axon analysis', 'locomotion_proxy');
save_figs = false;                                                % true: PNG (300 dpi) + vector PDF of every figure and the tables as CSV in fig_dir
fig_dir   = fullfile(repo_root, 'axon analysis', 'plots', 'spont_V1', session);

epoch_s   = [];           % [] = the 'spontaneous' block(s) of the file; or [t0 t1; ...] windows in session time (s)
roi_sel   = 'good';       % ROIs of the groups / mean: 'good' (signal-like, not motion-coupled, as V3) | 'signal' | 'all' | a list, e.g. [9 10 15 18 22 23 26]
loc_source = 'auto';      % 'auto' | 'ball' | 'proxy' | 'none'
loc_run_thr   = 20;       % mm/s: above this the animal counts as running (events)
loc_min_still_s = 3;      % a running onset needs this long below the threshold before it, and loc_min_run_s above it after it
loc_min_run_s   = 2;
dt_a      = 0.1;          % analysis grid (s)
corr_smooth_s = 0.5;      % smoothing of every signal before the correlations / the model (s)
corr_max_lag_s = 5;       % best lag searched within +- this (s); positive lag = dF/F follows the behaviour
min_shift_s    = 20;      % circular shifts used for the null: at least this long
enc_lags_s  = -1:0.5:3;   % encoding model: behaviour at t - lag predicts dF/F at t (positive = behaviour first)
enc_nfold   = 6;          % contiguous folds of the epoch
enc_margin_s = 5;         % training samples this close to the test fold are left out
eta_win_s   = [-3 6];     % event-triggered window
eta_base_s  = [-2 -0.5];  % baseline of the event-triggered averages
eta_amp_s   = [0 3];      % window of the response amplitude
whisk_thr   = 1.6;        % WhiskIndex (1 = rest level of the DLC face pipeline) above this = whisking
whisk_gap_s = 1.5;        % isolated whisking bout: whisking off this long before
whisk_min_s = 0.5;        % a whisking bout lasts at least this
sacc_thr_degs = 120;      % saccade: eye speed peak above this (deg/s, 2-frame difference) and a displacement >= sacc_min_deg (RB24_261023_1527: all horizontal, median 7 deg, no excess at blinks)
sacc_min_deg  = 4;
slow_ratio_thr = 4;  mot_r_thr = 0.15;  smooth_s = 0.5;  neuropil_factor = 0.5;   % ROI quality, as in V3
grp_r_thr   = 0.8;        % ROIs whose smoothed dF/F correlates above this are one group (average linkage)
n_boot      = 1000;       % bootstrap / null draws
rng(1);

col = struct('dff', [0 0 0], 'pupil', [0.85 0.15 0.15], 'whisk', [0.1 0.55 0.25], 'face', [0.9 0.55 0.05], 'loc', [0.1 0.35 0.8], 'eye', [0.5 0.3 0.7], 'mot', [0.45 0.45 0.45], 'pc1', [0.8 0.1 0.5]);
div_cmap = [linspace(0.1, 1, 128)' linspace(0.3, 1, 128)' linspace(0.8, 1, 128)'; ones(128, 1) linspace(1, 0.15, 128)' linspace(1, 0.15, 128)'];   % blue - white - red
newfig = @(nm, w, h) figure('Name', nm, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 w h], 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 8, ...
    'DefaultAxesBox', 'off', 'DefaultAxesTickDir', 'out', 'DefaultAxesLineWidth', 0.75, 'DefaultLineLineWidth', 0.75, 'DefaultTextFontName', 'Arial', 'DefaultLegendBox', 'off');
lst = @(v) ['[' strtrim(sprintf('%g ', v)) ']'];
zs  = @(x) (x - mean(x, 'omitnan')) ./ std(x, 'omitnan');

%% 1. Load, epoch, ROI quality
addpath(repo_root, fullfile(repo_root, 'analysis_functions'));
if save_figs && ~isfolder(fig_dir); mkdir(fig_dir); end
savef = @(fh, nm) save_fig(fh, fullfile(fig_dir, nm));
data_file = fullfile(agg_dir, [session '_preprocessed.mat']);
assert(isfile(data_file), 'Aggregated file not found: %s', data_file);
fmeta = whos('-file', data_file);  inFile = {fmeta.name};
want  = {'CaData', 'TimeCa', 'Stimuli', 'Triggers', 'Behav', 'LocomotionCal', 'LocomotionNonCal', 'VideoParams', 'AnimalInfo'};
fprintf('Loading %s ...\n', data_file);
data = load(data_file, want{ismember(want, inFile)});
nz = find(~cellfun(@isempty, {data.CaData.Ca_dFF}));  assert(~isempty(nz), 'CaData has no calcium traces');
ca = data.CaData(nz(1));
dFF = double(ca.Ca_dFF);  F = double(ca.Ca_F);  Fneu = double(ca.Ca_F_neuropil);
[nROI, nT] = size(dFF);
tCa = double(data.TimeCa(1, :));  fs = 1 / median(diff(tCa));
cent = ca.Ca_centroid_voxel;  npix = arrayfun(@(c) c.npix, ca.Ca_cellstats(:));
tcam = double(data.Triggers.TimeCamera(:))';  Bh = data.Behav;
px_per_mm = NaN;  if isfield(data, 'VideoParams') && isfield(data.VideoParams, 'Eye') && isfield(data.VideoParams.Eye, 'px_per_mm'); px_per_mm = data.VideoParams.Eye.px_per_mm; end
if isfield(data, 'AnimalInfo'); ai = data.AnimalInfo; fprintf('\n%s | %s | %s\n', ai.AnimalName, ai.Genotype, ai.Indicator.Type); end

% epoch: spontaneous block(s) with real frame times (a block of a stopped session can be all NaN)
S = data.Stimuli;
if isempty(epoch_s)
    for k = find(strcmp({S.type}, 'spontaneous'))
        tf = double(S(k).TimeStimulusFrame(:));
        if nnz(isfinite(tf)) > 100; epoch_s = [epoch_s; min(tf) max(tf)]; end %#ok<AGROW>
    end
    assert(~isempty(epoch_s), 'No spontaneous block with frame times in this file; set epoch_s = [t0 t1] (session time, s)');
end
epoch_s = sortrows(epoch_s);  nW = size(epoch_s, 1);
ep_dur = diff(epoch_s, 1, 2);
fprintf('%d ROIs, %d frames (%.0f s, %.2f Hz) | epoch: %d window(s) %s s, %.0f s in total\n', nROI, nT, tCa(end) - tCa(1), fs, nW, mat2str(round(epoch_s, 1)), sum(ep_dur));
in_ep = false(1, nT);  for w = 1:nW; in_ep = in_ep | (tCa >= epoch_s(w, 1) & tCa <= epoch_s(w, 2)); end
assert(nnz(diff(tCa) > 1.5 / fs & in_ep(2:end)) == 0, 'Imaging gaps inside the epoch');

% ROI quality as V3 section 5 / 6c (whole session), see the header of V3 for the definitions
Fc = F - neuropil_factor * Fneu;
w_sm = max(1, round(smooth_s * fs));  dsm = movmean(dFF, w_sm, 2);
med_s = median(dsm, 2);  dev = dsm - med_s;  neg = dev < 0;
sig_neg = max(sqrt(sum(dev .^ 2 .* neg, 2) ./ max(sum(neg, 2), 1)), eps('single'));
dd = diff(dFF, 1, 2);  hf_noise = 1.4826 * median(abs(dd - median(dd, 2)), 2) / sqrt(2);
slow_ratio = var(dsm, 0, 2) ./ (hf_noise .^ 2 / w_sm);
ev_pos = sum(diff([false(nROI, 1) dev > 3 * sig_neg], 1, 2) == 1, 2);  ev_neg = sum(diff([false(nROI, 1) dev < -3 * sig_neg], 1, 2) == 1, 2);
bad_F = mean(Fc <= 0, 2) > 0.01;
signal_like = ~bad_F & ~(slow_ratio < slow_ratio_thr | ev_pos < 2 * max(ev_neg, 1));
% motion: suite2p rigid shift (xoff, yoff) and non-rigid spread; ROI follows motion if |r| > circular-shift p99 and > mot_r_thr
have_motion = false;  rigid = nan(1, nT);  nonr = nan(1, nT);  mot_flag = false(nROI, 1);  motion_r = nan(nROI, 2);
fall_file = fullfile(s2p_dir, 'Fall.mat');
if isfile(fall_file)
    o = load(fall_file, 'ops');  xo = double(o.ops.xoff(:))';  yo = double(o.ops.yoff(:))';
    if numel(xo) == nT
        have_motion = true;
        rigid = hypot(xo - median(xo), yo - median(yo));
        nonr  = hypot(std(double(o.ops.xoff1), 0, 2), std(double(o.ops.yoff1), 0, 2))';
        mots = {rigid, nonr};
        dz_  = zs(double(dsm)');                                              % [nT x nROI], z-scored per ROI
        for m = 1:2
            mz_ = zs(movmean(mots{m}, w_sm));
            c_all = circ_xcorr(dz_, mz_');                                    % [nT x nROI] circular r at every lag (row 1 = lag 0)
            sh = (round(30 * fs):nT - round(30 * fs)) + 1;                    % circular shifts >= 30 s (all of them, no random draw)
            motion_r(:, m) = c_all(1, :)';
            thr99 = prctile(abs(c_all(sh, :)), 99, 1)';
            mot_flag = mot_flag | abs(motion_r(:, m)) > max(thr99, mot_r_thr);
        end
    else
        fprintf('suite2p ops has %d motion samples for %d frames: motion not used\n', numel(xo), nT);
    end
else
    fprintf('%s not found: motion-coupled ROIs cannot be flagged\n', fall_file);
end
good_roi = find(signal_like & ~mot_flag)';
switch lower(char(string(roi_sel)))
    case 'good';   sel_roi = good_roi;
    case 'signal'; sel_roi = find(signal_like)';
    case 'all';    sel_roi = 1:nROI;
    otherwise;     sel_roi = double(roi_sel(:))';
end
assert(~isempty(sel_roi), 'No ROI selected (roi_sel = %s)', char(string(roi_sel)));
sig_roi = find(signal_like)';
nSel = numel(sel_roi);
mot_txt = 'suite2p motion not available';  if have_motion; mot_txt = 'suite2p rigid / non-rigid'; end
fprintf('signal-like ROIs (%d): %s\nmotion-coupled (%s): %d ROIs, of them signal-like: %s\n''good'' ROIs: %s | used (roi_sel): %s\n', numel(sig_roi), lst(sig_roi), ...
    mot_txt, sum(mot_flag), lst(find(mot_flag & signal_like)'), lst(good_roi), lst(sel_roi));

%% 2. Behaviour on the analysis grid
% every signal is binned to dt_a inside the epoch windows (windows are concatenated; the joins are not special-cased) and smoothed corr_smooth_s
ta_e = [];  win_of = [];
for w = 1:nW
    e_ = epoch_s(w, 1):dt_a:epoch_s(w, 2);  ta_e = [ta_e, e_(1:end-1)]; win_of = [win_of, w * ones(1, numel(e_) - 1)]; %#ok<AGROW>
end
nA  = numel(ta_e);  ta = ta_e + dt_a / 2;                                   % bin centres
ed_w = cell(1, nW);  for w = 1:nW; ed_w{w} = epoch_s(w, 1):dt_a:epoch_s(w, 2); end
gridify = @(t, x) cell2mat(arrayfun(@(w) bin_mean(t, x, ed_w{w}), 1:nW, 'UniformOutput', false));   % row vector on the grid, NaN where nothing falls in a bin
fillx   = @(x) fillmissing(x, 'linear', 'EndValues', 'nearest');
smg     = @(x) movmean(x, max(1, round(corr_smooth_s / dt_a)));

% locomotion
loc_name = 'none';  loc_t = [];  loc_v = [];
use_ball = false;
if any(strcmp(loc_source, {'auto', 'ball'})) && isfield(data, 'LocomotionCal')
    ok_b = false(1, numel(data.LocomotionCal));
    for k = 1:numel(data.LocomotionCal)
        % complete log (one reading per ball trigger) AND its triggers cover the epoch; a block that was never shown (stopped session) can have a few triggers inside the epoch
        L_ = data.LocomotionNonCal(k);  nv = size(L_.Values, 1);  tt = double(L_.Triggers_sec(:));
        cov_ = (min(max(tt), epoch_s(end, 2)) - max(min(tt), epoch_s(1, 1))) / (epoch_s(end, 2) - epoch_s(1, 1));
        ok_b(k) = abs(nv - numel(tt)) <= 2 && nv > 1 && cov_ > 0.9;
    end
    if any(ok_b)
        use_ball = true;
        L = data.LocomotionCal(find(ok_b, 1));  n_ = min([numel(L.Triggers_sec), numel(L.Forward_mm), numel(L.SideRight_mm)]);
        tt = double(L.Triggers_sec(1:n_))';  dtl = median(diff(tt));
        inc = hypot(double(L.Forward_mm(1:n_))' - median(L.Forward_mm(1:n_)), double(L.SideRight_mm(1:n_))' - median(L.SideRight_mm(1:n_)));
        loc_t = tt;  loc_v = movmean(inc, max(1, round(0.5 / dtl))) / dtl;  loc_name = 'ball log';
    end
end
if ~use_ball && any(strcmp(loc_source, {'auto', 'proxy'}))
    pf = fullfile(proxy_dir, [session '_locproxy.mat']);
    if isfile(pf)
        pp = load(pf);  loc_t = double(pp.t(:))';  loc_v = double(pp.speed_mm_s(:))';  loc_name = 'video proxy';
    else
        fprintf('No usable ball log and no proxy file %s: locomotion not available\n', pf);
    end
end
have_loc = ~isempty(loc_t);

% behaviours (row vectors on the grid), z-scored inside the epoch after smoothing
beh_names = {'pupil', 'whisk', 'face', 'loc', 'eye speed', 'brain motion'};
beh_lab   = {'pupil area', 'whisking (log index)', 'face motion (log ME)', 'locomotion (log speed)', 'eye speed (log)', 'brain motion (log)'};
beh_col   = [col.pupil; col.whisk; col.face; col.loc; col.eye; col.mot];
Braw = nan(6, nA);
pa = double(Bh.PupilArea(:))';  if isfield(Bh, 'PupilBad'); pa(logical(Bh.PupilBad(:))') = NaN; end
if isfinite(px_per_mm); pa = pa / px_per_mm ^ 2; end
Braw(1, :) = gridify(tcam, pa);
Braw(2, :) = gridify(tcam, log(max(double(Bh.WhiskIndex(:))', 0.1)));
Braw(3, :) = gridify(tcam, log(max(double(Bh.FacePadME(:))', 0.01)));
if have_loc;  Braw(4, :) = interp1(loc_t, log1p(max(loc_v, 0)), ta, 'linear', NaN); end          % the proxy is sampled every 0.5 s: interpolated, not binned
% eye position (deg): speed over 2 camera frames from the median-filtered az / el; saccades are detected on it
cam_dt = median(diff(tcam));
az = movmedian(double(Bh.PupilAzEl(:, 1))', 3, 'omitnan');  el = movmedian(double(Bh.PupilAzEl(:, 2))', 3, 'omitnan');
bad_eye = false(size(az));  if isfield(Bh, 'PupilBad'); bad_eye = logical(Bh.PupilBad(:))'; end
az(bad_eye) = NaN;  el(bad_eye) = NaN;
eye_spd = hypot(circshift(az, -1) - circshift(az, 1), circshift(el, -1) - circshift(el, 1)) / (2 * cam_dt);  eye_spd([1 end]) = NaN;
Braw(5, :) = gridify(tcam, log1p(eye_spd));
if have_motion; Braw(6, :) = gridify(tCa, log1p(hypot([0 diff(xo)], [0 diff(yo)]) * 30)); end      % px per s
beh_ok = arrayfun(@(j) mean(isfinite(Braw(j, :))) > 0.7 && std(Braw(j, :), 'omitnan') > 0, 1:6);
fprintf('\nbehaviours on the %.1f s grid (%d bins): %s | locomotion: %s\n', dt_a, nA, strjoin(beh_names(beh_ok), ', '), loc_name);
if any(~beh_ok); fprintf('  not available: %s\n', strjoin(beh_names(~beh_ok), ', ')); end
B = nan(6, nA);
for j = find(beh_ok); B(j, :) = zs(smg(fillx(Braw(j, :)))); end
bidx = find(beh_ok);  nB = numel(bidx);  bn = beh_names(bidx);                  % behaviours that exist
% calcium on the same grid
dffA = nan(nROI, nA);  for r = 1:nROI; dffA(r, :) = gridify(tCa, dFF(r, :)); end
dffS = movmean(dffA, max(1, round(corr_smooth_s / dt_a)), 2, 'omitnan');
FcA = nan(nROI, nA);  FneuA = nan(nROI, nA);  for r = 1:nROI; FcA(r, :) = gridify(tCa, Fc(r, :)); FneuA(r, :) = gridify(tCa, Fneu(r, :)); end
rigA = gridify(tCa, rigid);
tE = (0:nA - 1) * dt_a;                                                      % epoch time (s) of the grid bins
Bc = B(bidx, :);                                                             % behaviours that exist [nB x nA], z-scored
Bfill = nan(6, nA);  for j = bidx; Bfill(j, :) = zs(fillx(Braw(j, :))); end  % the same without smoothing (event-triggered averages)
maxlag = round(corr_max_lag_s / dt_a);  minsh = round(min_shift_s / dt_a);  blk10 = round(10 / dt_a);

%% 3. Population structure: ROI groups, PCA, positions
Xs = zs(dffS(sig_roi, :)');                                                  % [nA x nsig] z-scored smoothed dF/F of the signal-like ROIs
[Uu, Ss, Vv] = svd(Xs, 'econ');  expl = diag(Ss) .^ 2 / sum(diag(Ss) .^ 2);
pc1 = (Uu(:, 1) * Ss(1, 1))';  ld1 = Vv(:, 1)';
if mean(ld1) < 0; pc1 = -pc1; ld1 = -ld1; end
part_ratio = sum(expl) ^ 2 / sum(expl .^ 2);
Csg = corrcoef(dffS(sig_roi, :)');
ord_sig = optimalleaforder(linkage(pdist(dffS(sig_roi, :), 'correlation'), 'average'), pdist(dffS(sig_roi, :), 'correlation'));
[gid, grp_rois] = roi_groups(dffS, sel_roi, grp_r_thr);  ng = numel(grp_rois);
groups = grp_rois;
off_ = Csg(triu(true(numel(sig_roi)), 1));
dist_px = squareform(pdist(cent(sig_roi, 1:2)));  dd_ = dist_px(triu(true(numel(sig_roi)), 1));
fprintf('\n--- 3. population structure (%d signal-like ROIs) ---\n', numel(sig_roi));
fprintf('PC1 explains %.0f%% of the variance (PC2 %.0f%%, PC3 %.0f%%); participation ratio %.1f of %d; %d of %d ROIs load on PC1 with the same sign\n', ...
    100 * expl(1), 100 * expl(2), 100 * expl(3), part_ratio, numel(sig_roi), sum(ld1 > 0), numel(sig_roi));
fprintf('pairwise r (0.5 s smoothed): median %.2f, %.0f%% of pairs > 0.5, %.0f%% > 0.8; far pairs (> 150 px) median r %.2f\n', median(off_), 100 * mean(off_ > 0.5), 100 * mean(off_ > 0.8), median(off_(dd_ > 150)));
fprintf('selected ROIs form %d group(s) at r > %g: %s\n', ng, grp_r_thr, strjoin(cellfun(@(g) lst(g), grp_rois, 'UniformOutput', false), '  '));
% units of the statistics
u_name = [{'PC1 (signal-like ROIs)'; 'mean of selected ROIs'}; arrayfun(@(g) sprintf('group %d: ROI %s', g, strjoin(string(grp_rois{g}), ',')), (1:ng)', 'UniformOutput', false)];
Uz = [zs(pc1); zs(mean(dffS(sel_roi, :), 1))];
for g = 1:ng; Uz(2 + g, :) = zs(mean(dffS(grp_rois{g}, :), 1)); end
nU = size(Uz, 1);
dff_sel = mean(dffA(sel_roi, :), 1);                                         % dF/F units, unsmoothed
r_pc1 = corr(pc1', dff_sel');
fprintf('PC1 of the signal-like ROIs vs the mean of the selected ROIs: r = %.2f\n', r_pc1);

%% 4. Correlation with behaviour: lag 0, best lag, partial r, circular-shift null
% positive lag = dF/F follows the behaviour. Partial: the OTHER behaviours (all in Bc, brain motion included) are regressed out of both signals at lag 0.
% p0 / p_best: rank among ALL admissible circular shifts of the behaviour (|r| at lag 0, or max |r| over the lag range), p_best is the one to read for 'best lag'.
rows = {};  ccf_store = cell(nU, nB);
for u = 1:nU
    for jj = 1:nB
        x = Uz(u, :);  y = Bc(jj, :);
        res = ccf_circ(x, y, maxlag, minsh);
        oth = setdiff(1:nB, jj);
        if isempty(oth); xr = x; yr = y; else; xr = resid_on(x', Bc(oth, :)')'; yr = resid_on(y', Bc(oth, :)')'; end
        resp_ = ccf_circ(xr, yr, maxlag, minsh);
        ci = boot_r(x, y, blk10, n_boot);
        ccf_store{u, jj} = res;
        rows(end + 1, :) = {u_name{u}, bn{jj}, res.r0, ci(1), ci(2), res.p0, res.r_best, res.lag_best * dt_a, res.p_best, resp_.r0, resp_.p0}; %#ok<SAGROW>
    end
end
spont_tbl = cell2table(rows, 'VariableNames', {'unit', 'behaviour', 'r_lag0', 'ci_lo', 'ci_hi', 'p_lag0', 'r_best', 'lag_s', 'p_best', 'r_partial', 'p_partial'});
fprintf('\n--- 4. correlation with behaviour (z-scored, %.1f s smoothing; partial = other behaviours regressed out) ---\n', corr_smooth_s);
for u = 1:min(nU, 2)
    fprintf('%s\n', u_name{u});
    disp(spont_tbl(strcmp(spont_tbl.unit, u_name{u}), 2:end));
end
% behaviour-behaviour
Cbb = corrcoef(Bc');
fprintf('behaviour-behaviour r:\n');  disp(array2table(Cbb, 'VariableNames', strrep(bn, ' ', '_'), 'RowNames', strrep(bn, ' ', '_')));

%% 5. Encoding model: held-out R^2 of the dF/F from lagged behaviours
% ridge (lambda chosen by an inner blocked CV), enc_nfold contiguous folds of the epoch with enc_margin_s purged around the test fold. Per unit:
% full model, each behaviour left out (unique contribution = R2 full - R2 without it), each behaviour alone (single). CIs: bootstrap over 10 s blocks of the held-out predictions.
lagb = round(enc_lags_s / dt_a);  nl = numel(lagb);
Xl = lag_matrix(Bc, lagb);  colb = repelem(1:nB, nl);
mg = round(enc_margin_s / dt_a);  lams = 10 .^ (0:0.5:5);
mod_names = [{'full'}, strcat('without_', strrep(bn, ' ', '_')), strcat('only_', strrep(bn, ' ', '_'))];
Yfull = nan(nU, nA);  enc_rows = {};  enc_R2 = nan(nU, numel(mod_names));  enc_ci = nan(nU, numel(mod_names), 2);  enc_dci = nan(nU, nB, 2);
for u = 1:nU
    y = Uz(u, :)';
    Yh = nan(nA, 1 + 2 * nB);
    Yh(:, 1) = cv_ridge(y, Xl, enc_nfold, mg, lams);  Yfull(u, :) = Yh(:, 1)';
    for jj = 1:nB; Yh(:, 1 + jj) = cv_ridge(y, Xl(:, colb ~= jj), enc_nfold, mg, lams); end
    for jj = 1:nB; Yh(:, 1 + nB + jj) = cv_ridge(y, Xl(:, colb == jj), enc_nfold, mg, lams); end
    rb = r2_boot(y, Yh, blk10, n_boot);
    enc_R2(u, :) = rb.r2;  enc_ci(u, :, :) = reshape(rb.ci', 1, [], 2);
    for jj = 1:nB; enc_dci(u, jj, :) = prctile(rb.boot(:, 1) - rb.boot(:, 1 + jj), [2.5 97.5]); end
    for m = 1:numel(mod_names); enc_rows(end + 1, :) = {u_name{u}, mod_names{m}, rb.r2(m), rb.ci(1, m), rb.ci(2, m)}; end %#ok<SAGROW>
end
enc_tbl = cell2table(enc_rows, 'VariableNames', {'unit', 'model', 'R2_heldout', 'ci_lo', 'ci_hi'});
fprintf('\n--- 5. encoding model (held-out R^2, %d folds, lags %g to %g s) ---\n', enc_nfold, enc_lags_s(1), enc_lags_s(end));
for u = 1:min(nU, 2)
    fprintf('%s: full R^2 %.2f [%.2f, %.2f]\n', u_name{u}, enc_R2(u, 1), enc_ci(u, 1, 1), enc_ci(u, 1, 2));
    for jj = 1:nB
        fprintf('   %-13s single %5.2f [%5.2f, %5.2f] | unique (full - without) %5.2f [%5.2f, %5.2f]\n', bn{jj}, enc_R2(u, 1 + nB + jj), enc_ci(u, 1 + nB + jj, 1), enc_ci(u, 1 + nB + jj, 2), ...
            enc_R2(u, 1) - enc_R2(u, 1 + jj), enc_dci(u, jj, 1), enc_dci(u, jj, 2));
    end
end

%% 6. Event-triggered averages
offs = round(eta_win_s(1) / dt_a):round(eta_win_s(2) / dt_a);  tt = offs * dt_a;
isb = tt >= eta_base_s(1) & tt <= eta_base_s(2);  isa = tt >= eta_amp_s(1) & tt <= eta_amp_s(2);
to_idx = @(t) arrayfun(@(x) ternary_idx(ta, x, dt_a), t);                     % absolute time -> grid index (NaN if not inside a window)
% the traces cut around events: dF/F of the mean of the selected ROIs (dF/F units), PC1 and each group (z), behaviours (z, unsmoothed)
sg_names = [{'mean of selected ROIs'}, {'PC1'}, cellfun(@(g) ['group ' char(strjoin(string(g), ','))], grp_rois, 'UniformOutput', false)];
Sg = [dff_sel; zs(pc1)];
for g = 1:ng; Sg(2 + g, :) = mean(dffA(grp_rois{g}, :), 1); end
Eb = Bfill(bidx, :);                                                         % behaviours [nB x nA]
% locomotion state
run = false(1, nA);  run_on = [];
if have_loc
    spg = movmedian(fillx(interp1(loc_t, loc_v, ta, 'linear', NaN)), max(1, round(1 / dt_a)));
    run = spg > loc_run_thr;
    nst = round(loc_min_still_s / dt_a);  nrn = round(loc_min_run_s / dt_a);
    cand = find(diff([run(1) run]) == 1);
    run_on = cand(arrayfun(@(i) i > nst && i + nrn - 1 <= nA && ~any(run(i - nst:i - 1)) && all(run(i:i + nrn - 1)) && win_of(i - nst) == win_of(i + nrn - 1), cand));
end
% whisking bouts (camera rate)
won = double(Bh.WhiskIndex(:))' > whisk_thr;  if isfield(Bh, 'WhiskArtifact'); won(logical(Bh.WhiskArtifact(:))') = false; end
won = movmedian(double(won), 5) > 0.5;
d_ = diff([0 won 0]);  w_st = find(d_ == 1);  w_en = find(d_ == -1) - 1;
kb = (w_en - w_st + 1) * cam_dt >= whisk_min_s;  w_st = w_st(kb);
gpf = round(whisk_gap_s / cam_dt);
iso = arrayfun(@(i) i > gpf && ~any(won(i - gpf:i - 1)), w_st);
w_idx = to_idx(tcam(w_st(iso)));  w_idx = w_idx(isfinite(w_idx));
near_run = false(size(w_idx));
if have_loc; near_run = arrayfun(@(i) any(run(max(1, i - round(2 / dt_a)):min(nA, i + round(2 / dt_a)))), w_idx); end
% saccades (camera rate)
cand = find(eye_spd > sacc_thr_degs);
cand = cand(cand > 4 & cand < numel(eye_spd) - 3);
cand = cand(arrayfun(@(c) eye_spd(c) == max(eye_spd(c - 3:c + 3)), cand));
dsp_ = hypot(az(cand + 3) - az(cand - 3), el(cand + 3) - el(cand - 3));
cand = cand(dsp_ >= sacc_min_deg);
cand = cand([true, diff(tcam(cand)) > 0.2]);
sac_t = tcam(cand);  s_idx = to_idx(sac_t);  keep_s = isfinite(s_idx);  s_idx = s_idx(keep_s);  cand = cand(keep_s);
sac_daz = az(cand + 3) - az(cand - 3);                                       % horizontal displacement (deg); the sign convention of the DLC azimuth is not tied to nasal / temporal here
fprintf('saccades: %d in the epoch, median amplitude %.1f deg, %.0f%% horizontal-dominated\n', numel(cand), median(hypot(sac_daz, el(cand + 3) - el(cand - 3))), 100 * mean(abs(sac_daz) > abs(el(cand + 3) - el(cand - 3))));
% events: struct array name / index list
ev = struct('name', {}, 'idx', {});
if have_loc; ev(end + 1) = struct('name', 'running onset', 'idx', run_on); end
ev(end + 1) = struct('name', 'whisk bout, still', 'idx', w_idx(~near_run));
if have_loc; ev(end + 1) = struct('name', 'whisk bout, running', 'idx', w_idx(near_run)); end
ev(end + 1) = struct('name', 'saccade', 'idx', s_idx);
ev(end + 1) = struct('name', 'saccade, az +', 'idx', s_idx(sac_daz > 0));
ev(end + 1) = struct('name', 'saccade, az -', 'idx', s_idx(sac_daz < 0));
% valid window centres for the null
vc = find(arrayfun(@(i) i + offs(1) >= 1 && i + offs(end) <= nA && win_of(i + offs(1)) == win_of(i + offs(end)), 1:nA));
fprintf('\n--- 6. event-triggered averages (%g to %g s, baseline %g to %g s, response = mean %g to %g s; null = %d circular shifts of the event set) ---\n', eta_win_s, eta_base_s, eta_amp_s, n_boot);
fprintf('events in the epoch: %s\n', strjoin(arrayfun(@(e) sprintf('%s %d', e.name, numel(e.idx)), ev, 'UniformOutput', false), ', '));
if ~isempty(w_idx) || ~isempty(cand); fprintf('(eye speed > %g deg/s peaks with >= %g deg displacement: %d saccades = %.1f per min; whisking bouts: %d isolated of %d)\n', sacc_thr_degs, sacc_min_deg, numel(cand), 60 * numel(cand) / sum(ep_dur), numel(w_idx), numel(w_st)); end
eta = struct('name', {}, 'n', {}, 'S', {}, 'B', {}, 'amp', {}, 'p', {}, 'lat', {}, 'lat_ci', {}, 'null_lo', {}, 'null_hi', {}, 'idx', {});
for e = 1:numel(ev)
    ie = ev(e).idx(:)';  ie = ie(ismember(ie, vc));
    if numel(ie) < 3; fprintf('  %-26s only %d usable events, skipped\n', ev(e).name, numel(ie)); continue; end
    [AS, AB] = cut_events(Sg, Eb, ie, offs);
    AS = AS - mean(AS(:, isb, :), 2);  AB = AB - mean(AB(:, isb, :), 2);
    amp = squeeze(mean(AS(:, isa, :), 2));                                    % [nrow x nev]
    nul = zeros(n_boot, size(Sg, 1));  nul_eta = zeros(n_boot, numel(tt));
    for b = 1:n_boot
        ir = [];
        for tries = 1:20                                                      % the event set shifted as a whole (keeps its clustering and the states it sits in)
            ir = mod(ie + randi([minsh, nA - minsh]) - 1, nA) + 1;  ir = ir(ismember(ir, vc));
            if numel(ir) >= max(3, round(0.6 * numel(ie))); break; end
            ir = [];
        end
        if isempty(ir); ir = vc(randi(numel(vc), 1, numel(ie))); end
        As_ = cut_events(Sg(1:2, :), zeros(0, nA), ir, offs);  As_ = As_ - mean(As_(:, isb, :), 2);
        nul(b, 1:2) = mean(squeeze(mean(As_(:, isa, :), 2)), 2)';
        nul_eta(b, :) = mean(squeeze(As_(1, :, :)), 2)';
    end
    [pp_, lat_, lci_] = deal(nan(1, 2));
    for r_ = 1:2
        pp_(r_) = (1 + sum(abs(nul(:, r_)) >= abs(mean(amp(r_, :))))) / (n_boot + 1);
        m_ = mean(squeeze(AS(r_, :, :)), 2)';  lat_(r_) = half_rise(tt, m_, eta_amp_s);
        lb = nan(1, 500);  for b = 1:500; ib = randi(numel(ie), 1, numel(ie)); lb(b) = half_rise(tt, mean(squeeze(AS(r_, :, ib)), 2)', eta_amp_s); end
        lci_(r_) = prctile(lb, 2.5);  lci_(r_ + 2) = prctile(lb, 97.5);
    end
    k_ = numel(eta) + 1;
    eta(k_) = struct('name', ev(e).name, 'n', numel(ie), 'S', AS, 'B', AB, 'amp', mean(amp, 2)', 'p', pp_, 'lat', lat_, 'lat_ci', lci_, ...
        'null_lo', prctile(nul_eta, 2.5, 1), 'null_hi', prctile(nul_eta, 97.5, 1), 'idx', ie);
    fprintf('  %-26s n = %3d | mean of selected ROIs: %+.3f dF/F (p %.3f, %.0f%% of events positive), half-rise %.2f s [%.2f, %.2f] | PC1: %+.2f SD (p %.3f), half-rise %.2f s\n', ev(e).name, numel(ie), ...
        mean(amp(1, :)), pp_(1), 100 * mean(amp(1, :) > 0), lat_(1), lci_(1), lci_(3), mean(amp(2, :)), pp_(2), lat_(2));
end
eta_tbl = table({eta.name}', [eta.n]', arrayfun(@(e) e.amp(1), eta)', arrayfun(@(e) e.p(1), eta)', arrayfun(@(e) e.lat(1), eta)', arrayfun(@(e) e.amp(2), eta)', arrayfun(@(e) e.p(2), eta)', ...
    'VariableNames', {'event', 'n', 'amp_sel_dFF', 'p_sel', 'half_rise_s', 'amp_pc1_SD', 'p_pc1'});
% behaviour around the large events of the shared signal: what precedes / accompanies a PC1 peak?
pk = islocalmax(movmean(Uz(1, :), round(1 / dt_a)), 'MinProminence', 1, 'MinSeparation', round(3 / dt_a)) & movmean(Uz(1, :), round(1 / dt_a)) > 1;
pk_idx = find(pk);  pwin = round(-8 / dt_a):round(4 / dt_a);  pt = pwin * dt_a;
vcp = find(arrayfun(@(i) i + pwin(1) >= 1 && i + pwin(end) <= nA && win_of(i + pwin(1)) == win_of(i + pwin(end)), 1:nA));
pk_idx = pk_idx(ismember(pk_idx, vcp));
pkB = [];  pk_amp = nan(1, nB);  pk_p = nan(1, nB);
if numel(pk_idx) >= 3
    pkB = cut_events(Eb, zeros(0, nA), pk_idx, pwin);  pkB = pkB - mean(pkB(:, pt >= -8 & pt <= -6, :), 2);
    ia = pt >= -2 & pt <= 0;  pk_amp = mean(squeeze(mean(pkB(:, ia, :), 2)), 2)';
    nl_ = zeros(n_boot, nB);
    for b = 1:n_boot; A_ = cut_events(Eb, zeros(0, nA), vcp(randi(numel(vcp), 1, numel(pk_idx))), pwin); A_ = A_ - mean(A_(:, pt >= -8 & pt <= -6, :), 2); nl_(b, :) = mean(squeeze(mean(A_(:, ia, :), 2)), 2)'; end
    pk_p = arrayfun(@(j) (1 + sum(abs(nl_(:, j)) >= abs(pk_amp(j)))) / (n_boot + 1), 1:nB);
    fprintf('behaviour in the 2 s before %d large PC1 events (z-units vs 8-6 s before): %s\n', numel(pk_idx), strjoin(arrayfun(@(j) sprintf('%s %+.2f (p %.3f)', bn{j}, pk_amp(j), pk_p(j)), 1:nB, 'UniformOutput', false), ', '));
end

%% 6b. States and the large events of the shared signal
fprintf('\n--- 6b. states, and what accompanies the large PC1 events ---\n');
won_g = gridify(tcam, double(won)) > 0.5;
st = ones(1, nA);  st(won_g) = 2;  if have_loc; st(run) = 3; end
st_lab = {'quiet', 'whisking, not running', 'running'};  nst3 = 2 + have_loc;
fprintf('time in states: %s\n', strjoin(arrayfun(@(k) sprintf('%s %.0f s (%.0f%%)', st_lab{k}, nnz(st == k) * dt_a, 100 * mean(st == k)), 1:nst3, 'UniformOutput', false), ', '));
bid = ceil((1:nA) / blk10);  nbk = max(bid);
st_rows = {};
for r_ = 1:size(Sg, 1)
    for k = 1:nst3
        m_ = st == k;  bs_ = zeros(n_boot, 1);
        for b = 1:n_boot; wgt = accumarray(randi(nbk, nbk, 1), 1, [nbk 1]);  w_ = wgt(bid)' .* m_;  bs_(b) = sum(w_ .* Sg(r_, :)) / sum(w_); end
        st_rows(end + 1, :) = {sg_names{r_}, st_lab{k}, nnz(m_) * dt_a, mean(Sg(r_, m_)), prctile(bs_, 2.5), prctile(bs_, 97.5)}; %#ok<SAGROW>
    end
end
state_tbl = cell2table(st_rows, 'VariableNames', {'unit', 'state', 'seconds', 'mean', 'ci_lo', 'ci_hi'});
fprintf('mean dF/F of the selected ROIs by state (95%% block-bootstrap interval): %s\n', strjoin(arrayfun(@(k) sprintf('%s %.2f [%.2f, %.2f]', st_lab{k}, state_tbl.mean(k), state_tbl.ci_lo(k), state_tbl.ci_hi(k)), 1:nst3, 'UniformOutput', false), ' | '));
fprintf('PC1 (z) by state: %s\n', strjoin(arrayfun(@(k) sprintf('%s %+.2f', st_lab{k}, state_tbl.mean(nst3 + k)), 1:nst3, 'UniformOutput', false), ' | '));
% the large events of PC1: what was the animal doing, and what did the held-out behaviour model predict?
pk_rows = {};
for i = 1:numel(pk_idx)
    ii = pk_idx(i);  w_ = max(1, ii - round(2 / dt_a)):min(nA, ii + round(2 / dt_a));
    lc_ = NaN;  if have_loc; lc_ = mean(spg(w_)); end
    pk_rows(end + 1, :) = {tE(ii), Uz(1, ii), Yfull(1, ii), exp(mean(Braw(2, w_), 'omitnan')), lc_, nnz(abs(tE(s_idx) - tE(ii)) <= 3), mean(exp(Braw(6, w_)) - 1, 'omitnan')}; %#ok<SAGROW>
end
pk_tbl = cell2table(pk_rows, 'VariableNames', {'t_epoch_s', 'PC1_z', 'predicted_z', 'whisk_idx', 'speed_mm_s', 'n_saccades_pm3s', 'shift_px_s'});
disp(pk_tbl);
if ~isempty(pk_idx)
    no_beh = pk_tbl.whisk_idx < whisk_thr & (~have_loc | pk_tbl.speed_mm_s < loc_run_thr);
    fprintf('of %d PC1 peaks, %d occur with no whisking and no running (+-2 s); the held-out behaviour model reaches %.0f%% of the peak height on average (median %.0f%%)\n', ...
        numel(pk_idx), nnz(no_beh), 100 * mean(pk_tbl.predicted_z ./ pk_tbl.PC1_z), 100 * median(pk_tbl.predicted_z ./ pk_tbl.PC1_z));
end
% kinetics of the peaks in the mean dF/F of the selected ROIs (0.2 s smoothing): calcium transients rise fast and decay slowly, a motion / gain change follows the movement
dsm2 = movmean(dff_sel, 2);
k_rise = nan(1, numel(pk_idx));  k_dec = nan(1, numel(pk_idx));  k_amp = nan(1, numel(pk_idx));
for i = 1:numel(pk_idx)
    ii = pk_idx(i);  lo = max(1, ii - round(1.5 / dt_a));  hi = min(nA, ii + round(1.5 / dt_a));
    [~, im] = max(dsm2(lo:hi));  im = lo + im - 1;
    b0 = max(1, im - round(6 / dt_a));  b1 = max(1, im - round(3 / dt_a));
    if b1 <= b0; continue; end
    base = median(dsm2(b0:b1));  A = dsm2(im) - base;  k_amp(i) = A;
    if A < 0.3; continue; end
    i10 = find(dsm2(1:im) <= base + 0.1 * A, 1, 'last');
    if isempty(i10); continue; end
    i90 = find(dsm2(i10:im) >= base + 0.9 * A, 1) + i10 - 1;
    if ~isempty(i90); k_rise(i) = (i90 - i10) * dt_a; end
    ih = find(dsm2(im:min(nA, im + round(8 / dt_a))) <= base + 0.5 * A, 1) + im - 1;
    if ~isempty(ih); k_dec(i) = (ih - im) * dt_a; end
end
fprintf('kinetics of %d PC1 peaks in the mean of the selected ROIs (amplitude median %.2f dF/F): 10-90%% rise %.2f s (median), time to half amplitude after the peak %.2f s (median; %d of %d did not decay to half within 8 s)\n', ...
    numel(pk_idx), median(k_amp, 'omitnan'), median(k_rise, 'omitnan'), median(k_dec, 'omitnan'), nnz(isnan(k_dec) & k_amp >= 0.3), nnz(k_amp >= 0.3));

%% 7. Controls: brain motion, neuropil, per-ROI coupling, what the behaviours do not explain
fprintf('\n--- 7. controls ---\n');
jm = find(strcmp(bn, 'brain motion'));
Mreg = [];   % three motion regressors: rigid shift speed, rigid displacement from the median position, non-rigid spread (0.5 s smoothed, z-scored)
if ~isempty(jm); Mreg = [Bc(jm, :); zs(smg(fillx(rigA))); zs(smg(fillx(gridify(tCa, nonr))))]'; end
% (a) fluorescence around running onsets: neuropil-subtracted ROI trace, neuropil, rigid shift, as % change of their median
Fc_rel = (mean(FcA(sel_roi, :), 1) ./ median(mean(FcA(sel_roi, :), 1), 'omitnan') - 1) * 100;
Fn_rel = (mean(FneuA(sel_roi, :), 1) ./ median(mean(FneuA(sel_roi, :), 1), 'omitnan') - 1) * 100;
ctl_sig = [Fc_rel; Fn_rel; rigA];  ctl_lab = {'F - 0.5 F_{neu} of selected ROIs (%)', 'neuropil F (%)', 'rigid shift (px)'};
ctl_eta = [];
if ~isempty(eta) && have_loc && any(strcmp({eta.name}, 'running onset'))
    ie = eta(strcmp({eta.name}, 'running onset')).idx;
    ctl_eta = cut_events(ctl_sig, zeros(0, nA), ie, offs);  ctl_eta = ctl_eta - mean(ctl_eta(:, isb, :), 2);
    fprintf('at running onset (n = %d): selected-ROI F %+.1f%%, neuropil F %+.1f%%, rigid shift %+.2f px (mean over %g-%g s after onset)\n', numel(ie), ...
        mean(mean(squeeze(ctl_eta(1, isa, :)), 1)), mean(mean(squeeze(ctl_eta(2, isa, :)), 1)), mean(mean(squeeze(ctl_eta(3, isa, :)), 1)), eta_amp_s);
end
% (a2) the same around the large PC1 events
ctl_pk = [];
if numel(pk_idx) >= 3
    ctl_pk = cut_events(ctl_sig, zeros(0, nA), pk_idx, pwin);  ctl_pk = ctl_pk - mean(ctl_pk(:, pt >= -8 & pt <= -6, :), 2);
    ia_ = pt >= -1 & pt <= 1;
    fprintf('at the %d PC1 peaks (+-1 s vs 8-6 s before): selected-ROI F %+.1f%%, neuropil F %+.1f%%, rigid shift %+.2f px\n', numel(pk_idx), mean(mean(squeeze(ctl_pk(1, ia_, :)), 1)), mean(mean(squeeze(ctl_pk(2, ia_, :)), 1)), mean(mean(squeeze(ctl_pk(3, ia_, :)), 1)));
end
% (b) per-ROI coupling: r(dF/F, behaviour) and r(dF/F, brain motion) for every ROI, and the partial r with the three motion regressors removed from both
roi_r = nan(nROI, nB);  roi_pr = nan(nROI, nB);
for r = 1:nROI
    for jj = 1:nB
        roi_r(r, jj) = corr(dffS(r, :)', Bc(jj, :)');
        if ~isempty(Mreg) && jj ~= jm; roi_pr(r, jj) = corr(resid_on(dffS(r, :)', Mreg), resid_on(Bc(jj, :)', Mreg)); end
    end
end
roi_tbl = array2table([(1:nROI)', signal_like, mot_flag, ismember(1:nROI, sel_roi)', roi_r, roi_pr], 'VariableNames', [{'roi', 'signal_like', 'motion_coupled', 'selected'}, strrep(bn, ' ', '_'), strcat('partial_', strrep(bn, ' ', '_'))]);
jl = find(strcmp(bn, 'loc'));
if ~isempty(jm) && ~isempty(jl)
    rr_ = corr(roi_r(signal_like, jl), roi_r(signal_like, jm), 'Type', 'Spearman');
    fprintf('across the %d signal-like ROIs, r(dF/F, locomotion) and r(dF/F, brain motion) correlate at rho = %.2f\n', nnz(signal_like), rr_);
    for q = {'loc', 'face', 'whisk'}
        jq = find(strcmp(bn, q{1}));  if isempty(jq); continue; end
        fm_ = signal_like & mot_flag;  fn_ = signal_like & ~mot_flag;
        fprintf('  r(dF/F, %s) raw -> with motion (3 regressors) regressed out: motion-coupled signal-like ROIs (%d) median %.2f -> %.2f (%d keep r > 0.2); other signal-like ROIs (%d) median %.2f -> %.2f\n', q{1}, nnz(fm_), ...
            median(roi_r(fm_, jq)), median(roi_pr(fm_, jq)), nnz(roi_pr(fm_, jq) > 0.2), nnz(fn_), median(roi_r(fn_, jq)), median(roi_pr(fn_, jq)));
    end
end
% (c) partial correlation with ONLY the brain motion regressors removed (the other behaviours untouched)
ctl_rows = {};
if ~isempty(jm)
    for u = 1:2
        for jj = setdiff(1:nB, jm)
            x = Uz(u, :);  y = Bc(jj, :);
            xr = resid_on(x', Mreg)';  yr = resid_on(y', Mreg)';
            r_raw = ccf_circ(x, y, maxlag, minsh);  r_par = ccf_circ(xr, yr, maxlag, minsh);
            ctl_rows(end + 1, :) = {u_name{u}, bn{jj}, r_raw.r0, r_par.r0, r_par.p0}; %#ok<SAGROW>
        end
    end
    ctl_tbl = cell2table(ctl_rows, 'VariableNames', {'unit', 'behaviour', 'r_lag0', 'r_partial_motion', 'p_partial_motion'});
    fprintf('r with brain motion (speed, displacement, non-rigid) regressed out (lag 0):\n');  disp(ctl_tbl);
end
% (d) activity the behaviours do not explain: out-of-fold residual of every signal-like ROI (smoothed dF/F, all lagged behaviours incl. motion)
res_sig = nan(numel(sig_roi), nA);  r2_roi = nan(1, numel(sig_roi));
for i = 1:numel(sig_roi)
    y = zs(dffS(sig_roi(i), :))';  yh = cv_ridge(y, Xl, enc_nfold, mg, lams);
    res_sig(i, :) = (y - yh)';  r2_roi(i) = 1 - sum((y - yh) .^ 2) / sum((y - mean(y)) .^ 2);
end
Cres = corrcoef(res_sig');
off_res = Cres(triu(true(numel(sig_roi)), 1));
sel_in_sig = find(ismember(sig_roi, sel_roi));
within = []; between = [];
if ng > 1 || numel(grp_rois{1}) > 1
    gof = zeros(1, numel(sig_roi));  for g = 1:ng; gof(ismember(sig_roi, grp_rois{g})) = g; end
    for a = 1:numel(sel_in_sig)
        for b = a + 1:numel(sel_in_sig)
            if gof(sel_in_sig(a)) == gof(sel_in_sig(b)); within(end + 1) = Cres(sel_in_sig(a), sel_in_sig(b)); %#ok<SAGROW>
            else; between(end + 1) = Cres(sel_in_sig(a), sel_in_sig(b)); %#ok<SAGROW>
            end
        end
    end
end
fprintf('held-out R^2 per signal-like ROI (behaviours incl. motion): median %.2f (range %.2f to %.2f)\n', median(r2_roi), min(r2_roi), max(r2_roi));
fprintf('pairwise r among signal-like ROIs: raw median %.2f -> residual median %.2f (%.0f%% of pairs still > 0.5); selected ROIs, same group %.2f / different groups %.2f (residual)\n', ...
    median(off_), median(off_res), 100 * mean(off_res > 0.5), mean(within, 'omitnan'), mean(between, 'omitnan'));
Xr = zs(res_sig');  [~, Sr, ~] = svd(Xr, 'econ');  expl_res = diag(Sr) .^ 2 / sum(diag(Sr) .^ 2);
fprintf('PC1 of the residuals explains %.0f%% (PC1 of the data: %.0f%%)\n', 100 * expl_res(1), 100 * expl(1));

%% 7b. The axial (z) motion concern
% Only one channel is saved and there is no z-stack, so axial displacement is not measured; suite2p sees x / y only. What the data can still say, with the x / y shift as
% the stand-in: (A) does motion carry information about dF/F beyond the behaviours? (B) behavioural state x displacement, matched; (C) the x / y sensitivity of each ROI
% measured in quiet frames, and what it predicts for a running onset; (D) the large PC1 bursts that occur with weak behaviour: do they come with a shift, and does their size follow it?
% (E) the time course after running stops; (F) the sign of the coupling across ALL suite2p ROIs, accepted and rejected (z motion moves structures on either side of the focal plane in
% opposite directions, a neuropil over-subtraction does the same for dim ROIs, calcium in synchronous axons does not). None of this excludes axial motion that is locked to arousal.
if ~isempty(jm)
    fprintf('\n--- 7b. axial-motion concern (x / y shift as stand-in; z is not observable) ---\n');
    rig_s  = zs(smg(fillx(rigA)));                                              % rigid displacement from the median position, z-scored
    rig_px = movmean(fillx(rigA), max(1, round(corr_smooth_s / dt_a)));         % the same in px
    nr_s   = zs(smg(fillx(gridify(tCa, nonr))));
    beh_nm = setdiff(1:nB, jm);
    act    = st >= 2;
    fprintf('(A0) rigid displacement (px, 0.5 s) median [p90]: quiet %.1f [%.1f], whisking %.1f [%.1f]', median(rig_px(st == 1)), prctile(rig_px(st == 1), 90), median(rig_px(st == 2)), prctile(rig_px(st == 2), 90));
    if have_loc; fprintf(', running %.1f [%.1f]', median(rig_px(st == 3)), prctile(rig_px(st == 3), 90)); end
    fprintf('\n');
    % (A) motion vs dF/F with the behaviours regressed out of both
    for u = 1:2
        for q = {'rigid displacement', 'rigid speed', 'non-rigid'}
            switch q{1}; case 'rigid displacement'; m = rig_s; case 'rigid speed'; m = Bc(jm, :); otherwise; m = nr_s; end
            xr = resid_on(Uz(u, :)', Bc(beh_nm, :)')';  yr = resid_on(m', Bc(beh_nm, :)')';
            r0 = ccf_circ(Uz(u, :), m, maxlag, minsh);  rp = ccf_circ(xr, yr, maxlag, minsh);
            fprintf('(A) %-24s vs %-18s r0 %+.2f -> behaviours regressed out %+.2f (p %.3f; best %+.2f at %+.1f s, p %.3f)\n', u_name{u}, q{1}, r0.r0, rp.r0, rp.p0, rp.r_best, rp.lag_best * dt_a, rp.p_best);
        end
    end
    % (B) matched design
    st_nm = {'quiet', 'active'};
    for k = 1:2
        m_ = (st == 1) * (k == 1) | act * (k == 2);  m_ = logical(m_);
        th = median(rig_px(m_));  lo_ = m_ & rig_px <= th;  hi_ = m_ & rig_px > th;
        fprintf('(B) %-6s low shift (%.1f px median, n = %d): %+.2f dF/F of the selected ROIs (PC1 %+.2f) | high shift (%.1f px, n = %d): %+.2f (PC1 %+.2f)\n', st_nm{k}, median(rig_px(lo_)), nnz(lo_), mean(dff_sel(lo_)), mean(Uz(1, lo_)), median(rig_px(hi_)), nnz(hi_), mean(dff_sel(hi_)), mean(Uz(1, hi_)));
    end
    lo_act = act & rig_px <= prctile(rig_px(act), 25);  hi_q = (st == 1) & rig_px >= prctile(rig_px(st == 1), 75);
    fprintf('(B) active frames with the lowest-quartile shift (%.1f px): %+.2f dF/F | quiet frames with the highest-quartile shift (%.1f px): %+.2f dF/F\n', median(rig_px(lo_act)), mean(dff_sel(lo_act)), median(rig_px(hi_q)), mean(dff_sel(hi_q)));
    % (C) x / y sensitivity in quiet frames (motion varies there without the behaviour) and the prediction for a running onset
    qm = st == 1;  Xq = [ones(nnz(qm), 1), rig_px(qm)'];
    b_sel = zeros(1, nSel);  for i = 1:nSel; bb = Xq \ dffS(sel_roi(i), qm)'; b_sel(i) = bb(2); end
    ie_on = find(strcmp({eta.name}, 'running onset'), 1);
    if ~isempty(ie_on) && ~isempty(ctl_eta)
        d_on = mean(mean(squeeze(ctl_eta(3, isa, :)), 1));
        fprintf('(C) quiet-frame slope (dF/F per px of displacement) per selected ROI: %s | the %.1f px shift at running onset predicts %+.2f dF/F (mean slope) vs observed %+.2f; per ROI predicted: %s\n', ...
            strjoin(arrayfun(@(i) sprintf('%d: %+.2f', sel_roi(i), b_sel(i)), 1:nSel, 'UniformOutput', false), ', '), d_on, mean(b_sel) * d_on, eta(ie_on).amp(1), mat2str(round(b_sel * d_on, 2)));
    end
    % (D) the large PC1 bursts: size vs the change in displacement
    dd_amp = nan(numel(pk_idx), 1);  dd_shift = nan(numel(pk_idx), 1);
    for i = 1:numel(pk_idx)
        ii = pk_idx(i);  w1 = ii + (round(-1 / dt_a):round(1 / dt_a));  w0 = ii + (round(-8 / dt_a):round(-6 / dt_a));
        dd_amp(i) = mean(dff_sel(w1)) - mean(dff_sel(w0));  dd_shift(i) = mean(rig_px(w1)) - mean(rig_px(w0));
    end
    if numel(pk_idx) >= 4
        weak = pk_tbl.whisk_idx < 3 & (~have_loc | pk_tbl.speed_mm_s < 50);
        fprintf('(D) %d PC1 bursts: dF/F rise vs change in displacement across bursts, Spearman rho %+.2f\n', numel(pk_idx), corr(dd_amp, dd_shift, 'Type', 'Spearman'));
        fprintf('(D) %d bursts with weak behaviour (whisk index < 3%s): t = %s s, dF/F rise %s, shift change %s px | the other bursts: rise median %.2f, shift change median %.2f px\n', nnz(weak), ternary_str(have_loc, ' and speed < 50 mm/s', ''), ...
            mat2str(round(pk_tbl.t_epoch_s(weak)'), 3), mat2str(dd_amp(weak)', 2), mat2str(dd_shift(weak)', 2), median(dd_amp(~weak)), median(dd_shift(~weak)));
    end
    % (E) after running stops: indicator decay vs an immediate return of the shift
    if have_loc
        run_off = find(diff([run(1) run]) == -1);  nr2 = round(loc_min_run_s / dt_a);  ns3 = round(loc_min_still_s / dt_a);
        run_off = run_off(arrayfun(@(i) i > nr2 && i + ns3 - 1 <= nA && all(run(i - nr2:i - 1)) && ~any(run(i:i + ns3 - 1)), run_off));  run_off = run_off(ismember(run_off, vc));
        if numel(run_off) >= 3
            AO = cut_events([dff_sel; rigA], zeros(0, nA), run_off, offs);  AO = AO - mean(AO(:, isb, :), 2);
            m1 = mean(squeeze(AO(1, :, :)), 2)';  m2 = mean(squeeze(AO(2, :, :)), 2)';
            fprintf('(E) %d running offsets: dF/F change 0-1 s %+.2f, 2-3 s %+.2f, 4-6 s %+.2f | rigid displacement change %+.2f, %+.2f, %+.2f px\n', numel(run_off), ...
                mean(m1(tt >= 0 & tt <= 1)), mean(m1(tt >= 2 & tt <= 3)), mean(m1(tt >= 4 & tt <= 6)), mean(m2(tt >= 0 & tt <= 1)), mean(m2(tt >= 2 & tt <= 3)), mean(m2(tt >= 4 & tt <= 6)));
        end
    end
    % (F) sign of r(relative F, face motion) across every suite2p ROI (needs Fall.mat, which also holds the rejected ROIs)
    jf = find(strcmp(bn, 'face'));
    if isfile(fall_file) && ~isempty(jf)
        o2 = load(fall_file, 'F', 'Fneu', 'iscell');
        Fca = double(o2.F) - neuropil_factor * double(o2.Fneu);  acc = o2.iscell(:, 1) > 0;  nAll = size(Fca, 1);
        if size(Fca, 2) == nT
            rel_ = (Fca - movmedian(Fca, 2 * round(15 * fs) + 1, 2)) ./ max(movmedian(Fca, 2 * round(15 * fs) + 1, 2), eps);
            bright = median(Fca, 2);  ra = nan(nAll, 1);  p_ra = nan(nAll, 1);
            for r = 1:nAll
                x = smg(fillx(gridify(tCa, rel_(r, :))));
                if bright(r) > 1 && std(x) > 0; q = ccf_circ(x, Bc(jf, :), maxlag, minsh); ra(r) = q.r0; p_ra(r) = q.p0; end
            end
            for grp = {'accepted', 'rejected'}
                m = isfinite(ra) & (acc == strcmp(grp{1}, 'accepted'));
                fprintf('(F) %-8s n = %3d: r(F_rel, face) > +0.2 in %d (p < 0.01: %d), < -0.2 in %d (p < 0.01: %d), median %+.2f\n', grp{1}, nnz(m), nnz(ra(m) > 0.2), nnz(ra(m) > 0.2 & p_ra(m) < 0.01), nnz(ra(m) < -0.2), nnz(ra(m) < -0.2 & p_ra(m) < 0.01), median(ra(m)));
            end
            q4 = prctile(bright(isfinite(ra)), [25 50 75]);  e4 = [-inf q4 inf];
            fprintf('(F) by brightness quartile (median Fc %s): fraction with r < -0.2: %s\n', mat2str(round(prctile(bright(isfinite(ra)), [12.5 37.5 62.5 87.5]))), ...
                mat2str(round(arrayfun(@(i) mean(ra(isfinite(ra) & bright >= e4(i) & bright < e4(i + 1)) < -0.2), 1:4), 2)));
            fprintf('(F) accepted ROIs with a negative r(dF/F, face) in the aggregated file: %s\n', lst(find(roi_r(:, strcmp(bn, 'face')) < 0)'));
        end
    end
end

%% 8. Figures
sel_lab = string(sig_roi(ord_sig));  star_ = ismember(sig_roi(ord_sig), sel_roi);  sel_lab(star_) = sel_lab(star_) + "*";
unit_col = [col.pc1; col.dff; lines(ng)];

% ---- 1 overview
fh = newfig('Overview', 18, 19);
tl = tiledlayout(fh, 10, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s | spontaneous epoch %.0f s | %d signal-like ROIs (* = selected, %d)', strrep(session, '_', '\_'), sum(ep_dur), numel(sig_roi), nSel), 'FontSize', 9, 'FontWeight', 'normal');
ax = gobjects(9, 1);
ax(1) = nexttile(tl, [3 1]);  imagesc(tE, 1:numel(sig_roi), dffS(sig_roi(ord_sig), :));  clim([-0.5 max(prctile(dffS(sig_roi, :), 99.5, 'all'), 1)]);  colormap(ax(1), turbo);
yticks(1:numel(sig_roi));  yticklabels(sel_lab);  ylabel('ROI (clustered)');  cb = colorbar;  cb.Label.String = 'dF/F (0.5 s)';  cb.FontSize = 7;
ax(2) = nexttile(tl);  hold on;  plot(tE, Uz(1, :), 'Color', col.pc1);  plot(tE, Uz(2, :), 'Color', col.dff);  ylabel('z');
text(0.01, 0.88, 'PC1 (signal-like)', 'Color', col.pc1, 'Units', 'normalized', 'FontSize', 7);  text(0.15, 0.88, 'mean of selected', 'Color', col.dff, 'Units', 'normalized', 'FontSize', 7);
ax(3) = nexttile(tl);  plot(tE, fillx(Braw(1, :)), 'Color', col.pupil);  ylabel('pupil (mm^2)');
ax(4) = nexttile(tl);  hold on;  plot(tE, exp(fillx(Braw(2, :))), 'Color', col.whisk);  yline(whisk_thr, ':');  ylabel('whisk idx');  set(gca, 'YScale', 'log');
ax(5) = nexttile(tl);  plot(tE, exp(fillx(Braw(3, :))), 'Color', col.face);  ylabel('face ME');  set(gca, 'YScale', 'log');
ax(6) = nexttile(tl);  hold on;
if have_loc; plot(tE, spg, 'Color', col.loc);  yline(loc_run_thr, ':');  plot(tE(run_on), 0 * run_on + max(spg) * 1.05, 'v', 'Color', col.loc, 'MarkerSize', 3, 'MarkerFaceColor', col.loc); else; text(0.5, 0.5, 'locomotion not available', 'Units', 'normalized', 'HorizontalAlignment', 'center'); end
ylabel('speed (mm/s)');  text(0.01, 0.88, loc_name, 'Units', 'normalized', 'FontSize', 6, 'Color', [0.4 0.4 0.4]);
ax(7) = nexttile(tl);  hold on;  plot(tE, gridify(tcam, eye_spd), 'Color', col.eye);  plot(tE(s_idx), 0 * s_idx + max(gridify(tcam, eye_spd)) * 1.05, 'v', 'Color', col.eye, 'MarkerSize', 3, 'MarkerFaceColor', col.eye);  ylabel('eye (deg/s)');
ax(8) = nexttile(tl);  if have_motion; plot(tE, exp(fillx(Braw(6, :))) - 1, 'Color', col.mot); end;  ylabel('shift (px/s)');  xlabel('epoch time (s)');
linkaxes(ax, 'x');  xlim(ax(1), [tE(1) tE(end)]);  set(ax(1:7), 'XTickLabel', []);
if save_figs; savef(fh, '01_overview'); end

% ---- 2 correlations
fh = newfig('Correlations', 18, 15);
tl = tiledlayout(fh, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'correlation with behaviour: CCF of PC1 (magenta) and the mean of the selected ROIs (black); grey = 99th pct of the circular-shift null', 'FontSize', 8, 'FontWeight', 'normal');
for jj = 1:min(nB, 6)
    nexttile;  hold on;  res1 = ccf_store{1, jj};  res2 = ccf_store{2, jj};  lg = res1.lags * dt_a;
    fill([lg fliplr(lg)], [res1.band * ones(1, numel(lg)) -res1.band * ones(1, numel(lg))], [0.85 0.85 0.85], 'EdgeColor', 'none');
    plot(lg, res1.r, 'Color', col.pc1);  plot(lg, res2.r, 'Color', col.dff);  xline(0, ':');  yline(0, ':');
    title(sprintf('%s: r0 %.2f / %.2f', bn{jj}, res1.r0, res2.r0), 'FontWeight', 'normal', 'FontSize', 7);  ylim([-0.3 0.9]);  xlim([-corr_max_lag_s corr_max_lag_s]);
    if jj > 3; xlabel('lag (s), + = dF/F follows'); end
end
nexttile;  imagesc(Cbb);  clim([-1 1]);  colormap(gca, div_cmap);  axis image;  xticks(1:nB);  xticklabels(bn);  xtickangle(45);  yticks(1:nB);  yticklabels(bn);
for a = 1:nB; for b = 1:nB; text(b, a, sprintf('%.1f', Cbb(a, b)), 'HorizontalAlignment', 'center', 'FontSize', 6); end; end
title('behaviour vs behaviour', 'FontWeight', 'normal', 'FontSize', 7);
nexttile;  hold on;
tb1 = spont_tbl(strcmp(spont_tbl.unit, u_name{1}), :);  tb2 = spont_tbl(strcmp(spont_tbl.unit, u_name{2}), :);
bar_v = [tb1.r_lag0, tb1.r_partial, tb2.r_lag0, tb2.r_partial];
hb = bar(bar_v, 'grouped');  hb(1).FaceColor = col.pc1;  hb(2).FaceColor = (col.pc1 + 1) / 2;  hb(3).FaceColor = col.dff;  hb(4).FaceColor = [0.6 0.6 0.6];
xticks(1:nB);  xticklabels(bn);  xtickangle(45);  ylabel('r at lag 0');  legend({'PC1 raw', 'PC1 partial', 'mean raw', 'mean partial'}, 'FontSize', 6, 'Location', 'northoutside', 'NumColumns', 2);
nexttile;  Rg = reshape(spont_tbl.r_lag0, nB, nU)';  imagesc(Rg);  clim([-0.8 0.8]);  colormap(gca, div_cmap);  colorbar;  xticks(1:nB);  xticklabels(bn);  xtickangle(45);  yticks(1:nU);
yticklabels(cellfun(@(s) s(1:min(end, 24)), u_name, 'UniformOutput', false));  title('r at lag 0, every unit', 'FontWeight', 'normal', 'FontSize', 7);
if save_figs; savef(fh, '02_correlations'); end

% ---- 3 encoding
fh = newfig('Encoding model', 22, 6.5);
tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;  hold on;
sv = enc_R2(1:2, 1 + nB + (1:nB))';  slo = sv - squeeze(enc_ci(1:2, 1 + nB + (1:nB), 1))';  shi = squeeze(enc_ci(1:2, 1 + nB + (1:nB), 2))' - sv;
hb = bar(sv, 'grouped');  hb(1).FaceColor = col.pc1;  hb(2).FaceColor = col.dff;
xe = cell2mat(arrayfun(@(h) h.XEndPoints', hb, 'UniformOutput', false));  errorbar(xe, sv, slo, shi, 'k.', 'CapSize', 3);
xticks(1:nB);  xticklabels(bn);  xtickangle(45);  ylabel('held-out R^2, behaviour alone (all lags)');  yline(0);
title(sprintf('single: full model R^2 PC1 %.2f, mean %.2f', enc_R2(1, 1), enc_R2(2, 1)), 'FontWeight', 'normal', 'FontSize', 7);
nexttile;  hold on;
uv = (enc_R2(1:2, 1) - enc_R2(1:2, 1 + (1:nB)))';  ulo = uv - squeeze(enc_dci(1:2, :, 1))';  uhi = squeeze(enc_dci(1:2, :, 2))' - uv;
hb = bar(uv, 'grouped');  hb(1).FaceColor = col.pc1;  hb(2).FaceColor = col.dff;
xe = cell2mat(arrayfun(@(h) h.XEndPoints', hb, 'UniformOutput', false));  errorbar(xe, uv, ulo, uhi, 'k.', 'CapSize', 3);
xticks(1:nB);  xticklabels(bn);  xtickangle(45);  ylabel('unique held-out R^2 (full - without)');  yline(0);
title('unique contribution', 'FontWeight', 'normal', 'FontSize', 7);
nexttile;  hold on;
rs_ = [1, 3:size(Sg, 1)];  sv2 = nan(numel(rs_), nst3);  slo2 = sv2;  shi2 = sv2;
for i = 1:numel(rs_)
    for k = 1:nst3; q = (rs_(i) - 1) * nst3 + k;  sv2(i, k) = state_tbl.mean(q);  slo2(i, k) = state_tbl.mean(q) - state_tbl.ci_lo(q);  shi2(i, k) = state_tbl.ci_hi(q) - state_tbl.mean(q); end
end
st_col = [0.8 0.8 0.8; 0.45 0.7 0.45; 0.2 0.4 0.8];
hb = bar(sv2, 'grouped');  for k = 1:nst3; hb(k).FaceColor = st_col(k, :); end
xe = cell2mat(arrayfun(@(h) h.XEndPoints', hb, 'UniformOutput', false));  errorbar(xe, sv2, slo2, shi2, 'k.', 'CapSize', 3);
xticks(1:numel(rs_));  xticklabels(strrep(sg_names(rs_), 'mean of selected ROIs', 'mean sel.'));  xtickangle(45);  ylabel('mean dF/F');  legend(st_lab(1:nst3), 'FontSize', 6, 'Location', 'northoutside', 'NumColumns', 3);
title('by behavioural state', 'FontWeight', 'normal', 'FontSize', 7);
if save_figs; savef(fh, '03_encoding'); end

% ---- 4 event-triggered averages
nE = numel(eta);
if nE > 0
    fh = newfig('Event-triggered averages', 3.2 * (nE + 1), 11);
    tl = tiledlayout(fh, 2, nE + 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    sm3 = @(x) movmean(x, 3);
    for e = 1:nE
        E = eta(e);  nexttile(e);  hold on;
        M1 = squeeze(E.S(1, :, :));
        fill([tt fliplr(tt)], [E.null_lo fliplr(E.null_hi)], [0.85 0.85 0.85], 'EdgeColor', 'none');
        for g = 1:ng; plot(tt, sm3(mean(squeeze(E.S(2 + g, :, :)), 2))', 'Color', [unit_col(2 + g, :) 0.6], 'LineWidth', 0.5); end
        mu = sm3(mean(M1, 2))';  se = sm3(std(M1, 0, 2) / sqrt(E.n))';
        fill([tt fliplr(tt)], [mu - se fliplr(mu + se)], col.dff, 'FaceAlpha', 0.2, 'EdgeColor', 'none');  plot(tt, mu, 'Color', col.dff, 'LineWidth', 1.5);
        xline(0, ':');  yline(0, ':');  title(sprintf('%s (n = %d)\n%+.2f dF/F, p %.3f', E.name, E.n, E.amp(1), E.p(1)), 'FontWeight', 'normal', 'FontSize', 6.5);
        if e == 1; ylabel('dF/F, mean of selected ROIs (groups thin)'); end
        nexttile(nE + 1 + e);  hold on;
        for jj = 1:nB; plot(tt, sm3(mean(squeeze(E.B(jj, :, :)), 2))', 'Color', beh_col(bidx(jj), :)); end
        xline(0, ':');  yline(0, ':');  xlabel('time from event (s)');  if e == 1; ylabel('behaviour (z, baseline subtracted)'); end
    end
    nexttile(nE + 1);  hold on;
    if ~isempty(pkB)
        for jj = 1:nB; plot(pt, mean(squeeze(pkB(jj, :, :)), 2), 'Color', beh_col(bidx(jj), :)); end
        xline(0, ':');  yline(0, ':');  title(sprintf('behaviour around %d PC1 peaks', numel(pk_idx)), 'FontWeight', 'normal', 'FontSize', 7);  xlabel('time from PC1 peak (s)');
        legend(bn, 'FontSize', 6, 'Location', 'southeast');
    end
    if save_figs; savef(fh, '04_event_triggered'); end
end

% ---- 5 controls
fh = newfig('Controls', 22, 11);
tl = tiledlayout(fh, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;  hold on;
if ~isempty(ctl_eta)
    for k = 1:2; plot(tt, movmean(mean(squeeze(ctl_eta(k, :, :)), 2), 3), 'LineWidth', 1.2); end
    xline(0, ':');  yline(0, ':');  legend({'ROI F - 0.5 F_{neu}', 'neuropil F'}, 'FontSize', 6, 'Location', 'northwest');  ylabel('% of median F');  xlabel('time from running onset (s)');  title(sprintf('raw fluorescence, %d running onsets', size(ctl_eta, 3)), 'FontWeight', 'normal', 'FontSize', 7);
end
nexttile;  hold on;
if ~isempty(ctl_pk)
    for k = 1:2; plot(pt, mean(squeeze(ctl_pk(k, :, :)), 2), 'LineWidth', 1.2); end
    xline(0, ':');  yline(0, ':');  ylabel('% of median F');  xlabel('time from PC1 peak (s)');  title(sprintf('raw fluorescence, %d PC1 peaks', size(ctl_pk, 3)), 'FontWeight', 'normal', 'FontSize', 7);
end
nexttile;  hold on;
if ~isempty(ctl_eta); plot(tt, movmean(mean(squeeze(ctl_eta(3, :, :)), 2), 3), 'Color', col.mot);  xline(0, ':');  yline(0, ':');  ylabel('rigid shift (px)');  xlabel('time from running onset (s)');  title('suite2p rigid shift', 'FontWeight', 'normal', 'FontSize', 7); end
nexttile;  hold on;
if ~isempty(jm) && ~isempty(jl)
    scatter(roi_r(~signal_like, jl), roi_pr(~signal_like, jl), 18, [0.7 0.7 0.7], 'filled');  scatter(roi_r(signal_like, jl), roi_pr(signal_like, jl), 22, [0.2 0.2 0.2], 'filled');
    scatter(roi_r(sel_roi, jl), roi_pr(sel_roi, jl), 30, col.pc1, 'filled');  plot([-0.3 0.9], [-0.3 0.9], 'k:');
    for r = 1:nROI; text(roi_r(r, jl), roi_pr(r, jl), sprintf(' %d', r), 'FontSize', 5); end
    xline(0, ':');  yline(0, ':');  xlabel('r(dF/F, locomotion)');  ylabel('r | brain motion');  title('every ROI (red = selected, grey = not signal-like)', 'FontWeight', 'normal', 'FontSize', 7);
end
nexttile;  hold on;
if ~isempty(jm)
    cr = ctl_tbl(strcmp(ctl_tbl.unit, u_name{1}), :);  cm = ctl_tbl(strcmp(ctl_tbl.unit, u_name{2}), :);
    hb = bar([cr.r_lag0 cr.r_partial_motion cm.r_lag0 cm.r_partial_motion], 'grouped');  hb(1).FaceColor = col.pc1;  hb(2).FaceColor = (col.pc1 + 1) / 2;  hb(3).FaceColor = col.dff;  hb(4).FaceColor = [0.6 0.6 0.6];
    xticks(1:height(cr));  xticklabels(cr.behaviour);  xtickangle(45);  ylabel('r at lag 0');  legend({'PC1', 'PC1 | motion', 'mean', 'mean | motion'}, 'FontSize', 6, 'Location', 'northoutside', 'NumColumns', 2);  title('brain motion regressed out', 'FontWeight', 'normal', 'FontSize', 7);
end
nexttile;  imagesc(Csg(ord_sig, ord_sig));  axis image;  clim([-1 1]);  colormap(gca, div_cmap);  title('ROI-ROI r, data', 'FontWeight', 'normal', 'FontSize', 7);  xticks([]);  yticks(1:numel(sig_roi));  yticklabels(sel_lab);
nexttile;  imagesc(Cres(ord_sig, ord_sig));  axis image;  clim([-1 1]);  colormap(gca, div_cmap);  colorbar;  title(sprintf('ROI-ROI r after behaviour (median %.2f)', median(off_res)), 'FontWeight', 'normal', 'FontSize', 7);  xticks([]);  yticks([]);
nexttile;  hold on;
if any(isfinite(k_rise)); scatter(k_rise, k_dec, 30, 'filled', 'MarkerFaceColor', col.dff);  xlabel('10-90% rise (s)');  ylabel('time to half amplitude after peak (s)');  title('PC1 peaks: kinetics', 'FontWeight', 'normal', 'FontSize', 7);  plot([0 3], [0 3], 'k:'); end
if save_figs; savef(fh, '05_controls'); end

% ---- 6 population
fh = newfig('Population', 18, 6.5);
tl = tiledlayout(fh, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;  hold on;  bar(100 * expl(1:min(8, numel(expl))), 'FaceColor', [0.5 0.5 0.5]);  xlabel('PC');  ylabel('% variance');  title(sprintf('PCA of %d signal-like ROIs', numel(sig_roi)), 'FontWeight', 'normal', 'FontSize', 7);
nexttile;  hold on;
scatter(cent(sig_roi, 2), cent(sig_roi, 1), 12 + 0.12 * npix(sig_roi), ld1, 'filled', 'MarkerEdgeColor', 'k');  colormap(gca, turbo);  cb = colorbar;  cb.Label.String = 'PC1 loading';  cb.FontSize = 7;
for i = 1:numel(sig_roi); text(cent(sig_roi(i), 2) + 6, cent(sig_roi(i), 1), sprintf('%d', sig_roi(i)), 'FontSize', 5); end
set(gca, 'YDir', 'reverse');  axis equal;  xlim([0 512]);  ylim([0 512]);  title('ROI centroids (full-frame px)', 'FontWeight', 'normal', 'FontSize', 7);
nexttile;  scatter(dd_, off_, 8, [0.3 0.3 0.3], 'filled');  xlabel('centroid distance (px)');  ylabel('r (0.5 s smoothed)');  yline(grp_r_thr, ':');  title('pairwise correlation vs distance', 'FontWeight', 'normal', 'FontSize', 7);
if save_figs; savef(fh, '06_population'); end

%% 9. Tables
if save_figs
    writetable(spont_tbl, fullfile(fig_dir, 'correlations.csv'));  writetable(enc_tbl, fullfile(fig_dir, 'encoding.csv'));  writetable(roi_tbl, fullfile(fig_dir, 'roi_behaviour_r.csv'));
    if ~isempty(eta); writetable(eta_tbl, fullfile(fig_dir, 'event_triggered.csv')); end
    if exist('ctl_tbl', 'var'); writetable(ctl_tbl, fullfile(fig_dir, 'motion_partial.csv')); end
    writetable(state_tbl, fullfile(fig_dir, 'states.csv'));  writetable(pk_tbl, fullfile(fig_dir, 'pc1_peaks.csv'));
end
fprintf('\naxon_analysis_spont_V1 done: %s, %d ROIs selected %s, locomotion = %s\n', session, nSel, lst(sel_roi), loc_name);

%% local functions
function save_fig(fh, base)
% PNG (300 dpi) and vector PDF of a figure
exportgraphics(fh, [base '.png'], 'Resolution', 300);
exportgraphics(fh, [base '.pdf'], 'ContentType', 'vector');
end

function xb = bin_mean(t, x, edges)
% mean of x inside each bin of edges (row vector, NaN where no sample falls into a bin)
ib = discretize(t(:), edges(:));  ok = ~isnan(ib) & isfinite(x(:));
nb = numel(edges) - 1;
s = accumarray(ib(ok), x(ok)', [nb 1], @sum, 0);  c = accumarray(ib(ok), 1, [nb 1]);
xb = (s ./ c)';  xb(c' == 0) = NaN;
end

function C = circ_xcorr(X, y)
% circular cross-correlation of the z-scored columns of X [n x k] with the z-scored column y [n x 1], FFT; C(m+1, :) = Pearson r between X(t+m) and y(t)
% (X follows y by m samples), m = 0 .. n-1. A circular shift of y by s only re-indexes it: r_s(m) = C(m+s+1).
n = size(X, 1);
C = real(ifft(fft(X) .* conj(fft(y)))) / (n - 1);
end

function res = ccf_circ(x, y, maxlag, minshift)
% cross-correlation of x (dF/F) with y (behaviour) at lags -maxlag..maxlag samples (positive = x follows y) and its circular-shift null:
%   p0     |r| at lag 0 vs |r| at lag 0 of the behaviour shifted by every s in minshift .. n - minshift (exact, no random draw)
%   p_best max |r| over the lag range vs the same maximum for every shifted behaviour
x = (x(:) - mean(x(:))) / std(x(:));  y = (y(:) - mean(y(:))) / std(y(:));
n = numel(x);  c = circ_xcorr(x, y);
lags = -maxlag:maxlag;  r = c(mod(lags, n) + 1)';
[~, ib] = max(abs(r));
sh = minshift:(n - minshift);
cs = c(sh + 1);
mx = max(abs(c(mod(sh(:) + lags, n) + 1)), [], 2);
res = struct('r', r, 'lags', lags, 'r0', r(lags == 0), 'p0', (1 + sum(abs(cs) >= abs(r(lags == 0)))) / (numel(cs) + 1), 'r_best', r(ib), 'lag_best', lags(ib), ...
    'p_best', (1 + sum(mx >= abs(r(ib)))) / (numel(mx) + 1), 'band', prctile(abs(cs), 99));
end

function e = resid_on(x, Z)
% residual of x [n x 1] after least squares on Z [n x k] and an intercept
A = [ones(size(Z, 1), 1) Z];
e = x - A * (A \ x);
end

function ci = boot_r(x, y, blk, nb)
% 95% interval of the Pearson r from a block bootstrap (blocks of blk samples)
n = numel(x);  nbk = ceil(n / blk);  rr = zeros(nb, 1);
for b = 1:nb
    pick = randi(nbk, 1, nbk);  idx = [];
    for k = pick; idx = [idx, (k - 1) * blk + 1:min(k * blk, n)]; end %#ok<AGROW>
    rr(b) = corr(x(idx)', y(idx)');
end
ci = prctile(rr, [2.5 97.5]);
end

function [gid, grp_rois] = roi_groups(dsm, sel_roi, thr)
% cluster the selected ROIs on the correlation of their smoothed dF/F (average linkage, cut at 1 - thr); group 1 = the largest group
nSel = numel(sel_roi);
ds = double(dsm(sel_roi, :));
if nSel >= 2; gid = cluster(linkage(pdist(ds, 'correlation'), 'average'), 'cutoff', 1 - thr, 'criterion', 'distance'); else; gid = ones(nSel, 1); end
gid = gid(:)';
gsz = accumarray(gid(:), 1)';  gmn = accumarray(gid(:), sel_roi(:), [], @min)';
[~, gord] = sort(-1000 * gsz + gmn);
[~, gid] = ismember(gid, gord);
grp_rois = arrayfun(@(g) sel_roi(gid == g), 1:max(gid), 'UniformOutput', false);
end

function X = lag_matrix(Bc, lagb)
% predictors [n x nB*nl]: behaviour j at t - lag (positive lag = behaviour earlier), edge samples repeated; columns behaviour-major
[nB, n] = size(Bc);  X = zeros(n, nB * numel(lagb));  c = 0;
for j = 1:nB
    for l = 1:numel(lagb)
        c = c + 1;  X(:, c) = Bc(j, min(max((1:n) - lagb(l), 1), n))';
    end
end
end

function yhat = cv_ridge(y, X, nfold, margin, lams)
% out-of-fold ridge predictions: nfold contiguous folds, `margin` samples around the test fold removed from the training set, lambda by an inner blocked 4-fold CV
n = numel(y);  e = round(linspace(0, n, nfold + 1));  yhat = nan(n, 1);
for f = 1:nfold
    te = e(f) + 1:e(f + 1);
    trn = setdiff(1:n, max(1, te(1) - margin):min(n, te(end) + margin));
    lam = pick_lambda(y(trn), X(trn, :), 4, margin, lams);
    [Xtr, mu, sd] = zscore_cols(X(trn, :));  ym = mean(y(trn));
    w = (Xtr' * Xtr + lam * eye(size(Xtr, 2))) \ (Xtr' * (y(trn) - ym));
    yhat(te) = ((X(te, :) - mu) ./ sd) * w + ym;
end
end

function lam = pick_lambda(y, X, k, margin, lams)
n = numel(y);  e = round(linspace(0, n, k + 1));  err = zeros(1, numel(lams));
for f = 1:k
    te = e(f) + 1:e(f + 1);  trn = setdiff(1:n, max(1, te(1) - margin):min(n, te(end) + margin));
    [Xtr, mu, sd] = zscore_cols(X(trn, :));  ym = mean(y(trn));  Xte = (X(te, :) - mu) ./ sd;
    XtX = Xtr' * Xtr;  Xty = Xtr' * (y(trn) - ym);
    for l = 1:numel(lams); w = (XtX + lams(l) * eye(size(XtX, 1))) \ Xty; err(l) = err(l) + sum((y(te) - (Xte * w + ym)) .^ 2); end
end
[~, il] = min(err);  lam = lams(il);
end

function [Z, mu, sd] = zscore_cols(X)
mu = mean(X, 1);  sd = std(X, 0, 1);  sd(sd == 0) = 1;  Z = (X - mu) ./ sd;
end

function out = r2_boot(y, Yh, blk, nb)
% held-out R^2 of each column of Yh and its 95% interval (bootstrap over blocks of blk samples; the draws are returned for differences between models)
n = numel(y);  nbk = ceil(n / blk);
out.r2 = 1 - sum((y - Yh) .^ 2, 1) / sum((y - mean(y)) .^ 2);
bt = zeros(nb, size(Yh, 2));
for b = 1:nb
    pick = randi(nbk, 1, nbk);  idx = [];
    for k = pick; idx = [idx, (k - 1) * blk + 1:min(k * blk, n)]; end %#ok<AGROW>
    yy = y(idx);  bt(b, :) = 1 - sum((yy - Yh(idx, :)) .^ 2, 1) / sum((yy - mean(yy)) .^ 2);
end
out.ci = prctile(bt, [2.5 97.5], 1);  out.boot = bt;
end

function [A1, A2] = cut_events(X1, X2, idx, offs)
% traces cut around events: X1 [r1 x n], X2 [r2 x n] (may have 0 rows), idx event indices, offs sample offsets -> A [rows x numel(offs) x numel(idx)]
ii = (idx(:) + offs(:)')';
A1 = reshape(X1(:, ii), size(X1, 1), numel(offs), numel(idx));
A2 = reshape(X2(:, ii), size(X2, 1), numel(offs), numel(idx));
end

function lat = half_rise(tt, m, amp_win)
% time at which the mean event-triggered trace m first reaches half of its peak (peak of the sign of its mean response in 0..amp_win(2)+1 s), linear interpolation; NaN if none
s = sign(mean(m(tt >= amp_win(1) & tt <= amp_win(2))));  if s == 0; s = 1; end
mm = s * m;  pk = max(mm(tt >= 0 & tt <= amp_win(2) + 1));
lat = NaN;
if ~(pk > 0); return; end
i = find(tt >= -1 & mm >= 0.5 * pk, 1);
if isempty(i) || i == 1; return; end
lat = interp1(mm([i - 1 i]), tt([i - 1 i]), 0.5 * pk);
end

function i = ternary_idx(ta, x, dt)
% index of the grid bin (centres ta) that contains the time x; NaN if x is not inside the epoch windows
[d, i] = min(abs(ta - x));
if d > dt; i = NaN; end
end

function t = ternary_str(c, a, b)
if c; t = a; else; t = b; end
end
