%% axon_analysis_behaviour_V1
% How much of the axon activity can be explained by behaviour? (after Collins, Francis, Emanuel & McCormick 2023, eLife 12:e81826, Fig. 3-4: ridge regression of
% LC-NA / BF-ACh axon dF/F on pupil, whisker-pad / face motion and walking; variance explained per behaviour; shared signal across axon segments before / after
% removing the behaviour.) For sessions aggregated with newAggregation/save_preprocessed.m; the whole session is used (every stimulus block, the grey windows
% around the looming and the corridor), not one epoch.
%
% Model: dF/F(t) of a unit = ridge regression on the behaviours at lags enc_lags_s (positive lag = behaviour first), 1 s smoothing, dt_a grid. Held-out R^2:
% enc_nfold contiguous folds of the session with enc_margin_s purged around the test fold, lambda from an inner blocked CV, so a R^2 is always a prediction of
% time that the model has not seen. Units: the mean of the selected ROIs, each ROI group (ROIs with r > grp_r_thr are one axon) and every single ROI.
%   alone     R^2 with one behaviour only (all its lags)                   -> what each behaviour could explain by itself
%   unique    R^2(full) - R^2(without that behaviour)                      -> what only this behaviour explains (the behaviours are collinear)
%   shared    R^2(full) - sum of the unique parts                          -> explained, but not attributable to one behaviour
%   null      the same model with the dF/F circularly shifted by >= null_min_shift_s relative to the behaviour (the paper's session permutation)
% Further: R^2 per epoch type (the pooled model, within-epoch variance), correlation between ROIs before / after removing the behaviour (paper Fig. 4F),
% control with the suite2p brain motion, and, if the ball log is usable in part of the session (the open-loop blocks are truncated, only the corridor is complete),
% the same model restricted to those blocks with locomotion as a behaviour (model B).
%
% Statistics: one session = one animal. 95% intervals = block bootstrap (boot_blk_s) of the held-out predictions. Exploratory.
% Figures: 1 = summary (example, variance partition, per-behaviour, per ROI), 2 = details and controls. Outputs in the workspace: enc_tbl, roi_tbl, epoch_tbl,
% sel_roi, groups, R2 (struct with the model fits).

%% 0. Settings (edit only this block)
session   = 'AnimalRB24_261023_1538';
agg_dir   = 'Z:\Group Members\Rima\Aggregated';
s2p_dir   = fullfile('Z:\Group Members\Rima\DATA_2P', session, 'suite2p', 'plane0');          % Fall.mat: rigid shift for the brain-motion control
repo_root = 'C:\Users\rbondarenko\projects\2P_Scripts';
save_figs = false;                                                % true: PNG (300 dpi) + vector PDF of both figures and the tables as CSV in fig_dir
fig_dir   = fullfile(repo_root, 'axon analysis', 'plots', 'behaviour_V1', session);

roi_sel   = 'good';       % ROIs of the mean / groups: 'good' (signal-like, not motion-coupled, as V3) | 'signal' | 'all' | a list, e.g. [9 18 20 21 23 24 25 27]
beh_use   = {'PupilArea', 'WhiskIndex', 'FacePadME', 'EyeLidOpening'};   % camera behaviours (fields of Behav) of the main model
log_beh   = {'WhiskIndex', 'FacePadME'};                          % heavy-tailed indices: log transform (as axon_analysis_spont_V1)
use_loc   = true;         % ball speed: in the main model if it covers >= loc_cov_main of the analysed time, else only in model B (blocks with a usable ball log)
loc_cov_main = 0.7;
min_loc_s    = 300;       % model B needs at least this much usable ball time
dt_a      = 0.2;          % analysis grid (s)
smooth_s  = 1;            % smoothing of dF/F and behaviour (s); the paper low-passes at 1 Hz
max_gap_s = 1;            % gaps of the behaviour up to this long (blinks, lost fits) are interpolated, longer ones are left out of the analysis
enc_lags_s   = -1:0.5:3;  % behaviour at t - lag predicts dF/F at t
enc_nfold    = 10;        % contiguous folds of the session (model B: enc_nfold_B)
enc_nfold_B  = 5;
enc_margin_s = 5;         % purge around the test fold
lams      = 10 .^ (0:0.5:6);
n_boot    = 1000;         % block bootstrap draws
boot_blk_s = 30;
n_null    = 30;           % shifted-dF/F nulls for the mean and the groups
null_min_shift_s = 120;
ex_win_s  = 150;          % length of the example stretch (figure 1A)
corr_max_lag_s = 3;       % cross-correlation range (figure 2A)
max_groups_shown = 6;
slow_ratio_thr = 4;  mot_r_thr = 0.15;  smooth_roi_s = 0.5;  neuropil_factor = 0.5;   % ROI quality, as in V3
grp_r_thr = 0.8;
rng(1);

col = struct('dff', [0 0 0], 'pupil', [0.85 0.15 0.15], 'whisk', [0.1 0.55 0.25], 'face', [0.9 0.55 0.05], 'eye', [0.5 0.3 0.7], 'loc', [0.1 0.35 0.8], 'mot', [0.45 0.45 0.45], 'pred', [0.8 0.1 0.5]);
beh_cmap = struct('PupilArea', col.pupil, 'WhiskIndex', col.whisk, 'FacePadME', col.face, 'EyeLidOpening', col.eye, 'LocSpeed', col.loc);
beh_lab  = struct('PupilArea', 'pupil', 'WhiskIndex', 'whisking', 'FacePadME', 'face motion', 'EyeLidOpening', 'eyelid', 'LocSpeed', 'locomotion');
beh_short = struct('PupilArea', 'pupil', 'WhiskIndex', 'whisk', 'FacePadME', 'face', 'EyeLidOpening', 'eyelid', 'LocSpeed', 'loc');
gcol = [0.15 0.15 0.15; 0.0 0.45 0.45; 0.65 0.4 0.0; 0.45 0.45 0.8; 0.6 0.6 0.2; 0.55 0.2 0.2; 0.3 0.6 0.3; 0.6 0.3 0.6];
newfig = @(nm, w, h) figure('Name', nm, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 w h], 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 8, ...
    'DefaultAxesBox', 'off', 'DefaultAxesTickDir', 'out', 'DefaultAxesLineWidth', 0.75, 'DefaultLineLineWidth', 0.75, 'DefaultTextFontName', 'Arial', 'DefaultLegendBox', 'off');
lst = @(v) ['[' strtrim(sprintf('%g ', v)) ']'];
zs  = @(x) (x - mean(x, 'omitnan')) ./ std(x, 'omitnan');
letter = @(ax, s) panel_letter(ax, s);
ep_short = @(n) strrep(strrep(strrep(strrep(strrep(n, 'sparse_local_global_flashes', 'flashes'), 'full_field_flash', 'chirp'), 'moving_bar', 'bar sweep'), 'alleyCL', 'corridor'), '_', ' ');

%% 1. Load, ROI quality (as V3 / axon_analysis_spont_V1)
addpath(repo_root, fullfile(repo_root, 'analysis_functions'));
if save_figs && ~isfolder(fig_dir); mkdir(fig_dir); end
savef = @(fh, nm) save_fig(fh, fullfile(fig_dir, nm));
data_file = fullfile(agg_dir, [session '_preprocessed.mat']);
assert(isfile(data_file), 'Aggregated file not found: %s', data_file);
fmeta = whos('-file', data_file);  inFile = {fmeta.name};
want  = {'CaData', 'TimeCa', 'Stimuli', 'Triggers', 'Behav', 'LocomotionCal', 'LocomotionNonCal', 'AnimalInfo'};
fprintf('Loading %s ...\n', data_file);
data = load(data_file, want{ismember(want, inFile)});
assert(isfield(data, 'Behav') && ~isempty(data.Triggers.TimeCamera), 'No behaviour in this file');
nz = find(~cellfun(@isempty, {data.CaData.Ca_dFF}));  assert(~isempty(nz), 'CaData has no calcium traces');
ca = data.CaData(nz(1));
dFF = double(ca.Ca_dFF);  F = double(ca.Ca_F);  Fneu = double(ca.Ca_F_neuropil);
[nROI, nT] = size(dFF);
tCa = double(data.TimeCa(1, :));  fs = 1 / median(diff(tCa));
tcam = double(data.Triggers.TimeCamera(:))';  Bh = data.Behav;
if isfield(data, 'AnimalInfo'); ai = data.AnimalInfo; fprintf('\n%s | %s | %s\n', ai.AnimalName, ai.Genotype, ai.Indicator.Type); end

Fc = F - neuropil_factor * Fneu;
w_sm = max(1, round(smooth_roi_s * fs));  dsm = movmean(dFF, w_sm, 2);
med_s = median(dsm, 2);  dev = dsm - med_s;  neg = dev < 0;
sig_neg = max(sqrt(sum(dev .^ 2 .* neg, 2) ./ max(sum(neg, 2), 1)), eps('single'));
dd = diff(dFF, 1, 2);  hf_noise = 1.4826 * median(abs(dd - median(dd, 2)), 2) / sqrt(2);
slow_ratio = var(dsm, 0, 2) ./ (hf_noise .^ 2 / w_sm);
ev_pos = sum(diff([false(nROI, 1) dev > 3 * sig_neg], 1, 2) == 1, 2);  ev_neg = sum(diff([false(nROI, 1) dev < -3 * sig_neg], 1, 2) == 1, 2);
bad_F = mean(Fc <= 0, 2) > 0.01;
signal_like = ~bad_F & ~(slow_ratio < slow_ratio_thr | ev_pos < 2 * max(ev_neg, 1));
have_motion = false;  rigid = nan(1, nT);  mot_flag = false(nROI, 1);  xo = nan(1, nT);  yo = nan(1, nT);
fall_file = fullfile(s2p_dir, 'Fall.mat');
if isfile(fall_file)
    o = load(fall_file, 'ops');  xo = double(o.ops.xoff(:))';  yo = double(o.ops.yoff(:))';
    if numel(xo) == nT
        have_motion = true;
        rigid = hypot(xo - median(xo), yo - median(yo));
        nonr  = hypot(std(double(o.ops.xoff1), 0, 2), std(double(o.ops.yoff1), 0, 2))';
        mots = {rigid, nonr};
        dz_  = zs(double(dsm)');
        for m = 1:2
            mz_ = zs(movmean(mots{m}, w_sm));
            c_all = circ_xcorr(dz_, mz_');
            sh = (round(30 * fs):nT - round(30 * fs)) + 1;
            thr99 = prctile(abs(c_all(sh, :)), 99, 1)';
            mot_flag = mot_flag | abs(c_all(1, :)') > max(thr99, mot_r_thr);
        end
    else
        fprintf('suite2p ops has %d motion samples for %d frames: motion not used\n', numel(xo), nT);
        xo = nan(1, nT);  yo = nan(1, nT);
    end
else
    fprintf('%s not found: motion-coupled ROIs cannot be flagged, no brain-motion control\n', fall_file);
end
good_roi = find(signal_like & ~mot_flag)';  sig_roi = find(signal_like)';
switch lower(char(string(roi_sel)))
    case 'good';   sel_roi = good_roi;
    case 'signal'; sel_roi = sig_roi;
    case 'all';    sel_roi = 1:nROI;
    otherwise;     sel_roi = double(roi_sel(:))';
end
assert(~isempty(sel_roi), 'No ROI selected (roi_sel = %s; signal-like %s, motion-coupled %s)', char(string(roi_sel)), lst(sig_roi), lst(find(mot_flag)'));
nSel = numel(sel_roi);
fprintf('%d ROIs, %d frames (%.0f s, %.1f Hz) | signal-like %s | ''good'' %s | used %s\n', nROI, nT, tCa(end) - tCa(1), fs, lst(sig_roi), lst(good_roi), lst(sel_roi));

%% 2. Blocks, epochs and the analysis grid
S = data.Stimuli;  nBk = numel(S);
b_type = {S.type};  b_t0 = nan(1, nBk);  b_t1 = nan(1, nBk);  b_pre = zeros(1, nBk);  b_post = zeros(1, nBk);  b_dir = repmat({''}, 1, nBk);
for k = 1:nBk
    tf = double(S(k).TimeStimulusFrame(:));
    if nnz(isfinite(tf)) > 100; b_t0(k) = min(tf);  b_t1(k) = max(tf); end                % blocks of a stopped session have no frame times
    sp_ = struct();  if isfield(S, 'specParams') && isstruct(S(k).specParams); sp_ = S(k).specParams; end
    if isfield(sp_, 'pre_grey_t') && isscalar(sp_.pre_grey_t) && isfinite(sp_.pre_grey_t); b_pre(k) = sp_.pre_grey_t; end
    if isfield(sp_, 'post_grey_t') && isscalar(sp_.post_grey_t) && isfinite(sp_.post_grey_t); b_post(k) = sp_.post_grey_t; end
    if strcmp(b_type{k}, 'looming') && isfield(sp_, 'loom_maxexpand') && ~isempty(sp_.loom_maxexpand) && sp_.loom_maxexpand(1) < 0; b_dir{k} = '_receding'; end
end
live = find(isfinite(b_t0));
assert(~isempty(live), 'No stimulus block with frame times');
edges = tCa(1):dt_a:tCa(end);  nG = numel(edges) - 1;  tg = edges(1:end-1) + dt_a / 2;
ep_lbl = repmat({''}, 1, nG);                                                  % epoch label of every grid bin ('' = between blocks)
for k = live
    ia = tg >= b_t0(k) - b_pre(k) & tg <= b_t1(k) + b_post(k);
    ib = tg >= b_t0(k) & tg <= b_t1(k);
    ep_lbl(ia) = {'grey'};  ep_lbl(ib) = {[b_type{k} b_dir{k}]};
end
ep_names = unique(ep_lbl(~cellfun(@isempty, ep_lbl)), 'stable');
ep_id = zeros(1, nG);  for e = 1:numel(ep_names); ep_id(strcmp(ep_lbl, ep_names{e})) = e; end
nEp = numel(ep_names);
ecol = repmat([0.35 0.35 0.35; 0.6 0.6 0.6], ceil(nEp / 2), 1);                % epochs: alternating greys (the colours are for the behaviours)
fprintf('%d blocks with frame times: %s | epochs: %s\n', numel(live), strjoin(b_type(live), ', '), strjoin(ep_names, ', '));

%% 3. Behaviour on the grid
w_g = max(1, round(smooth_s / dt_a));
bname = {};  Bg = zeros(0, nG);                                                % z-scored, smoothed behaviours (NaN where missing)
for i = 1:numel(beh_use)
    if ~isfield(Bh, beh_use{i}); fprintf('Behav.%s not in the file, skipped\n', beh_use{i}); continue; end
    x = double(Bh.(beh_use{i})(:))';
    if strcmp(beh_use{i}, 'PupilArea') && isfield(Bh, 'PupilBad'); x(logical(Bh.PupilBad(:))') = NaN; end
    if ismember(beh_use{i}, log_beh)
        floor_ = 0.01;  if strcmp(beh_use{i}, 'WhiskIndex'); floor_ = 0.1; end                % 1 = rest level of the whisking index
        x = log(max(x, floor_));
    end
    x = fillmissing(x, 'linear', 'SamplePoints', tcam, 'MaxGap', max_gap_s);
    xg = bin_mean(tcam, x, edges);
    xs = movmean(xg, w_g, 'omitnan');  xs(isnan(xg)) = NaN;
    if mean(isfinite(xs)) < 0.3 || std(xs, 'omitnan') == 0; fprintf('Behav.%s: too little valid data, skipped\n', beh_use{i}); continue; end
    bname{end + 1} = beh_use{i};  Bg(end + 1, :) = xs; %#ok<SAGROW>
end
nBeh = numel(bname);  assert(nBeh >= 2, 'Fewer than two behaviours available');

% locomotion: ball speed of the blocks whose log has one reading per ball trigger (otherwise it cannot be aligned)
loc_g = nan(1, nG);  loc_blk = [];
if use_loc && isfield(data, 'LocomotionCal') && isfield(data, 'LocomotionNonCal')
    for k = 1:min(numel(data.LocomotionCal), numel(data.LocomotionNonCal))
        nv = size(data.LocomotionNonCal(k).Values, 1);  tt = double(data.LocomotionNonCal(k).Triggers_sec(:));
        if abs(nv - numel(tt)) > 2 || nv < 2; continue; end
        L = data.LocomotionCal(k);  n_ = min([numel(L.Triggers_sec), numel(L.Forward_mm), numel(L.SideRight_mm)]);
        tt = double(L.Triggers_sec(1:n_))';  dtl = median(diff(tt));
        inc = hypot(double(L.Forward_mm(1:n_))' - median(L.Forward_mm(1:n_)), double(L.SideRight_mm(1:n_))' - median(L.SideRight_mm(1:n_)));
        sp = movmean(inc, max(1, round(0.5 / dtl))) / dtl;
        xg = bin_mean(tt, log1p(max(sp, 0)), edges);
        xs = movmean(xg, w_g, 'omitnan');  xs(isnan(xg)) = NaN;
        ii = isfinite(xs);  loc_g(ii) = xs(ii);  loc_blk(end + 1) = k; %#ok<SAGROW>
    end
end
have_loc = any(isfinite(loc_g));
fprintf('behaviours: %s | ball log usable in blocks %s (%.0f s)\n', strjoin(bname, ', '), lst(loc_blk), sum(isfinite(loc_g)) * dt_a);

% z-score every behaviour over the bins where all camera behaviours exist, then the lagged design matrix (NaN outside the data: lagged rows that
% reach into a pause or over the edge of a block are dropped, not padded)
cam_ok = all(isfinite(Bg), 1) & ep_id > 0;
Bz = (Bg - mean(Bg(:, cam_ok), 2)) ./ std(Bg(:, cam_ok), 0, 2);
lagb = round(enc_lags_s / dt_a);  nl = numel(lagb);
Xcam = lag_matrix(Bz, lagb);
in_main = false;
if have_loc
    lz = (loc_g - mean(loc_g(cam_ok), 'omitnan')) / std(loc_g(cam_ok), 'omitnan');
    Xloc = lag_matrix(lz, lagb);
end

% imaging on the grid
dffG = nan(nROI, nG);  for r = 1:nROI; dffG(r, :) = bin_mean(tCa, dFF(r, :), edges); end
dffS = movmean(dffG, w_g, 2, 'omitnan');
ok_y = all(isfinite(dffS), 1);

vr = find(all(isfinite(Xcam), 2)' & ok_y & ep_id > 0);                         % analysed bins (main model)
if have_loc && use_loc
    loc_cov = mean(all(isfinite(Xloc(vr, :)), 2));
    fprintf('ball speed covers %.0f%% of the analysed time\n', 100 * loc_cov);
    if loc_cov >= loc_cov_main
        in_main = true;  Bz = [Bz; lz];  bname{end + 1} = 'LocSpeed';  nBeh = nBeh + 1;  Xcam = lag_matrix(Bz, lagb);
        vr = find(all(isfinite(Xcam), 2)' & ok_y & ep_id > 0);
    end
end
nV = numel(vr);
assert(nV * dt_a > 120, 'Fewer than 120 s of analysable data (%d bins)', nV);
fprintf('analysed: %d bins = %.0f s of %.0f s (%.0f%%), %d behaviours x %d lags = %d predictors\n', nV, nV * dt_a, nG * dt_a, 100 * nV / nG, nBeh, nl, nBeh * nl);
colb = repelem(1:nBeh, nl);
bsh = cellfun(@(b) beh_short.(b), bname, 'UniformOutput', false);
bcol = cell2mat(cellfun(@(b) beh_cmap.(b), bname(:), 'UniformOutput', false));
Xv = zscore_cols(Xcam(vr, :));                                                 % standardised predictors of the analysed bins

%% 4. Units
[gid, grp_rois] = roi_groups(dffS(:, vr), sel_roi, grp_r_thr);  ng = numel(grp_rois);
groups = grp_rois;
u_name = [{'mean of selected ROIs'}, arrayfun(@(g) sprintf('group %d: ROI %s', g, strjoin(string(grp_rois{g}), ',')), 1:ng, 'UniformOutput', false)];
u_short = [{'mean sel.'}, arrayfun(@(g) sprintf('g%d (n=%d)', g, numel(grp_rois{g})), 1:ng, 'UniformOutput', false)];
UG = [mean(dffS(sel_roi, :), 1); cell2mat(arrayfun(@(g) mean(dffS(grp_rois{g}, :), 1), (1:ng)', 'UniformOutput', false))];       % unit series on the full grid [nU x nG]
nU = size(UG, 1);
fprintf('%d selected ROIs form %d group(s) at r > %g: %s\n', nSel, ng, grp_r_thr, strjoin(cellfun(@(g) lst(g), grp_rois, 'UniformOutput', false), '  '));

%% 5. Encoding models (held-out R^2)
mg = round(enc_margin_s / dt_a);  blk = round(boot_blk_s / dt_a);
mod_names = [{'full'}, strcat('without_', bname), strcat('only_', bname)];  nM = numel(mod_names);
fprintf('\n--- 5. ridge model, %d folds, lags %g to %g s (positive = behaviour first), %.0f s purged around the test fold ---\n', enc_nfold, enc_lags_s(1), enc_lags_s(end), enc_margin_s);
Yh = nan(nV, nM, nU);  Yo = nan(nV, nU);  fold_id = [];
R2 = nan(nU, nM);  CI = nan(nU, nM, 2);  dCI = nan(nU, nBeh, 2);  enc_rows = {};
tic;
for u = 1:nU
    y = zs(UG(u, vr))';  Yo(:, u) = y;
    [Yh(:, 1, u), fold_id] = cv_ridge_v(y, Xv, vr, enc_nfold, mg, lams);
    for j = 1:nBeh; Yh(:, 1 + j, u) = cv_ridge_v(y, Xv(:, colb ~= j), vr, enc_nfold, mg, lams); end
    for j = 1:nBeh; Yh(:, 1 + nBeh + j, u) = cv_ridge_v(y, Xv(:, colb == j), vr, enc_nfold, mg, lams); end
    rb = r2_boot(y, squeeze(Yh(:, :, u)), blk, n_boot);
    R2(u, :) = rb.r2;  CI(u, :, :) = reshape(rb.ci', 1, [], 2);
    for j = 1:nBeh; dCI(u, j, :) = prctile(rb.boot(:, 1) - rb.boot(:, 1 + j), [2.5 97.5]); end
    for m = 1:nM; enc_rows(end + 1, :) = {u_name{u}, mod_names{m}, rb.r2(m), rb.ci(1, m), rb.ci(2, m)}; end %#ok<SAGROW>
end
enc_tbl = cell2table(enc_rows, 'VariableNames', {'unit', 'model', 'R2_heldout', 'ci_lo', 'ci_hi'});
fprintf('%.0f s\n', toc);
uniq = R2(:, 1) - R2(:, 2:1 + nBeh);  alone = R2(:, 2 + nBeh:end);
for u = 1:nU
    fprintf('%s: R^2 full %.2f [%.2f, %.2f]\n', u_name{u}, R2(u, 1), CI(u, 1, 1), CI(u, 1, 2));
    for j = 1:nBeh
        fprintf('   %-13s alone %5.2f [%5.2f, %5.2f] | unique %5.2f [%5.2f, %5.2f]\n', bname{j}, alone(u, j), CI(u, 1 + nBeh + j, 1), CI(u, 1 + nBeh + j, 2), uniq(u, j), dCI(u, j, 1), dCI(u, j, 2));
    end
end

% null: dF/F shifted against the behaviour
null_R2 = nan(n_null, nU);
tic;
for it = 1:n_null
    sh = round(null_min_shift_s / dt_a) + randi(nV - 2 * round(null_min_shift_s / dt_a));      % shift in the order of the analysed bins (no pauses / gaps are shifted into the data)
    for u = 1:nU
        y = circshift(zs(UG(u, vr))', sh);
        yh = cv_ridge_v(y, Xv, vr, enc_nfold, mg, lams);
        null_R2(it, u) = 1 - sum((y - yh) .^ 2) / sum((y - mean(y)) .^ 2);
    end
end
fprintf('null (%d shifts of the dF/F, %.0f s): mean unit R^2 %.3f +- %.3f (95th percentile %.3f), observed %.2f\n', n_null, toc, mean(null_R2(:, 1)), std(null_R2(:, 1)), prctile(null_R2(:, 1), 95), R2(1, 1));

% every ROI: full model (per-ROI R^2, residuals for the shared-signal analysis)
roi_R2 = nan(nROI, 1);  Yh_roi = nan(nV, nROI);  Yo_roi = nan(nV, nROI);
tic;
for r = 1:nROI
    y = zs(dffS(r, vr))';  Yo_roi(:, r) = y;
    Yh_roi(:, r) = cv_ridge_v(y, Xv, vr, enc_nfold, mg, lams);
    roi_R2(r) = 1 - sum((y - Yh_roi(:, r)) .^ 2) / sum((y - mean(y)) .^ 2);
end
fprintf('per ROI (%d ROIs, %.0f s): selected median R^2 %.2f [%.2f-%.2f], other signal-like %.2f, flat / neuropil-limited %.2f\n', nROI, toc, median(roi_R2(sel_roi)), min(roi_R2(sel_roi)), max(roi_R2(sel_roi)), ...
    median(roi_R2(setdiff(sig_roi, sel_roi))), median(roi_R2(~signal_like)));
grp_of = zeros(nROI, 1);  for g = 1:ng; grp_of(grp_rois{g}) = g; end
ceil_r2 = nan(nROI, 1);                                                          % variance an ROI shares with the other ROIs of its group (no behaviour involved)
for g = 1:ng
    if numel(grp_rois{g}) < 2; continue; end
    for r = grp_rois{g}
        oth = setdiff(grp_rois{g}, r);  ceil_r2(r) = corr(Yo_roi(:, r), mean(Yo_roi(:, oth), 2)) ^ 2;
    end
end
roi_tbl = table((1:nROI)', signal_like, mot_flag, ismember((1:nROI)', sel_roi), grp_of, roi_R2, ceil_r2, 'VariableNames', {'roi', 'signal_like', 'motion_coupled', 'selected', 'group', 'R2_full', 'R2_shared_with_group'});

%% 6. Epochs, brain-motion control, locomotion model B, cross-correlation
% R^2 per epoch type from the out-of-fold predictions of the pooled model (SST around the mean of the epoch: how much of the fluctuations inside the epoch)
ep_of = ep_id(vr);  epoch_rows = {};
ep_R2 = nan(nEp, 1);  ep_ci = nan(nEp, 2);  ep_sec = zeros(nEp, 1);  ep_sd = nan(nEp, 1);
for e = 1:nEp
    ii = find(ep_of == e);  ep_sec(e) = numel(ii) * dt_a;
    if numel(ii) * dt_a < 20; continue; end
    rb = r2_boot(Yo(ii, 1), Yh(ii, 1, 1), blk, n_boot);
    ep_R2(e) = rb.r2;  ep_ci(e, :) = rb.ci(:)';  ep_sd(e) = std(UG(1, vr(ii)));
    epoch_rows(end + 1, :) = {ep_names{e}, ep_sec(e), rb.r2, rb.ci(1), rb.ci(2), ep_sd(e)}; %#ok<SAGROW>
end
epoch_tbl = cell2table(epoch_rows, 'VariableNames', {'epoch', 'seconds', 'R2_heldout', 'ci_lo', 'ci_hi', 'sd_dFF'});
fprintf('\n--- 6. R^2 of the mean unit per epoch ---\n');  disp(epoch_tbl);
assert(any(isfinite(ep_R2)), 'No epoch with at least 20 s of analysed data');

% brain motion (suite2p rigid shift: displacement and speed) added to the model
R2_mot = nan(1, 3);                                                             % camera model | + brain motion | brain motion alone
if have_motion
    dispm = bin_mean(tCa, rigid, edges);  spd = bin_mean(tCa, log1p(hypot([0 diff(xo)], [0 diff(yo)]) * fs), edges);
    Mz = zs([movmean(dispm, w_g, 'omitnan'); movmean(spd, w_g, 'omitnan')]')';
    Xm = zscore_cols(lag_matrix(Mz, lagb));  Xm = Xm(vr, :);
    if all(isfinite(Xm(:)))
        y = Yo(:, 1);
        yh1 = cv_ridge_v(y, [Xv, Xm], vr, enc_nfold, mg, lams);  yh2 = cv_ridge_v(y, Xm, vr, enc_nfold, mg, lams);
        sst = sum((y - mean(y)) .^ 2);
        R2_mot = [R2(1, 1), 1 - sum((y - yh1) .^ 2) / sst, 1 - sum((y - yh2) .^ 2) / sst];
        fprintf('brain motion: behaviour only %.2f, behaviour + brain motion %.2f, brain motion alone %.2f\n', R2_mot);
    end
end

% model B: the blocks with a usable ball log, camera behaviours + locomotion
B = struct('ok', false);
if have_loc && use_loc && ~in_main
    XB = [Xcam, Xloc];  vB = find(all(isfinite(XB), 2)' & ok_y & ep_id > 0);
    if numel(vB) * dt_a >= min_loc_s
        B.ok = true;  B.vr = vB;  B.names = [bname, {'LocSpeed'}];  nBB = numel(B.names);  B.col = repelem(1:nBB, nl);
        XBz = zscore_cols(XB(vB, :));
        B.R2 = nan(nU, 1 + 2 * nBB);  B.CI = nan(nU, 1 + 2 * nBB, 2);  B.dCI = nan(nU, nBB, 2);
        for u = 1:nU
            y = zs(UG(u, vB))';  yh = nan(numel(vB), 1 + 2 * nBB);
            yh(:, 1) = cv_ridge_v(y, XBz, vB, enc_nfold_B, mg, lams);
            for j = 1:nBB; yh(:, 1 + j) = cv_ridge_v(y, XBz(:, B.col ~= j), vB, enc_nfold_B, mg, lams); end
            for j = 1:nBB; yh(:, 1 + nBB + j) = cv_ridge_v(y, XBz(:, B.col == j), vB, enc_nfold_B, mg, lams); end
            rb = r2_boot(y, yh, blk, n_boot);
            B.R2(u, :) = rb.r2;  B.CI(u, :, :) = reshape(rb.ci', 1, [], 2);
            for j = 1:nBB; B.dCI(u, j, :) = prctile(rb.boot(:, 1) - rb.boot(:, 1 + j), [2.5 97.5]); end
        end
        B.uniq = B.R2(:, 1) - B.R2(:, 2:1 + nBB);  B.alone = B.R2(:, 2 + nBB:end);
        fprintf('model B (ball usable, %.0f s, blocks %s): %s R^2 full %.2f; loc alone %.2f [%.2f, %.2f], unique %.2f [%.2f, %.2f]\n', numel(vB) * dt_a, lst(loc_blk), u_name{1}, B.R2(1, 1), ...
            B.alone(1, end), B.CI(1, end, 1), B.CI(1, end, 2), B.uniq(1, end), B.dCI(1, end, 1), B.dCI(1, end, 2));
    else
        fprintf('model B skipped: ball log usable for %.0f s (< %g s)\n', numel(vB) * dt_a, min_loc_s);
    end
end

% cross-correlation of the mean unit with each behaviour (positive lag = dF/F follows the behaviour), and its null band (dF/F shifted)
L = round(corr_max_lag_s / dt_a);  lg = -L:L;
Bv_grid = Bz;  Bv_grid(:, ~(cam_ok & ok_y)) = NaN;
xc_obs = xc_valid(UG(1, :), Bv_grid, lg);
xc_nul = zeros(n_null, nBeh, numel(lg));
for it = 1:n_null
    sh = null_min_shift_s / dt_a + randi(nG - 2 * round(null_min_shift_s / dt_a));
    xc_nul(it, :, :) = reshape(xc_valid(circshift(UG(1, :), sh, 2), Bv_grid, lg), 1, nBeh, []);
end

% kernel: weights of the full model fitted on all analysed bins (per SD of behaviour), descriptive
y1 = Yo(:, 1);  lam1 = pick_lambda_v(y1, Xv, vr, 4, mg, lams);
mu1 = mean(Xv, 1);  w1 = ((Xv - mu1)' * (Xv - mu1) + lam1 * eye(size(Xv, 2))) \ ((Xv - mu1)' * (y1 - mean(y1)));
kern = reshape(w1, nl, nBeh)';

% correlation between ROIs before / after removing the behaviour (held-out residuals)
sg = unique([sig_roi, sel_roi]);  nS = numel(sg);
C_raw = corr(Yo_roi(:, sg));  C_res = corr(Yo_roi(:, sg) - Yh_roi(:, sg));
[pi_, pj_] = find(triu(true(nS), 1));  pr_raw = C_raw(sub2ind([nS nS], pi_, pj_));  pr_res = C_res(sub2ind([nS nS], pi_, pj_));
same_ax = grp_of(sg(pi_)) > 0 & grp_of(sg(pi_)) == grp_of(sg(pj_));
fprintf('\nROI pairs (%d signal-like ROIs, %d pairs): median r %.2f -> %.2f after removing the behaviour; pairs inside one group %.2f -> %.2f; other pairs %.2f -> %.2f\n', nS, numel(pr_raw), ...
    median(pr_raw), median(pr_res), median(pr_raw(same_ax)), median(pr_res(same_ax)), median(pr_raw(~same_ax)), median(pr_res(~same_ax)));
beh_r = corrcoef(Bz(:, cam_ok)', 'Rows', 'pairwise');

%% 7. Figure 1: how much of the activity is explained by behaviour
W = 18;  x0 = 2.7;  wA = W - x0 - 0.5;
hObs = 2.6;  hBeh = 0.8;  gB = 0.12;  hStrip = 0.22;
hA = 0.6 + hObs + 0.1 + hStrip + 0.1 + nBeh * (hBeh + gB);                    % legend, measured / predicted, epoch strip, behaviours
hRow = 6.2;  H = 0.9 + hA + 0.3 + hRow;
fh1 = newfig('Behaviour explains axon activity', W, H);
annotation(fh1, 'textbox', [0, 1 - 0.65 / H, 1, 0.55 / H], 'String', sprintf('%s | held-out variance of the axon dF/F explained by behaviour (%d ROIs, %.0f s analysed)', session, nSel, nV * dt_a), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 9, 'Interpreter', 'none');
% A: example stretch
nW = round(ex_win_s / dt_a);  starts = 1:round(30 / dt_a):nG - nW;  isv = false(1, nG);  isv(vr) = true;
pos_in_vr = zeros(1, nG);  pos_in_vr(vr) = 1:nV;
cand = starts(arrayfun(@(s) mean(isv(s:s + nW - 1)) > 0.98, starts));
if isempty(cand); [~, ib] = max(movsum(double(isv), nW));  cand = max(1, ib - nW + 1); end
loc_r2 = nan(size(cand));
for i = 1:numel(cand)
    ii = pos_in_vr(cand(i):cand(i) + nW - 1);  ii = ii(ii > 0);
    loc_r2(i) = 1 - sum((Yo(ii, 1) - Yh(ii, 1, 1)) .^ 2) / sum((Yo(ii, 1) - mean(Yo(ii, 1))) .^ 2);
end
[~, ic] = min(abs(loc_r2 - median(loc_r2)));  s0 = cand(ic);  ws = s0:s0 + nW - 1;  tw = (0:nW - 1) * dt_a;
yA = H - 0.9 - 0.6;
yA = yA - hObs;  axo = axes('Units', 'centimeters', 'Position', [x0 yA wA hObs]);  hold(axo, 'on');
okw = isv(ws);
yo_w = nan(1, nW);  yh_w = nan(1, nW);  yo_w(okw) = Yo(pos_in_vr(ws(okw)), 1);  yh_w(okw) = Yh(pos_in_vr(ws(okw)), 1, 1);
plot(axo, tw, yo_w, 'Color', col.dff, 'LineWidth', 1);  plot(axo, tw, yh_w, 'Color', col.pred, 'LineWidth', 1.2);
xlim(axo, [0 ex_win_s]);  set(axo, 'XTick', []);  ylabel(axo, 'dF/F (z)');
legend(axo, {'measured, mean of selected ROIs', sprintf('predicted from behaviour (held out), R^2 = %.2f in this stretch', loc_r2(ic))}, 'FontSize', 6.5, 'Location', 'northoutside', 'NumColumns', 2);
letter(axo, 'A');
yA = yA - 0.1 - hStrip;  axs = axes('Units', 'centimeters', 'Position', [x0 yA wA hStrip]);  hold(axs, 'on');          % epoch strip: which stimulus / state the stretch is in
r0_ = 1;
while r0_ <= nW
    e_ = ep_id(ws(r0_));  r1_ = r0_;  while r1_ < nW && ep_id(ws(r1_ + 1)) == e_; r1_ = r1_ + 1; end
    if e_ > 0
        patch(axs, tw([r0_ r1_ r1_ r0_]), [0 0 1 1], ecol(e_, :), 'EdgeColor', 'none', 'FaceAlpha', 0.7);
        if (r1_ - r0_) * dt_a > 15; text(axs, mean(tw([r0_ r1_])), 0.5, strrep(ep_names{e_}, '_', ' '), 'HorizontalAlignment', 'center', 'FontSize', 6, 'Color', 'w'); end
    end
    r0_ = r1_ + 1;
end
xlim(axs, [0 ex_win_s]);  ylim(axs, [0 1]);  axis(axs, 'off');
yA = yA - 0.1;
for j = 1:nBeh
    yA = yA - hBeh;  axb = axes('Units', 'centimeters', 'Position', [x0 yA wA hBeh]);  hold(axb, 'on');
    plot(axb, tw, Bz(j, ws), 'Color', bcol(j, :), 'LineWidth', 0.9);
    yl_ = [min(Bz(j, ws)) max(Bz(j, ws))];  yl_ = yl_ + [-0.15 0.15] * max(diff(yl_), 0.5);                 % each trace autoscaled to the stretch (z-scored over the session)
    xlim(axb, [0 ex_win_s]);  ylim(axb, yl_);  yticks(axb, []);  ylabel(axb, beh_lab.(bname{j}), 'FontSize', 7, 'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    if j < nBeh; set(axb, 'XTick', []); else; xlabel(axb, sprintf('time in the stretch (s); session time %.0f s', tg(ws(1)))); end
    yA = yA - gB;
end
% B: variance partition
yR = 0;  axH = hRow - 3.2;                                                      % the three axes of the second row share bottom (ticks / labels below) and top (titles / legend above)
axB = axes('Units', 'centimeters', 'Position', [x0 yR + 1.3 4.4 axH]);  hold(axB, 'on');
nShow = min(nU, 1 + max_groups_shown);
M = zeros(nShow, nBeh + 2);
for u = 1:nShow
    rf = max(R2(u, 1), 0);  un = max(uniq(u, :), 0);  un = un * min(1, rf / max(sum(un), eps));
    M(u, :) = [un, max(rf - sum(un), 0), 1 - rf];
end
hb = barh(axB, 1:nShow, 100 * M, 'stacked', 'BarWidth', 0.72);
for j = 1:nBeh; set(hb(j), 'FaceColor', bcol(j, :), 'EdgeColor', 'none'); end
set(hb(nBeh + 1), 'FaceColor', [0.72 0.72 0.72], 'EdgeColor', 'none');  set(hb(nBeh + 2), 'FaceColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.5);
set(axB, 'YDir', 'reverse', 'YTick', 1:nShow, 'YTickLabel', u_short(1:nShow), 'TickLabelInterpreter', 'none');  xlim(axB, [0 100]);  ylim(axB, [0.4 nShow + 0.6]);
xlabel(axB, 'variance of dF/F (%)');
for u = 1:nShow; text(axB, 100 * min(R2(u, 1), 1) + 1.5, u, sprintf('%.0f%%', 100 * R2(u, 1)), 'FontSize', 7, 'VerticalAlignment', 'middle', 'Color', [0.2 0.2 0.2]); end
lg_h = [hb(1:nBeh + 1)];  lg_s = [cellfun(@(b) ['unique: ' b], bsh, 'UniformOutput', false), {'shared by several'}];
legend(axB, lg_h, lg_s, 'FontSize', 6, 'Location', 'northoutside', 'NumColumns', 2);
title(axB, 'held-out variance explained (white = not)', 'FontSize', 7, 'FontWeight', 'normal');  letter(axB, 'B');
% C: each behaviour alone
axC = axes('Units', 'centimeters', 'Position', [8.7 yR + 1.3 3.6 axH]);  hold(axC, 'on');
vals = [alone(1, :), R2(1, 1)];  lo = [reshape(CI(1, 1 + nBeh + (1:nBeh), 1), 1, []), CI(1, 1, 1)];  hi = [reshape(CI(1, 1 + nBeh + (1:nBeh), 2), 1, []), CI(1, 1, 2)];
for j = 1:nBeh + 1
    if j <= nBeh; cc = bcol(j, :); else; cc = [0.25 0.25 0.25]; end
    bar(axC, j, vals(j), 0.7, 'FaceColor', cc, 'EdgeColor', 'none');
end
errorbar(axC, 1:nBeh + 1, vals, vals - lo, hi - vals, 'k.', 'CapSize', 3, 'MarkerSize', 1);
yline(axC, prctile(null_R2(:, 1), 95), '--', 'Color', [0.4 0.4 0.4]);  yline(axC, 0, '-');
xticks(axC, 1:nBeh + 1);  xticklabels(axC, [bsh, {'all'}]);  xtickangle(axC, 45);  xlim(axC, [0.4 nBeh + 1.6]);  ylabel(axC, 'held-out R^2');
title(axC, {'each behaviour alone', 'dashed: null (95th pct)'}, 'FontSize', 7, 'FontWeight', 'normal');  letter(axC, 'C');
% D: every ROI
axD = axes('Units', 'centimeters', 'Position', [13.8 yR + 1.3 W - 13.8 - 0.4 axH]);  hold(axD, 'on');
ord = [sel_roi(:); setdiff(sig_roi(:), sel_roi(:)); find(~signal_like)];
ord = ord(:)';
for i = 1:numel(ord)
    r = ord(i);
    if ismember(r, sel_roi); c_ = gcol(mod(grp_of(r) - 1, size(gcol, 1)) + 1, :); plot(axD, i, roi_R2(r), 'o', 'MarkerFaceColor', c_, 'MarkerEdgeColor', 'none', 'MarkerSize', 4.5);
    elseif signal_like(r); plot(axD, i, roi_R2(r), 'o', 'MarkerFaceColor', [0.65 0.65 0.65], 'MarkerEdgeColor', 'none', 'MarkerSize', 4.5);
    else; plot(axD, i, roi_R2(r), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.55 0.55 0.55], 'MarkerSize', 4.5);
    end
    if isfinite(ceil_r2(r)); plot(axD, i + [-0.35 0.35], ceil_r2(r) * [1 1], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1); end
end
yline(axD, 0, '-');  yline(axD, prctile(null_R2(:, 1), 95), '--', 'Color', [0.4 0.4 0.4]);
xlim(axD, [0 numel(ord) + 1]);  ylim(axD, [min(-0.1, min(roi_R2) - 0.02), 1]);  ylabel(axD, 'held-out R^2');  xlabel(axD, 'ROI');
tk = 1:max(1, ceil(numel(ord) / 20)):numel(ord);  xticks(axD, tk);  xticklabels(axD, string(ord(tk)));  xtickangle(axD, 90);  set(axD, 'FontSize', 6);
title(axD, {'every ROI, full model', 'colour = group, open = flat', 'dash = r^2 shared with its group'}, 'FontSize', 6, 'FontWeight', 'normal');  letter(axD, 'D');
if save_figs; savef(fh1, '01_behaviour_explains_activity'); end

%% 8. Figure 2: details and controls
fh2 = newfig('Behaviour model: details and controls', 18, 17);
tl = tiledlayout(fh2, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s | details and controls', session), 'FontSize', 9, 'FontWeight', 'normal', 'Interpreter', 'none');
% a: cross-correlation
a1 = nexttile(tl);  hold(a1, 'on');
band = prctile(reshape(xc_nul, [], 1), [2.5 97.5]);                                % null of every behaviour and lag pooled
fill(a1, [lg * dt_a, fliplr(lg * dt_a)], [band(2) * ones(1, numel(lg)), band(1) * ones(1, numel(lg))], [0.85 0.85 0.85], 'EdgeColor', 'none');
hl = gobjects(1, nBeh);
for j = 1:nBeh; hl(j) = plot(a1, lg * dt_a, xc_obs(j, :), 'Color', bcol(j, :), 'LineWidth', 1.3); end
xline(a1, 0, ':');  yline(a1, 0, ':');  xlabel(a1, 'lag (s), > 0: dF/F follows behaviour');  ylabel(a1, 'r');  xlim(a1, lg([1 end]) * dt_a);
yl1 = ylim(a1);  ylim(a1, [yl1(1), yl1(2) + 0.35 * diff(yl1)]);                    % headroom for the legend
legend(a1, hl, bsh, 'FontSize', 6, 'Location', 'northeast');
title(a1, {'cross-correlation, mean of selected ROIs', 'grey: 95% of the shifted-dF/F null'}, 'FontSize', 7, 'FontWeight', 'normal');  letter(a1, 'a');
% b: predicted vs observed
a2 = nexttile(tl);
yh1_ = Yh(:, 1, 1);
lim_ = prctile([yh1_; Yo(:, 1)], [0.5 99.5]);  ed_ = linspace(lim_(1), lim_(2), 41);  ct_ = (ed_(1:end-1) + ed_(2:end)) / 2;
Nh_ = histcounts2(yh1_, Yo(:, 1), ed_, ed_);
imagesc(a2, ct_, ct_, log10(1 + Nh_'));  set(a2, 'YDir', 'normal');  colormap(a2, flipud(gray(256)));  hold(a2, 'on');
plot(a2, lim_, lim_, 'r:');  xlim(a2, lim_);  ylim(a2, lim_);  axis(a2, 'square');
xlabel(a2, 'predicted, held out (z)');  ylabel(a2, 'measured (z)');
title(a2, {'held-out prediction, mean unit', sprintf('R^2 = %.2f, r = %.2f', R2(1, 1), corr(yh1_, Yo(:, 1)))}, 'FontSize', 7, 'FontWeight', 'normal');  letter(a2, 'b');
% c: kernels
a3 = nexttile(tl);  hold(a3, 'on');
for j = 1:nBeh; plot(a3, enc_lags_s, kern(j, :), '-o', 'Color', bcol(j, :), 'MarkerSize', 3, 'MarkerFaceColor', bcol(j, :)); end
xline(a3, 0, ':');  yline(a3, 0, ':');  xlabel(a3, 'behaviour leads dF/F by (s)');  ylabel(a3, 'weight (z dF/F per SD)');
title(a3, 'model weights (fit on all analysed time)', 'FontSize', 7, 'FontWeight', 'normal');  letter(a3, 'c');
% d: per epoch
a4 = nexttile(tl);  hold(a4, 'on');
ev = find(isfinite(ep_R2))';
for i = 1:numel(ev)
    bar(a4, i, ep_R2(ev(i)), 0.7, 'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none');
    errorbar(a4, i, ep_R2(ev(i)), ep_R2(ev(i)) - ep_ci(ev(i), 1), ep_ci(ev(i), 2) - ep_R2(ev(i)), 'k.', 'CapSize', 3, 'MarkerSize', 1);
end
yline(a4, 0, '-');  yline(a4, R2(1, 1), '--', 'Color', col.pred);  xlim(a4, [0.4 numel(ev) + 0.6]);
xticks(a4, 1:numel(ev));  xticklabels(a4, arrayfun(@(e) sprintf('%s (%.0f s)', ep_short(ep_names{e}), ep_sec(e)), ev, 'UniformOutput', false));  xtickangle(a4, 45);  set(a4, 'FontSize', 6.5);
ylabel(a4, 'held-out R^2');  title(a4, {'per epoch, pooled model', 'dashed: whole session'}, 'FontSize', 7, 'FontWeight', 'normal');  letter(a4, 'd');
% e: shared signal across ROIs
a5 = nexttile(tl);  hold(a5, 'on');
scatter(a5, pr_raw(~same_ax), pr_res(~same_ax), 8, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.5);
scatter(a5, pr_raw(same_ax), pr_res(same_ax), 14, [0 0 0], 'filled');
lo5 = max(-1, floor(10 * min([0; pr_raw; pr_res])) / 10);
plot(a5, [lo5 1], [lo5 1], 'k:');  yline(a5, 0, ':');  xline(a5, 0, ':');  xlim(a5, [lo5 1]);  ylim(a5, [lo5 1]);  axis(a5, 'square');
xlabel(a5, 'ROI-ROI r of dF/F');  ylabel(a5, 'r after removing behaviour');
title(a5, {sprintf('median r %.2f -> %.2f', median(pr_raw), median(pr_res)), 'black = same group'}, 'FontSize', 7, 'FontWeight', 'normal');  letter(a5, 'e');
% f: model B
a6 = nexttile(tl);  hold(a6, 'on');
if B.ok
    nBB = numel(B.names);  vals = [B.alone(1, :), B.R2(1, 1)];  lo = [reshape(B.CI(1, 1 + nBB + (1:nBB), 1), 1, []), B.CI(1, 1, 1)];  hi = [reshape(B.CI(1, 1 + nBB + (1:nBB), 2), 1, []), B.CI(1, 1, 2)];
    un = B.uniq(1, :);  uc = reshape(B.dCI(1, :, :), nBB, 2);
    for j = 1:nBB
        cc = beh_cmap.(B.names{j});
        bar(a6, j - 0.18, vals(j), 0.34, 'FaceColor', cc, 'EdgeColor', 'none');
        bar(a6, j + 0.18, un(j), 0.34, 'FaceColor', 'w', 'EdgeColor', cc, 'LineWidth', 1);
        errorbar(a6, j - 0.18, vals(j), vals(j) - lo(j), hi(j) - vals(j), 'k.', 'CapSize', 2, 'MarkerSize', 1);
        errorbar(a6, j + 0.18, un(j), un(j) - uc(j, 1), uc(j, 2) - un(j), 'k.', 'CapSize', 2, 'MarkerSize', 1);
    end
    bar(a6, nBB + 1, vals(end), 0.5, 'FaceColor', [0.25 0.25 0.25], 'EdgeColor', 'none');
    errorbar(a6, nBB + 1, vals(end), vals(end) - lo(end), hi(end) - vals(end), 'k.', 'CapSize', 2, 'MarkerSize', 1);
    yline(a6, 0, '-');  xticks(a6, 1:nBB + 1);  xticklabels(a6, [cellfun(@(b) beh_short.(b), B.names, 'UniformOutput', false), {'all'}]);  xtickangle(a6, 45);  xlim(a6, [0.4 nBB + 1.6]);
    ylabel(a6, 'held-out R^2');  title(a6, {sprintf('ball log only (%.0f s) + locomotion', numel(B.vr) * dt_a), 'filled: alone, open: unique'}, 'FontSize', 7, 'FontWeight', 'normal');
elseif in_main
    axis(a6, 'off');  text(a6, 0.5, 0.5, 'locomotion is part of the main model', 'HorizontalAlignment', 'center', 'FontSize', 7);
else
    axis(a6, 'off');  text(a6, 0.5, 0.5, {'locomotion not available:', 'ball log unusable or too short', 'in this session'}, 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
end
letter(a6, 'f');
% g: brain motion control
a7 = nexttile(tl);  hold(a7, 'on');
if all(isfinite(R2_mot))
    bar(a7, 1:3, R2_mot, 0.6, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
    for i = 1:3; text(a7, i, R2_mot(i), sprintf('%.2f', R2_mot(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7); end
    yline(a7, 0, '-');  xticks(a7, 1:3);  xticklabels(a7, {'behaviour', '+ brain motion', 'brain motion alone'});  xtickangle(a7, 30);  xlim(a7, [0.4 3.6]);
    ylabel(a7, 'held-out R^2');  title(a7, 'control: suite2p rigid shift', 'FontSize', 7, 'FontWeight', 'normal');
else
    axis(a7, 'off');  text(a7, 0.5, 0.5, 'no suite2p motion for this session', 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
end
letter(a7, 'g');
% h: null
a8 = nexttile(tl);  hold(a8, 'on');
histogram(a8, null_R2(:, 1), 'BinWidth', 0.005, 'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none');
xline(a8, R2(1, 1), 'Color', col.pred, 'LineWidth', 2);
xlim(a8, [-0.05, max(0.2, 1.15 * R2(1, 1))]);
xlabel(a8, 'held-out R^2');  ylabel(a8, 'draws');  title(a8, {sprintf('null: dF/F shifted by >= %g s', null_min_shift_s), 'grey: shifted, red: observed'}, 'FontSize', 7, 'FontWeight', 'normal');  letter(a8, 'h');
% i: numbers
a9 = nexttile(tl);  axis(a9, 'off');
txt = {sprintf('%d ROIs in %d group(s)', nSel, ng), sprintf('%.0f s analysed (%.0f%% of session)', nV * dt_a, 100 * nV / nG), ...
    sprintf('predictors: %s', strjoin(bsh, ', ')), sprintf('lags %g to %g s', enc_lags_s(1), enc_lags_s(end)), sprintf('%d folds, %g s purged', enc_nfold, enc_margin_s), ...
    sprintf('R^2 mean unit %.2f [%.2f, %.2f]', R2(1, 1), CI(1, 1, 1), CI(1, 1, 2)), sprintf('null %.3f +- %.3f', mean(null_R2(:, 1)), std(null_R2(:, 1))), ...
    'behaviour-behaviour |r|:', sprintf('  median %.2f, max %.2f', median(abs(beh_r(triu(true(nBeh), 1)))), max(abs(beh_r(triu(true(nBeh), 1)))))};
text(a9, 0, 1, txt, 'VerticalAlignment', 'top', 'FontSize', 7, 'Interpreter', 'none');
if save_figs
    savef(fh2, '02_details_and_controls');
    writetable(enc_tbl, fullfile(fig_dir, 'encoding.csv'));  writetable(roi_tbl, fullfile(fig_dir, 'roi_R2.csv'));  writetable(epoch_tbl, fullfile(fig_dir, 'epoch_R2.csv'));
end
fprintf('\naxon_analysis_behaviour_V1 done: %s, %d ROIs %s, full model R^2 (mean unit) %.2f [%.2f, %.2f]\n', session, nSel, lst(sel_roi), R2(1, 1), CI(1, 1, 1), CI(1, 1, 2));

%% local functions
function panel_letter(ax, s)
% bold panel letter 0.9 cm left of and 0.3 cm above the top-left corner of the axes
u = ax.Units;  ax.Units = 'centimeters';  p = ax.Position;  ax.Units = u;
text(ax, -0.9 / p(3), 1 + 0.3 / p(4), s, 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Clipping', 'off');
end

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
n = size(X, 1);
C = real(ifft(fft(X) .* conj(fft(y)))) / (n - 1);
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
% predictors [n x nB*nl]: behaviour j at t - lag (positive lag = behaviour earlier), NaN where t - lag is outside the data; columns behaviour-major
[nB, n] = size(Bc);  X = nan(n, nB * numel(lagb));  c = 0;
for j = 1:nB
    for l = 1:numel(lagb)
        c = c + 1;  src = (1:n) - lagb(l);  ok = src >= 1 & src <= n;
        X(ok, c) = Bc(j, src(ok))';
    end
end
end

function [Z, mu, sd] = zscore_cols(X)
mu = mean(X, 1, 'omitnan');  sd = std(X, 0, 1, 'omitnan');  sd(sd == 0 | isnan(sd)) = 1;  Z = (X - mu) ./ sd;
end

function [yhat, fold] = cv_ridge_v(y, X, tidx, nfold, margin, lams)
% out-of-fold ridge predictions. Rows are the analysed time bins in time order (tidx = their bin numbers): nfold contiguous folds of equal size, the training set
% leaves out every bin within `margin` bins of the test fold, lambda from an inner blocked 4-fold CV. X columns are already standardised.
n = numel(y);  e = round(linspace(0, n, nfold + 1));  yhat = nan(n, 1);  fold = zeros(n, 1);  p = size(X, 2);
for f = 1:nfold
    te = e(f) + 1:e(f + 1);  fold(te) = f;
    trn = find(tidx < tidx(te(1)) - margin | tidx > tidx(te(end)) + margin);
    lam = pick_lambda_v(y(trn), X(trn, :), tidx(trn), 4, margin, lams);
    mu = mean(X(trn, :), 1);  Xc = X(trn, :) - mu;  ym = mean(y(trn));
    w = (Xc' * Xc + lam * eye(p)) \ (Xc' * (y(trn) - ym));
    yhat(te) = (X(te, :) - mu) * w + ym;
end
end

function lam = pick_lambda_v(y, X, tidx, k, margin, lams)
n = numel(y);  e = round(linspace(0, n, k + 1));  err = zeros(1, numel(lams));
for f = 1:k
    te = e(f) + 1:e(f + 1);  trn = find(tidx < tidx(te(1)) - margin | tidx > tidx(te(end)) + margin);
    mu = mean(X(trn, :), 1);  Xc = X(trn, :) - mu;  ym = mean(y(trn));
    [V, D] = eig(Xc' * Xc);  d = diag(D);  Vb = V' * (Xc' * (y(trn) - ym));
    Xte = (X(te, :) - mu) * V;
    for l = 1:numel(lams); err(l) = err(l) + sum((y(te) - (Xte * (Vb ./ (d + lams(l))) + ym)) .^ 2); end
end
[~, il] = min(err);  lam = lams(il);
end

function out = r2_boot(y, Yh, blk, nb)
% held-out R^2 of each column of Yh and its 95% interval (bootstrap over blocks of blk consecutive rows; the draws are returned for differences between models)
n = numel(y);  nbk = ceil(n / blk);
out.r2 = 1 - sum((y - Yh) .^ 2, 1) / sum((y - mean(y)) .^ 2);
bt = zeros(nb, size(Yh, 2));
for b = 1:nb
    pick = randi(nbk, 1, nbk);  idx = reshape((pick - 1) * blk + (1:blk)', [], 1);  idx = idx(idx <= n);
    yy = y(idx);  bt(b, :) = 1 - sum((yy - Yh(idx, :)) .^ 2, 1) / sum((yy - mean(yy)) .^ 2);
end
out.ci = prctile(bt, [2.5 97.5], 1);  out.boot = bt;
end

function R = xc_valid(x, Y, lags)
% Pearson r between x(t) [1 x n] and each row of Y [nB x n] at t - lag, over the bins where both exist (positive lag = x follows Y)
nB = size(Y, 1);  n = numel(x);  R = nan(nB, numel(lags));
for li = 1:numel(lags)
    l = lags(li);  src = (1:n) - l;  ok = src >= 1 & src <= n;
    xi = x(ok);
    for j = 1:nB
        yj = Y(j, src(ok));  g = isfinite(yj) & isfinite(xi);
        if nnz(g) > 10; R(j, li) = corr(xi(g)', yj(g)'); end
    end
end
end
