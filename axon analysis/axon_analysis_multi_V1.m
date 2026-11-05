%% axon_analysis_multi_V1
% Exploratory multi-recording average: the same stimulus in several aggregated axon sessions (V3-type sensory sessions and task sessions), one figure.
% Per session: ROIs -> signal-like ROIs (V3 criterion) -> ROI groups (one axon = several ROIs, r > grp_r_thr) -> mean of the group means = one trace per
% session. That trace and the behaviour trace (beh_name, z-scored over the session) are aligned to the onsets of every stimulus type present, baseline-
% subtracted, and averaged over trials, then over sessions. ROIs are NOT matched between sessions: the average is over fields of view, not over axons.
% Test: the evoked amplitude (mean in awin minus mean in bwin) against onsets shifted at random inside their own block (same trial spacing, wrong phase);
% for several sessions the null is the mean over sessions of independent shifts. A response that goes with an equally large behaviour change is arousal.
% Figure: row 1/2 = trial-averaged dF/F / behaviour (thin = session +- SEM over its trials, black = mean over sessions), row 3/4 = evoked amplitude per
% session (dots +- SEM over trials, black bar = mean over sessions, grey = 95% of the null). Not done: motion-coupled ROI exclusion, hit/miss split of the task.

%% 0. Settings
agg_dir   = 'Z:\Group Members\Rima\Aggregated';
repo_root = 'C:\Users\rbondarenko\projects\2P_Scripts';
sess = {'AnimalRB24_261021_1630', 'signal';      % {session, ROIs}: 'signal' (signal-like) | 'all' | a list such as [21 23 24 25 27] (V3 'good_roi' of that session)
        'AnimalRB24_261023_1538', 'signal';      % (AnimalRB24_261023_1527 is left out: interrupted after 365 s, blocks 3-10 have no frame times)
        'AnimalRB24_261030_1744', 'signal'};
min_trials = 5;           % a session contributes to a stimulus type only with at least this many trials
beh_name  = 'WhiskIndex';   % behaviour shown next to dF/F: WhiskIndex | PupilArea | FacePadME
save_figs = false;          % true: PNG (300 dpi) to fig_dir
fig_dir   = fullfile(repo_root, 'axon analysis', 'plots', 'multi_V1');

% window (s from onset) shown, baseline, evoked amplitude. As V3, except looming amplitude 0-4 s (disc 0-1 s) and the task windows of task_analysis_V1
keys = {'flash', 'fff', 'loom_expand', 'loom_recede', 'bar', 'square'};
labs = {'local/global flashes', 'full-field flash', 'looming, expanding', 'looming, receding', 'moving bar', 'square + reward'};
win  = struct('flash', [-1 4], 'fff', [-2 32], 'loom_expand', [-2 8], 'loom_recede', [-2 8], 'bar', [-1 6], 'square', [-5 12]);
bwin = struct('flash', [-1 0], 'fff', [0 2],   'loom_expand', [-2 0], 'loom_recede', [-2 0], 'bar', [-1 0], 'square', [-3 -0.5]);
awin = struct('flash', [0 2],  'fff', [2 32],  'loom_expand', [0 4],  'loom_recede', [0 4],  'bar', [0.5 4], 'square', [0 3]);
fs_plot = 30;             % samples/s of the common time axis
slow_ratio_thr = 4;       % signal-like: var(0.5 s smoothed dFF) >= thr x white-noise expectation (V3 section 5)
neuropil_factor = 0.5;  smooth_s = 0.5;  max_interp_s = 0.5;  grp_r_thr = 0.8;
n_null = 500;  rng(1);
addpath(repo_root, fullfile(repo_root, 'analysis_functions'));
sc = lines(size(sess, 1));  slab = regexprep(sess(:, 1), '^.*_\d{2}(\d{4})_(\d{4})$', '$1-$2');

%% 1. Per session: ROI traces, behaviour, trials of every stimulus type
nS = size(sess, 1);  R = cell(nS, 1);
for i = 1:nS
    f = fullfile(agg_dir, [sess{i, 1} '_preprocessed.mat']);  m = whos('-file', f);  inf_ = {m.name};
    want = {'CaData', 'TimeCa', 'Stimuli', 'Triggers', 'Behav'};  D = load(f, want{ismember(want, inf_)});
    ca = D.CaData(find(~cellfun(@isempty, {D.CaData.Ca_dFF}), 1));
    dFF = double(ca.Ca_dFF);  Fc = double(ca.Ca_F) - neuropil_factor * double(ca.Ca_F_neuropil);  tCa = double(D.TimeCa(1, :));
    [nROI, nT] = size(dFF);  fs = 1 / median(diff(tCa));  w_sm = max(1, round(smooth_s * fs));

    % signal-like ROIs, same criterion as V3 section 5 (slow variance vs white-noise expectation, upward events, Fc > 0)
    dsm = movmean(dFF, w_sm, 2);  dev = dsm - median(dsm, 2);  neg = dev < 0;
    sig_neg = max(sqrt(sum(dev .^ 2 .* neg, 2) ./ max(sum(neg, 2), 1)), eps('single'));
    dd = diff(dFF, 1, 2);  hf = 1.4826 * median(abs(dd - median(dd, 2)), 2) / sqrt(2);
    ev_pos = sum(diff([false(nROI, 1) dev > 3 * sig_neg], 1, 2) == 1, 2);  ev_neg = sum(diff([false(nROI, 1) dev < -3 * sig_neg], 1, 2) == 1, 2);
    signal_like = mean(Fc <= 0, 2) <= 0.01 & var(dsm, 0, 2) ./ (hf .^ 2 / w_sm) >= slow_ratio_thr & ev_pos >= 2 * max(ev_neg, 1);
    if ischar(sess{i, 2}) && strcmp(sess{i, 2}, 'all'); sel = 1:nROI; elseif ischar(sess{i, 2}); sel = find(signal_like)'; else; sel = sess{i, 2}(:)'; end
    if isempty(sel); warning('%s: no ROI selected, session skipped', sess{i, 1}); continue; end
    grp = roi_groups(dsm(sel, :), grp_r_thr);
    x_ca = mean(cell2mat(cellfun(@(g) mean(dFF(sel(g), :), 1), grp(:), 'UniformOutput', false)), 1);        % mean of the group means

    x_wh = nan(1, nT);                                                                % behaviour on the imaging grid (NaN where the camera is farther than max_interp_s)
    if isfield(D, 'Behav') && isfield(D.Behav, beh_name) && numel(D.Behav.(beh_name)) == numel(D.Triggers.TimeCamera)
        b = double(D.Behav.(beh_name)(:))';  tcam = double(D.Triggers.TimeCamera(:))';  ok = isfinite(b);  tv = tcam(ok);
        x_wh = interp1(tv, b(ok), tCa, 'linear', NaN);  ix = interp1(tv, 1:numel(tv), tCa, 'nearest', 'extrap');  x_wh(abs(tv(ix) - tCa) > max_interp_s) = NaN;
        x_wh = (x_wh - mean(x_wh, 'omitnan')) / std(x_wh, 'omitnan');
    end

    % onsets of every stimulus type (V3 section 8) and the time span of the block each belongs to (for the shifted-onset null)
    on = struct();  sp = struct();
    for k = 1:numel(D.Stimuli)
        s = D.Stimuli(k);  span = [];
        switch s.type
            case 'sparse_local_global_flashes'; key = 'flash'; o = get_flash_onsets(s);
            case 'full_field_flash';            key = 'fff';   o = first_frames(s);
            case 'moving_bar';                  key = 'bar';   o = first_frames(s);
            case 'looming'
                key = 'loom_expand';  o = first_frames(s);
                if isfield(s, 'specParams') && isstruct(s.specParams) && isfield(s.specParams, 'loom_maxexpand') && ~isempty(s.specParams.loom_maxexpand) && s.specParams.loom_maxexpand(1) < 0; key = 'loom_recede'; end
            case 'training'                                                           % task: onset = first valid square frame (rows wait1 | stim | wait2, as task_analysis_V1)
                key = 'square';  Fi = s.Frameinfo;  n_w1 = size(Fi.frame_count_wait1, 1);  n_st = size(Fi.frame_count_stim, 1);
                FT = s.TimeStimulusFrame;  FV = s.FrameValid & isfinite(FT);  o = nan(size(FT, 2), 1);
                for t = 1:size(FT, 2); j = find(FV(n_w1 + (1:n_st), t), 1); if ~isempty(j); o(t) = FT(n_w1 + j, t); end; end
                o = o(isfinite(o));  span = [min(FT(FV)) max(FT(FV))];
            otherwise; continue
        end
        o = o(isfinite(o));  if isempty(o); continue; end                              % interrupted blocks have NaN frame times
        if isempty(span); span = [min(s.TimeStimulusFrame(:)) max(s.TimeStimulusFrame(:))]; end
        if ~isfield(on, key); on.(key) = zeros(0, 1);  sp.(key) = zeros(0, 2); end
        on.(key) = [on.(key); o(:)];  sp.(key) = [sp.(key); repmat(span, numel(o), 1)];
    end
    for key = fieldnames(on)'; if numel(on.(key{1})) < min_trials; on = rmfield(on, key{1}); end; end

    r_ = struct();  cum = {cum_of(x_ca), cum_of(x_wh)};  X = {x_ca, x_wh};
    for q = 1:numel(keys)
        key = keys{q};  if ~isfield(on, key); continue; end
        w = win.(key);  n_com = round(diff(w) * fs_plot) + 1;  tcom = linspace(w(1), w(2), n_com);  rq = struct('tcom', tcom, 'n', numel(on.(key)), 'C', {cell(1, 2)}, 'amp', {cell(1, 2)}, 'nul', {cell(1, 2)});
        for c = 1:2
            C = reshape(align_trials_to_onsets(on.(key), X{c}, tCa, -w(1), w(2), n_com), n_com, []);                      % [time x trials]
            rq.C{c} = C - mean(C(tcom >= bwin.(key)(1) & tcom < bwin.(key)(2), :), 1, 'omitnan');
            rq.amp{c} = amp_at(cum{c}, tCa, on.(key), bwin.(key), awin.(key));
            rq.nul{c} = arrayfun(@(b) mean(amp_at(cum{c}, tCa, shift_in_blocks(on.(key), sp.(key)), bwin.(key), awin.(key)), 'omitnan'), 1:n_null)';
        end
        r_.(key) = rq;
    end
    R{i} = r_;
    fprintf('%s: %d/%d ROIs signal-like, %d selected [%s] in %d group(s); stimuli: %s\n', slab{i}, sum(signal_like), nROI, numel(sel), strtrim(sprintf('%d ', sel)), numel(grp), ...
        strjoin(cellfun(@(k) sprintf('%s x%d', k, r_.(k).n), fieldnames(r_)', 'UniformOutput', false), ', '));
end
present = keys(cellfun(@(k) any(cellfun(@(r) isfield(r, k), R)), keys));  nP = numel(present);

%% 2. Average over sessions: statistics
rows = cell(0, 6);  z_ = @(v) (v - mean(v, 'omitnan')) / std(v, 'omitnan');
catv = @(c) vertcat(c{:});  cath = @(c) horzcat(c{:});
nulm = @(key, c, si) mean(cath(cellfun(@(r) r.(key).nul{c}, R(si), 'UniformOutput', false)), 2);      % null of the mean over sessions [n_null x 1]
for q = 1:nP
    key = present{q};  si = find(cellfun(@(r) isfield(r, key), R))';  st = cell(1, 2);
    for c = 1:2
        A = cellfun(@(r) mean(r.(key).amp{c}, 'omitnan'), R(si))';  nul = nulm(key, c, si);
        st{c} = sprintf('%+.3f (null %+.3f..%+.3f), p=%.3f', mean(A, 'omitnan'), prctile(nul, 2.5), prctile(nul, 97.5), (1 + sum(abs(nul) >= abs(mean(A, 'omitnan')))) / (n_null + 1));
    end
    zc = catv(cellfun(@(r) z_(r.(key).amp{1}), R(si), 'UniformOutput', false));  zw = catv(cellfun(@(r) z_(r.(key).amp{2}), R(si), 'UniformOutput', false));
    rows(end+1, :) = {labs{strcmp(keys, key)}, numel(si), sum(cellfun(@(r) r.(key).n, R(si))), st{1}, st{2}, corr(zc, zw, 'Type', 'Spearman', 'Rows', 'complete')}; %#ok<SAGROW>
end
tbl = cell2table(rows, 'VariableNames', {'stimulus', 'sessions', 'trials', 'dFF_amp', ['beh_amp_' beh_name], 'rho_trials_dFF_vs_beh'});
disp(tbl);

%% 3. Figure
FW = 3.7 * nP + 1.6;  FH = 15;
fh = figure('Name', 'multi-recording average', 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 FW FH], 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 8, ...
    'DefaultAxesBox', 'off', 'DefaultAxesTickDir', 'out', 'DefaultAxesLineWidth', 0.75, 'DefaultLineLineWidth', 0.75, 'DefaultLegendBox', 'off');
tl = tiledlayout(fh, 4, nP, 'TileSpacing', 'compact', 'Padding', 'compact');
ylab = {'\DeltadF/F', ['\Delta ' beh_name ' (SD)'], 'dF/F amplitude', [beh_name ' amplitude (SD)']};
for q = 1:nP
    key = present{q};  si = find(cellfun(@(r) isfield(r, key), R))';  tcom = R{si(1)}.(key).tcom;
    for c = 1:2                                                                        % traces
        ax = nexttile(tl, (c - 1) * nP + q);  hold(ax, 'on');  M = nan(numel(tcom), numel(si));  if q == nP && c == 1; ax_leg = ax; end
        for j = 1:numel(si)
            C = R{si(j)}.(key).C{c};  mu = mean(C, 2, 'omitnan');  se = std(C, 0, 2, 'omitnan') / sqrt(size(C, 2));  M(:, j) = mu;  ok = isfinite(mu) & isfinite(se);
            fill(ax, [tcom(ok) fliplr(tcom(ok))], [(mu(ok) + se(ok))' fliplr((mu(ok) - se(ok))')], sc(si(j), :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');  plot(ax, tcom, mu, 'Color', sc(si(j), :), 'LineWidth', 0.75);
        end
        if numel(si) > 1; plot(ax, tcom, mean(M, 2, 'omitnan'), 'k', 'LineWidth', 1.6); end
        xline(ax, 0, ':');  yline(ax, 0, ':');  xlim(ax, tcom([1 end]));
        if q == 1; ylabel(ax, ylab{c}); end
        if c == 1; title(ax, {labs{strcmp(keys, key)}, sprintf('%d session%s, %d trials', numel(si), repmat('s', 1, numel(si) > 1), sum(cellfun(@(r) r.(key).n, R(si))))}, 'FontSize', 8, 'FontWeight', 'normal'); else; xlabel(ax, 'time from onset (s)'); end
    end
    for c = 1:2                                                                        % evoked amplitude vs null
        ax = nexttile(tl, (c + 1) * nP + q);  hold(ax, 'on');
        nul = nulm(key, c, si);  nb = numel(si);
        patch(ax, [0.4 nb + 0.6 nb + 0.6 0.4], prctile(nul, [2.5 2.5 97.5 97.5]), [0.8 0.8 0.8], 'EdgeColor', 'none');
        A = zeros(nb, 1);
        for j = 1:nb
            a = R{si(j)}.(key).amp{c};  A(j) = mean(a, 'omitnan');
            errorbar(ax, j, A(j), std(a, 'omitnan') / sqrt(sum(isfinite(a))), 'o', 'Color', sc(si(j), :), 'MarkerFaceColor', sc(si(j), :), 'MarkerSize', 4, 'CapSize', 3);
        end
        plot(ax, [0.5 nb + 0.5], mean(A, 'omitnan') * [1 1], 'k-', 'LineWidth', 1.6);  yline(ax, 0, ':');  xlim(ax, [0.4 nb + 0.6]);
        xticks(ax, 1:nb);  xticklabels(ax, slab(si));  xtickangle(ax, 90);
        title(ax, sprintf('p = %.3f', (1 + sum(abs(nul) >= abs(mean(A, 'omitnan')))) / (n_null + 1)), 'FontSize', 7, 'FontWeight', 'normal');
        if q == 1; ylabel(ax, ylab{c + 2}); end
    end
end
h = arrayfun(@(i) plot(ax_leg, nan, nan, '-', 'Color', sc(i, :), 'LineWidth', 1.5), 1:nS);                                 % legend in the last panel of row 1 (empty top left)
legend(ax_leg, h, slab, 'FontSize', 6, 'Location', 'northwest');
title(tl, sprintf('multi-recording average | ROIs of a session averaged (group means) | grey = 95%% of onsets shifted inside the block (n = %d)', n_null), 'FontSize', 8, 'FontWeight', 'normal');
if save_figs; if ~isfolder(fig_dir); mkdir(fig_dir); end; exportgraphics(fh, fullfile(fig_dir, 'multi_average.png'), 'Resolution', 300); end

%% local functions
function o = first_frames(s)
% trial onsets of full-field flash / looming / moving bar blocks: first frame of every trial (V3 section 8)
sfc = s.Frameinfo.subtrial_frame_count(~isnan(s.Frameinfo.frame_count));
o = s.TimeStimulusFrame(find(sfc == 1));
o = o(:);
end

function grp = roi_groups(ds, thr)
% cluster the rows of ds on their correlation (average linkage, cut at 1 - thr); cell array of row indices per group
if size(ds, 1) < 2; grp = {1}; return; end
g = cluster(linkage(pdist(double(ds), 'correlation'), 'average'), 'cutoff', 1 - thr, 'criterion', 'distance');
grp = arrayfun(@(k) find(g == k)', 1:max(g), 'UniformOutput', false);
end

function c = cum_of(x)
% cumulative sum / count of the finite samples of x, with a leading 0 (window means in amp_at)
ok = isfinite(x);  x(~ok) = 0;
c = struct('cs', [0 cumsum(x)], 'cn', [0 cumsum(double(ok))]);
end

function a = amp_at(c, tCa, ons, bw, aw)
% mean of the signal in aw minus its mean in bw (s from each onset); NaN where a window leaves the recording or holds no finite sample
a = win_mean(c, tCa, ons, aw) - win_mean(c, tCa, ons, bw);
a = a(:);
end

function m = win_mean(c, tCa, ons, w)
i1 = interp1(tCa, 1:numel(tCa), ons + w(1), 'nearest', NaN);  i2 = interp1(tCa, 1:numel(tCa), ons + w(2), 'nearest', NaN);
bad = isnan(i1) | isnan(i2);  i1(bad) = 1;  i2(bad) = 1;
m = (c.cs(i2 + 1) - c.cs(i1)) ./ (c.cn(i2 + 1) - c.cn(i1));
m(bad) = NaN;
end

function o2 = shift_in_blocks(ons, sp)
% circularly shift the onsets inside their own block span sp [n x 2] by a random offset (keeps the trial spacing, destroys the timing against the data)
o2 = ons;  [ub, ~, g] = unique(sp, 'rows');
for k = 1:size(ub, 1)
    ii = g == k;  dur = ub(k, 2) - ub(k, 1);  d = (0.2 + 0.6 * rand) * dur;
    o2(ii) = ub(k, 1) + mod(ons(ii) - ub(k, 1) + d, dur);
end
end
