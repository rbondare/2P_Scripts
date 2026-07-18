%% Looming: full-window individual-cell analysis -- 260312 yohimbine pair
%
% Shows what happens THROUGHOUT the loom at the individual-cell level:
%   - Heatmaps (trial-median, row z-scored): pre-loom activity, during-loom
%     response, post-loom inhibition -- sorted by drug pre-loom activity so
%     cells most active before the loom in the drug condition rise to the top,
%     with the SAME row order kept in the baseline panel for direct comparison.
%   - Per-trial spaghetti: each of the 5 trials as a thin line + bold mean,
%     showing how consistent the pattern is trial to trial.
%   - Population mean ± SEM overlay: clearest quantitative comparison of
%     the three temporal phases (pre / during / post) between conditions.
%
% Window is [-2, +10s] -- standard extract_looming_trials only goes to +5s
% (stimulus is 7s, so that never even shows the loom offset). Here we use
% align_trials_to_onsets directly so we can see the full loom and what
% happens afterward.

clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

%% ====================== CONFIG ======================
baseline_file   = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat";
drug_file       = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat";
roi_match_file  = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0312_V1.mat";
selected_plane_idx = 1;

pre_window_sec  = 2;    % s before loom onset
post_window_sec = 10;   % s after loom onset  (loom offset at +7s, post-loom visible +7..+10s)
n_com           = 400;

c_base = [0.5  0.5  0.5];       % grey  – baseline
c_drug = [0.40 0.65 0.40];      % green – drug (yohimbine)

outdir = fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', ...
    'AnimalRB19_260312_1115_vs_AnimalRB19_260312_1243');
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ====================== LOAD ======================
fprintf('Loading baseline...\n');  B = load(baseline_file);
fprintf('Loading drug...\n');      D = load(drug_file);

% Plane selection
centZ_b = B.CaData(1).Ca_centroid_voxel(:,3); up_b = unique(centZ_b);
sel_b   = find(centZ_b == up_b(selected_plane_idx));
dff_base    = B.CaData(1).Ca_dFF(sel_b, :);
TimeCa_base = B.TimeCa(selected_plane_idx, :);

centZ_d = D.CaData(1).Ca_centroid_voxel(:,3); up_d = unique(centZ_d);
sel_d   = find(centZ_d == up_d(selected_plane_idx));
dff_drug    = D.CaData(1).Ca_dFF(sel_d, :);
TimeCa_drug = D.TimeCa(selected_plane_idx, :);

% Matched ROIs. allSessionMapping stores per-plane-LOCAL suite2p ROI indices
% (1..cellCount for that plane), not global indices into the multi-plane
% aggregated array -- confirmed against roiMatchPub.m's generation code,
% which assigns iCell (1..cellCount within a single session) directly into
% the mapping. These local indices ARE the row indices into sel_b/sel_d.
M = load(roi_match_file);
bl = M.roiMatchData.allSessionMapping(:,1);
dl = M.roiMatchData.allSessionMapping(:,2);
valid = (bl>=1 & bl<=numel(sel_b)) & (dl>=1 & dl<=numel(sel_d));
bl = bl(valid);  dl = dl(valid);
n_matched = numel(bl);
fprintf('%d matched cells in plane %d\n', n_matched, selected_plane_idx);

%% ====================== LOOM ONSETS ======================
loom_idx_b = find(strcmp({B.Stimuli.type}, 'looming'), 1);
sb = B.Stimuli(loom_idx_b);
PT_b = B.Triggers.TimeProjector;
onsets_b = PT_b(PT_b >= sb.TimeStimulusFrame(1) & PT_b <= sb.TimeStimulusFrame(end))';

loom_idx_d = find(strcmp({D.Stimuli.type}, 'looming'), 1);
sd = D.Stimuli(loom_idx_d);
PT_d = D.Triggers.TimeProjector;
onsets_d = PT_d(PT_d >= sd.TimeStimulusFrame(1) & PT_d <= sd.TimeStimulusFrame(end))';

trial_t  = sb.stimulus_trial_t;   % 7s -- loom offset at t=+7s relative to onset
n_trials = min(numel(onsets_b), numel(onsets_d));
fprintf('%d baseline trials, %d drug trials -- using first %d from each\n', ...
    numel(onsets_b), numel(onsets_d), n_trials);

%% ====================== ALIGN TRIALS (wide window) ======================
t_com = linspace(-pre_window_sec, post_window_sec, n_com);

% Full ROI extraction, then restrict to matched cells
Ca_b_all = align_trials_to_onsets(onsets_b, dff_base, TimeCa_base, pre_window_sec, post_window_sec, n_com);
Ca_d_all = align_trials_to_onsets(onsets_d, dff_drug,  TimeCa_drug,  pre_window_sec, post_window_sec, n_com);
% [n_rois_plane x n_com x n_trials_per_condition]

Ca_b = Ca_b_all(bl, :, 1:n_trials);  % [n_matched x n_com x n_trials]
Ca_d = Ca_d_all(dl, :, 1:n_trials);

%% ====================== HEATMAP PREP ======================
pre_idx = t_com < 0;

% Trial-median per cell, then full-trace row z-score for display
med_b = median(Ca_b, 3, 'omitnan');   % [n_matched x n_com]
med_d = median(Ca_d, 3, 'omitnan');

Zh_b = row_zscore(med_b);
Zh_d = row_zscore(med_d);

% Sort rows by drug pre-loom raw dF/F (mean across pre window and trials)
pre_raw_d = mean(reshape(Ca_d(:, pre_idx, :), n_matched, [], n_trials), [2 3], 'omitnan');
[~, ord]  = sort(pre_raw_d, 'descend');

Zh_b_s = Zh_b(ord, :);
Zh_d_s = Zh_d(ord, :);

%% ====================== POPULATION TRACES ======================
% Per-trial population mean (raw dF/F, mean across matched cells): [n_trials x n_com]
pop_b_trials = reshape(mean(Ca_b, 1, 'omitnan'), n_com, n_trials)';
pop_d_trials = reshape(mean(Ca_d, 1, 'omitnan'), n_com, n_trials)';

% Mean ± SEM across matched cells (using per-cell trial-median to remove
% trial noise before the cross-cell average)
pop_mean_b = mean(med_b, 1, 'omitnan');
pop_sem_b  = std(med_b, 0, 1, 'omitnan') / sqrt(n_matched);
pop_mean_d = mean(med_d, 1, 'omitnan');
pop_sem_d  = std(med_d, 0, 1, 'omitnan') / sqrt(n_matched);

%% ====================== FIGURE ======================
cmap      = diverging_cmap();
clim_hm   = max(prctile(abs([Zh_b_s(:); Zh_d_s(:)]), 99), 0.1) * [-1 1];
trial_clrs = lines(n_trials);

fig = figure('Name', 'Looming_Individual_Cells_0312', 'NumberTitle', 'off', ...
    'Position', [40 40 1400 980]);

% ---- Row 1: heatmaps ----
ax1 = subplot(3, 2, 1);
imagesc(t_com, 1:n_matched, Zh_b_s, clim_hm);
colormap(ax1, cmap);
hold on;
xline(0,       'k--', 'LineWidth', 1.5);
xline(trial_t, 'k:',  'LineWidth', 1.5);
add_loom_shade(trial_t);
hold off;
xlabel('Time from loom onset (s)');
ylabel(sprintf('Cell # (n=%d, sorted by drug pre-loom activity)', n_matched));
title('Baseline  –  trial-median, row z-scored');
cb = colorbar; cb.Label.String = 'z-score'; cb.TickDirection = 'out';
set(ax1, 'YDir', 'normal');
pos1 = get(ax1, 'Position');

ax2 = subplot(3, 2, 2);
imagesc(t_com, 1:n_matched, Zh_d_s, clim_hm);
colormap(ax2, cmap);
hold on;
xline(0,       'k--', 'LineWidth', 1.5);
xline(trial_t, 'k:',  'LineWidth', 1.5);
add_loom_shade(trial_t);
hold off;
xlabel('Time from loom onset (s)');
title('Drug (yohimbine)  –  same row order');
colorbar; set(ax2, 'YDir', 'normal');
pos2 = get(ax2, 'Position');
set(ax2, 'Position', pos2);
set(ax1, 'Position', [pos1(1) pos2(2) pos2(3) pos2(4)]);

% ---- Row 2: per-trial spaghetti (raw dF/F) ----
trial_labels = arrayfun(@(k) sprintf('Trial %d', k), 1:n_trials, 'UniformOutput', false);

ax3 = subplot(3, 2, 3);
hold on;
ht_b = gobjects(1, n_trials);
for k = 1:n_trials
    ht_b(k) = plot(t_com, pop_b_trials(k, :), 'Color', [trial_clrs(k,:) 0.45], 'LineWidth', 1);
end
hm_b = plot(t_com, mean(pop_b_trials, 1, 'omitnan'), 'Color', c_base, 'LineWidth', 2.5);
suppress_from_legend(xline(0,       'k--', 'LineWidth', 1.2));
suppress_from_legend(xline(trial_t, 'k:',  'LineWidth', 1.2));
suppress_from_legend(yline(0, 'Color', [0.75 0.75 0.75]));
hold off;
xlabel('Time from loom onset (s)'); ylabel('Population mean dF/F');
title('Baseline  –  individual trials (bold = mean)');
xlim([-pre_window_sec post_window_sec]);
legend([ht_b, hm_b], [trial_labels, {'Mean'}], 'Location', 'northeast', 'FontSize', 7);

ax4 = subplot(3, 2, 4);
hold on;
ht_d = gobjects(1, n_trials);
for k = 1:n_trials
    ht_d(k) = plot(t_com, pop_d_trials(k, :), 'Color', [trial_clrs(k,:) 0.45], 'LineWidth', 1);
end
hm_d = plot(t_com, mean(pop_d_trials, 1, 'omitnan'), 'Color', c_drug, 'LineWidth', 2.5);
suppress_from_legend(xline(0,       'k--', 'LineWidth', 1.2));
suppress_from_legend(xline(trial_t, 'k:',  'LineWidth', 1.2));
suppress_from_legend(yline(0, 'Color', [0.75 0.75 0.75]));
hold off;
xlabel('Time from loom onset (s)');
title('Drug (yohimbine)  –  individual trials (bold = mean)');
xlim([-pre_window_sec post_window_sec]);
legend([ht_d, hm_d], [trial_labels, {'Mean'}], 'Location', 'northeast', 'FontSize', 7);

% ---- Row 3: mean ± SEM overlay ----
ax5 = subplot(3, 1, 3);
hold on;
suppress_from_legend(fill([t_com fliplr(t_com)], [pop_mean_b+pop_sem_b, fliplr(pop_mean_b-pop_sem_b)], ...
    c_base, 'FaceAlpha', 0.25, 'EdgeColor', 'none'));
hb = plot(t_com, pop_mean_b, 'Color', c_base, 'LineWidth', 2.5);

suppress_from_legend(fill([t_com fliplr(t_com)], [pop_mean_d+pop_sem_d, fliplr(pop_mean_d-pop_sem_d)], ...
    c_drug, 'FaceAlpha', 0.25, 'EdgeColor', 'none'));
hd = plot(t_com, pop_mean_d, 'Color', c_drug, 'LineWidth', 2.5);

suppress_from_legend(xline(0,       'k--', 'LineWidth', 1.2));
suppress_from_legend(xline(trial_t, 'k:',  'LineWidth', 1.2));
suppress_from_legend(yline(0, 'Color', [0.75 0.75 0.75]));
hold off;

% Phase labels just below the top of the axis
yl = ylim;
text(-pre_window_sec/2, yl(2)*0.97, 'Pre-loom', 'HorizontalAlignment', 'center', ...
    'FontSize', 8.5, 'Color', [0.45 0.45 0.45]);
text(trial_t/2, yl(2)*0.97, 'During loom', 'HorizontalAlignment', 'center', ...
    'FontSize', 8.5, 'Color', [0.45 0.45 0.45]);
text(trial_t + (post_window_sec - trial_t)/2, yl(2)*0.97, 'Post-loom', ...
    'HorizontalAlignment', 'center', 'FontSize', 8.5, 'Color', [0.45 0.45 0.45]);

xlabel('Time from loom onset (s)');
ylabel(sprintf('Population mean dF/F \\pm SEM  (n=%d matched cells)', n_matched));
title('Trial-median per cell, then mean ± SEM across cells');
legend([hb, hd], {'Baseline', 'Drug (yohimbine)'}, 'Location', 'best');
xlim([-pre_window_sec post_window_sec]);

sgtitle(sprintf('AnimalRB19 260312  |  looming full window  |  1115 baseline vs 1243 yohimbine'), ...
    'FontWeight', 'bold', 'Interpreter', 'none', 'FontSize', 12);

% ---- Save ----
saveas(fig, fullfile(outdir, 'looming_full_window_individual_cells.png'));
fprintf('Saved → %s\n', fullfile(outdir, 'looming_full_window_individual_cells.png'));

%% ====================== LOCAL FUNCTIONS ======================
function Z = row_zscore(X)
    m = mean(X, 2, 'omitnan');
    s = std(X, 0, 2, 'omitnan');
    Z = (X - m) ./ s;
end

function suppress_from_legend(h)
    h.HandleVisibility = 'off';
end

function add_loom_shade(trial_t)
    yl = ylim; xl = xlim;
    x0 = max(0, xl(1)); x1 = min(trial_t, xl(2));
    h = patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.95 0.95 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h.HandleVisibility = 'off';
end

function cmap = diverging_cmap(n)
    if nargin < 1; n = 256; end
    half = floor(n/2);
    neg = [linspace(0.05,1,half)', linspace(0.05,1,half)', ones(half,1)];
    pos = [ones(n-half,1), linspace(1,0.05,n-half)', linspace(1,0.05,n-half)'];
    cmap = [neg; pos];
end
