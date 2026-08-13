%% TEMP: quick unmatched looming check -- 260319 baseline (1459) vs 1606
%
% Throwaway sanity-check script. No ROI matching available/possible for
% this recording (planes don't align between sessions at all), so this
% just pools ALL cells across ALL planes for each recording independently
% (no cell-to-cell correspondence) and looks at the population loom
% response, individual trials + mean/median, baseline vs 1606. Purpose:
% decide whether this recording is even worth pursuing proper ROI matching
% for.
%
% CaData(1).Ca_dFF is already iscell-filtered -- verified directly against
% the raw suite2p iscell.npy for both recordings: summed iscell=true across
% the 3 functional planes (flyback plane3 excluded) matches the aggregated
% row count exactly (baseline 2347+2436+2482=7265; drug 2421+2333+2192=6946).
% So there's no separate "non-cell junk" to filter out here -- what's pooled
% below is already suite2p's classified cells. Median (not just mean) across
% cells is shown per trial so a small number of high-amplitude/noisy cells
% dominating the mean is directly visible as a mean/median split.

clear; clc; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260319_1459_preprocessed.mat";
drug_file     = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260319_1606_preprocessed.mat";

c_base = [0.3 0.3 0.3];   % grey
c_drug = [0.55 0.35 0.75]; % purple -- unassigned drug condition for this check

%% ====================== LOAD (all planes pooled, unmatched) ======================
B = load(baseline_file);
D = load(drug_file);

dff_base = B.CaData(1).Ca_dFF;   % all ROIs, all planes
dff_drug = D.CaData(1).Ca_dFF;
fprintf('Baseline (1459): %d ROIs across all planes\n', size(dff_base, 1));
fprintf('Drug (1606):     %d ROIs across all planes\n', size(dff_drug, 1));

[Ca_base, t_com] = extract_looming_trials(B.Stimuli, dff_base, B.TimeCa, B.Triggers);
[Ca_drug, ~]      = extract_looming_trials(D.Stimuli, dff_drug, D.TimeCa, D.Triggers);

if isempty(Ca_base) || isempty(Ca_drug)
    error('No looming stimulus found in one of these recordings.');
end

%% ====================== POPULATION PER-TRIAL TRACES ======================
% Mean across ALL cells (unmatched) per trial: [n_trials x n_com]. Individual
% trial lines use this (cross-cell mean per trial, same as before).
pop_base_trials = squeeze(mean(Ca_base, 1, 'omitnan'))';
pop_drug_trials = squeeze(mean(Ca_drug, 1, 'omitnan'))';
n_trials_base = size(pop_base_trials, 1);
n_trials_drug = size(pop_drug_trials, 1);
fprintf('%d baseline trials, %d drug trials\n', n_trials_base, n_trials_drug);

% Cross-cell MEDIAN per trial (separate from the mean above) -- if the mean
% and median summary lines diverge, a small number of high-amplitude/noisy
% cells are dominating the mean rather than a genuine population response.
pop_base_trials_med = squeeze(median(Ca_base, 1, 'omitnan'))';
pop_drug_trials_med = squeeze(median(Ca_drug, 1, 'omitnan'))';

global_ymin = min([pop_base_trials(:); pop_drug_trials(:); pop_base_trials_med(:); pop_drug_trials_med(:)]);
global_ymax = max([pop_base_trials(:); pop_drug_trials(:); pop_base_trials_med(:); pop_drug_trials_med(:)]);
pad = 0.05 * (global_ymax - global_ymin);
global_ylim = [global_ymin - pad, global_ymax + pad];

%% ====================== FIGURE ======================
% Trial lines colored light -> dark = trial 1 -> N (within each condition's
% own hue) so any trend of peak size with loom repetition number is visible.
trial_clrs_base = trial_gradient(c_base, n_trials_base);
trial_clrs_drug = trial_gradient(c_drug, n_trials_drug);

fig = figure('Name', 'Temp_Looming_0319_UnmatchedCheck', 'NumberTitle', 'off', ...
    'Position', [80 80 1100 450]);

ax_b = subplot(1, 2, 1);
hold on;
for k = 1:n_trials_base
    plot(t_com, pop_base_trials(k, :), 'Color', trial_clrs_base(k,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end
plot(t_com, mean(pop_base_trials, 1, 'omitnan'), 'Color', c_base, 'LineWidth', 2.5, 'DisplayName', 'Mean (across cells)');
plot(t_com, mean(pop_base_trials_med, 1, 'omitnan'), '--', 'Color', c_base, 'LineWidth', 2, 'DisplayName', 'Median (across cells)');
xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
yline(0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
hold off;
xlim([t_com(1), t_com(end)]); ylim(global_ylim);
xlabel('Time from onset (s)'); ylabel('dF/F');
title(sprintf('260319 baseline (1459) -- all planes, unmatched (n=%d cells, %d trials)', size(dff_base,1), n_trials_base));
subtitle('light = trial 1, dark = trial N', 'FontSize', 8, 'FontAngle', 'italic');
legend('Location', 'northeast', 'FontSize', 8);

ax_d = subplot(1, 2, 2);
hold on;
for k = 1:n_trials_drug
    plot(t_com, pop_drug_trials(k, :), 'Color', trial_clrs_drug(k,:), 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end
plot(t_com, mean(pop_drug_trials, 1, 'omitnan'), 'Color', c_drug, 'LineWidth', 2.5, 'DisplayName', 'Mean (across cells)');
plot(t_com, mean(pop_drug_trials_med, 1, 'omitnan'), '--', 'Color', c_drug, 'LineWidth', 2, 'DisplayName', 'Median (across cells)');
xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
yline(0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
hold off;
xlim([t_com(1), t_com(end)]); ylim(global_ylim);
xlabel('Time from onset (s)');
title(sprintf('260319 1606 -- all planes, unmatched (n=%d cells, %d trials)', size(dff_drug,1), n_trials_drug));
subtitle('light = trial 1, dark = trial N', 'FontSize', 8, 'FontAngle', 'italic');
legend('Location', 'northeast', 'FontSize', 8);

sgtitle('Looming, unmatched population (all planes pooled) -- 260319 1459 vs 1606', 'FontWeight', 'bold');

%% ====================== LOCAL FUNCTIONS ======================
function colors = trial_gradient(base_color, n)
    % Light tint of base_color (trial 1) linearly graded to base_color
    % itself (trial n), so line darkness tracks trial/loom-repetition order.
    if n <= 1
        colors = base_color;
        return;
    end
    light_color = base_color + (1 - base_color) * 0.8;  % near-white tint
    frac = linspace(0, 1, n)';
    colors = light_color + (base_color - light_color) .* frac;
end
