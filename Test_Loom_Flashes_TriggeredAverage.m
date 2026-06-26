%% Quick validation: triggered-average response for looming and flashes
% Baseline recording only. Confirms that (1) looming's per-trial red-frame
% onsets and (2) flashes' reconstructed onsets (via reconstructFlashes.m)
% both line up with real calcium transients -- a rise shortly after t=0 in
% both means timing is correct; a flat/noisy average means it isn't.
%
% TimeStimulusFrame is now directly camera-clock-correct end-to-end
% (RedFrameSynchronized=true for every stimulus block, including looming)
% so no manual anchor-correction against TimeProjector is needed for
% either stimulus anymore -- see project_aggregation_pipeline_fixes memory.

addpath(genpath(fullfile(pwd, 'violin_plot_utils')));

%% ====================== CONFIGURATION ======================
baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat";
selected_plane_idx = 1;

%% ====================== LOAD DATA ======================
fprintf('Loading baseline: %s\n', baseline_file);
B = load(baseline_file);

base_dff_full = B.CaData(1).Ca_dFF;
centroid = B.CaData(1).Ca_centroid_voxel;
centroidZ = centroid(:, 3);
unique_planes = unique(centroidZ);
selected_plane = unique_planes(selected_plane_idx);
selected_roi_idx = find(centroidZ == selected_plane);
base_dff = base_dff_full(selected_roi_idx, :);

%% ====================== LOOMING ======================
loom_idx = find(strcmp({B.Stimuli.type}, 'looming'), 1);
loom_s = B.Stimuli(loom_idx);
if ~loom_s.RedFrameSynchronized
    warning('Looming: RedFrameSynchronized=false -- timing fell back to PC clock, expect an offset');
end

PT = B.Triggers.TimeProjector;
% Inclusive bounds: TimeStimulusFrame(1) now lands exactly on the trial-1
% pulse, so a strict > would silently drop it.
loom_onsets = PT(PT >= loom_s.TimeStimulusFrame(1) & PT <= loom_s.TimeStimulusFrame(end))';
fprintf('Looming: %d onsets found (expected %d trials)\n', numel(loom_onsets), loom_s.trials);

[Ca_loom, t_com_loom] = align_trials_to_onsets(loom_onsets, base_dff, B.TimeCa, selected_plane, 2, 5, 100);
plot_triggered_average(t_com_loom, Ca_loom, 'Looming', [0.2 0.6 1]);

%% ====================== FLASHES ======================
flash_idx = find(strcmp({B.Stimuli.type}, 'sparse_local_global_flashes'), 1);
flash_s = B.Stimuli(flash_idx);
if ~flash_s.RedFrameSynchronized
    warning('Flashes: RedFrameSynchronized=false -- timing fell back to PC clock, expect an offset');
end

[flash_onsets_rel, ~, ~] = reconstructFlashes(flash_s);
flash_onsets_abs = flash_s.TimeStimulusFrame(1) + flash_onsets_rel;
fprintf('Flashes: %d reconstructed onsets\n', numel(flash_onsets_abs));

[Ca_flash, t_com_flash] = align_trials_to_onsets(flash_onsets_abs, base_dff, B.TimeCa, selected_plane, 1, 4, 100);
plot_triggered_average(t_com_flash, Ca_flash, 'Sparse local/global flashes', [1 0.5 0.2]);

%% ============= HELPER FUNCTIONS =============
function [Ca_trials, t_com] = align_trials_to_onsets(onsets, dff, TimeCa, plane, pre_window_sec, post_window_sec, n_com)
    t_com     = linspace(-pre_window_sec, post_window_sec, n_com);
    n_rois    = size(dff, 1);
    n_trials  = numel(onsets);
    Ca_trials = nan(n_rois, n_com, n_trials);
    time_vec  = TimeCa(plane, :);

    for k = 1:n_trials
        idx = find(time_vec > onsets(k) - pre_window_sec & time_vec < onsets(k) + post_window_sec);
        if numel(idx) < 2; continue; end
        t_rel = time_vec(idx) - onsets(k);
        Ca_trials(:, :, k) = interp1(t_rel(:), dff(:, idx)', t_com(:), 'linear', 'extrap')';
    end
end

function plot_triggered_average(t_com, Ca_trials, fig_name, mean_color)
    n_trials      = size(Ca_trials, 3);
    pop_per_trial = squeeze(mean(Ca_trials, 1, 'omitnan'));   % [n_com x n_trials]
    grand_mean    = mean(pop_per_trial, 2, 'omitnan');
    sem_trace     = std(pop_per_trial, 0, 2, 'omitnan') / sqrt(n_trials);

    figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [100 100 900 500], 'Visible', 'off');
    plot(t_com, pop_per_trial, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); hold on;
    x_sh = [t_com, fliplr(t_com)];
    y_sh = [(grand_mean + sem_trace); flipud(grand_mean - sem_trace)]';
    fill(x_sh, y_sh, mean_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, grand_mean, '-', 'Color', mean_color, 'LineWidth', 2.5);
    xline(0, 'r--', 'LineWidth', 1);
    xlabel('Time from onset (s)', 'FontSize', 11);
    ylabel('dF/F (population mean)', 'FontSize', 11);
    title(sprintf('%s -- triggered average (n = %d trials)', fig_name, n_trials), 'FontSize', 12, 'Interpreter', 'none');
    set(gca, 'Box', 'off'); xlim([t_com(1), t_com(end)]);

    safe_name = regexprep(fig_name, '[^a-zA-Z0-9]+', '_');
    outfile = fullfile(fileparts(mfilename('fullpath')), sprintf('test_%s.png', safe_name));
    saveas(gcf, outfile);
    fprintf('Saved %s\n', outfile);
end
