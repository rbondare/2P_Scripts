%% Flexible heatmap explorer -- one script, pick your own options, re-run
%
% Replaces having to ask for a new one-off script every time you want to
% look at the data a different way. Edit the CONFIGURATION block below
% and re-run (or run section-by-section, Ctrl+Enter) -- nothing else
% needs to change.
%
% view_mode:
%   'entire'       -> the continuous whole-block trace, real elapsed time,
%                     NOT averaged across trials/repeats/directions, gaps
%                     between trials INCLUDED. The "raw, as recorded" view.
%   'trial_blocks' -> each trial/event/subtrial's own aligned snippet
%                     concatenated back-to-back, in order, with the
%                     (often random, uninformative) gaps BETWEEN them
%                     skipped entirely -- only the actual stimulus-relevant
%                     windows, no dead time. Not available for
%                     spontaneous/checkers2 (no discrete trial structure).
%   'averaged'     -> each cell's trial/repeat-averaged response:
%                       looming/flashes:        median across trials, aligned to onset
%                       grating/moving_bar:      each cell's own preferred-direction snippet
%                       spontaneous/checkers2:   same as 'entire' (no trial structure
%                                                to average over -- there's only 1 block)
%   'per_trial'    -> rows are TRIALS, not cells: each row is the
%                     population mean (across the chosen cell_set) for
%                     ONE trial/event/subtrial, so you can see how
%                     different each individual repeat is from the others
%                     BEFORE collapsing them together in 'averaged' mode.
%                     Not available for spontaneous/checkers2.
%
% Self-contained: does its own extraction directly from the lower-level
% analysis_functions (extract_looming_trials, extract_direction_trials,
% etc.) rather than depending on the production pipeline's analyze_*_stimulus
% wrapper functions, so this script can show ANY cell subset (matched or
% all) for ANY stimulus type without needing pipeline changes.

clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

%% ====================== CONFIGURATION ======================
baseline_file   = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1444_preprocessed.mat";
drug_file       = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1609_preprocessed.mat";
roi_match_file  = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_1444_V1.mat";
selected_plane_idx = 1; 

%%
stimulus_type   = 'looming';    % 'looming' | 'sparse_local_global_flashes' | 'checkers2' | 'grating' | 'moving_bar' | 'spontaneous'
view_mode       = 'trial_blocks'; % 'entire' | 'trial_blocks' | 'averaged' | 'per_trial'  (see notes above)
cell_set        = 'matched';    % 'matched' | 'all'
n_cells_to_plot = 500;           % how many rows to show (top-N by peak activity); set Inf to show all
sort_by         = 'baseline';   % 'baseline' | 'drug' | 'none'
metric          = 'zscore';        % 'raw' | 'zscore'

%% ====================== LOAD DATA ======================
% Self-contained -- doesn't depend on Preprocess_Stimulus_Data.m's
% current state, so changing baseline_file/drug_file/roi_match_file above
% is the only thing you need to touch to switch recording pairs.
pd = load_pair(baseline_file, drug_file, roi_match_file, selected_plane_idx);

if strcmp(cell_set, 'matched') && ~pd.matched_rois_available
    error('cell_set=''matched'' requested but no matched ROIs are available for this pair.');
end
if strcmp(cell_set, 'matched')
    base_idx = pd.match_idx_local_base; drug_idx = pd.match_idx_local_drug;
else
    base_idx = 1:size(pd.dff_base, 1); drug_idx = 1:size(pd.dff_drug, 1);
end

%% ====================== EXTRACT ======================
boundaries = []; block_labels = {};
effective_view_mode = view_mode;
if ismember(stimulus_type, {'spontaneous', 'checkers2'}) && ~strcmp(view_mode, 'entire')
    fprintf('  Note: %s has no discrete trial structure -- using ''entire'' regardless of view_mode\n', stimulus_type);
    effective_view_mode = 'entire';
end

switch effective_view_mode
    case 'entire'
        [t_axis, resp_base, resp_drug] = extract_entire(stimulus_type, pd.Stimuli_base, pd.dff_base, pd.TimeCa_base, ...
            pd.Stimuli_drug, pd.dff_drug, pd.TimeCa_drug, base_idx, drug_idx);
    case 'trial_blocks'
        [t_axis, resp_base, resp_drug, boundaries, block_labels] = extract_trial_blocks(stimulus_type, ...
            pd.Stimuli_base, pd.dff_base, pd.TimeCa_base, pd.Triggers_base, ...
            pd.Stimuli_drug, pd.dff_drug, pd.TimeCa_drug, pd.Triggers_drug, base_idx, drug_idx);
    case 'averaged'
        [t_axis, resp_base, resp_drug] = extract_averaged(stimulus_type, pd.Stimuli_base, pd.dff_base, pd.TimeCa_base, pd.Triggers_base, ...
            pd.Stimuli_drug, pd.dff_drug, pd.TimeCa_drug, pd.Triggers_drug, base_idx, drug_idx);
    case 'per_trial'
        [t_axis, resp_base, resp_drug] = extract_per_trial_population(stimulus_type, ...
            pd.Stimuli_base, pd.dff_base, pd.TimeCa_base, pd.Triggers_base, ...
            pd.Stimuli_drug, pd.dff_drug, pd.TimeCa_drug, pd.Triggers_drug, base_idx, drug_idx);
    otherwise
        error('Unknown view_mode ''%s''.', view_mode);
end
row_label = 'Cell';
if strcmp(effective_view_mode, 'per_trial'); row_label = 'Trial'; end

if strcmp(metric, 'zscore')
    [resp_base, resp_drug] = quick_zscore(resp_base, resp_drug, t_axis);
    cbar_label = 'z-score'; center_zero = true;
else
    cbar_label = 'dF/F'; center_zero = false;
end

%% ====================== PLOT ======================
fig_title = sprintf('%s_%s_%s_%s_top%s', stimulus_type, view_mode, cell_set, metric, ...
    num2str(min(n_cells_to_plot, size(resp_base,1))));
plot_explorer_heatmap(t_axis, resp_base, resp_drug, n_cells_to_plot, sort_by, center_zero, cbar_label, fig_title, ...
    fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', '_explorer'), boundaries, block_labels);

%% ====================== LOCAL FUNCTIONS ======================
function pair = load_pair(baseline_file, drug_file, roi_match_file, selected_plane_idx)
    B = load(baseline_file);
    base_dff_full = B.CaData(1).Ca_dFF;
    centroidZ = B.CaData(1).Ca_centroid_voxel(:, 3);
    unique_planes = unique(centroidZ);
    selected_roi_idx = find(centroidZ == unique_planes(selected_plane_idx));
    base_dff = base_dff_full(selected_roi_idx, :);

    D = load(drug_file);
    drug_dff_full = D.CaData(1).Ca_dFF;
    drug_centroidZ = D.CaData(1).Ca_centroid_voxel(:, 3);
    unique_planes_drug = unique(drug_centroidZ);
    drug_selected_roi_idx = find(drug_centroidZ == unique_planes_drug(selected_plane_idx));
    drug_dff = drug_dff_full(drug_selected_roi_idx, :);

    matched_rois_available = false;
    base_match_idx_local = []; drug_match_idx_local = [];
    if ~isempty(roi_match_file) && isfile(roi_match_file)
        M = load(roi_match_file);
        base_match_idx_global = M.roiMatchData.allSessionMapping(:, 1);
        drug_match_idx_global = M.roiMatchData.allSessionMapping(:, 2);
        [~, base_match_idx_local] = ismember(base_match_idx_global, selected_roi_idx);
        [~, drug_match_idx_local] = ismember(drug_match_idx_global, drug_selected_roi_idx);
        valid_matches = (base_match_idx_local > 0) & (drug_match_idx_local > 0);
        base_match_idx_local = base_match_idx_local(valid_matches);
        drug_match_idx_local = drug_match_idx_local(valid_matches);
        matched_rois_available = ~isempty(base_match_idx_local);
    end

    pair.Stimuli_base = B.Stimuli;
    pair.Stimuli_drug = D.Stimuli;
    pair.dff_base = base_dff;
    pair.dff_drug = drug_dff;
    pair.TimeCa_base = B.TimeCa;
    pair.TimeCa_drug = D.TimeCa;
    pair.Triggers_base = B.Triggers;
    pair.Triggers_drug = D.Triggers;
    pair.match_idx_local_base = base_match_idx_local;
    pair.match_idx_local_drug = drug_match_idx_local;
    pair.matched_rois_available = matched_rois_available;
end

function [t_axis, resp_base, resp_drug] = extract_entire(stype, base_Stimuli, base_dff, base_TimeCa, ...
        drug_Stimuli, drug_dff, drug_TimeCa, base_idx, drug_idx)
    resp_base_cell = extract_full_stimulus_responses(base_Stimuli, base_dff, stype, base_TimeCa);
    resp_drug_cell = extract_full_stimulus_responses(drug_Stimuli, drug_dff, stype, drug_TimeCa);
    if isempty(resp_base_cell) || isempty(resp_drug_cell)
        error('%s missing in baseline or drug.', stype);
    end
    block_base = resp_base_cell{1}; block_drug = resp_drug_cell{1};
    n_com = 300;
    dt_base = median(diff(base_TimeCa(1, :)), 'omitnan'); dur_base = (size(block_base, 2) - 1) * dt_base;
    dt_drug = median(diff(drug_TimeCa(1, :)), 'omitnan'); dur_drug = (size(block_drug, 2) - 1) * dt_drug;
    % Baseline and drug can genuinely have different real durations (e.g.
    % looming's inter-trial waits are randomized per session) -- use
    % whichever is SHORTER so both conditions' query points always stay
    % within their own actually-recorded range. Using dur_base alone would
    % ask interp1 for points beyond drug's real data when drug is shorter,
    % which returns NaN -- and imagesc renders NaN using the colormap's
    % first color (deep blue here), making "no data" look like "extreme
    % negative." Truncating to the shared/overlapping duration avoids
    % fabricating or mis-coloring anything.
    dur_common = min(dur_base, dur_drug);
    if abs(dur_base - dur_drug) > 1
        fprintf('  Note: baseline (%.1fs) and drug (%.1fs) block durations differ -- showing the shared %.1fs\n', ...
            dur_base, dur_drug, dur_common);
    end
    t_axis = linspace(0, dur_common, n_com);
    src_base = linspace(0, dur_base, size(block_base, 2));
    src_drug = linspace(0, dur_drug, size(block_drug, 2));
    resp_base = interp1(src_base, block_base(base_idx, :)', t_axis)';
    resp_drug = interp1(src_drug, block_drug(drug_idx, :)', t_axis)';
end

function [t_axis, resp_base, resp_drug, boundaries, block_labels] = extract_trial_blocks(stype, ...
        base_Stimuli, base_dff, base_TimeCa, base_Triggers, drug_Stimuli, drug_dff, drug_TimeCa, drug_Triggers, base_idx, drug_idx)
    % Concatenates each trial/event/subtrial's own aligned snippet
    % end-to-end, in natural order, skipping the gaps BETWEEN them
    % entirely -- avoids 'entire' mode's baseline/drug duration-mismatch
    % risk for stimuli with randomized inter-trial waits (looming), since
    % trial COUNT is what's compared here, not real elapsed time.
    switch stype
        case 'looming'
            [Ca_base, t_snippet] = extract_looming_trials(base_Stimuli, base_dff, base_TimeCa, base_Triggers);
            [Ca_drug, ~] = extract_looming_trials(drug_Stimuli, drug_dff, drug_TimeCa, drug_Triggers);
        case 'sparse_local_global_flashes'
            pre_window_sec = 1; post_window_sec = 4; n_com = 100;
            base_onsets = get_flash_onsets(base_Stimuli);
            drug_onsets = get_flash_onsets(drug_Stimuli);
            Ca_base = align_trials_to_onsets(base_onsets, base_dff, base_TimeCa, pre_window_sec, post_window_sec, n_com);
            Ca_drug = align_trials_to_onsets(drug_onsets, drug_dff, drug_TimeCa, pre_window_sec, post_window_sec, n_com);
            t_snippet = linspace(-pre_window_sec, post_window_sec, n_com);
        case {'grating', 'moving_bar'}
            if strcmp(stype, 'grating'); orientation_field = 'grating_orientations'; else; orientation_field = 'bar_orientations'; end
            [Ca_base, t_snippet, ~] = extract_direction_trials(base_Stimuli, base_dff, base_TimeCa, stype, orientation_field, 1, 6, 100);
            [Ca_drug, ~, ~] = extract_direction_trials(drug_Stimuli, drug_dff, drug_TimeCa, stype, orientation_field, 1, 6, 100);
        otherwise
            error('trial_blocks mode is not defined for ''%s'' (no discrete trial/event structure).', stype);
    end
    n_trials_base = size(Ca_base, 3); n_trials_drug = size(Ca_drug, 3);
    n_trials = min(n_trials_base, n_trials_drug);
    if n_trials_base ~= n_trials_drug
        fprintf('  Note: baseline has %d trials/events, drug has %d -- showing the first %d from each\n', ...
            n_trials_base, n_trials_drug, n_trials);
    end
    n_com = numel(t_snippet);
    resp_base = reshape(Ca_base(base_idx, :, 1:n_trials), numel(base_idx), n_com * n_trials);
    resp_drug = reshape(Ca_drug(drug_idx, :, 1:n_trials), numel(drug_idx), n_com * n_trials);
    t_axis = 1:(n_com * n_trials);
    boundaries = (0:n_trials-1) * n_com + 1;
    block_labels = arrayfun(@(k) sprintf('Trial %d', k), 1:n_trials, 'UniformOutput', false);
end

function [t_axis, resp_base, resp_drug] = extract_averaged(stype, base_Stimuli, base_dff, base_TimeCa, base_Triggers, ...
        drug_Stimuli, drug_dff, drug_TimeCa, drug_Triggers, base_idx, drug_idx)
    switch stype
        case 'looming'
            [Ca_base, t_axis] = extract_looming_trials(base_Stimuli, base_dff, base_TimeCa, base_Triggers);
            [Ca_drug, ~] = extract_looming_trials(drug_Stimuli, drug_dff, drug_TimeCa, drug_Triggers);
            resp_base = median(Ca_base(base_idx, :, :), 3, 'omitnan');
            resp_drug = median(Ca_drug(drug_idx, :, :), 3, 'omitnan');

        case 'sparse_local_global_flashes'
            pre_window_sec = 1; post_window_sec = 4; n_com = 100;
            base_onsets = get_flash_onsets(base_Stimuli);
            drug_onsets = get_flash_onsets(drug_Stimuli);
            Ca_base = align_trials_to_onsets(base_onsets, base_dff, base_TimeCa, pre_window_sec, post_window_sec, n_com);
            Ca_drug = align_trials_to_onsets(drug_onsets, drug_dff, drug_TimeCa, pre_window_sec, post_window_sec, n_com);
            t_axis = linspace(-pre_window_sec, post_window_sec, n_com);
            resp_base = median(Ca_base(base_idx, :, :), 3, 'omitnan');
            resp_drug = median(Ca_drug(drug_idx, :, :), 3, 'omitnan');

        case {'grating', 'moving_bar'}
            if strcmp(stype, 'grating'); orientation_field = 'grating_orientations'; else; orientation_field = 'bar_orientations'; end
            pre_window_sec = 1; post_window_sec = 6; n_com = 100;
            [Ca_base, t_axis, dir_base] = extract_direction_trials(base_Stimuli, base_dff, base_TimeCa, ...
                stype, orientation_field, pre_window_sec, post_window_sec, n_com);
            [Ca_drug, ~, dir_drug] = extract_direction_trials(drug_Stimuli, drug_dff, drug_TimeCa, ...
                stype, orientation_field, pre_window_sec, post_window_sec, n_com);
            [Ca_by_dir_base, dirs_base] = chunk_by_direction(Ca_base(base_idx, :, :), dir_base);
            [Ca_by_dir_drug, dirs_drug] = chunk_by_direction(Ca_drug(drug_idx, :, :), dir_drug);
            common_dirs = intersect(dirs_base, dirs_drug);
            [~, ib] = ismember(common_dirs, dirs_base); [~, idr] = ismember(common_dirs, dirs_drug);
            Ca_by_dir_base = Ca_by_dir_base(:, :, ib); Ca_by_dir_drug = Ca_by_dir_drug(:, :, idr);
            peak_idx = t_axis >= 0.5 & t_axis <= min(4, post_window_sec);
            amp_base = squeeze(mean(Ca_by_dir_base(:, peak_idx, :), 2, 'omitnan'));
            [~, pref] = max(amp_base, [], 2);
            n = size(Ca_by_dir_base, 1);
            resp_base = nan(n, numel(t_axis)); resp_drug = nan(n, numel(t_axis));
            for c = 1:n
                resp_base(c, :) = Ca_by_dir_base(c, :, pref(c));
                resp_drug(c, :) = Ca_by_dir_drug(c, :, pref(c));
            end

        otherwise
            error('No averaged-mode extraction defined for stimulus type ''%s''.', stype);
    end
end

function [Zb, Zd] = quick_zscore(resp_base, resp_drug, t_axis)
    % Pre-onset baseline z-score if there's a meaningful pre-onset window
    % (t_axis has negative values, i.e. 'averaged' mode for a triggered
    % stimulus); otherwise population-relative z-score (no pre-onset
    % period exists in 'entire'/blockwise mode).
    if any(t_axis < 0)
        pre = t_axis < 0;
        bm = mean(resp_base(:, pre), 2, 'omitnan'); bs = std(resp_base(:, pre), 0, 2, 'omitnan');
        dm = mean(resp_drug(:, pre), 2, 'omitnan'); ds = std(resp_drug(:, pre), 0, 2, 'omitnan');
        Zb = (resp_base - bm) ./ bs;
        Zd = (resp_drug - dm) ./ ds;
    else
        cm_b = mean(resp_base, 2, 'omitnan'); pm_b = mean(cm_b, 'omitnan'); ps_b = std(cm_b, 'omitnan');
        cm_d = mean(resp_drug, 2, 'omitnan'); pm_d = mean(cm_d, 'omitnan'); ps_d = std(cm_d, 'omitnan');
        Zb = (resp_base - pm_b) / ps_b;
        Zd = (resp_drug - pm_d) / ps_d;
    end
end

function plot_explorer_heatmap(t_axis, resp_base, resp_drug, n_cells_to_plot, sort_by, center_zero, cbar_label, fig_title, outdir, boundaries, block_labels)
    if ~exist(outdir, 'dir'); mkdir(outdir); end
    switch sort_by
        case 'baseline'
            key = max(resp_base, [], 2); ylab_extra = 'sorted by baseline peak';
        case 'drug'
            key = max(resp_drug, [], 2); ylab_extra = 'sorted by drug peak';
        otherwise
            key = []; ylab_extra = 'unsorted';
    end
    if isempty(key)
        order = 1:size(resp_base, 1);
    else
        key(isnan(key)) = -Inf;
        [~, order] = sort(key, 'descend');
    end
    n_show = min(n_cells_to_plot, numel(order));
    order = order(1:n_show);
    Zb = resp_base(order, :);
    Zd = resp_drug(order, :);

    if center_zero
        clim = prctile([Zb(:); Zd(:)], [1 99]);
        clim = max(abs(clim)) * [-1 1];
        if any(~isfinite(clim)) || clim(1) == clim(2); clim = [-1 1]; end
        colormap_to_use = diverging_colormap_local();
        onset_color = 'k--';
    else
        clim = [0, prctile([Zb(:); Zd(:)], 99)];
        if any(~isfinite(clim)) || clim(1) == clim(2); clim = [0 1]; end
        colormap_to_use = flipud(gray);
        onset_color = 'r--';
    end

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 1100 650]);
    colormap(colormap_to_use);

    % AlphaData makes any NaN cell transparent (shows the figure's plain
    % background) instead of imagesc's default of coloring NaN as the
    % colormap's first entry -- which can look like a real extreme value
    % (e.g. solid blue on the diverging map) when it's actually just
    % missing data (out-of-range time points, std=0 cells, etc).
    has_blocks = nargin >= 11 && ~isempty(boundaries);

    ax1 = subplot(1, 2, 1);
    h1 = imagesc(t_axis, 1:size(Zb, 1), Zb, clim);
    set(h1, 'AlphaData', ~isnan(Zb));
    hold on;
    if has_blocks
        for b = boundaries(2:end); xline(b - 0.5, onset_color, 'LineWidth', 1); end
        set(gca, 'XTick', boundaries, 'XTickLabel', block_labels); xtickangle(45);
    else
        xline(0, onset_color, 'LineWidth', 1);
    end
    hold off;
    xlabel('Time / position'); ylabel(sprintf('Cell # (%s)', ylab_extra)); title('Baseline');
    pos1 = get(ax1, 'Position');

    ax2 = subplot(1, 2, 2);
    h2 = imagesc(t_axis, 1:size(Zd, 1), Zd, clim);
    set(h2, 'AlphaData', ~isnan(Zd));
    hold on;
    if has_blocks
        for b = boundaries(2:end); xline(b - 0.5, onset_color, 'LineWidth', 1); end
        set(gca, 'XTick', boundaries, 'XTickLabel', block_labels); xtickangle(45);
    else
        xline(0, onset_color, 'LineWidth', 1);
    end
    hold off;
    xlabel('Time / position'); title('Drug');
    pos2 = get(ax2, 'Position');
    cb = colorbar; cb.Label.String = cbar_label;
    set(ax2, 'Position', pos2);
    set(ax1, 'Position', [pos1(1) pos2(2) pos2(3) pos2(4)]);

    sgtitle(sprintf('%s (showing %d cells)', fig_title, size(Zb,1)), 'FontWeight', 'bold', 'Interpreter', 'none');

    safe_name = regexprep(fig_title, '[^a-zA-Z0-9]+', '_');
    outfile = fullfile(outdir, [safe_name '.png']);
    saveas(gcf, outfile);
    fprintf('Saved %s\n', outfile);
end

function cmap = diverging_colormap_local(n)
    if nargin < 1; n = 256; end
    half = floor(n/2);
    neg = [linspace(0.05,1,half)', linspace(0.05,1,half)', ones(half,1)];
    pos = [ones(n-half,1), linspace(1,0.05,n-half)', linspace(1,0.05,n-half)'];
    cmap = [neg; pos];
end
