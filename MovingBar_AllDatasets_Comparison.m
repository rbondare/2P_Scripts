%% Moving_bar comparison across all 3 recording pairs
%
% One-off comparison script (not part of the 6-stimulus production
% pipeline) -- for each of the 3 baseline/drug pairs collected so far,
% produces:
%   1) The continuous whole-block heatmap (all 8 directions x 3 repeats
%      concatenated in real chronological order, NOT collapsed to each
%      cell's preferred direction) -- same resp_base_block/resp_drug_block
%      fields analyze_direction_stimulus.m already computes for the main
%      pipeline's "_heatmap_block_*" figures.
%   2) For the 20 highest-SNR matched cells (selected from the continuous
%      block, SNR = peak(baseline trace)/std(baseline trace)): the
%      PER-DIRECTION averaged dF/F (mean +/- SEM across just those 20
%      cells, one subplot per direction), baseline vs drug overlaid.
%      Directions are aligned via the same common_dirs intersection
%      analyze_direction_stimulus.m uses -- direction X in baseline is
%      guaranteed to be matched against the SAME physical angle X in drug,
%      not just the same array position.
%
% Self-contained: loads all 3 pairs itself (doesn't depend on
% Preprocess_Stimulus_Data.m's current config), via a load_pair local
% function so repeated calls can't clobber each other (MATLAB local
% functions have their own workspace).

clear; clc; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

%% ====================== CONFIGURATION: all 3 pairs ======================
pairs = struct( ...
    'name', {}, 'baseline_file', {}, 'drug_file', {}, 'roi_match_file', {});

pairs(1).name = 'AnimalRB19_260312_1115_vs_1243_yohimbine';
pairs(1).baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat";
pairs(1).drug_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat";
pairs(1).roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0312_V1.mat";

pairs(2).name = 'AnimalRB19_260305_1327_vs_1409_guanfacine';
pairs(2).baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat";
pairs(2).drug_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat";
pairs(2).roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_V2.mat";

pairs(3).name = 'AnimalRB19_260320_1444_vs_1609';
pairs(3).baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1444_preprocessed.mat";
pairs(3).drug_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1609_preprocessed.mat";
pairs(3).roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_1444_V1.mat";

selected_plane_idx = 1;
n_top_cells = 10;
n_heatmap_cells = 200;   % how many cells to show in the heatmap (top-N by baseline peak, after sorting); set to Inf to show all matched cells

outdir = fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', '_movingbar_all_datasets_comparison');
if ~exist(outdir, 'dir'); mkdir(outdir); end

%% ====================== PER-PAIR ANALYSIS ======================
for p = 1:numel(pairs)
    fprintf('\n========== %s ==========\n', pairs(p).name);
    pd = load_pair(pairs(p).baseline_file, pairs(p).drug_file, pairs(p).roi_match_file, selected_plane_idx);

    result = analyze_moving_bar_stimulus(pd.Stimuli_base, pd.dff_base, pd.TimeCa_base, ...
        pd.Stimuli_drug, pd.dff_drug, pd.TimeCa_drug, ...
        pd.match_idx_local_base, pd.match_idx_local_drug, pd.matched_rois_available);

    if isempty(result.resp_base_block)
        fprintf('  (no block data -- skipping)\n');
        continue;
    end

    t_block = result.t_com_block;
    resp_base_block = result.resp_base_block;
    resp_drug_block = result.resp_drug_block;
    n_matched = size(resp_base_block, 1);
    fprintf('  %d matched cells, block duration %.1fs\n', n_matched, t_block(end));

    %% Heatmap: all directions/repeats concatenated, sorted by baseline peak
    plot_block_heatmap(t_block, resp_base_block, resp_drug_block, n_heatmap_cells, ...
        sprintf('%s_moving_bar_heatmap_block', pairs(p).name), outdir);

    %% 20 highest-SNR cells: baseline vs drug overlay
    peak_base = max(resp_base_block, [], 2);
    std_base = std(resp_base_block, 0, 2, 'omitnan');
    snr = peak_base ./ std_base;
    snr(std_base == 0 | isnan(snr)) = -Inf;
    [sorted_snr, order] = sort(snr, 'descend');
    n_show = min(n_top_cells, n_matched);
    top_cells = order(1:n_show);

    fprintf('  Top %d SNR values: %s\n', n_show, mat2str(round(sorted_snr(1:n_show), 1)));

    %% Per-direction averaged dF/F for those same top-N cells, baseline vs drug
    pre_window_sec = 1; post_window_sec = 6; n_com = 100;
    [Ca_base, t_dir, dir_base] = extract_direction_trials(pd.Stimuli_base, pd.dff_base, pd.TimeCa_base, ...
        'moving_bar', 'bar_orientations', pre_window_sec, post_window_sec, n_com);
    [Ca_drug, ~, dir_drug] = extract_direction_trials(pd.Stimuli_drug, pd.dff_drug, pd.TimeCa_drug, ...
        'moving_bar', 'bar_orientations', pre_window_sec, post_window_sec, n_com);

    [Ca_by_dir_base, dirs_base] = chunk_by_direction(Ca_base(pd.match_idx_local_base, :, :), dir_base);
    [Ca_by_dir_drug, dirs_drug] = chunk_by_direction(Ca_drug(pd.match_idx_local_drug, :, :), dir_drug);

    % Align directions: same physical angle in both, not just same array
    % position (same logic as analyze_direction_stimulus.m).
    common_dirs = intersect(dirs_base, dirs_drug);
    if numel(common_dirs) < numel(union(dirs_base, dirs_drug))
        warning('%s: baseline/drug direction sets differ -- using %d common directions', pairs(p).name, numel(common_dirs));
    end
    [~, ib] = ismember(common_dirs, dirs_base);
    [~, idr] = ismember(common_dirs, dirs_drug);
    Ca_by_dir_base = Ca_by_dir_base(top_cells, :, ib);   % [n_show x n_com x n_dir], top cells only
    Ca_by_dir_drug = Ca_by_dir_drug(top_cells, :, idr);

    %% One figure PER CELL: all directions for that single cell, baseline vs drug
    for k = 1:n_show
        cell_base = squeeze(Ca_by_dir_base(k, :, :));   % [n_com x n_dir]
        cell_drug = squeeze(Ca_by_dir_drug(k, :, :));
        plot_single_cell_direction_grid(t_dir, common_dirs, cell_base, cell_drug, top_cells(k), sorted_snr(k), ...
            sprintf('%s_moving_bar_cell%d_SNR%.1f_by_direction', pairs(p).name, top_cells(k), sorted_snr(k)), outdir);
    end
end

fprintf('\nAll figures saved under: %s\n', outdir);

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

    M = load(roi_match_file);
    base_match_idx_global = M.roiMatchData.allSessionMapping(:, 1);
    drug_match_idx_global = M.roiMatchData.allSessionMapping(:, 2);
    [~, base_match_idx_local] = ismember(base_match_idx_global, selected_roi_idx);
    [~, drug_match_idx_local] = ismember(drug_match_idx_global, drug_selected_roi_idx);
    valid_matches = (base_match_idx_local > 0) & (drug_match_idx_local > 0);
    base_match_idx_local = base_match_idx_local(valid_matches);
    drug_match_idx_local = drug_match_idx_local(valid_matches);
    matched_rois_available = ~isempty(base_match_idx_local);

    pair.Stimuli_base = B.Stimuli;
    pair.Stimuli_drug = D.Stimuli;
    pair.dff_base = base_dff;
    pair.dff_drug = drug_dff;
    pair.TimeCa_base = B.TimeCa;
    pair.TimeCa_drug = D.TimeCa;
    pair.match_idx_local_base = base_match_idx_local;
    pair.match_idx_local_drug = drug_match_idx_local;
    pair.matched_rois_available = matched_rois_available;
end

function plot_block_heatmap(t_block, resp_base, resp_drug, n_cells_show, fig_title, outdir)
    base_max = max(resp_base, [], 2); base_max(isnan(base_max)) = -Inf;
    [~, order] = sort(base_max, 'descend');
    n_cells_show = min(n_cells_show, numel(order));
    order = order(1:n_cells_show);
    Zb = resp_base(order, :);
    Zd = resp_drug(order, :);
    clim = [0, prctile([Zb(:); Zd(:)], 99)];

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 1100 650]);
    colormap(flipud(gray));

    ax1 = subplot(1,2,1);
    imagesc(t_block, 1:size(Zb,1), Zb, clim);
    hold on; xline(0, 'r--', 'LineWidth', 1); hold off;
    xlabel('Time (s)'); ylabel(sprintf('Matched cell # (top %d, sorted by baseline peak)', n_cells_show)); title('Baseline');
    pos1 = get(ax1, 'Position');

    ax2 = subplot(1,2,2);
    imagesc(t_block, 1:size(Zd,1), Zd, clim);
    hold on; xline(0, 'r--', 'LineWidth', 1); hold off;
    xlabel('Time (s)'); title('Drug');
    pos2 = get(ax2, 'Position');
    cb = colorbar; cb.Label.String = 'dF/F';
    set(ax2, 'Position', pos2);
    set(ax1, 'Position', [pos1(1) pos2(2) pos2(3) pos2(4)]);

    sgtitle(sprintf('%s (n=%d matched cells)', fig_title, size(Zb,1)), 'FontWeight', 'bold', 'Interpreter', 'none');
    save_fig_local(fig_title, outdir);
end

function plot_single_cell_direction_grid(t_dir, common_dirs, cell_base, cell_drug, cell_id, snr_val, fig_title, outdir)
    % One subplot per direction, for ONE matched cell -- its own actual
    % trace (no averaging across cells), baseline vs drug overlaid.
    % cell_base/cell_drug: [n_com x n_dir].
    n_dir = numel(common_dirs);
    grid_cols = ceil(sqrt(n_dir));
    grid_rows = ceil(n_dir / grid_cols);
    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [50 50 min(280*grid_cols, 1800) min(260*grid_rows, 1000)]);
    for d = 1:n_dir
        subplot(grid_rows, grid_cols, d);
        hold on;
        plot(t_dir, cell_base(:, d), '-', 'Color', [0.2 0.6 1], 'LineWidth', 1.8, 'DisplayName', 'Baseline');
        plot(t_dir, cell_drug(:, d), '-', 'Color', [1 0.5 0.2], 'LineWidth', 1.8, 'DisplayName', 'Drug');
        xline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        title(sprintf('%g deg', common_dirs(d)), 'FontSize', 10);
        set(gca, 'Box', 'off', 'FontSize', 8); xlim([t_dir(1), t_dir(end)]);
        if d == 1; legend('Location', 'best', 'FontSize', 7); end
        hold off;
    end
    sgtitle(sprintf('%s (cell %d, SNR=%.1f)', fig_title, cell_id, snr_val), 'FontWeight', 'bold', 'Interpreter', 'none');
    save_fig_local(fig_title, outdir);
end

function save_fig_local(fig_title, outdir)
    safe_name = regexprep(fig_title, '[^a-zA-Z0-9]+', '_');
    outfile = fullfile(outdir, [safe_name '.png']);
    saveas(gcf, outfile);
    fprintf('  Saved %s\n', outfile);
end
