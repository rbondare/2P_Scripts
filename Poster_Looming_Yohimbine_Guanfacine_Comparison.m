%% Poster figure: pooled baseline vs yohimbine vs guanfacine looming response
%
% One-off poster-prep script (not part of the 6-stimulus production
% pipeline in Plot_Stimulus_Results.m) -- same category as
% Test_MovingBar_Chunking_Walkthrough.m.
%
% Four recording pairs, same animal (AnimalRB19), two per drug condition:
%   Guanfacine (lowers NA, looming response INCREASES): 260320 (1444/1609), 260305 (1327/1409)
%   Yohimbine  (raises NA,  looming response DECREASES): 260312 (1115/1243), 260310 (1246/1341)
% All 4 baselines are pooled into one reference population (same animal's
% untreated state across different days); the 2 drug recordings per
% condition are pooled into that condition's population.
%
% All pairs' looming stimulus_trial_t = 7.0s, so extract_looming_trials.m
% produces an IDENTICAL t_com = linspace(-2, 5, 100) for all of them --
% their resp_base/resp_drug arrays share the same column structure and can
% be pooled directly with no interpolation (checked with an assert below).
%
% Plane indexing: each roi_match_file's allSessionMapping stores per-plane
% LOCAL suite2p ROI indices, and which plane that numbering refers to is
% NOT guaranteed to be the same index on both sides of a pair (confirmed:
% the 260310 pair needs plane 3 for the baseline file and plane 1 for the
% drug file) or across pairs -- so load_pair auto-detects the correct
% plane per file by matching against the roi_match_file's own recorded
% cellCount, rather than trusting a manually-specified plane index.

clear; clc; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

%% ====================== CONFIGURATION ======================
pairs = struct('baseline_file', {}, 'drug_file', {}, 'roi_match_file', {}, 'condition', {});

pairs(end+1) = struct( ...
    'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1444_preprocessed.mat", ...
    'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1609_preprocessed.mat", ...
    'roi_match_file', "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane2_0320_1444_1609_V1.mat", ...
    'condition', 'guanfacine');

pairs(end+1) = struct( ...
    'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat", ...
    'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat", ...
    'roi_match_file', "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0305_1327_1409_V1.mat", ...
    'condition', 'guanfacine');

pairs(end+1) = struct( ...
    'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat", ...
    'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat", ...
    'roi_match_file', "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane1_0312_1115_1243_V1.mat", ...
    'condition', 'yohimbine');

pairs(end+1) = struct( ...
    'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260310_1246_preprocessed.mat", ...
    'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260310_1341_preprocessed.mat", ...
    'roi_match_file', "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane1_0310_1246_1341_V1.mat", ...
    'condition', 'yohimbine');

poster_outdir = fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', '_poster_figures');
if ~exist(poster_outdir, 'dir'); mkdir(poster_outdir); end

%% ====================== LOAD & ANALYZE ALL PAIRS ======================
resp_base_all = {};
resp_yoh_all = {};
resp_gua_all = {};
t_com = [];

for iP = 1:numel(pairs)
    p = pairs(iP);
    fprintf('\n========== LOADING PAIR %d/%d (%s) ==========\n', iP, numel(pairs), p.condition);
    pd = load_pair(p.baseline_file, p.drug_file, p.roi_match_file);
    result = analyze_looming_stimulus(pd.Stimuli_base, pd.dff_base, pd.TimeCa_base, pd.Triggers_base, ...
        pd.Stimuli_drug, pd.dff_drug, pd.TimeCa_drug, pd.Triggers_drug, ...
        pd.match_idx_local_base, pd.match_idx_local_drug, pd.matched_rois_available);
    fprintf('  %d matched cells\n', size(result.resp_base, 1));

    if isempty(t_com)
        t_com = result.t_com;
    else
        assert(isequal(t_com, result.t_com), 'time axes differ across pairs -- pooling invalid');
    end

    resp_base_all{end+1} = result.resp_base; %#ok<AGROW>
    switch p.condition
        case 'guanfacine'
            resp_gua_all{end+1} = result.resp_drug; %#ok<AGROW>
        case 'yohimbine'
            resp_yoh_all{end+1} = result.resp_drug; %#ok<AGROW>
        otherwise
            error('Unknown condition label: %s', p.condition);
    end
end

%% ====================== POOL ======================
resp_base_pooled = vertcat(resp_base_all{:});
resp_yohimbine_pooled = vertcat(resp_yoh_all{:});
resp_guanfacine_pooled = vertcat(resp_gua_all{:});
n_pooled = size(resp_base_pooled, 1);
n_yohimbine = size(resp_yohimbine_pooled, 1);
n_guanfacine = size(resp_guanfacine_pooled, 1);
fprintf('\nn_pooled=%d, n_yohimbine=%d, n_guanfacine=%d\n', n_pooled, n_yohimbine, n_guanfacine);

%% ====================== PLOT ======================
plot_transient_3cond(t_com, resp_base_pooled, resp_yohimbine_pooled, resp_guanfacine_pooled, ...
    'Looming_Yohimbine_Guanfacine_Comparison', poster_outdir, n_pooled, n_yohimbine, n_guanfacine);

%% ====================== LOCAL FUNCTIONS ======================
function pair = load_pair(baseline_file, drug_file, roi_match_file)
    B = load(baseline_file);
    base_dff_full = B.CaData(1).Ca_dFF;
    base_centroidZ = B.CaData(1).Ca_centroid_voxel(:, 3);

    D = load(drug_file);
    drug_dff_full = D.CaData(1).Ca_dFF;
    drug_centroidZ = D.CaData(1).Ca_centroid_voxel(:, 3);

    M = load(roi_match_file);
    % allSessionMapping stores per-plane-LOCAL suite2p ROI indices (1..cellCount
    % for that plane), not global indices into the multi-plane aggregated array
    % -- confirmed against roiMatchPub.m's generation code, which assigns iCell
    % (1..cellCount within a single session) directly into the mapping. Find
    % which plane in each file that local numbering refers to by matching the
    % match file's own recorded per-session cellCount (robust to plane-index
    % drift between sessions, unlike trusting a shared selected_plane_idx).
    base_cell_count = M.roiMatchData.rois{1}.cellCount;
    drug_cell_count = M.roiMatchData.rois{2}.cellCount;
    selected_roi_idx = find_plane_by_cellcount(base_centroidZ, base_cell_count, roi_match_file, 'baseline');
    drug_selected_roi_idx = find_plane_by_cellcount(drug_centroidZ, drug_cell_count, roi_match_file, 'drug');

    base_dff = base_dff_full(selected_roi_idx, :);
    drug_dff = drug_dff_full(drug_selected_roi_idx, :);

    base_match_idx_local = M.roiMatchData.allSessionMapping(:, 1);
    drug_match_idx_local = M.roiMatchData.allSessionMapping(:, 2);
    valid_matches = (base_match_idx_local >= 1 & base_match_idx_local <= numel(selected_roi_idx)) & ...
        (drug_match_idx_local >= 1 & drug_match_idx_local <= numel(drug_selected_roi_idx));
    base_match_idx_local = base_match_idx_local(valid_matches);
    drug_match_idx_local = drug_match_idx_local(valid_matches);
    matched_rois_available = ~isempty(base_match_idx_local);

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

function roi_idx = find_plane_by_cellcount(centroidZ, target_count, roi_match_file, side_label)
    unique_planes = unique(centroidZ);
    counts = arrayfun(@(z) sum(centroidZ == z), unique_planes);
    match = find(counts == target_count);
    if isempty(match)
        error('load_pair:noPlaneMatch', ...
            ['No plane has %d ROIs to match the %s side of %s.\nPer-plane ROI counts found: %s'], ...
            target_count, side_label, roi_match_file, mat2str(counts'));
    elseif numel(match) > 1
        error('load_pair:ambiguousPlaneMatch', ...
            'Multiple planes (%s) have %d ROIs -- cannot uniquely determine the %s plane for %s.', ...
            mat2str(match'), target_count, side_label, roi_match_file);
    end
    roi_idx = find(centroidZ == unique_planes(match));
end

function plot_transient_3cond(t_com, resp_base, resp_yohimbine, resp_guanfacine, fig_title, outdir, n_pooled, n_yohimbine, n_guanfacine)
    % mean/std/sem computed ACROSS CELLS (dim=1: one row per matched cell,
    % resp_base is already the pooled array) at each of the 100 timepoints
    % -- standard population mean+-SEM convention, same formula as
    % Plot_Stimulus_Results.m's plot_transient. SEM = std / sqrt(n_cells),
    % so a larger pooled n mechanically produces a tighter band -- not a
    % separate effect, just sqrt(n) at work. Diagnostic print below shows
    % the actual numbers at each condition's own peak so this can be
    % checked directly rather than just trusted.
    mean_base = mean(resp_base, 1, 'omitnan');
    sem_base  = std(resp_base, 0, 1, 'omitnan') / sqrt(size(resp_base, 1));
    mean_yoh = mean(resp_yohimbine, 1, 'omitnan');
    sem_yoh  = std(resp_yohimbine, 0, 1, 'omitnan') / sqrt(size(resp_yohimbine, 1));
    mean_gua = mean(resp_guanfacine, 1, 'omitnan');
    sem_gua  = std(resp_guanfacine, 0, 1, 'omitnan') / sqrt(size(resp_guanfacine, 1));

    fprintf('\n  Peak diagnostics (mean +/- SEM at each condition''s own peak timepoint):\n');
    print_peak_stats('Baseline (pooled)', t_com, mean_base, sem_base, resp_base);
    print_peak_stats('Yohimbine', t_com, mean_yoh, sem_yoh, resp_yohimbine);
    print_peak_stats('Guanfacine', t_com, mean_gua, sem_gua, resp_guanfacine);

    color_base = [0.3 0.3 0.3];     % grey
    color_yoh  = [0.20 0.45 0.80];  % cool blue -- yohimbine DECREASES the loom response
    color_gua  = [0.85 0.25 0.15];  % warm red -- guanfacine INCREASES the loom response

    figure('Name', fig_title, 'NumberTitle', 'off', 'Position', [100 100 900 500]);
    hold on;
    x_sh = [t_com, fliplr(t_com)];

    fill(x_sh, [(mean_base+sem_base), fliplr(mean_base-sem_base)], color_base, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_base, '-', 'Color', color_base, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Baseline (n=%d cells)', n_pooled));

    fill(x_sh, [(mean_yoh+sem_yoh), fliplr(mean_yoh-sem_yoh)], color_yoh, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_yoh, '-', 'Color', color_yoh, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Yohimbine (n=%d cells)', n_yohimbine));

    fill(x_sh, [(mean_gua+sem_gua), fliplr(mean_gua-sem_gua)], color_gua, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_gua, '-', 'Color', color_gua, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Guanfacine (n=%d cells)', n_guanfacine));

    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time from onset (s)', 'FontSize', 11);
    ylabel('dF/F (mean \pm SEM)', 'FontSize', 11);
    title(sprintf('Looming response: noradrenergic modulation (baseline n=%d, yohimbine n=%d, guanfacine n=%d)', ...
        n_pooled, n_yohimbine, n_guanfacine), 'FontSize', 12);
    legend('Location', 'best'); set(gca, 'Box', 'off'); xlim([t_com(1), t_com(end)]);
    save_fig_local(fig_title, outdir);
end

function print_peak_stats(label, t_com, mean_trace, sem_trace, resp)
    [pk_val, pk_idx] = max(mean_trace);
    n_cells = size(resp, 1);
    raw_std = std(resp(:, pk_idx), 0, 1, 'omitnan');
    fprintf('    %-20s n=%-5d  t=%.2fs  mean=%.4f  std(across cells)=%.4f  SEM=std/sqrt(n)=%.4f\n', ...
        label, n_cells, t_com(pk_idx), pk_val, raw_std, sem_trace(pk_idx));
end

function save_fig_local(fig_title, outdir)
    safe_name = regexprep(fig_title, '[^a-zA-Z0-9]+', '_');
    outfile = fullfile(outdir, [safe_name '.png']);
    saveas(gcf, outfile);
    fprintf('  Saved %s\n', outfile);
end
