%% Poster figure (diagnostic): per-recording baseline vs drug, individual trials
%
% Same recording pairs, same cell-matching logic as
% Poster_Looming_Yohimbine_Guanfacine_Comparison.m, but instead of
% collapsing everything down to one pooled 3-condition summary, this shows
% each recording pair on its own row: baseline subplot | drug subplot, each
% with every individual trial as a thin line (population mean across
% matched cells, per trial) plus the across-trial mean as a bold line.
% Purpose: sanity-check what's actually driving the pooled comparison --
% trial-to-trial consistency, obvious outlier trials, etc. -- before
% trusting the pooled figure. All panels share one y-axis for direct
% visual comparison across recordings/conditions.
%
% Unlike analyze_looming_stimulus.m (which uses the per-cell trial MEDIAN),
% here each trial's line is the mean across matched cells for that single
% trial, so trial count and per-trial variability are directly visible.
%
% Multiple planes per recording: list more than one roi_match_file under
% a pair's roi_match_files to pool matched cells across planes for that
% recording (e.g. if separate ROI-matching was run per plane). Comment out
% any entry to drop that plane -- each is resolved independently via its
% own cellCount, including its own TimeCa row, so order doesn't matter.

clear; clc; close all;

addpath(fullfile(fileparts(mfilename('fullpath')), 'analysis_functions'));

%% ====================== CONFIGURATION ======================
pairs = struct('baseline_file', {}, 'drug_file', {}, 'roi_match_files', {}, 'condition', {}, 'label', {});

roi_match_files_1 = {
    "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane2_0320_1444_1609_V1.mat"
    "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0320_1444_1609_V1.mat"
};
pairs(end+1) = struct( ...
    'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1444_preprocessed.mat", ...
    'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1609_preprocessed.mat", ...
    'roi_match_files', {roi_match_files_1}, ...
    'condition', 'guanfacine', 'label', '260320');

% --- 260303 DISABLED: baseline file AnimalRB19_260303_1546_preprocessed.mat is
% incompletely aggregated (only AnimalInfo/SI_Info/TimeCa present; missing
% CaData/Stimuli/Triggers). Re-run aggregation for this recording, then
% uncomment. (Drug file _1705_ is fine; only the baseline is broken.)
% roi_match_files_2 = {
%     "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane1_0303_1546_1705_V1.mat"
% };
% pairs(end+1) = struct( ...
%     'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260303_1546_preprocessed.mat", ...
%     'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260303_1705_preprocessed.mat", ...
%     'roi_match_files', {roi_match_files_2}, ...
%     'condition', 'guanfacine', 'label', '260303');

roi_match_files_3 = {
    "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0305_1327_1409_V1.mat"
    % "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_planeX_0305_1327_1409_V1.mat"
};
pairs(end+1) = struct( ...
    'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat", ...
    'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat", ...
    'roi_match_files', {roi_match_files_3}, ...
    'condition', 'guanfacine', 'label', '260305');

roi_match_files_4 = {
    "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane1_0312_1115_1243_V1.mat"
    "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0312_1115_1243_V1.mat"
};
pairs(end+1) = struct( ...
    'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat", ...
    'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat", ...
    'roi_match_files', {roi_match_files_4}, ...
    'condition', 'yohimbine', 'label', '260312');

roi_match_files_5 = {
    "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane1_0310_1246_1341_V1.mat"
    % "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_planeX_0310_1246_1341_V1.mat"
};
pairs(end+1) = struct( ...
    'baseline_file', "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260310_1246_preprocessed.mat", ...
    'drug_file',     "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260310_1341_preprocessed.mat", ...
    'roi_match_files', {roi_match_files_5}, ...
    'condition', 'yohimbine', 'label', '260310');

poster_outdir = fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', '_poster_figures');
if ~exist(poster_outdir, 'dir'); mkdir(poster_outdir); end

c_base = [0.3 0.3 0.3];  % grey -- baseline mean line

%% ====================== LOAD & EXTRACT ALL PAIRS ======================
n_pairs = numel(pairs);
pair_data = cell(1, n_pairs);
global_ymin = Inf;
global_ymax = -Inf;

for iP = 1:n_pairs
    p = pairs(iP);
    fprintf('\n========== PAIR %d/%d (%s, %s) ==========\n', iP, n_pairs, p.label, p.condition);
    [Ca_base, Ca_drug, t_com, n_matched] = extract_pair_matched_trials( ...
        p.baseline_file, p.drug_file, p.roi_match_files);
    fprintf('  %d matched cells (pooled across %d plane(s))\n', n_matched, numel(p.roi_match_files));

    if isempty(Ca_base) || isempty(Ca_drug) || n_matched == 0
        fprintf('  Skipping (missing looming data or no matched cells)\n');
        pair_data{iP} = [];
        continue;
    end

    % Population trace per trial (mean across matched cells): [n_trials x n_com]
    pop_base_trials = squeeze(mean(Ca_base, 1, 'omitnan'))';
    pop_drug_trials = squeeze(mean(Ca_drug, 1, 'omitnan'))';

    pair_data{iP} = struct('t_com', t_com, 'n_matched', n_matched, ...
        'pop_base_trials', pop_base_trials, 'pop_drug_trials', pop_drug_trials);
    global_ymin = min([global_ymin, pop_base_trials(:)', pop_drug_trials(:)']);
    global_ymax = max([global_ymax, pop_base_trials(:)', pop_drug_trials(:)']);
end

pad = 0.05 * (global_ymax - global_ymin);
global_ylim = [global_ymin - pad, global_ymax + pad];

%% ====================== FIGURE ======================
fig = figure('Name', 'Looming_PerTrial_BaselineVsDrug', 'NumberTitle', 'off', ...
    'Position', [40 20 1300 300 * n_pairs]);
all_axes = gobjects(0);

for iP = 1:n_pairs
    p = pairs(iP);
    pd = pair_data{iP};
    if isempty(pd); continue; end

    t_com = pd.t_com;
    pop_base_trials = pd.pop_base_trials;
    pop_drug_trials = pd.pop_drug_trials;
    n_trials_base = size(pop_base_trials, 1);
    n_trials_drug = size(pop_drug_trials, 1);

    % Trial lines colored light -> dark = trial 1 -> N (within each
    % subplot's own hue) so any trend of peak size with loom repetition
    % number is visible.
    c_drug = strcmp(p.condition, 'guanfacine') * [0.85 0.25 0.15] + ...
             strcmp(p.condition, 'yohimbine') * [0.20 0.45 0.80];
    trial_clrs_base = trial_gradient(c_base, n_trials_base);
    trial_clrs_drug = trial_gradient(c_drug, n_trials_drug);

    % ---- Baseline subplot ----
    ax_b = subplot(n_pairs, 2, 2*iP - 1);
    hold on;
    for k = 1:n_trials_base
        plot(t_com, pop_base_trials(k, :), 'Color', trial_clrs_base(k,:), 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
    plot(t_com, mean(pop_base_trials, 1, 'omitnan'), 'Color', c_base, 'LineWidth', 2.5, ...
        'DisplayName', 'Mean');
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    hold off;
    xlim([t_com(1), t_com(end)]);
    ylabel('dF/F');
    title(sprintf('%s baseline  (n=%d cells, %d trials)', p.label, pd.n_matched, n_trials_base));
    if iP == n_pairs; xlabel('Time from onset (s)'); end
    if iP == 1; legend('Location', 'northeast', 'FontSize', 8); end
    all_axes(end+1) = ax_b; %#ok<AGROW>

    % ---- Drug subplot ----
    ax_d = subplot(n_pairs, 2, 2*iP);
    hold on;
    for k = 1:n_trials_drug
        plot(t_com, pop_drug_trials(k, :), 'Color', trial_clrs_drug(k,:), 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
    plot(t_com, mean(pop_drug_trials, 1, 'omitnan'), 'Color', c_drug, 'LineWidth', 2.5, ...
        'DisplayName', 'Mean');
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    hold off;
    xlim([t_com(1), t_com(end)]);
    title(sprintf('%s drug (%s)  (n=%d cells, %d trials)', p.label, p.condition, pd.n_matched, n_trials_drug));
    if iP == n_pairs; xlabel('Time from onset (s)'); end
    if iP == 1; legend('Location', 'northeast', 'FontSize', 8); end
    all_axes(end+1) = ax_d; %#ok<AGROW>
end

set(all_axes, 'YLim', global_ylim);

sgtitle('Looming: baseline vs drug, individual trials per recording', 'FontWeight', 'bold', 'FontSize', 13);

outfile = fullfile(poster_outdir, 'Looming_PerTrial_BaselineVsDrug.png');
saveas(fig, outfile);
fprintf('\nSaved %s\n', outfile);

%% ====================== HEATMAP FIGURE ======================
% Same per-trial data as the line figure above, drawn as trial x time
% heatmaps: baseline | drug, one row per recording pair. Rows are trials
% ordered first (top) -> last (bottom); columns are time from onset. Each
% pixel is the population mean (across matched cells) for that single trial,
% so trial-to-trial evolution and any baseline-vs-drug shift read at a
% glance. All panels share one colour scale (diverging, centred on
% dF/F = 0) so trials, conditions and recordings are directly comparable.

% Global symmetric colour limits from the 1st/99th percentile of every
% value (robust to a few outlier pixels), centred on zero so white always
% marks dF/F = 0 rather than the middle of the data range.
all_vals = [];
for iP = 1:n_pairs
    pd = pair_data{iP};
    if isempty(pd); continue; end
    all_vals = [all_vals; pd.pop_base_trials(:); pd.pop_drug_trials(:)]; %#ok<AGROW>
end
cmax = max(abs(prctile(all_vals, [1 99])));
if isempty(cmax) || ~(cmax > 0); cmax = max(abs(all_vals)); end  % degenerate fallback
global_clim = [-cmax, cmax];

fig_hm = figure('Name', 'Looming_PerTrial_BaselineVsDrug_Heatmap', 'NumberTitle', 'off', ...
    'Position', [40 20 1300 300 * n_pairs]);
hm_axes = gobjects(0);

for iP = 1:n_pairs
    p = pairs(iP);
    pd = pair_data{iP};
    if isempty(pd); continue; end

    t_com = pd.t_com;
    pop_base_trials = pd.pop_base_trials;
    pop_drug_trials = pd.pop_drug_trials;
    n_trials_base = size(pop_base_trials, 1);
    n_trials_drug = size(pop_drug_trials, 1);

    % ---- Baseline heatmap ----
    ax_b = subplot(n_pairs, 2, 2*iP - 1);
    imagesc(t_com, 1:n_trials_base, pop_base_trials);
    set(ax_b, 'YDir', 'reverse');   % trial 1 on top, trial N at the bottom
    clim(global_clim);
    hold on; xline(0, 'k--', 'LineWidth', 1); hold off;
    ylabel('Trial (1 = first)');
    if n_trials_base <= 20; yticks(1:n_trials_base); end
    title(sprintf('%s baseline  (n=%d cells, %d trials)', p.label, pd.n_matched, n_trials_base));
    if iP == n_pairs; xlabel('Time from onset (s)'); end
    hm_axes(end+1) = ax_b; %#ok<AGROW>

    % ---- Drug heatmap ----
    ax_d = subplot(n_pairs, 2, 2*iP);
    imagesc(t_com, 1:n_trials_drug, pop_drug_trials);
    set(ax_d, 'YDir', 'reverse');
    clim(global_clim);
    hold on; xline(0, 'k--', 'LineWidth', 1); hold off;
    ylabel('Trial (1 = first)');
    if n_trials_drug <= 20; yticks(1:n_trials_drug); end
    title(sprintf('%s drug (%s)  (n=%d cells, %d trials)', p.label, p.condition, pd.n_matched, n_trials_drug));
    if iP == n_pairs; xlabel('Time from onset (s)'); end
    hm_axes(end+1) = ax_d; %#ok<AGROW>
end

colormap(fig_hm, diverging_colormap());

% One shared colourbar on the right margin. colorbar() shrinks its host
% axis, so capture and restore that axis's Position to keep every panel the
% same width (matching the pipeline's equal-panel-width convention).
if ~isempty(hm_axes)
    host = hm_axes(end);
    pos_host = get(host, 'Position');
    cb = colorbar(host);
    set(host, 'Position', pos_host);
    cb.Position = [0.935 0.11 0.012 0.8];
    cb.Label.String = '\DeltaF/F (population mean)';
end

sgtitle('Looming: baseline vs drug, per-trial heatmaps (trial 1 = top)', 'FontWeight', 'bold', 'FontSize', 13);

outfile_hm = fullfile(poster_outdir, 'Looming_PerTrial_BaselineVsDrug_Heatmap.png');
saveas(fig_hm, outfile_hm);
fprintf('\nSaved %s\n', outfile_hm);

%% ====================== LOCAL FUNCTIONS ======================
function [Ca_base_all, Ca_drug_all, t_com, n_matched] = extract_pair_matched_trials(baseline_file, drug_file, roi_match_files)
    Ca_base_all = [];
    Ca_drug_all = [];
    t_com       = [];
    n_matched   = 0;

    B = load(baseline_file);
    D = load(drug_file);

    % Guard against incompletely-preprocessed aggregated files (e.g. a run
    % that only saved AnimalInfo/SI_Info/TimeCa before the calcium/stimulus
    % pipeline finished). Skip the whole pair with a clear warning rather
    % than crashing the entire figure on the first missing field.
    required = {'CaData', 'TimeCa', 'Stimuli', 'Triggers'};
    missB = required(~isfield(B, required));
    missD = required(~isfield(D, required));
    if ~isempty(missB) || ~isempty(missD)
        if ~isempty(missB)
            warning('extract_pair:incompleteFile', ...
                'Incomplete baseline file (missing: %s) -- skipping pair:\n    %s', ...
                strjoin(missB, ', '), baseline_file);
        end
        if ~isempty(missD)
            warning('extract_pair:incompleteFile', ...
                'Incomplete drug file (missing: %s) -- skipping pair:\n    %s', ...
                strjoin(missD, ', '), drug_file);
        end
        return;
    end

    base_dff_full = B.CaData(1).Ca_dFF;
    base_centroidZ = B.CaData(1).Ca_centroid_voxel(:, 3);
    drug_dff_full = D.CaData(1).Ca_dFF;
    drug_centroidZ = D.CaData(1).Ca_centroid_voxel(:, 3);

    for iM = 1:numel(roi_match_files)
        roi_match_file = roi_match_files{iM};
        M = load(roi_match_file);
        % allSessionMapping stores per-plane-LOCAL suite2p ROI indices (1..cellCount
        % for that plane), not global indices into the multi-plane aggregated array
        % -- confirmed against roiMatchPub.m's generation code, which assigns iCell
        % (1..cellCount within a single session) directly into the mapping. Find
        % which plane in each file that local numbering refers to by matching the
        % match file's own recorded per-session cellCount (robust to plane-index
        % drift between sessions and across different roi_match_files/planes).
        base_cell_count = M.roiMatchData.rois{1}.cellCount;
        drug_cell_count = M.roiMatchData.rois{2}.cellCount;
        [base_roi_idx, base_plane_idx] = find_plane_by_cellcount(base_centroidZ, base_cell_count, roi_match_file, 'baseline');
        [drug_roi_idx, drug_plane_idx] = find_plane_by_cellcount(drug_centroidZ, drug_cell_count, roi_match_file, 'drug');

        base_match_idx_local = M.roiMatchData.allSessionMapping(:, 1);
        drug_match_idx_local = M.roiMatchData.allSessionMapping(:, 2);
        valid_matches = (base_match_idx_local >= 1 & base_match_idx_local <= numel(base_roi_idx)) & ...
            (drug_match_idx_local >= 1 & drug_match_idx_local <= numel(drug_roi_idx));
        base_match_idx_local = base_match_idx_local(valid_matches);
        drug_match_idx_local = drug_match_idx_local(valid_matches);
        if isempty(base_match_idx_local); continue; end

        base_dff_plane = base_dff_full(base_roi_idx, :);
        drug_dff_plane = drug_dff_full(drug_roi_idx, :);
        TimeCa_base_plane = B.TimeCa(base_plane_idx, :);
        TimeCa_drug_plane = D.TimeCa(drug_plane_idx, :);

        [Ca_base_plane, t_com_plane] = extract_looming_trials(B.Stimuli, base_dff_plane, TimeCa_base_plane, B.Triggers);
        [Ca_drug_plane, ~]           = extract_looming_trials(D.Stimuli, drug_dff_plane, TimeCa_drug_plane, D.Triggers);
        if isempty(Ca_base_plane) || isempty(Ca_drug_plane); continue; end

        if isempty(t_com)
            t_com = t_com_plane;
        else
            assert(isequal(t_com, t_com_plane), 'time axes differ across planes -- pooling invalid');
        end

        Ca_base_all = cat(1, Ca_base_all, Ca_base_plane(base_match_idx_local, :, :));
        Ca_drug_all = cat(1, Ca_drug_all, Ca_drug_plane(drug_match_idx_local, :, :));
    end

    n_matched = size(Ca_base_all, 1);
end

function [roi_idx, plane_idx] = find_plane_by_cellcount(centroidZ, target_count, roi_match_file, side_label)
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
    plane_idx = match;
    roi_idx = find(centroidZ == unique_planes(match));
end

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

function cmap = diverging_colormap()
    % Blue -> white -> red diverging colormap with white at the centre, so
    % (with a clim symmetric about 0) white always marks dF/F = 0, negative
    % (suppression) reads blue and positive (activation) reads red. Matches
    % the diverging-map convention used by the pipeline's z-scored heatmaps.
    n    = 256;
    half = floor(n / 2);
    blue  = [0.23 0.30 0.75];
    red   = [0.80 0.10 0.15];
    white = [1 1 1];
    top = [linspace(blue(1),  white(1), half)', ...
           linspace(blue(2),  white(2), half)', ...
           linspace(blue(3),  white(3), half)'];
    bot = [linspace(white(1), red(1),   n - half)', ...
           linspace(white(2), red(2),   n - half)', ...
           linspace(white(3), red(3),   n - half)'];
    cmap = [top; bot];
end
