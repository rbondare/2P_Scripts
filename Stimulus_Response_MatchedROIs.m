%% MATCHED ROI STIMULUS RESPONSE COMPARISON
%  V1 of the script to compare baseline vs drug using ROI matches from ROIMatchPub.
% Inputs:
%   1) Two preprocessed files (baseline and drug)
%   2) ROIMatch file containing roiMatchData.allSessionMapping
% Assumptions:
%   - allSessionMapping columns correspond to session order used in ROIMatchPub
%   - Mapping indices are "valid cell" indices within a single plane
%   - CaData(...).Ca_centroid_voxel(:,3) stores plane index (1-based)

clc;
close all;
clear;

%% USER SETTINGS
% Preprocessed recordings
baseline_file = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat";
drug_file     = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat";

% ROI match file saved by ROIMatchPub (contains roiMatchData)
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_V2.mat";

% Session columns in roiMatchData.allSessionMapping
% Use 1 for baseline if it was first experiment loaded in ROIMatchPub
baseline_session_col = 1;
drug_session_col = 2;

% Plane in suite2p convention (plane0, plane1, ...)
selected_plane_s2p = 0;

% Calcium type index (usually 1 in your current pipeline)
ca_type = 1;

% Optional stimulus filter. Keep {} to use all stimulus entries.
% Example: {'moving_bar', 'grating'}
stimulus_types = {'spontaneous'};

% Heatmap frame selection (edit to choose which frames to display)
plot_frame_start = 1;       % starting frame index to display
plot_num_frames  = 1000;    % number of frames to display (choose 1000)

% Color limits for heatmaps
heatmap_clim = [0 5];       % colorbar axis limits

%% LOAD DATA
B = load(baseline_file);
D = load(drug_file);
M = load(roi_match_file);

if ~isfield(M, 'roiMatchData') || ~isfield(M.roiMatchData, 'allSessionMapping')
    error('roi_match_file must contain roiMatchData.allSessionMapping');
end

required_fields = {'CaData', 'TimeCa', 'Stimuli'};
for i = 1:numel(required_fields)
    f = required_fields{i};
    if ~isfield(B, f)
        error('Baseline file missing field: %s', f);
    end
    if ~isfield(D, f)
        error('Drug file missing field: %s', f);
    end
end

if numel(B.CaData) < ca_type || numel(D.CaData) < ca_type
    error('ca_type=%d not found in one or both files.', ca_type);
end

if isempty(B.CaData(ca_type).Ca_dFF) || isempty(D.CaData(ca_type).Ca_dFF)
    error('CaData(%d).Ca_dFF is empty in one or both files.', ca_type);
end

if isempty(B.CaData(ca_type).Ca_centroid_voxel) || isempty(D.CaData(ca_type).Ca_centroid_voxel)
    error('Ca_centroid_voxel is required for plane selection but is missing/empty.');
end

all_mapping = M.roiMatchData.allSessionMapping;
if baseline_session_col > size(all_mapping, 2) || drug_session_col > size(all_mapping, 2)
    error('Session column exceeds allSessionMapping width (%d).', size(all_mapping, 2));
end

%% PREPARE PLANE-SPECIFIC ROW MAPS
selected_plane_internal = selected_plane_s2p + 1; % suite2p plane0 -> internal 1

base_centroid = B.CaData(ca_type).Ca_centroid_voxel;
drug_centroid = D.CaData(ca_type).Ca_centroid_voxel;

if size(base_centroid, 2) < 3 || size(drug_centroid, 2) < 3
    error('Ca_centroid_voxel must have at least 3 columns [x y plane].');
end

base_plane_rows = find(base_centroid(:, 3) == selected_plane_internal);
drug_plane_rows = find(drug_centroid(:, 3) == selected_plane_internal);

if isempty(base_plane_rows) || isempty(drug_plane_rows)
    error('Selected plane not found in one or both recordings.');
end

%% EXTRACT MATCHED IDS FOR BASELINE/DRUG
base_ids = all_mapping(:, baseline_session_col);
drug_ids = all_mapping(:, drug_session_col);

valid_pairs = base_ids > 0 & drug_ids > 0;
base_ids = base_ids(valid_pairs);
drug_ids = drug_ids(valid_pairs);

if isempty(base_ids)
    error('No matched ROI pairs found for selected session columns.');
end

% Keep only pairs that are valid for selected plane row counts.
in_range = base_ids <= numel(base_plane_rows) & drug_ids <= numel(drug_plane_rows);
base_ids = base_ids(in_range);
drug_ids = drug_ids(in_range);

if isempty(base_ids)
    error('No matched ROI pairs remain after filtering to selected plane.');
end

base_row_idx = base_plane_rows(base_ids);
drug_row_idx = drug_plane_rows(drug_ids);

%% GET MATCHED DFF TRACES
base_dff = B.CaData(ca_type).Ca_dFF(base_row_idx, :);
drug_dff = D.CaData(ca_type).Ca_dFF(drug_row_idx, :);

base_time = get_time_vector(B.TimeCa);
drug_time = get_time_vector(D.TimeCa);

fprintf('Matched ROI pairs used: %d\n', size(base_dff, 1));
fprintf('Selected plane: suite2p plane%d\n', selected_plane_s2p);
fprintf('Baseline time range: %.2f - %.2f s (%d frames)\n', base_time(1), base_time(end), numel(base_time));
fprintf('Drug time range:     %.2f - %.2f s (%d frames)\n', drug_time(1), drug_time(end), numel(drug_time));

% Determine frame indices to plot (use same window for both recordings)
nFramesBase = size(base_dff, 2);
nFramesDrug = size(drug_dff, 2);
common_max = min(nFramesBase, nFramesDrug);
plot_frame_end = min(plot_frame_start + plot_num_frames - 1, common_max);
plot_idx = plot_frame_start:plot_frame_end;

fprintf('Plotting frames %d:%d (max %d)\n', plot_frame_start, plot_frame_end, common_max);

%% VALIDATE STIMULUS DATA
% Check if Stimuli contain valid timing and type fields
if ~isempty(B.Stimuli)
    has_type = false;
    has_timing = false;
    for i = 1:min(3, numel(B.Stimuli))
        if isfield(B.Stimuli(i), 'type') && ~isempty(B.Stimuli(i).type)
            has_type = true;
        end
        if isfield(B.Stimuli(i), 'TimeStimulusFrame') && ~isempty(B.Stimuli(i).TimeStimulusFrame)
            has_timing = true;
        end
    end
    if has_type
        fprintf('Stimulus data validated: contains type field.\n');
    end
    if has_timing
        fprintf('Stimulus data contains TimeStimulusFrame (frame-based timing).\n');
    end
else
    warning('No stimulus data in baseline file.');
end

%% PLOT 1: MEAN TRACES WITH SEM SHADING
figure('Name', 'Matched ROI Mean Traces with SEM', 'NumberTitle', 'off', 'Position', [100 100 1200 600]);

subplot(1,2,1);
    trace_base = mean(base_dff, 1);
    sem_base = std(base_dff, [], 1) / sqrt(size(base_dff, 1));
    fill([base_time fliplr(base_time)], [trace_base + sem_base fliplr(trace_base - sem_base)], ...
        'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on;
    plot(base_time, trace_base, 'r-', 'LineWidth', 2, 'Label', 'Mean');
    xlabel('Time (s)'); ylabel('Mean dF/F');
    title(sprintf('Baseline: Mean ± SEM (%d matched ROIs)', size(base_dff, 1)));
    legend;
    grid on;

subplot(1,2,2);
    trace_drug = mean(drug_dff, 1);
    sem_drug = std(drug_dff, [], 1) / sqrt(size(drug_dff, 1));
    fill([drug_time fliplr(drug_time)], [trace_drug + sem_drug fliplr(trace_drug - sem_drug)], ...
        'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on;
    plot(drug_time, trace_drug, 'b-', 'LineWidth', 2, 'Label', 'Mean');
    xlabel('Time (s)'); ylabel('Mean dF/F');
    title(sprintf('Drug: Mean ± SEM (%d matched ROIs)', size(drug_dff, 1)));
    legend;
    grid on;

%% PLOT 2: INDIVIDUAL ROI TRACES (subset for clarity)
% Show ~20 randomly selected individual traces overlaid
n_show = min(20, size(base_dff, 1));
show_idx = randperm(size(base_dff, 1), n_show);

figure('Name', 'Individual ROI Traces (Sample)', 'NumberTitle', 'off', 'Position', [100 100 1200 600]);

subplot(1,2,1);
    hold on;
    for i = show_idx
        plot(base_time, base_dff(i, :), 'r-', 'Alpha', 0.4, 'LineWidth', 0.5);
    end
    plot(base_time, trace_base, 'k-', 'LineWidth', 2, 'Label', 'Population mean');
    xlabel('Time (s)'); ylabel('dF/F');
    title(sprintf('Baseline: %d of %d ROIs shown', n_show, size(base_dff, 1)));
    legend;
    grid on;

subplot(1,2,2);
    hold on;
    for i = show_idx
        plot(drug_time, drug_dff(i, :), 'b-', 'Alpha', 0.4, 'LineWidth', 0.5);
    end
    plot(drug_time, trace_drug, 'k-', 'LineWidth', 2, 'Label', 'Population mean');
    xlabel('Time (s)'); ylabel('dF/F');
    title(sprintf('Drug: %d of %d ROIs shown', n_show, size(drug_dff, 1)));
    legend;
    grid on;

%% PLOT 3: HEATMAPS (sorted by response magnitude)
% Sort by mean response for better visualization
[~, sort_idx_base] = sort(base_metrics.mean_stim, 'descend');
[~, sort_idx_drug] = sort(drug_metrics.mean_stim, 'descend');

figure('Name', 'Heatmaps Sorted by Response', 'NumberTitle', 'off', 'Position', [100 100 1400 600]);

subplot(1,2,1);
    imagesc(base_dff(sort_idx_base, plot_idx));
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
    c = colorbar;
    caxis(heatmap_clim);
    ylabel(c, 'dF/F');
    xlabel('Frame');
    ylabel('ROI (sorted by response)');
    title(sprintf('Baseline Heatmap (frames %d-%d, n=%d)', plot_frame_start, plot_frame_end, size(base_dff, 1)));

subplot(1,2,2);
    imagesc(drug_dff(sort_idx_drug, plot_idx));
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
    c = colorbar;
    caxis(heatmap_clim);
    ylabel(c, 'dF/F');
    xlabel('Frame');
    ylabel('ROI (sorted by response)');
    title(sprintf('Drug Heatmap (frames %d-%d, n=%d)', plot_frame_start, plot_frame_end, size(drug_dff, 1)));

%% PLOT 4: RESPONSE AMPLITUDE DISTRIBUTIONS
figure('Name', 'Response Distributions', 'NumberTitle', 'off', 'Position', [100 100 1400 800]);

% Mean stimulus response
subplot(2,3,1);
    valid_mean = ~isnan(base_metrics.mean_stim) & ~isnan(drug_metrics.mean_stim);
    boxplot([base_metrics.mean_stim(valid_mean)', drug_metrics.mean_stim(valid_mean)'], ...
        'Labels', {'Baseline', 'Drug'});
    ylabel('Mean dF/F during stimulus');
    title('Response Amplitude');
    grid on;

subplot(2,3,2);
    hold on;
    histogram(base_metrics.mean_stim(valid_mean), 15, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'k');
    histogram(drug_metrics.mean_stim(valid_mean), 15, 'FaceColor', 'b', 'FaceAlpha', 0.6, 'EdgeColor', 'k');
    xlabel('Mean dF/F');
    ylabel('Number of ROIs');
    legend('Baseline', 'Drug');
    title('Distribution of Mean Response');
    grid on;

% Evoked response (stimulus - baseline)
subplot(2,3,3);
    valid_delta = ~isnan(base_metrics.delta_stim) & ~isnan(drug_metrics.delta_stim);
    violin_data = [base_metrics.delta_stim(valid_delta)', drug_metrics.delta_stim(valid_delta)'];
    plot_violin(violin_data, {'Baseline', 'Drug'});
    ylabel('Evoked Response (ΔdF/F)');
    title('Evoked Response Magnitude');
    grid on;

% Response reliability (% responsive ROIs)
subplot(2,3,4);
    resp_threshold = 0.5; % dF/F threshold for "responsive"
    base_resp_pct = 100 * sum(base_metrics.mean_stim > resp_threshold) / size(base_dff, 1);
    drug_resp_pct = 100 * sum(drug_metrics.mean_stim > resp_threshold) / size(drug_dff, 1);
    bar([base_resp_pct, drug_resp_pct], 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'LineWidth', 2);
    set(gca, 'XTickLabel', {'Baseline', 'Drug'});
    ylabel('% Responsive ROIs');
    title(sprintf('Responsiveness (threshold = %.1f dF/F)', resp_threshold));
    ylim([0 100]);
    grid on;

% Cumulative distribution
subplot(2,3,5);
    [sorted_base, idx_base] = sort(base_metrics.mean_stim(valid_mean));
    [sorted_drug, idx_drug] = sort(drug_metrics.mean_stim(valid_mean));
    cdf_x = linspace(0, max([sorted_base; sorted_drug]), 100);
    plot(cdf_x, ksdensity(sorted_base, cdf_x, 'cumulative'), 'r-', 'LineWidth', 2, 'Label', 'Baseline');
    hold on;
    plot(cdf_x, ksdensity(sorted_drug, cdf_x, 'cumulative'), 'b-', 'LineWidth', 2, 'Label', 'Drug');
    xlabel('Mean dF/F');
    ylabel('Cumulative Probability');
    title('Cumulative Distribution of Responses');
    legend;
    grid on;

% Modulation index (drug vs baseline)
subplot(2,3,6);
    mod_idx = (drug_metrics.mean_stim(valid_mean) - base_metrics.mean_stim(valid_mean)) ./ ...
              (abs(drug_metrics.mean_stim(valid_mean)) + abs(base_metrics.mean_stim(valid_mean)) + eps);
    histogram(mod_idx, 20, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
    xlabel('Modulation Index [(Drug-Base)/(|Drug|+|Base|)]');
    ylabel('Number of ROIs');
    title('Response Modulation by Drug');
    hold on;
    yl = ylim;
    plot([0 0], yl, 'k--', 'LineWidth', 2, 'Label', 'No modulation');
    hold off;
    grid on;

%% PLOT 5: PAIRED RESPONSE COMPARISON SCATTER
figure('Name', 'Paired Response Scatter', 'NumberTitle', 'off', 'Position', [100 100 1200 500]);

subplot(1,2,1);
    valid_idx = ~isnan(base_metrics.mean_stim) & ~isnan(drug_metrics.mean_stim);
    scatter(base_metrics.mean_stim(valid_idx), drug_metrics.mean_stim(valid_idx), ...
        50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
    hold on;
    lims = [min([base_metrics.mean_stim(valid_idx); drug_metrics.mean_stim(valid_idx)]), ...
            max([base_metrics.mean_stim(valid_idx); drug_metrics.mean_stim(valid_idx)])];
    plot(lims, lims, 'k--', 'LineWidth', 1, 'Label', 'No change');
    xlabel('Baseline Mean dF/F');
    ylabel('Drug Mean dF/F');
    title(sprintf('Paired ROI Response (n=%d)', sum(valid_idx)));
    legend;
    axis equal;
    grid on;

subplot(1,2,2);
    valid_delta = ~isnan(base_metrics.delta_stim) & ~isnan(drug_metrics.delta_stim);
    scatter(base_metrics.delta_stim(valid_delta), drug_metrics.delta_stim(valid_delta), ...
        50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
    hold on;
    lims = [min([base_metrics.delta_stim(valid_delta); drug_metrics.delta_stim(valid_delta)]), ...
            max([base_metrics.delta_stim(valid_delta); drug_metrics.delta_stim(valid_delta)])];
    plot(lims, lims, 'k--', 'LineWidth', 1, 'Label', 'No change');
    xlabel('Baseline Evoked Response (ΔdF/F)');
    ylabel('Drug Evoked Response (ΔdF/F)');
    title(sprintf('Paired ROI Evoked Response (n=%d)', sum(valid_delta)));
    legend;
    axis equal;
    grid on;

%% RESPONSE METRICS 
% Uses mean dF/F during stimulus minus pre-stimulus baseline.
fprintf('\n--- Computing response metrics for BASELINE ---\\n');
base_metrics = compute_stimulus_metrics(base_dff, base_time, B.Stimuli, stimulus_types);
fprintf('--- Computing response metrics for DRUG ---\\n');
drug_metrics = compute_stimulus_metrics(drug_dff, drug_time, D.Stimuli, stimulus_types);

% Validation: check if metrics contain valid (non-NaN) values
base_valid = ~isnan(base_metrics.mean_stim);
drug_valid = ~isnan(drug_metrics.mean_stim);
n_valid = sum(base_valid & drug_valid);
fprintf('\\nValid paired responses: %d / %d ROIs\\n', n_valid, size(base_dff, 1));

if n_valid < size(base_dff, 1) * 0.5
    warning('Less than 50%% of ROIs have valid responses. Check stimulus timing alignment.');
end

%% COMPREHENSIVE STATISTICS
fprintf('\n');
fprintf('====================================================================\n');
fprintf('STIMULUS RESPONSE ANALYSIS: Baseline vs Drug Manipulation\n');
fprintf('====================================================================\n');
fprintf('Analysis Date: %s\n', datetime('now'));
fprintf('Matched ROI Pairs: %d\n', size(base_dff, 1));
fprintf('Selected Plane: suite2p plane%d\n', selected_plane_s2p);
fprintf('Stimulus Type(s): %s\n', strjoin(stimulus_types, ', '));
fprintf('Baseline frames: %d | Drug frames: %d\n', size(base_dff, 2), size(drug_dff, 2));

% Response metrics summary
valid_mean = ~isnan(base_metrics.mean_stim) & ~isnan(drug_metrics.mean_stim);
valid_delta = ~isnan(base_metrics.delta_stim) & ~isnan(drug_metrics.delta_stim);
n_valid_mean = sum(valid_mean);
n_valid_delta = sum(valid_delta);

fprintf('\n--- RESPONSE AMPLITUDE ---\n');
fprintf('Baseline Mean dF/F:  %.4f ± %.4f (mean ± std, n=%d)\n', ...
    mean(base_metrics.mean_stim(valid_mean)), std(base_metrics.mean_stim(valid_mean)), n_valid_mean);
fprintf('Drug Mean dF/F:      %.4f ± %.4f (mean ± std, n=%d)\n', ...
    mean(drug_metrics.mean_stim(valid_mean)), std(drug_metrics.mean_stim(valid_mean)), n_valid_mean);
fprintf('Median difference: %.4f dF/F\n', ...
    median(drug_metrics.mean_stim(valid_mean) - base_metrics.mean_stim(valid_mean)));

fprintf('\n--- EVOKED RESPONSE (Stimulus - Baseline) ---\n');
fprintf('Baseline Evoked:     %.4f ± %.4f (mean ± std, n=%d)\n', ...
    mean(base_metrics.delta_stim(valid_delta)), std(base_metrics.delta_stim(valid_delta)), n_valid_delta);
fprintf('Drug Evoked:         %.4f ± %.4f (mean ± std, n=%d)\n', ...
    mean(drug_metrics.delta_stim(valid_delta)), std(drug_metrics.delta_stim(valid_delta)), n_valid_delta);
fprintf('Median difference: %.4f dF/F\n', ...
    median(drug_metrics.delta_stim(valid_delta) - base_metrics.delta_stim(valid_delta)));

% Responsiveness
resp_threshold = 0.5;
base_resp = sum(base_metrics.mean_stim > resp_threshold) / size(base_dff, 1) * 100;
drug_resp = sum(drug_metrics.mean_stim > resp_threshold) / size(drug_dff, 1) * 100;
fprintf('\n--- RESPONSIVENESS (threshold: %.1f dF/F) ---\n', resp_threshold);
fprintf('Baseline: %.1f%% responsive (%d / %d ROIs)\n', base_resp, ...
    sum(base_metrics.mean_stim > resp_threshold), size(base_dff, 1));
fprintf('Drug:     %.1f%% responsive (%d / %d ROIs)\n', drug_resp, ...
    sum(drug_metrics.mean_stim > resp_threshold), size(drug_dff, 1));

% Statistical tests
fprintf('\n--- STATISTICAL TESTS (Wilcoxon Signed-Rank) ---\n');
if n_valid_mean < 3
    fprintf('WARNING: Fewer than 3 valid ROI pairs. Skipping statistical tests.\n');
else
    [p_mean, ~, stats_mean] = signrank(base_metrics.mean_stim(valid_mean), drug_metrics.mean_stim(valid_mean));
    fprintf('Mean stimulus response:\n');
    fprintf('  p-value: %.4g (significant at α=0.05: %s)\n', p_mean, iif(p_mean < 0.05, 'YES', 'NO'));
    fprintf('  Signed-rank statistic: %d\n', stats_mean.signedrank);
    fprintf('  Median effect: %.4f dF/F\n', ...
        median(drug_metrics.mean_stim(valid_mean)) - median(base_metrics.mean_stim(valid_mean)));
end

if n_valid_delta < 3
    fprintf('WARNING: Fewer than 3 valid ROI pairs for delta. Skipping delta test.\n');
else
    [p_delta, ~, stats_delta] = signrank(base_metrics.delta_stim(valid_delta), drug_metrics.delta_stim(valid_delta));
    fprintf('\nEvoked response (delta):\n');
    fprintf('  p-value: %.4g (significant at α=0.05: %s)\n', p_delta, iif(p_delta < 0.05, 'YES', 'NO'));
    fprintf('  Signed-rank statistic: %d\n', stats_delta.signedrank);
    fprintf('  Median effect: %.4f dF/F\n', ...
        median(drug_metrics.delta_stim(valid_delta)) - median(base_metrics.delta_stim(valid_delta)));
end

% Effect sizes (percent change)
fprintf('\n--- EFFECT SIZES ---\n');
pct_change_mean = 100 * (mean(drug_metrics.mean_stim(valid_mean)) - mean(base_metrics.mean_stim(valid_mean))) / ...
                       (abs(mean(base_metrics.mean_stim(valid_mean))) + eps);
pct_change_delta = 100 * (mean(drug_metrics.delta_stim(valid_delta)) - mean(base_metrics.delta_stim(valid_delta))) / ...
                        (abs(mean(base_metrics.delta_stim(valid_delta))) + eps);
fprintf('Mean response change: %+.1f%%\n', pct_change_mean);
fprintf('Evoked response change: %+.1f%%\n', pct_change_delta);

% ROI-level modulation
mod_idx = (drug_metrics.mean_stim(valid_mean) - base_metrics.mean_stim(valid_mean)) ./ ...
          (abs(drug_metrics.mean_stim(valid_mean)) + abs(base_metrics.mean_stim(valid_mean)) + eps);
n_increased = sum(mod_idx > 0.1); % >10% increase
n_decreased = sum(mod_idx < -0.1); % >10% decrease
n_unchanged = sum(abs(mod_idx) <= 0.1); % <10% change

fprintf('\n--- ROI-LEVEL MODULATION ---\n');
fprintf('Increased response (>10%%):  %d ROIs (%.1f%%)\n', n_increased, 100*n_increased/numel(mod_idx));
fprintf('Decreased response (<-10%%): %d ROIs (%.1f%%)\n', n_decreased, 100*n_decreased/numel(mod_idx));
fprintf('Unchanged (±10%%):           %d ROIs (%.1f%%)\n', n_unchanged, 100*n_unchanged/numel(mod_idx));

fprintf('====================================================================\n\n');

% Optional outputs in workspace
matched_results = struct();
matched_results.selected_plane_s2p = selected_plane_s2p;
matched_results.ca_type = ca_type;
matched_results.num_matched = size(base_dff, 1);
matched_results.base_row_idx = base_row_idx;
matched_results.drug_row_idx = drug_row_idx;
matched_results.base_metrics = base_metrics;
matched_results.drug_metrics = drug_metrics;
matched_results.base_dff = base_dff;
matched_results.drug_dff = drug_dff;
matched_results.base_time = base_time;
matched_results.drug_time = drug_time;
matched_results.stimulus_types = stimulus_types;


%% ------------------------- LOCAL FUNCTIONS -------------------------

function result = iif(condition, true_val, false_val)
% Inline if function for concise ternary operations
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function plot_violin(data, labels)
% Simple violin plot visualization using boxplot with overlaid distributions
% data: matrix where each column is a group (groups as columns)
% labels: cell array of group labels
    
    if nargin < 2
        labels = cellfun(@(x) sprintf('Group%d', x), num2cell(1:size(data,2)), 'UniformOutput', false);
    end
    
    % Create boxplot
    boxplot(data, 'Labels', labels, 'Widths', 0.5);
    
    % Overlay individual points with jitter
    hold on;
    for g = 1:size(data, 2)
        group_data = data(~isnan(data(:,g)), g);
        n = numel(group_data);
        % Add small random jitter to x position
        x_jitter = g + (rand(n,1)-0.5)*0.15;
        plot(x_jitter, group_data, 'o', 'MarkerSize', 4, 'MarkerFaceAlpha', 0.5, ...
            'MarkerEdgeColor', 'k', 'LineStyle', 'none');
    end
    hold off;
    
end

function t = get_time_vector(TimeCa)
% Handle either 1xN vector or 2xN format robustly.
    if isvector(TimeCa)
        t = TimeCa(:)';
    elseif size(TimeCa, 1) >= 2
        t = TimeCa(2, :);
    else
        t = TimeCa(1, :);
    end
end

function metrics = compute_stimulus_metrics(dff, t, Stimuli, stimulus_types)
% Compute simple per-ROI response metrics from stimulus windows.
% mean_stim: mean dF/F during all stimulus windows
% delta_stim: mean(dF/F during stim - dF/F during immediately preceding window)

    if isempty(Stimuli)
        warning('Stimuli is empty. Returning NaN metrics.');
        metrics = struct('mean_stim', nan(size(dff,1),1), 'delta_stim', nan(size(dff,1),1));
        return;
    end

    stim_idx_list = collect_stimulus_indices(t, Stimuli, stimulus_types);
    pre_idx_list = collect_prestim_indices(t, Stimuli, stimulus_types);

    % If no windows found for requested stimulus_types, try falling back to
    % using all stimuli and report available stimulus types to the user.
    if isempty(stim_idx_list)
        if ~isempty(stimulus_types)
            types = {};
            for k = 1:numel(Stimuli)
                if isfield(Stimuli(k), 'type') && ~isempty(Stimuli(k).type)
                    types{end+1} = Stimuli(k).type; %#ok<AGROW>
                end
            end
            types = unique(types);
            if isempty(types)
                warning('No valid stimulus windows found. Stimuli contain no ''type'' fields. Returning NaN metrics.');
                metrics = struct('mean_stim', nan(size(dff,1),1), 'delta_stim', nan(size(dff,1),1));
                return;
            else
                fprintf('WARNING: No stimulus windows found for type(s): %s\n', strjoin(stimulus_types, ', '));
                fprintf('Available stimulus types: %s\n', strjoin(types, ', '));
                fprintf('Retrying with all stimuli as fallback.\n');
                stim_idx_list = collect_stimulus_indices(t, Stimuli, {});
                pre_idx_list = collect_prestim_indices(t, Stimuli, {});
            end
        else
            warning('No valid stimulus windows found. Returning NaN metrics.');
            metrics = struct('mean_stim', nan(size(dff,1),1), 'delta_stim', nan(size(dff,1),1));
            return;
        end
    end

    if isempty(stim_idx_list)
        warning('No valid stimulus windows found after fallback. Returning NaN metrics.');
        metrics = struct('mean_stim', nan(size(dff,1),1), 'delta_stim', nan(size(dff,1),1));
        return;
    end

    % Concatenate all windows for a compact essential metric.
    stim_idx = unique([stim_idx_list{:}]);
    mean_stim = mean(dff(:, stim_idx), 2, 'omitnan');
    
    fprintf('Used %d frames from stimulus windows (out of %d total frames).\\n', numel(stim_idx), size(dff, 2));

    if isempty(pre_idx_list)
        delta_stim = nan(size(mean_stim));
        fprintf('No pre-stimulus baseline available. Delta response set to NaN.\\n');
    else
        pre_idx = unique([pre_idx_list{:}]);
        mean_pre = mean(dff(:, pre_idx), 2, 'omitnan');
        delta_stim = mean_stim - mean_pre;
        fprintf('Used %d frames from pre-stimulus baseline.\\n', numel(pre_idx));
    end

    metrics = struct('mean_stim', mean_stim, 'delta_stim', delta_stim);
end

function stim_idx_list = collect_stimulus_indices(t, Stimuli, stimulus_types)
    % Collect frame indices corresponding to stimulus windows.
    % Assumes TimeStimulusFrame contains frame indices (integers).
    stim_idx_list = {};
    for i = 1:numel(Stimuli)
        if ~isfield(Stimuli(i), 'TimeStimulusFrame') || isempty(Stimuli(i).TimeStimulusFrame)
            continue;
        end
        if ~should_keep_stimulus(Stimuli(i), stimulus_types)
            continue;
        end

        % TimeStimulusFrame is expected to be a vector of frame indices
        stim_frames = Stimuli(i).TimeStimulusFrame;
        stim_frames = unique(round(stim_frames)); % Ensure integer frame indices
        stim_frames = stim_frames(stim_frames >= 1 & stim_frames <= numel(t));
        
        if ~isempty(stim_frames)
            stim_idx_list{end+1} = stim_frames; %#ok<AGROW>
        end
    end
end

function pre_idx_list = collect_prestim_indices(t, Stimuli, stimulus_types)
    % Collect frame indices for pre-stimulus baseline (duration-matched window before stimulus).
    % Assumes TimeStimulusFrame contains frame indices.
    pre_idx_list = {};
    for i = 1:numel(Stimuli)
        if ~isfield(Stimuli(i), 'TimeStimulusFrame') || isempty(Stimuli(i).TimeStimulusFrame)
            continue;
        end
        if ~should_keep_stimulus(Stimuli(i), stimulus_types)
            continue;
        end

        % Extract stimulus window frames
        stim_frames = Stimuli(i).TimeStimulusFrame;
        stim_frames = unique(round(stim_frames));
        stim_frames = stim_frames(stim_frames >= 1 & stim_frames <= numel(t));
        
        if isempty(stim_frames)
            continue;
        end

        st_frame = min(stim_frames);
        en_frame = max(stim_frames);
        dur = en_frame - st_frame + 1;
        
        % Pre-stimulus window: same duration, ending just before stimulus
        pre_st_frame = max(1, st_frame - dur);
        pre_en_frame = st_frame - 1;
        
        if pre_en_frame >= pre_st_frame
            pre_frames = pre_st_frame:pre_en_frame;
            pre_idx_list{end+1} = pre_frames; %#ok<AGROW>
        end
    end
end

function keep = should_keep_stimulus(stim, stimulus_types)
    if isempty(stimulus_types)
        keep = true;
        return;
    end

    if ~isfield(stim, 'type') || isempty(stim.type)
        keep = false;
        return;
    end

    keep = any(strcmp(stim.type, stimulus_types));
end
