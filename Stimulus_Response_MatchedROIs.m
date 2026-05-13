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
stimulus_types = {};

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

%% PLOT 1: AVERAGE TRACE ACROSS MATCHED ROIS
figure('Name', 'Matched ROI Mean Traces', 'NumberTitle', 'off');
subplot(2,1,1);
plot(base_time, mean(base_dff, 1), 'r', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Mean dF/F');
title(sprintf('Baseline mean trace (%d matched ROIs)', size(base_dff, 1)));
grid on;

subplot(2,1,2);
plot(drug_time, mean(drug_dff, 1), 'b', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Mean dF/F');
title(sprintf('Drug mean trace (%d matched ROIs)', size(drug_dff, 1)));
grid on;

%% PLOT 2: HEATMAPS OF MATCHED ROIS
figure('Name', 'Matched ROI Heatmaps', 'NumberTitle', 'off');
subplot(1,2,1);
    imagesc(base_dff);
    axis tight;
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
    colorbar;
xlabel('Frame'); ylabel('Matched ROI');
title('Baseline matched ROIs');

subplot(1,2,2);
    imagesc(drug_dff);
    axis tight;
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
    colorbar;
xlabel('Frame'); ylabel('Matched ROI');
title('Drug matched ROIs');

%% RESPONSE METRICS (ESSENTIAL)
% Uses mean dF/F during stimulus minus pre-stimulus baseline.
base_metrics = compute_stimulus_metrics(base_dff, base_time, B.Stimuli, stimulus_types);
drug_metrics = compute_stimulus_metrics(drug_dff, drug_time, D.Stimuli, stimulus_types);

%% PLOT 3: PAIRED RESPONSE COMPARISON
figure('Name', 'Matched ROI Response Comparison', 'NumberTitle', 'off');
subplot(1,2,1);
plot([1 2], [base_metrics.mean_stim, drug_metrics.mean_stim]', '-o', 'Color', [0.7 0.7 0.7]);
hold on;
plot([1 2], [mean(base_metrics.mean_stim) mean(drug_metrics.mean_stim)], 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Baseline', 'Drug'});
ylabel('Mean dF/F during stimulus');
title('Per-ROI stimulus mean (paired)');
grid on;

subplot(1,2,2);
plot([1 2], [base_metrics.delta_stim, drug_metrics.delta_stim]', '-o', 'Color', [0.7 0.7 0.7]);
hold on;
plot([1 2], [mean(base_metrics.delta_stim) mean(drug_metrics.delta_stim)], 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
set(gca, 'XTick', [1 2], 'XTickLabel', {'Baseline', 'Drug'});
ylabel('Stimulus - pre-stimulus dF/F');
title('Per-ROI evoked delta (paired)');
grid on;

%% BASIC STATS OUTPUT
[p_mean, ~, stats_mean] = signrank(base_metrics.mean_stim, drug_metrics.mean_stim);
[p_delta, ~, stats_delta] = signrank(base_metrics.delta_stim, drug_metrics.delta_stim);

fprintf('\n=== Response Summary (Matched ROIs) ===\n');
fprintf('Baseline mean(stim): %.4f\n', mean(base_metrics.mean_stim));
fprintf('Drug mean(stim):     %.4f\n', mean(drug_metrics.mean_stim));
fprintf('signrank p (mean_stim): %.3g, signed-rank=%g\n', p_mean, stats_mean.signedrank);

fprintf('Baseline mean(delta): %.4f\n', mean(base_metrics.delta_stim));
fprintf('Drug mean(delta):     %.4f\n', mean(drug_metrics.delta_stim));
fprintf('signrank p (delta): %.3g, signed-rank=%g\n', p_delta, stats_delta.signedrank);

% Optional outputs in workspace
matched_results = struct();
matched_results.selected_plane_s2p = selected_plane_s2p;
matched_results.ca_type = ca_type;
matched_results.num_matched = size(base_dff, 1);
matched_results.base_row_idx = base_row_idx;
matched_results.drug_row_idx = drug_row_idx;
matched_results.base_metrics = base_metrics;
matched_results.drug_metrics = drug_metrics;


%% ------------------------- LOCAL FUNCTIONS -------------------------
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

    if isempty(stim_idx_list)
        warning('No valid stimulus windows found. Returning NaN metrics.');
        metrics = struct('mean_stim', nan(size(dff,1),1), 'delta_stim', nan(size(dff,1),1));
        return;
    end

    % Concatenate all windows for a compact essential metric.
    stim_idx = unique([stim_idx_list{:}]);
    mean_stim = mean(dff(:, stim_idx), 2, 'omitnan');

    if isempty(pre_idx_list)
        delta_stim = nan(size(mean_stim));
    else
        pre_idx = unique([pre_idx_list{:}]);
        mean_pre = mean(dff(:, pre_idx), 2, 'omitnan');
        delta_stim = mean_stim - mean_pre;
    end

    metrics = struct('mean_stim', mean_stim, 'delta_stim', delta_stim);
end

function stim_idx_list = collect_stimulus_indices(t, Stimuli, stimulus_types)
    stim_idx_list = {};
    for i = 1:numel(Stimuli)
        if ~isfield(Stimuli(i), 'TimeStimulusFrame') || isempty(Stimuli(i).TimeStimulusFrame)
            continue;
        end
        if ~should_keep_stimulus(Stimuli(i), stimulus_types)
            continue;
        end

        st = Stimuli(i).TimeStimulusFrame(1);
        en = Stimuli(i).TimeStimulusFrame(end);
        idx = find(t >= st & t <= en);
        if ~isempty(idx)
            stim_idx_list{end+1} = idx; %#ok<AGROW>
        end
    end
end

function pre_idx_list = collect_prestim_indices(t, Stimuli, stimulus_types)
    pre_idx_list = {};
    for i = 1:numel(Stimuli)
        if ~isfield(Stimuli(i), 'TimeStimulusFrame') || isempty(Stimuli(i).TimeStimulusFrame)
            continue;
        end
        if ~should_keep_stimulus(Stimuli(i), stimulus_types)
            continue;
        end

        st = Stimuli(i).TimeStimulusFrame(1);
        en = Stimuli(i).TimeStimulusFrame(end);
        dur = en - st;
        if dur <= 0
            continue;
        end

        pre_st = st - dur;
        pre_en = st;
        idx = find(t >= pre_st & t < pre_en);
        if ~isempty(idx)
            pre_idx_list{end+1} = idx; %#ok<AGROW>
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
