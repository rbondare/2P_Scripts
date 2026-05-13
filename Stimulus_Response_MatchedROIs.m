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
stimulus_types = {'grating'};

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

%% PLOT 1: AVERAGE TRACE ACROSS MATCHED ROIS
figure('Name', 'Matched ROI Mean Traces', 'NumberTitle', 'off');
subplot(2,1,1);
    trace_base = mean(base_dff, 1);
    plot(base_time, trace_base, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Mean dF/F');
    title(sprintf('Baseline mean trace (%d matched ROIs, %d frames)', size(base_dff, 1), size(base_dff, 2)));
    grid on;
    ylim([min(trace_base) - 0.5, max(trace_base) + 0.5]);

subplot(2,1,2);
    trace_drug = mean(drug_dff, 1);
    plot(drug_time, trace_drug, 'b', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Mean dF/F');
    title(sprintf('Drug mean trace (%d matched ROIs, %d frames)', size(drug_dff, 1), size(drug_dff, 2)));
    grid on;
    ylim([min(trace_drug) - 0.5, max(trace_drug) + 0.5]);

%% PLOT 2: HEATMAPS OF MATCHED ROIS
figure('Name', 'Matched ROI Heatmaps', 'NumberTitle', 'off');
subplot(1,2,1);
    imagesc(base_dff(:, plot_idx));
    axis tight;
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
    c = colorbar;
    caxis(heatmap_clim);
    ylabel(c, 'dF/F');
    xlabel('Frame (from frame ' + string(plot_frame_start) + ')');
    ylabel('Matched ROI');
    title(sprintf('Baseline matched ROIs (frames %d-%d)', plot_frame_start, plot_frame_end));

subplot(1,2,2);
    imagesc(drug_dff(:, plot_idx));
    axis tight;
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
    c = colorbar;
    caxis(heatmap_clim);
    ylabel(c, 'dF/F');
    xlabel('Frame (from frame ' + string(plot_frame_start) + ')');
    ylabel('Matched ROI');
    title(sprintf('Drug matched ROIs (frames %d-%d)', plot_frame_start, plot_frame_end));

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

%% PLOT 3: PAIRED RESPONSE COMPARISON
% Only plot valid (non-NaN) data
valid_idx = ~isnan(base_metrics.mean_stim) & ~isnan(drug_metrics.mean_stim);

figure('Name', 'Matched ROI Response Comparison', 'NumberTitle', 'off');
subplot(1,2,1);
    if sum(valid_idx) > 0
        plot([1 2], [base_metrics.mean_stim(valid_idx), drug_metrics.mean_stim(valid_idx)]', '-o', 'Color', [0.7 0.7 0.7]);
        hold on;
        plot([1 2], [mean(base_metrics.mean_stim(valid_idx)), mean(drug_metrics.mean_stim(valid_idx))], ...
            'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    else
        warning('No valid data for mean stimulus response plot.');
    end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Baseline', 'Drug'});
    ylabel('Mean dF/F during stimulus');
    title(sprintf('Per-ROI stimulus mean (n=%d)', sum(valid_idx)));
    grid on;

subplot(1,2,2);
    valid_delta_idx = ~isnan(base_metrics.delta_stim) & ~isnan(drug_metrics.delta_stim);
    if sum(valid_delta_idx) > 0
        plot([1 2], [base_metrics.delta_stim(valid_delta_idx), drug_metrics.delta_stim(valid_delta_idx)]', '-o', 'Color', [0.7 0.7 0.7]);
        hold on;
        plot([1 2], [mean(base_metrics.delta_stim(valid_delta_idx)), mean(drug_metrics.delta_stim(valid_delta_idx))], ...
            'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    else
        warning('No valid data for delta response plot.');
    end
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Baseline', 'Drug'});
    ylabel('Stimulus - pre-stimulus dF/F');
    title(sprintf('Per-ROI evoked delta (n=%d)', sum(valid_delta_idx)));
    grid on;

%% BASIC STATS OUTPUT
% Only compute statistics on valid (non-NaN) paired data
valid_pairs = ~isnan(base_metrics.mean_stim) & ~isnan(drug_metrics.mean_stim);
if sum(valid_pairs) < 3
    fprintf('WARNING: Fewer than 3 valid paired ROIs for statistics. Skipping statistical tests.\\n');
else
    [p_mean, ~, stats_mean] = signrank(base_metrics.mean_stim(valid_pairs), drug_metrics.mean_stim(valid_pairs));
    [p_delta, ~, stats_delta] = signrank(base_metrics.delta_stim(valid_pairs), drug_metrics.delta_stim(valid_pairs));
    
    fprintf('\\n=== Response Summary (Matched ROIs, n=%d valid pairs) ===\\n', sum(valid_pairs));
    fprintf('Baseline mean(stim): %.4f ± %.4f\\n', mean(base_metrics.mean_stim(valid_pairs), 'omitnan'), ...
            std(base_metrics.mean_stim(valid_pairs), 'omitnan'));
    fprintf('Drug mean(stim):     %.4f ± %.4f\\n', mean(drug_metrics.mean_stim(valid_pairs), 'omitnan'), ...
            std(drug_metrics.mean_stim(valid_pairs), 'omitnan'));
    fprintf('signrank p (mean_stim): %.3g, signed-rank=%g\\n', p_mean, stats_mean.signedrank);
    
    fprintf('\\nBaseline mean(delta): %.4f ± %.4f\\n', mean(base_metrics.delta_stim(valid_pairs), 'omitnan'), ...
            std(base_metrics.delta_stim(valid_pairs), 'omitnan'));
    fprintf('Drug mean(delta):     %.4f ± %.4f\\n', mean(drug_metrics.delta_stim(valid_pairs), 'omitnan'), ...
            std(drug_metrics.delta_stim(valid_pairs), 'omitnan'));
    fprintf('signrank p (delta): %.3g, signed-rank=%g\\n', p_delta, stats_delta.signedrank);
end

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
