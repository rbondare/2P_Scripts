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

%% STIMULUS ALIGNMENT SETTINGS
% Stimulus type to focus on for detailed response analysis
selected_stimulus_type = 'spontaneous';  % e.g., 'spontaneous', 'grating', 'moving_bar'

% Frame selection within stimulus window 
stim_frame_start = 1;      % start frame within stimulus window (e.g., 1 to start from beginning)
stim_frame_count = 1000;    % number of frames to show 

% Plotting options
n_top_snr = 8;             % number of top SNR ROIs to display 
heatmap_clim = [0 5];      % colorbar axis limits for heatmaps

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

%% STIMULUS ALIGNMENT AND SNR CALCULATION
fprintf('\n--- ALIGNING STIMULUS RESPONSES ---\n');

% Extract frame ranges for stimulus presentations of selected type
base_stim_frame_ranges = {};
drug_stim_frame_ranges = {};

fprintf('Searching for "%s" stimulus presentations...\n', selected_stimulus_type);

% Get TimeStimulusFrame for each stimulus of selected type
for i = 1:numel(B.Stimuli)
    if isfield(B.Stimuli(i), 'type') && strcmp(B.Stimuli(i).type, selected_stimulus_type)
        if isfield(B.Stimuli(i), 'TimeStimulusFrame') && ~isempty(B.Stimuli(i).TimeStimulusFrame)
            stim_frames = B.Stimuli(i).TimeStimulusFrame;
            st_frame = min(stim_frames(:));
            en_frame = max(stim_frames(:));
            base_stim_frame_ranges{end+1} = [st_frame, en_frame];
        end
    end
end

for i = 1:numel(D.Stimuli)
    if isfield(D.Stimuli(i), 'type') && strcmp(D.Stimuli(i).type, selected_stimulus_type)
        if isfield(D.Stimuli(i), 'TimeStimulusFrame') && ~isempty(D.Stimuli(i).TimeStimulusFrame)
            stim_frames = D.Stimuli(i).TimeStimulusFrame;
            st_frame = min(stim_frames(:));
            en_frame = max(stim_frames(:));
            drug_stim_frame_ranges{end+1} = [st_frame, en_frame];
        end
    end
end

% Ensure we have matching lengths
n_matched = min(numel(base_stim_frame_ranges), numel(drug_stim_frame_ranges));
if n_matched == 0
    error('No matching stimulus presentations found. Check stimulus type: "%s"', selected_stimulus_type);
end

fprintf('Extracted %d frame ranges for each condition.\n', n_matched);

% Calculate SNR for each ROI across matched stimulus presentations
base_snr = zeros(size(base_dff, 1), n_matched);
drug_snr = zeros(size(drug_dff, 1), n_matched);

fprintf('Calculating SNR across stimulus presentations...\n');
for i = 1:n_matched
    % Baseline
    if ~isempty(base_stim_frame_ranges{i})
        st_idx = base_stim_frame_ranges{i}(1);
        en_idx = base_stim_frame_ranges{i}(2);
        if en_idx <= size(base_dff, 2) && st_idx >= 1
            base_window = base_dff(:, st_idx:en_idx);
            base_snr(:, i) = calculate_snr_window(base_window);
        end
    end
    
    % Drug
    if ~isempty(drug_stim_frame_ranges{i})
        st_idx = drug_stim_frame_ranges{i}(1);
        en_idx = drug_stim_frame_ranges{i}(2);
        if en_idx <= size(drug_dff, 2) && st_idx >= 1
            drug_window = drug_dff(:, st_idx:en_idx);
            drug_snr(:, i) = calculate_snr_window(drug_window);
        end
    end
end

% Average SNR across stimulus presentations for each ROI
mean_snr_base = mean(base_snr, 2, 'omitnan');
mean_snr_drug = mean(drug_snr, 2, 'omitnan');

% Get top SNR ROI indices
[~, base_top_idx] = sort(mean_snr_base, 'descend');
[~, drug_top_idx] = sort(mean_snr_drug, 'descend');
top_idx_display = unique([base_top_idx(1:min(n_top_snr, numel(base_top_idx))); ...
                          drug_top_idx(1:min(n_top_snr, numel(drug_top_idx)))]);
top_idx_display = top_idx_display(1:min(n_top_snr, numel(top_idx_display)));

fprintf('Top %d SNR ROIs identified for visualization.\n', numel(top_idx_display));

%% PLOT 1: TOP SNR ROI RESPONSES 
% Display individual ROI responses during stimulus window
figure('Name', sprintf('High SNR dFF Responses - %s', selected_stimulus_type), ...
       'NumberTitle', 'off', 'Position', [100 100 1200 1000]);

n_plots = numel(top_idx_display);
n_cols = 2;  % baseline left, drug right
n_rows = ceil(n_plots / 1);  % one row per ROI

for plot_num = 1:n_plots
    roi_idx = top_idx_display(plot_num);
    
    % ===== BASELINE RESPONSES (Left) =====
    subplot(n_rows, 2, 2*plot_num - 1);
    
    hold on;
    baseline_traces = [];
    
    for stim_idx = 1:n_matched
        if ~isempty(base_stim_frame_ranges{stim_idx})
            st_frame = base_stim_frame_ranges{stim_idx}(1);
            en_frame = base_stim_frame_ranges{stim_idx}(2);
            
            % Extract selected frames from this stimulus window
            disp_st = st_frame + stim_frame_start - 1;
            disp_en = min(st_frame + stim_frame_start + stim_frame_count - 1, en_frame);
            
            if disp_en > disp_st && disp_st >= 1 && disp_en <= size(base_dff, 2)
                trace = base_dff(roi_idx, disp_st:disp_en);
                baseline_traces(:, stim_idx) = trace';
            end
        end
    end
    
    % Plot all traces with lighter colors
    if ~isempty(baseline_traces)
        x_frames = 1:size(baseline_traces, 1);
        for s = 1:size(baseline_traces, 2)
            trace_s = baseline_traces(:, s);
            % Only plot if trace is valid (not all NaN or empty)
            if ~all(isnan(trace_s))
                h = plot(x_frames, trace_s, 'r-', 'LineWidth', 0.8);
                h.Color = [1 0.4 0.4];  % Light red
            end
        end
        % Overlay mean
        mean_trace = mean(baseline_traces, 2, 'omitnan');
        plot(x_frames, mean_trace, 'r-', 'LineWidth', 2.5);
    end
    
    hold off;
    set(gca, 'YLim', heatmap_clim, 'LineWidth', 1.5, 'FontSize', 10);
    ylabel('dF/F', 'FontSize', 10);
    title(sprintf('ROI %d - Baseline (SNR=%.2f)', roi_idx, mean_snr_base(roi_idx)), ...
          'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    box on;
    
    % ===== DRUG RESPONSES (Right) =====
    subplot(n_rows, 2, 2*plot_num);
    
    hold on;
    drug_traces = [];
    
    for stim_idx = 1:n_matched
        if ~isempty(drug_stim_frame_ranges{stim_idx})
            st_frame = drug_stim_frame_ranges{stim_idx}(1);
            en_frame = drug_stim_frame_ranges{stim_idx}(2);
            
            % Extract selected frames from this stimulus window
            disp_st = st_frame + stim_frame_start - 1;
            disp_en = min(st_frame + stim_frame_start + stim_frame_count - 1, en_frame);
            
            if disp_en > disp_st && disp_st >= 1 && disp_en <= size(drug_dff, 2)
                trace = drug_dff(roi_idx, disp_st:disp_en);
                drug_traces(:, stim_idx) = trace';
            end
        end
    end
    
    % Plot all traces with lighter colors
    if ~isempty(drug_traces)
        x_frames = 1:size(drug_traces, 1);
        for s = 1:size(drug_traces, 2)
            trace_s = drug_traces(:, s);
            % Only plot if trace is valid (not all NaN or empty)
            if ~all(isnan(trace_s))
                h = plot(x_frames, trace_s, 'b-', 'LineWidth', 0.8);
                h.Color = [0.4 0.4 1];  % Light blue
            end
        end
        % Overlay mean
        mean_trace = mean(drug_traces, 2, 'omitnan');
        plot(x_frames, mean_trace, 'b-', 'LineWidth', 2.5);
    end
    
    hold off;
    set(gca, 'YLim', heatmap_clim, 'LineWidth', 1.5, 'FontSize', 10);
    ylabel('dF/F', 'FontSize', 10);
    title(sprintf('ROI %d - Drug (SNR=%.2f)', roi_idx, mean_snr_drug(roi_idx)), ...
          'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    box on;
end

sgtitle(sprintf('High SNR dFF Responses - %s (frames %d-%d)', ...
    selected_stimulus_type, stim_frame_start, stim_frame_start + stim_frame_count - 1), ...
    'FontSize', 13, 'FontWeight', 'bold');

%% PLOT 2: HEATMAP SORTED BY SNR
figure('Name', sprintf('SNR-Sorted Heatmap - %s', selected_stimulus_type), ...
       'NumberTitle', 'off', 'Position', [100 100 1300 600]);

% Build heatmap matrices from selected frames within each stimulus window
base_heatmap = [];
drug_heatmap = [];

for stim_idx = 1:n_matched
    % Baseline heatmap
    if ~isempty(base_stim_frame_ranges{stim_idx})
        st_frame = base_stim_frame_ranges{stim_idx}(1);
        en_frame = base_stim_frame_ranges{stim_idx}(2);
        
        disp_st = st_frame + stim_frame_start - 1;
        disp_en = min(st_frame + stim_frame_start + stim_frame_count - 1, en_frame);
        
        if disp_en > disp_st && disp_st >= 1 && disp_en <= size(base_dff, 2)
            base_heatmap = [base_heatmap, base_dff(:, disp_st:disp_en)];
        end
    end
    
    % Drug heatmap
    if ~isempty(drug_stim_frame_ranges{stim_idx})
        st_frame = drug_stim_frame_ranges{stim_idx}(1);
        en_frame = drug_stim_frame_ranges{stim_idx}(2);
        
        disp_st = st_frame + stim_frame_start - 1;
        disp_en = min(st_frame + stim_frame_start + stim_frame_count - 1, en_frame);
        
        if disp_en > disp_st && disp_st >= 1 && disp_en <= size(drug_dff, 2)
            drug_heatmap = [drug_heatmap, drug_dff(:, disp_st:disp_en)];
        end
    end
end

% Sort by SNR (descending)
[~, sort_idx_base] = sort(mean_snr_base, 'descend');
[~, sort_idx_drug] = sort(mean_snr_drug, 'descend');

% Baseline heatmap
subplot(1, 2, 1);
if ~isempty(base_heatmap)
    imagesc(base_heatmap(sort_idx_base, :));
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
end
c = colorbar;
caxis(heatmap_clim);
ylabel(c, 'dF/F');
xlabel('Frame');
ylabel('ROI (sorted by SNR)');
title(sprintf('Baseline - %s', selected_stimulus_type), 'FontSize', 11, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1.5, 'FontSize', 10);

% Drug heatmap
subplot(1, 2, 2);
if ~isempty(drug_heatmap)
    imagesc(drug_heatmap(sort_idx_drug, :));
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
end
c = colorbar;
caxis(heatmap_clim);
ylabel(c, 'dF/F');
xlabel('Frame');
ylabel('ROI (sorted by SNR)');
title(sprintf('Drug - %s', selected_stimulus_type), 'FontSize', 11, 'FontWeight', 'bold');
set(gca, 'LineWidth', 1.5, 'FontSize', 10);


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

function snr_roi = calculate_snr_window(response_window)
% Calculate signal-to-noise ratio for a response window.
% Assumes the entire window is stimulus response (no pre-stimulus period).
%
% response_window: matrix of dFF responses (nROIs x nFrames)
% snr_roi: SNR for each ROI (nROIs x 1)
    
    % Signal: mean response across entire window
    signal = mean(response_window, 2);
    
    % Noise: std dev of response across entire window
    noise = std(response_window, [], 2);
    
    % SNR = signal / noise
    snr_roi = signal ./ (noise + eps);  % Add eps to avoid division by zero
end

