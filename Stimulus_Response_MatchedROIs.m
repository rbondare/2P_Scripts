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

% Time window around stimulus (in seconds)
time_before_stim = 0.5;  % seconds before stimulus start
time_after_stim = 1.0;   % seconds after stimulus end

% Plotting options
n_top_snr = 12;          % number of top SNR ROIs to display (e.g., 12 for 2 rows x 6 cols)
heatmap_clim = [0 5];    % colorbar axis limits for heatmaps

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

% Find stimulus presentations of selected type in both recordings
base_stim_frames = find_stimulus_presentations(B.Stimuli, selected_stimulus_type);
drug_stim_frames = find_stimulus_presentations(D.Stimuli, selected_stimulus_type);

fprintf('Found %d presentations of "%s" in baseline.\n', numel(base_stim_frames), selected_stimulus_type);
fprintf('Found %d presentations of "%s" in drug.\n', numel(drug_stim_frames), selected_stimulus_type);

% Convert time window to frame counts
base_frame_before = round(time_before_stim * (numel(base_time) / base_time(end)));
base_frame_after = round(time_after_stim * (numel(base_time) / base_time(end)));
drug_frame_before = round(time_before_stim * (numel(drug_time) / drug_time(end)));
drug_frame_after = round(time_after_stim * (numel(drug_time) / drug_time(end)));

fprintf('Baseline: %.2f sec before (~%d frames), %.2f sec after (~%d frames)\n', ...
    time_before_stim, base_frame_before, time_after_stim, base_frame_after);
fprintf('Drug:     %.2f sec before (~%d frames), %.2f sec after (~%d frames)\n', ...
    time_before_stim, drug_frame_before, time_after_stim, drug_frame_after);

% Extract response windows and calculate SNR
base_responses = [];
drug_responses = [];
base_snr = [];
drug_snr = [];

n_matched = min(numel(base_stim_frames), numel(drug_stim_frames));
for i = 1:n_matched
    % Extract baseline response window
    base_st = max(1, base_stim_frames(i) - base_frame_before);
    base_en = min(size(base_dff, 2), base_stim_frames(i) + base_frame_after);
    if base_en > base_st
        base_window = base_dff(:, base_st:base_en);
        base_responses(:, i) = mean(base_window, 2);
        base_snr(:, i) = calculate_snr(base_window, base_stim_frames(i) - base_st + 1);
    end
    
    % Extract drug response window
    drug_st = max(1, drug_stim_frames(i) - drug_frame_before);
    drug_en = min(size(drug_dff, 2), drug_stim_frames(i) + drug_frame_after);
    if drug_en > drug_st
        drug_window = drug_dff(:, drug_st:drug_en);
        drug_responses(:, i) = mean(drug_window, 2);
        drug_snr(:, i) = calculate_snr(drug_window, drug_stim_frames(i) - drug_st + 1);
    end
end

% Average SNR across stimulus presentations for each ROI
mean_snr_base = mean(base_snr, 2, 'omitnan');
mean_snr_drug = mean(drug_snr, 2, 'omitnan');

fprintf('Processed %d matched stimulus presentations.\n', n_matched);

% Get top SNR ROI indices
[~, base_top_idx] = sort(mean_snr_base, 'descend');
[~, drug_top_idx] = sort(mean_snr_drug, 'descend');
top_idx_display = unique([base_top_idx(1:min(n_top_snr, numel(base_top_idx))); ...
                          drug_top_idx(1:min(n_top_snr, numel(drug_top_idx)))]);
top_idx_display = top_idx_display(1:min(n_top_snr, numel(top_idx_display)));

fprintf('Top %d SNR ROIs identified for visualization.\n', numel(top_idx_display));

%% PLOT 1: HIGHEST SNR RESPONSES (Baseline vs Drug)
figure('Name', sprintf('Top SNR Responses - %s', selected_stimulus_type), ...
       'NumberTitle', 'off', 'Position', [100 100 1400 800]);

n_plots = numel(top_idx_display);
n_cols = 6;
n_rows = ceil(n_plots / n_cols);

for plot_idx = 1:n_plots
    roi_idx = top_idx_display(plot_idx);
    
    % Baseline response (left subplot)
    ax_base = subplot(n_rows, 2*n_cols, 2*plot_idx - 1);
    hold on;
    for stim_idx = 1:min(3, n_matched)  % Show up to 3 stimulus presentations
        base_st = max(1, base_stim_frames(stim_idx) - base_frame_before);
        base_en = min(size(base_dff, 2), base_stim_frames(stim_idx) + base_frame_after);
        if base_en > base_st
            window_data = base_dff(roi_idx, base_st:base_en);
            x_axis = 1:numel(window_data);
            plot(x_axis, window_data, 'r-', 'Alpha', 0.5);
        end
    end
    hold off;
    ylabel('dF/F');
    title(sprintf('ROI %d Baseline (SNR=%.2f)', roi_idx, mean_snr_base(roi_idx)));
    set(gca, 'YLim', heatmap_clim);
    grid on;
    
    % Drug response (right subplot)
    ax_drug = subplot(n_rows, 2*n_cols, 2*plot_idx);
    hold on;
    for stim_idx = 1:min(3, n_matched)  % Show up to 3 stimulus presentations
        drug_st = max(1, drug_stim_frames(stim_idx) - drug_frame_before);
        drug_en = min(size(drug_dff, 2), drug_stim_frames(stim_idx) + drug_frame_after);
        if drug_en > drug_st
            window_data = drug_dff(roi_idx, drug_st:drug_en);
            x_axis = 1:numel(window_data);
            plot(x_axis, window_data, 'b-', 'Alpha', 0.5);
        end
    end
    hold off;
    ylabel('dF/F');
    title(sprintf('ROI %d Drug (SNR=%.2f)', roi_idx, mean_snr_drug(roi_idx)));
    set(gca, 'YLim', heatmap_clim);
    grid on;
end

sgtitle(sprintf('Highest SNR Responses - Stimulus Type: %s', selected_stimulus_type));

%% PLOT 2: HEATMAP SORTED BY SNR
figure('Name', sprintf('SNR-Sorted Heatmap - %s', selected_stimulus_type), ...
       'NumberTitle', 'off', 'Position', [100 100 1300 600]);

% Build heatmap matrices from all stimulus presentations
base_heatmap = [];
drug_heatmap = [];

for stim_idx = 1:n_matched
    base_st = max(1, base_stim_frames(stim_idx) - base_frame_before);
    base_en = min(size(base_dff, 2), base_stim_frames(stim_idx) + base_frame_after);
    if base_en > base_st
        base_heatmap = [base_heatmap, base_dff(:, base_st:base_en)];
    end
    
    drug_st = max(1, drug_stim_frames(stim_idx) - drug_frame_before);
    drug_en = min(size(drug_dff, 2), drug_stim_frames(stim_idx) + drug_frame_after);
    if drug_en > drug_st
        drug_heatmap = [drug_heatmap, drug_dff(:, drug_st:drug_en)];
    end
end

% Sort by SNR (descending)
[~, sort_idx_base] = sort(mean_snr_base, 'descend');
[~, sort_idx_drug] = sort(mean_snr_drug, 'descend');

% Baseline heatmap
subplot(1, 2, 1);
imagesc(base_heatmap(sort_idx_base, :));
colormap(flipud(gray));
set(gca, 'YDir', 'reverse');
c = colorbar;
caxis(heatmap_clim);
ylabel(c, 'dF/F');
xlabel('Frame');
ylabel('ROI (sorted by SNR)');
title(sprintf('Baseline - %s\n(n=%d ROIs, %d stim presentations)', ...
    selected_stimulus_type, size(base_dff, 1), n_matched));

% Drug heatmap
subplot(1, 2, 2);
imagesc(drug_heatmap(sort_idx_drug, :));
colormap(flipud(gray));
set(gca, 'YDir', 'reverse');
c = colorbar;
caxis(heatmap_clim);
ylabel(c, 'dF/F');
xlabel('Frame');
ylabel('ROI (sorted by SNR)');
title(sprintf('Drug - %s\n(n=%d ROIs, %d stim presentations)', ...
    selected_stimulus_type, size(drug_dff, 1), n_matched));


%% ------------------------- LOCAL FUNCTIONS -------------------------

function result = iif(condition, true_val, false_val)
% Inline if function for concise ternary operations
    if condition
        result = true_val;
    else
        result = false_val;
    end
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

function stim_frames = find_stimulus_presentations(Stimuli, stim_type)
% Find all frame indices where a particular stimulus type is presented.
% Returns a vector of frame indices (one per stimulus presentation).
%
% stim_type: stimulus type to search for (string, e.g., 'spontaneous', 'grating')
% stim_frames: vector of frame indices where stimulus starts
    
    stim_frames = [];
    
    if isempty(Stimuli)
        return;
    end
    
    for i = 1:numel(Stimuli)
        if ~isfield(Stimuli(i), 'type') || isempty(Stimuli(i).type)
            continue;
        end
        
        if ~strcmp(Stimuli(i).type, stim_type)
            continue;
        end
        
        if ~isfield(Stimuli(i), 'TimeStimulusFrame') || isempty(Stimuli(i).TimeStimulusFrame)
            continue;
        end
        
        % Get the first frame of this stimulus presentation
        stim_frame_window = Stimuli(i).TimeStimulusFrame;
        first_frame = min(round(stim_frame_window));
        stim_frames = [stim_frames, first_frame]; %#ok<AGROW>
    end
    
    % Return as column vector
    stim_frames = stim_frames(:);
end

function snr_roi = calculate_snr(response_window, stim_start_idx)
% Calculate signal-to-noise ratio for a response window.
% Assumes stimulus starts at stim_start_idx.
%
% response_window: matrix of dFF responses (nROIs x nFrames)
% stim_start_idx: frame index where stimulus starts (1-based index into window)
% snr_roi: SNR for each ROI (nROIs x 1)
    
    [n_rois, n_frames] = size(response_window);
    snr_roi = zeros(n_rois, 1);
    
    if stim_start_idx < 1 || stim_start_idx > n_frames
        snr_roi = nan(n_rois, 1);
        return;
    end
    
    % Signal: mean response during stimulus (from stim_start to end of window)
    signal = mean(response_window(:, stim_start_idx:end), 2);
    
    % Noise: std dev of response before stimulus
    if stim_start_idx > 1
        noise = std(response_window(:, 1:stim_start_idx-1), [], 2);
    else
        noise = std(response_window, [], 2);
    end
    
    % SNR = signal / noise
    snr_roi = signal ./ (noise + eps);  % Add eps to avoid division by zero
end

