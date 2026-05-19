%% Analysis for 2P data with ROI MATCHED OPTION and multi-recording comparison
%
% This script creates master response matrices for all neurons across stimuli,
% generates population and matched-ROI heatmaps, and provides stimulus-specific analysis 

%
% STRUCTURE:
%   1. Load and organize data
%   2. Extract master response matrices (all neurons, full stimulus duration)
%   3. Generate population heatmaps 
%   4. Generate matched ROI heatmaps 
%   5. Stimulus-specific analysis 

clear; clc; close all;

%% ====================== CONFIGURATION ======================

% Data files
baseline_file = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat";
drug_file     = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat";
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_V2.mat";

% Analysis parameters
ca_type = 1;                  % 1=FVDff, 2=deconvolved, 3=F
selected_plane = 1;           % Plane index (1-based)
stimulus_types_to_analyze = {'spontaneous', 'grating', 'moving_bar', 'sparse_local_global_flashes', 'checkers2'};

% Visualization
max_neurons_display = 500;    % For heatmap readability
caxis_lim = [0, 5];           % dF/F range for heatmaps
colormap_type = 'flipud(gray)';


%% ====================== LOAD DATA ======================

fprintf('\n========== DATA LOADING ==========\n');

fprintf('Loading baseline: %s\n', baseline_file);
B = load(baseline_file);

% Extract calcium data and organize by plane
n_planes = length(B.CaData);
dff_plane_base = cell(n_planes, 1);  % dFF for each plane
selected_rois_by_plane = cell(n_planes, 1);  % ROI indices for each plane

% Process each plane
for plane_idx = 1:n_planes
    if ~isempty(B.CaData(plane_idx).Ca_dFF)
        dff_plane_base{plane_idx} = B.CaData(plane_idx).Ca_dFF;
        centroid = B.CaData(plane_idx).Ca_centroid_voxel;
        centroidZ = centroid(:, 3);
        selected_rois_by_plane{plane_idx} = find(centroidZ == plane_idx);
    else
        dff_plane_base{plane_idx} = [];
        selected_rois_by_plane{plane_idx} = [];
    end
end

% For backward compatibility: use selected_plane
base_dff = B.CaData(selected_plane).Ca_dFF;
centroid = B.CaData(selected_plane).Ca_centroid_voxel;
centroidX = centroid(:,1);
centroidY = centroid(:,2);
centroidZ = centroid(:,3);

selected_roi = find(centroidZ == selected_plane);

% Create full dFF (all ROIs from selected plane)
dff_plane = dff_plane_base{selected_plane};
selected_rois = selected_rois_by_plane{selected_plane};

% Create combined dFF for all valid planes (if needed for cross-plane analysis)
dff_all = [];  % Will contain concatenated dFF from all planes
n_rois_all = 0;
for plane_idx = 1:n_planes
    if ~isempty(dff_plane_base{plane_idx})
        dff_all = [dff_all; dff_plane_base{plane_idx}];
        n_rois_all = n_rois_all + size(dff_plane_base{plane_idx}, 1);
    end
end

n_rois_base = size(base_dff, 1);
n_frames_base = size(base_dff, 2);





fprintf('Loading drug: %s\n', drug_file);
D = load(drug_file);

% Extract drug calcium data and organize by plane
n_planes_drug = length(D.CaData);
dff_plane_drug = cell(n_planes_drug, 1);  % dFF for each plane
drug_rois_by_plane = cell(n_planes_drug, 1);  % ROI indices for each plane

% Process each plane for drug data
for plane_idx = 1:n_planes_drug
    if ~isempty(D.CaData(plane_idx).Ca_dFF)
        dff_plane_drug{plane_idx} = D.CaData(plane_idx).Ca_dFF;
        centroid = D.CaData(plane_idx).Ca_centroid_voxel;
        centroidZ = centroid(:, 3);
        drug_rois_by_plane{plane_idx} = find(centroidZ == plane_idx);
    else
        dff_plane_drug{plane_idx} = [];
        drug_rois_by_plane{plane_idx} = [];
    end
end

% For selected plane
drug_dff = get_calcium_data(D.CaData(selected_plane), ca_type);
drug_dff_plane = dff_plane_drug{selected_plane};
drug_selected_rois = drug_rois_by_plane{selected_plane};

% Create combined dFF for all valid planes (drug, if needed for cross-plane analysis)
drug_dff_all = [];  % Will contain concatenated dFF from all drug planes
n_rois_drug_all = 0;
for plane_idx = 1:n_planes_drug
    if ~isempty(dff_plane_drug{plane_idx})
        drug_dff_all = [drug_dff_all; dff_plane_drug{plane_idx}];
        n_rois_drug_all = n_rois_drug_all + size(dff_plane_drug{plane_idx}, 1);
    end
end

n_rois_drug = size(drug_dff, 1);
n_frames_drug = size(drug_dff, 2);

fprintf('Baseline: %d ROIs × %d frames\n', n_rois_base, n_frames_base);
fprintf('Drug:     %d ROIs × %d frames\n', n_rois_drug, n_frames_drug);

% Load ROI matching if available
matched_rois_available = false;  % Default: no matched ROIs
base_match_idx = [];
drug_match_idx = [];
n_matched = 0;

if ~isempty(roi_match_file) && isfile(roi_match_file)
    fprintf('Loading ROI matches: %s\n', roi_match_file);
    M = load(roi_match_file);
    if isfield(M.roiMatchData, 'allSessionMapping')
        base_match_idx = M.roiMatchData.allSessionMapping(:, 1);
        drug_match_idx = M.roiMatchData.allSessionMapping(:, 2);
        n_matched = length(base_match_idx);
        matched_rois_available = true;
        fprintf('Found %d matched ROI pairs\n', n_matched);
    else
        fprintf('  Warning: allSessionMapping not found in ROI match file\n');
    end
else
    if ~isempty(roi_match_file)
        fprintf('  Warning: ROI match file not found: %s\n', roi_match_file);
    end
end

% Validate stimulus types
available_stim_types = unique({B.Stimuli(:).type});
fprintf('\nAvailable stimulus types: %s\n', strjoin(available_stim_types, ', '));

%% ====================== EXTRACT MASTER RESPONSE MATRICES ======================

fprintf('\n========== EXTRACTING MASTER RESPONSE MATRICES ==========\n');

% Create master structs: one entry per stimulus type
% Each contains: baseline_responses (neurons × time × presentations)
%                drug_responses (neurons × time × presentations)
%                frame_ranges, properties, etc.

master_data = struct();

for stim_idx = 1:length(stimulus_types_to_analyze)
    stim_type = stimulus_types_to_analyze{stim_idx};
    
    % Check if this stimulus exists
    stim_exists_base = any(strcmp({B.Stimuli(:).type}, stim_type));
    stim_exists_drug = any(strcmp({D.Stimuli(:).type}, stim_type));
    
    if ~stim_exists_base || ~stim_exists_drug
        fprintf('Skipping "%s" (not in both baseline and drug)\n', stim_type);
        continue;
    end
    
    fprintf('\nProcessing "%s"...\n', stim_type);
    
    % Extract full stimulus responses for baseline
    [base_responses, base_frame_ranges, base_time_ranges, base_properties] = extract_full_stimulus_responses(...
        B.Stimuli, base_dff, stim_type, B.TimeCa);
    
    % Extract full stimulus responses for drug
    [drug_responses, drug_frame_ranges, drug_time_ranges, drug_properties] = extract_full_stimulus_responses(...
        D.Stimuli, drug_dff, stim_type, D.TimeCa);
    
    % Store in master struct
    field_name = matlab.lang.makeValidName(stim_type);
    master_data.(field_name).baseline_responses = base_responses;
    master_data.(field_name).drug_responses = drug_responses;
    master_data.(field_name).baseline_frame_ranges = base_frame_ranges;
    master_data.(field_name).drug_frame_ranges = drug_frame_ranges;
    master_data.(field_name).baseline_time_ranges = base_time_ranges;  % NEW: time in seconds
    master_data.(field_name).drug_time_ranges = drug_time_ranges;      % NEW: time in seconds
    master_data.(field_name).baseline_properties = base_properties;
    master_data.(field_name).drug_properties = drug_properties;
    master_data.(field_name).stimulus_type = stim_type;
    
    n_pres = length(base_responses);
    fprintf('  Extracted: %d neurons × variable time × %d presentations\n', ...
        size(base_responses{1}, 1), n_pres);
    fprintf('  Baseline timing: %.2f - %.2f sec\n', min(base_time_ranges(:,1)), max(base_time_ranges(:,2)));
    fprintf('  Drug timing: %.2f - %.2f sec\n', min(drug_time_ranges(:,1)), max(drug_time_ranges(:,2)));
end

fprintf('\nMaster data matrices created successfully.\n');

%% ====================== GENERATE POPULATION HEATMAPS ======================

fprintf('\n========== GENERATING POPULATION HEATMAPS ==========\n');

stim_fields = fieldnames(master_data);

for field_idx = 1:length(stim_fields)
    field_name = stim_fields{field_idx};
    data = master_data.(field_name);
    stim_type = data.stimulus_type;
    
    fprintf('Generating heatmap for "%s"\n', stim_type);
    
    % Create figure: left=baseline all neurons, right=drug all neurons
    fig = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
        'Name', sprintf('Population Heatmap: %s', stim_type));
    
    % BASELINE - All neurons ranked by max activity
    subplot(1, 2, 1);
    baseline_heatmap = average_stimulus_presentations(data.baseline_responses);
    [~, sort_idx_base] = sort(max(baseline_heatmap, [], 2), 'descend');
    
    % Limit display to max_neurons_display for clarity
    display_idx_base = sort_idx_base(1:min(max_neurons_display, length(sort_idx_base)));
    imagesc(baseline_heatmap(display_idx_base, :));
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(caxis_lim);
    colorbar;
    xlabel('Time (frames)');
    ylabel('ROI (ranked by max activity)');
    title(sprintf('Baseline - %s\n(%d/%d neurons shown)', stim_type, ...
        length(display_idx_base), size(baseline_heatmap, 1)), 'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % DRUG - All neurons ranked by max activity
    subplot(1, 2, 2);
    drug_heatmap = average_stimulus_presentations(data.drug_responses);
    [~, sort_idx_drug] = sort(max(drug_heatmap, [], 2), 'descend');
    
    display_idx_drug = sort_idx_drug(1:min(max_neurons_display, length(sort_idx_drug)));
    imagesc(drug_heatmap(display_idx_drug, :));
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(caxis_lim);
    c = colorbar;
    ylabel(c, 'dF/F');
    xlabel('Time (frames)');
    ylabel('ROI (ranked by max activity)');
    title(sprintf('Drug - %s\n(%d/%d neurons shown)', stim_type, ...
        length(display_idx_drug), size(drug_heatmap, 1)), 'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle(sprintf('Population Heatmaps: %s', stim_type), 'FontSize', 12, 'FontWeight', 'bold');
    
    % Store sort indices for later use
    data.population_sort_idx_baseline = sort_idx_base;
    data.population_sort_idx_drug = sort_idx_drug;
    master_data.(field_name) = data;
end

%% ====================== GENERATE MATCHED ROI HEATMAPS ======================

if matched_rois_available
    fprintf('\n========== GENERATING MATCHED ROI HEATMAPS ==========\n');
    
    for field_idx = 1:length(stim_fields)
        field_name = stim_fields{field_idx};
        data = master_data.(field_name);
        stim_type = data.stimulus_type;
        
        fprintf('Generating matched ROI heatmap for "%s"\n', stim_type);
        
        % Convert cell arrays to 3D matrices first (for matched ROI extraction)
        baseline_3d = cell_responses_to_3d(data.baseline_responses);
        drug_3d = cell_responses_to_3d(data.drug_responses);
        
        if isempty(baseline_3d) || isempty(drug_3d)
            fprintf('  (Skipping - insufficient data)\n');
            continue;
        end
        
        % Extract matched ROI responses
        baseline_matched = baseline_3d(base_match_idx, :, :);
        drug_matched = drug_3d(drug_match_idx, :, :);
        
        % Average across presentations
        baseline_matched_avg = squeeze(mean(baseline_matched, 3));
        drug_matched_avg = squeeze(mean(drug_matched, 3));
        
        % Ensure 2D
        if isvector(baseline_matched_avg)
            baseline_matched_avg = baseline_matched_avg(:)';
        end
        if isvector(drug_matched_avg)
            drug_matched_avg = drug_matched_avg(:)';
        end
        
        % Rank baseline by max activity
        [~, sort_idx] = sort(max(baseline_matched_avg, [], 2), 'descend');
        
        % Create figure
        fig = figure('Position', [100 100 1400 600], 'NumberTitle', 'off', ...
            'Name', sprintf('Matched ROI Heatmap: %s', stim_type));
        
        % BASELINE MATCHED - Ranked by activity
        subplot(1, 2, 1);
        imagesc(baseline_matched_avg(sort_idx, :));
        eval(['colormap(', colormap_type, ');']);
        set(gca, 'YDir', 'reverse');
        caxis(caxis_lim);
        colorbar;
        xlabel('Time (frames)');
        ylabel('Matched ROI (ranked by baseline activity)');
        title(sprintf('Baseline - %s\n(%d matched ROIs)', stim_type, n_matched), 'FontWeight', 'bold');
        set(gca, 'LineWidth', 1.5, 'FontSize', 10);
        
        % DRUG MATCHED - SAME ORDER AS BASELINE
        subplot(1, 2, 2);
        imagesc(drug_matched_avg(sort_idx, :));
        eval(['colormap(', colormap_type, ');']);
        set(gca, 'YDir', 'reverse');
        caxis(caxis_lim);
        c = colorbar;
        ylabel(c, 'dF/F');
        xlabel('Time (frames)');
        ylabel('Matched ROI');
        title(sprintf('Drug - %s\n(same ROI order)', stim_type), 'FontWeight', 'bold');
        set(gca, 'LineWidth', 1.5, 'FontSize', 10);
        
        sgtitle(sprintf('Matched ROI Heatmaps: %s', stim_type), 'FontSize', 12, 'FontWeight', 'bold');
        
        % Store for analysis
        data.matched_baseline_avg = baseline_matched_avg;
        data.matched_drug_avg = drug_matched_avg;
        data.matched_sort_idx = sort_idx;
        master_data.(field_name) = data;
    end
end

%% ====================== STIMULUS-SPECIFIC ANALYSIS ======================

fprintf('\n========== STIMULUS-SPECIFIC ANALYSIS ==========\n');

% For each stimulus type, perform specialized analysis

for field_idx = 1:length(stim_fields)
    field_name = stim_fields{field_idx};
    data = master_data.(field_name);
    stim_type = data.stimulus_type;
    
    fprintf('\nAnalyzing "%s"...\n', stim_type);
    
    switch stim_type
        case 'grating'
            analyze_grating_stimulus(B.Stimuli, D.Stimuli, base_dff, drug_dff, ...
                data, base_match_idx, drug_match_idx, matched_rois_available);
            
        case 'moving_bar'
            analyze_moving_bar_stimulus(B.Stimuli, D.Stimuli, base_dff, drug_dff, ...
                data, base_match_idx, drug_match_idx, matched_rois_available);
            
        case 'full_field_flash'
            analyze_full_field_flash_stimulus(B.Stimuli, D.Stimuli, base_dff, drug_dff, ...
                data, base_match_idx, drug_match_idx, matched_rois_available);
            
        case 'spontaneous'
            % Spontaneous activity: no specific trigger, but can analyze statistics
            fprintf('  (Spontaneous activity - no stimulus-triggered analysis)\n');
            
        otherwise
            fprintf('  (No specific analysis for this stimulus type)\n');
    end
end

fprintf('\n========== ANALYSIS COMPLETE ==========\n');
fprintf('Master data saved in "master_data" struct.\n');
fprintf('Use master_data.stimulus_field.baseline/drug_responses for further analysis.\n');

%% ====================== LOCAL FUNCTIONS ======================

function dff = get_calcium_data(ca_data, ca_type)
    % Extract calcium data by type
    switch ca_type
        case 1
            dff = ca_data.Ca_dFF;
        case 2
            dff = ca_data.Ca_deconvolved;
        case 3
            dff = ca_data.Ca_F;
        otherwise
            error('Unknown calcium type: %d', ca_type);
    end
end

function [responses, frame_ranges, time_ranges, properties] = extract_full_stimulus_responses(Stimuli, dff, stim_type, TimeCa)
    % Extract full dF/F responses for all presentations of a stimulus type
    % 
    % MATCHED TO USER'S PROVEN METHOD:
    % - Uses TimeStimulusFrame(1) as start time
    % - Calculates end time: stim_start + (stimulus_trial_t * trials)
    % - Matches against TimeCa(1,:) using strict > and < (excludes boundaries)
    % - Creates trial-chunked organization
    %
    % Inputs:
    %   Stimuli: stimulus array from preprocessed file
    %   dff: calcium data (n_neurons × n_frames)
    %   stim_type: stimulus type to extract (string)
    %   TimeCa: time matrix from preprocessed file (2×n_frames, row 1 = time)
    %
    % Returns:
    %   responses: cell array, each element is (n_neurons × stim_duration_frames)
    %   frame_ranges: (n_presentations × 2) [frame_start, frame_end]
    %   time_ranges: (n_presentations × 2) [time_start, time_end]
    %   properties: struct array with stimulus metadata
    
    n_neurons = size(dff, 1);
    n_frames_max = size(dff, 2);
    
    % Extract time vector from TimeCa (row 1 contains time values)
    if size(TimeCa, 1) < 1 || ~ismatrix(TimeCa)
        error('TimeCa must be at least 1×n_frames (preferably 2×n_frames)');
    end
    if isvector(TimeCa)
        time_vector = TimeCa(:)';
    else
        time_vector = TimeCa(1, :);  % First row contains time values
    end
    
    if length(time_vector) ~= n_frames_max
        error('TimeCa row 1 length (%d) must match dff columns (%d)', length(time_vector), n_frames_max);
    end
    
    responses = {};
    frame_ranges = [];
    time_ranges = [];
    properties = {};
    
    pres_count = 0;
    
    for i = 1:numel(Stimuli)
        if ~isfield(Stimuli(i), 'type') || ~strcmp(Stimuli(i).type, stim_type)
            continue;
        end
        
        if ~isfield(Stimuli(i), 'TimeStimulusFrame') || isempty(Stimuli(i).TimeStimulusFrame)
            continue;
        end
        
        % === MATCHED TO USER'S PROVEN METHOD ===
        % Extract start time from FIRST element of TimeStimulusFrame
        stim_time_values = Stimuli(i).TimeStimulusFrame;
        time_start = stim_time_values(1);  % Use FIRST value as start
        
        % Calculate end time explicitly: start + (duration × num_trials)
        % This matches user's: stim_end = stim_start + stim_total_time
        if isfield(Stimuli(i), 'stimulus_trial_t') && isfield(Stimuli(i), 'trials')
            stim_total_time = Stimuli(i).stimulus_trial_t * Stimuli(i).trials;
            time_end = time_start + stim_total_time;
        else
            % Fallback if metadata missing
            time_end = max(stim_time_values(:));
        end
        
        % Find frame indices using STRICT inequalities (exclude boundaries)
        % Matches user's: find(TimeCa(1,:) > stim_start & TimeCa(1,:) < stim_end)
        frame_idx = find(time_vector > time_start & time_vector < time_end);
        
        % Validate
        if isempty(frame_idx) || length(frame_idx) < 2
            continue;
        end
        
        pres_count = pres_count + 1;
        
        % Extract response for all neurons using found frame indices
        st_frame = frame_idx(1);
        en_frame = frame_idx(end);
        response_matrix = dff(:, st_frame:en_frame);
        responses{pres_count, 1} = response_matrix;
        
        % Store frame indices (for reference)
        frame_ranges = [frame_ranges; st_frame, en_frame];
        
        % Store actual time ranges (matched to user's approach)
        time_ranges = [time_ranges; time_start, time_end];
        
        % Store properties
        props = struct();
        props.stimulus_index = i;
        props.frame_start = st_frame;
        props.frame_end = en_frame;
        props.time_start_sec = time_start;
        props.time_end_sec = time_end;
        props.duration_frames = en_frame - st_frame + 1;
        props.duration_sec = time_end - time_start;
        props.n_frames_matched = length(frame_idx);
        
        % Extract metadata
        if isfield(Stimuli(i), 'stimulus_trial_t')
            props.stimulus_trial_t = Stimuli(i).stimulus_trial_t;
        end
        if isfield(Stimuli(i), 'trials')
            props.trials = Stimuli(i).trials;
        end
        if isfield(Stimuli(i), 'specParams') && isstruct(Stimuli(i).specParams)
            props.spec_params = Stimuli(i).specParams;
        end
        
        properties{pres_count, 1} = props;
    end
end

function avg_response = average_stimulus_presentations(response_cell)
    % Average response across presentations (handles variable time lengths)
    % Input: response_cell - cell array where each element is (n_neurons × time)
    %
    % For simplicity, interpolate all to same length then average
    
    if isempty(response_cell)
        avg_response = [];
        return;
    end
    
    n_neurons = size(response_cell{1}, 1);
    
    % Find maximum time length
    time_lengths = cellfun(@(x) size(x, 2), response_cell);
    max_time = max(time_lengths);
    
    % Interpolate all to max length
    interpolated = zeros(n_neurons, max_time, length(response_cell));
    
    for i = 1:length(response_cell)
        original_time = size(response_cell{i}, 2);
        if original_time < max_time
            % Interpolate
            x_orig = linspace(0, 1, original_time);
            x_new = linspace(0, 1, max_time);
            interpolated(:, :, i) = interp1(x_orig, response_cell{i}', x_new)';
        else
            interpolated(:, :, i) = response_cell{i}(:, 1:max_time);
        end
    end
    
    % Average across presentations
    avg_response = mean(interpolated, 3);
end

function response_3d = cell_responses_to_3d(response_cell)
    % Convert cell array of responses to 3D matrix
    % Interpolates all presentations to common time length
    % 
    % Input: response_cell - cell array, each element is (n_neurons × time_i)
    % Output: response_3d - (n_neurons × max_time × n_presentations)
    
    if isempty(response_cell)
        response_3d = [];
        return;
    end
    
    n_neurons = size(response_cell{1}, 1);
    n_presentations = length(response_cell);
    
    % Find maximum time length
    time_lengths = cellfun(@(x) size(x, 2), response_cell);
    max_time = max(time_lengths);
    
    % Initialize 3D array
    response_3d = zeros(n_neurons, max_time, n_presentations);
    
    % Fill 3D array with interpolated responses
    for i = 1:n_presentations
        original_time = size(response_cell{i}, 2);
        
        if original_time == max_time
            % Already correct length
            response_3d(:, :, i) = response_cell{i};
        else
            % Interpolate to max_time
            x_orig = linspace(0, 1, original_time);
            x_new = linspace(0, 1, max_time);
            response_3d(:, :, i) = interp1(x_orig, response_cell{i}', x_new)';
        end
    end
end

function analyze_grating_stimulus(B_Stimuli, D_Stimuli, base_dff, drug_dff, ...
    data, base_match_idx, drug_match_idx, matched_available)
    % Analyze grating stimulus: find preferred directions, average responses
    
    fprintf('    - Extracting direction preferences\n');
    
    % Find all directions
    directions = {};
    for i = 1:length(B_Stimuli)
        if isfield(B_Stimuli(i), 'type') && strcmp(B_Stimuli(i).type, 'grating')
            if isfield(B_Stimuli(i), 'specParams') && isfield(B_Stimuli(i).specParams, 'direction')
                dir_val = B_Stimuli(i).specParams.direction;
                directions{end+1} = dir_val;
            end
        end
    end
    
    if isempty(directions)
        fprintf('      (Could not extract direction information)\n');
        return;
    end
    
    unique_dirs = unique(directions);
    fprintf('    - Found %d unique directions: %s\n', length(unique_dirs), sprintf('%.0f° ', [unique_dirs{:}]));
    
    % For each direction, plot stimulus-triggered average at onset
    % This would require syncing to stimulus onset - implementation depends on
    % exact requirements for grating onset detection
    fprintf('    - Direction-tuning analysis: available but requires stimulus onset detection\n');
end

function analyze_moving_bar_stimulus(B_Stimuli, D_Stimuli, base_dff, drug_dff, ...
    data, base_match_idx, drug_match_idx, matched_available)
    % Analyze moving bar stimulus: extract stimulus-triggered responses at onset
    
    fprintf('    - Extracting moving bar onset responses\n');
    
    % Extract onset responses (first ~100 frames of each stimulus)
    onset_window = 100;
    
    baseline_onsets = {};
    drug_onsets = {};
    
    for i = 1:length(data.baseline_responses)
        onset = data.baseline_responses{i}(:, 1:min(onset_window, size(data.baseline_responses{i}, 2)));
        baseline_onsets{i} = onset;
    end
    
    for i = 1:length(data.drug_responses)
        onset = data.drug_responses{i}(:, 1:min(onset_window, size(data.drug_responses{i}, 2)));
        drug_onsets{i} = onset;
    end
    
    fprintf('    - Stimulus-triggered average: extracted %d baseline and %d drug presentations\n', ...
        length(baseline_onsets), length(drug_onsets));
end

function analyze_full_field_flash_stimulus(B_Stimuli, D_Stimuli, base_dff, drug_dff, ...
    data, base_match_idx, drug_match_idx, matched_available)
    % Analyze full-field flash: plot individual trial responses + averages
    
    fprintf('    - Analyzing full-field flash responses\n');
    
    % Create stimulus-triggered average plot
    fig = figure('Position', [100 100 1200 700], 'NumberTitle', 'off', ...
        'Name', 'Full-Field Flash Analysis');
    
    % Subplot 1: Baseline individual trials + mean
    subplot(2, 2, 1);
    baseline_responses = data.baseline_responses;
    colors_base = repmat([1 0.4 0.4], length(baseline_responses), 1);
    
    max_len = max(cellfun(@(x) size(x, 2), baseline_responses));
    
    for trial = 1:min(length(baseline_responses), 5)  % Show up to 5 trials
        resp = baseline_responses{trial};
        % Get top 10 ROIs by max activity
        [~, top_idx] = sort(max(resp, [], 2), 'descend');
        top_resp = resp(top_idx(1:min(10, size(resp, 1))), :);
        mean_resp = mean(top_resp, 1);
        plot(mean_resp, 'Color', colors_base(trial, :) * 0.7, 'LineWidth', 1.5);
        hold on;
    end
    xlabel('Time (frames)');
    ylabel('dF/F');
    title('Baseline - Individual Flash Responses (top 10 ROIs)');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    grid on;
    
    % Subplot 2: Drug individual trials + mean
    subplot(2, 2, 2);
    drug_responses = data.drug_responses;
    colors_drug = repmat([0.4 0.4 1], length(drug_responses), 1);
    
    for trial = 1:min(length(drug_responses), 5)
        resp = drug_responses{trial};
        [~, top_idx] = sort(max(resp, [], 2), 'descend');
        top_resp = resp(top_idx(1:min(10, size(resp, 1))), :);
        mean_resp = mean(top_resp, 1);
        plot(mean_resp, 'Color', colors_drug(trial, :) * 0.7, 'LineWidth', 1.5);
        hold on;
    end
    xlabel('Time (frames)');
    ylabel('dF/F');
    title('Drug - Individual Flash Responses (top 10 ROIs)');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    grid on;
    
    % Subplot 3: Baseline heatmap
    subplot(2, 2, 3);
    baseline_avg = average_stimulus_presentations(baseline_responses);
    [~, sort_idx] = sort(max(baseline_avg, [], 2), 'descend');
    display_idx = sort_idx(1:min(50, length(sort_idx)));
    imagesc(baseline_avg(display_idx, :));
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
    caxis([0 5]);
    colorbar;
    xlabel('Time (frames)');
    ylabel('ROI (ranked by activity)');
    title('Baseline - Flash Response Heatmap');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % Subplot 4: Drug heatmap
    subplot(2, 2, 4);
    drug_avg = average_stimulus_presentations(drug_responses);
    [~, sort_idx] = sort(max(drug_avg, [], 2), 'descend');
    display_idx = sort_idx(1:min(50, length(sort_idx)));
    imagesc(drug_avg(display_idx, :));
    colormap(flipud(gray));
    set(gca, 'YDir', 'reverse');
    caxis([0 5]);
    colorbar;
    xlabel('Time (frames)');
    ylabel('ROI (ranked by activity)');
    title('Drug - Flash Response Heatmap');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle('Full-Field Flash: Stimulus-Triggered Responses', 'FontSize', 12, 'FontWeight', 'bold');
end
