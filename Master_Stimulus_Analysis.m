%% Analysis for 2P data with ROI MATCHED OPTION and multi-recording comparison
%
% This script creates master response matrices for all neurons across stimuli,
% generates population and matched-ROI heatmaps, analyzes stimulus-specific
% responses, and creates activity distribution histograms and violin plots.
%
% STRUCTURE:
%   1. Load and organize data
%   2. Extract master response matrices (all neurons, full stimulus duration)
%   3. Generate 4-level progressive heatmaps (per stimulus)
%   4. Generate activity distribution histograms (all neurons and matched neurons)
%   5. Generate violin plots (all neurons and matched neurons)
%   6. Stimulus-specific analysis 

clear; clc; close all;

% Add violin plot utilities to path
addpath(genpath(fullfile(pwd, 'violin_plot_utils')));

%% ====================== CONFIGURATION ======================

% Data files
%baseline_file = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat";
baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat";
%drug_file = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat";
drug_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat";
%roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_V2.mat";
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_0312_V1.mat";

% Analysis parameters
ca_type = 1;                  % 1=FVDff, 2=deconvolved, 3=F
selected_plane_idx = 1;           % Plane index (1-based)
stimulus_types_to_analyze = {'spontaneous','sparse_local_global_flashes', 'checkers2', 'grating', 'looming'};

% Visualization
max_neurons_display = 1000;    % For heatmap readability
caxis_lim = [0, 5];           % dF/F range for heatmaps
colormap_type = 'flipud(gray)';


%% ====================== LOAD DATA ======================

fprintf('\n========== DATA LOADING ==========\n');

fprintf('Loading baseline: %s\n', baseline_file);
B = load(baseline_file);

% Load baseline Ca_dFF and centroid data
base_dff_full = B.CaData(1).Ca_dFF;  % ALL ROIs (all planes combined)
centroid = B.CaData(1).Ca_centroid_voxel;
centroidZ = centroid(:, 3);

% Determine number of planes from unique Z values in centroid
unique_planes = unique(centroidZ);
n_planes = length(unique_planes);

% Create separate dFF matrix for each plane
base_dff_plane = cell(n_planes, 1);
for plane_idx = 1:n_planes
    plane_z_value = unique_planes(plane_idx);
    roi_indices = find(centroidZ == plane_z_value);
    base_dff_plane{plane_idx} = base_dff_full(roi_indices, :);
end

% Get data for selected plane
selected_roi_idx = find(centroidZ == unique_planes(selected_plane_idx));
if isempty(selected_roi_idx)
    error('ERROR: Selected plane %d not found in baseline data', unique_planes(selected_plane_idx));
end
base_dff = base_dff_plane{selected_plane_idx};
n_rois_selected_plane = size(base_dff, 1);

n_rois_base_full = size(base_dff_full, 1);  % Total ROIs across all planes
n_frames_base = size(base_dff_full, 2);


% Verification: check that subsetting preserves calcium values correctly
test_roi_global = selected_roi_idx(1);  % First ROI in selected plane
test_roi_local = 1;  % This should be row 1 in base_dff_plane
val_global = base_dff_full(test_roi_global, 1);  % First timepoint
val_local = base_dff_plane{selected_plane_idx}(test_roi_local, 1);
fprintf('Verification check: base_dff_full(%d, 1) = %.4f, base_dff_plane{%d}(%d, 1) = %.4f\n', ...
    test_roi_global, val_global, selected_plane_idx, test_roi_local, val_local);
if abs(val_global - val_local) < 1e-10
    fprintf('Values match - subsetting is correct!\n');
else
    fprintf('ERROR: Values do not match - subsetting logic is broken!\n');
end


fprintf('Loading drug: %s\n', drug_file);
D = load(drug_file);

% Load drug Ca_dFF and centroid data
drug_dff_full = D.CaData(1).Ca_dFF;  % ALL ROIs (all planes combined)
drug_centroid = D.CaData(1).Ca_centroid_voxel;
drug_centroidZ = drug_centroid(:, 3);

% Determine number of planes from unique Z values in drug centroid
unique_planes_drug = unique(drug_centroidZ);
n_planes_drug = length(unique_planes_drug);

% Create separate dFF matrix for each plane
drug_dff_plane = cell(n_planes_drug, 1);
for plane_idx = 1:n_planes_drug
    plane_z_value = unique_planes_drug(plane_idx);
    roi_indices = find(drug_centroidZ == plane_z_value);
    drug_dff_plane{plane_idx} = drug_dff_full(roi_indices, :);
end

% Get data for selected plane
if selected_plane_idx > length(unique_planes_drug)
    error('ERROR: Selected plane %d exceeds number of planes in drug data (%d)', selected_plane_idx, n_planes_drug);
end
drug_selected_roi_idx = find(drug_centroidZ == unique_planes_drug(selected_plane_idx));
if isempty(drug_selected_roi_idx)
    error('ERROR: Selected plane %d not found in drug data', unique_planes_drug(selected_plane_idx));
end
drug_dff = drug_dff_plane{selected_plane_idx};
n_rois_drug_selected_plane = size(drug_dff, 1);

n_rois_drug_full = size(drug_dff_full, 1);  % Total ROIs across all planes
n_frames_drug = size(drug_dff_full, 2);

fprintf('Drug plane details:\n');
for plane_idx = 1:n_planes_drug
    plane_z = unique_planes_drug(plane_idx);
    n_rois_plane = size(drug_dff_plane{plane_idx}, 1);
    fprintf('  Plane %d: %d ROIs\n', plane_z, n_rois_plane);
end
fprintf('Drug selected plane %d: %d ROIs (global indices: %d to %d)\n', ...
    selected_plane_idx, n_rois_drug_selected_plane, min(drug_selected_roi_idx), max(drug_selected_roi_idx));

% Load ROI matching if available
matched_rois_available = false;  % Default: no matched ROIs
base_match_idx_global = [];      % Global ROI indices (from file)
drug_match_idx_global = [];      % Global ROI indices (from file)
base_match_idx_local = [];       % Local row indices in base_dff
drug_match_idx_local = [];       % Local row indices in drug_dff
n_matched = 0;

if ~isempty(roi_match_file) && isfile(roi_match_file)
    fprintf('Loading ROI matches: %s\n', roi_match_file);
    M = load(roi_match_file);
    if isfield(M.roiMatchData, 'allSessionMapping')
        base_match_idx_global = M.roiMatchData.allSessionMapping(:, 1);
        drug_match_idx_global = M.roiMatchData.allSessionMapping(:, 2);
        fprintf('Found %d matched ROI pairs (global indices)\n', length(base_match_idx_global));
        
        % Convert global ROI indices to LOCAL indices within selected plane
        % For baseline: find which row index in base_dff corresponds to each global index
        [~, base_match_idx_local] = ismember(base_match_idx_global, selected_roi_idx);
        [~, drug_match_idx_local] = ismember(drug_match_idx_global, drug_selected_roi_idx);
        
        % Filter to keep only matches where BOTH neurons are in the selected plane
        valid_matches = (base_match_idx_local > 0) & (drug_match_idx_local > 0);
        base_match_idx_local = base_match_idx_local(valid_matches);
        drug_match_idx_local = drug_match_idx_local(valid_matches);
        base_match_idx_global = base_match_idx_global(valid_matches);
        drug_match_idx_global = drug_match_idx_global(valid_matches);
        
        n_matched = length(base_match_idx_local);
        
        if n_matched > 0
            matched_rois_available = true;
            fprintf('\n--- MATCHED ROI ANALYSIS ---\n');
            fprintf('Baseline selected plane (Z=%d): %d local neuron indices found\n', ...
                selected_plane_idx, n_rois_selected_plane);
            fprintf('Drug selected plane (Z=%d): %d local neuron indices found\n', ...
                selected_plane_idx, n_rois_drug_selected_plane);
            fprintf('Matched neuron pairs in selected plane: %d\n', n_matched);
            fprintf('Local baseline match indices: %d to %d\n', min(base_match_idx_local), max(base_match_idx_local));
            fprintf('Local drug match indices: %d to %d\n', min(drug_match_idx_local), max(drug_match_idx_local));
            fprintf('Ready for matched neuron analysis\n');
            
            % Verify conversion: check first matched pair
            if n_matched > 0
                base_global = base_match_idx_global(1);
                base_local = base_match_idx_local(1);
                drug_global = drug_match_idx_global(1);
                drug_local = drug_match_idx_local(1);
                fprintf('\nVerification (first matched pair):\n');
                fprintf('  Baseline: global ROI %d -> local row %d in base_dff(%d, :)\n', ...
                    base_global, base_local, size(base_dff, 1));
                fprintf('Drug: global ROI %d -> local row %d in drug_dff(%d, :)\n', ...
                    drug_global, drug_local, size(drug_dff, 1));
                
                % Check that the global index is actually in selected_roi_idx
                if selected_roi_idx(base_local) == base_global
                    fprintf('Baseline conversion verified\n');
                else
                    fprintf('ERROR: Baseline conversion mismatch!\n');
                end
                if drug_selected_roi_idx(drug_local) == drug_global
                    fprintf('Drug conversion verified\n');
                else
                    fprintf('Drug conversion mismatch!\n');
                end
            end
        else
            fprintf('WARNING: No matched ROIs found in selected plane %d!\n', selected_plane_idx);
            fprintf('All %d matched pairs are either not in baseline plane or drug plane\n', ...
                length(base_match_idx_global));
        end
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
    
    % CRITICAL: Pass FULL dataset (base_dff_full) with global ROI indices to maintain numbering
    % This ensures global ROI identities are preserved for matching later
    fprintf('  Baseline: Extracting from base_dff_full (%d×%d) using selected_roi_idx (%d ROIs)\n', ...
        size(base_dff_full, 1), size(base_dff_full, 2), length(selected_roi_idx));
    fprintf('  Global indices: %d to %d\n', min(selected_roi_idx), max(selected_roi_idx));
    
    % Extract stimulus responses for selected plane only
    % Pass full dataset and roi_indices to maintain global ROI numbering
    [base_responses, base_frame_ranges, base_time_ranges, base_properties] = extract_full_stimulus_responses(...
        B.Stimuli, base_dff_full, stim_type, B.TimeCa, selected_roi_idx);
    
    fprintf('  Drug: Extracting from drug_dff_full (%d×%d) using drug_selected_roi_idx (%d ROIs)\n', ...
        size(drug_dff_full, 1), size(drug_dff_full, 2), length(drug_selected_roi_idx));
    fprintf('  Global indices: %d to %d\n', min(drug_selected_roi_idx), max(drug_selected_roi_idx));
    
    [drug_responses, drug_frame_ranges, drug_time_ranges, drug_properties] = extract_full_stimulus_responses(...
        D.Stimuli, drug_dff_full, stim_type, D.TimeCa, drug_selected_roi_idx);
    
    % Verify extraction
    if ~isempty(base_responses) && ~isempty(drug_responses)
        fprintf('Extraction successful:\n');
        fprintf('Baseline responses: %d presentations, %d neurons per presentation\n', ...
            length(base_responses), size(base_responses{1}, 1));
        fprintf('Drug responses: %d presentations, %d neurons per presentation\n', ...
            length(drug_responses), size(drug_responses{1}, 1));
    else
        fprintf('WARNING: Empty responses!\n');
    end
    
    % Store in master struct
    field_name = matlab.lang.makeValidName(stim_type);
    master_data.(field_name).baseline_responses = base_responses;
    master_data.(field_name).drug_responses = drug_responses;
    master_data.(field_name).baseline_frame_ranges = base_frame_ranges;
    master_data.(field_name).drug_frame_ranges = drug_frame_ranges;
    master_data.(field_name).baseline_time_ranges = base_time_ranges; 
    master_data.(field_name).drug_time_ranges = drug_time_ranges;      
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

%% ====================== HEATMAPS======================

fprintf('\n========== LEVEL 1: Full dFF Heatmaps (all neurons, all timepoints) ==========\n');

% ===== Axis limits for Level 1 (adjust as needed) =====
L1_caxis_lim = [0, 5];  % Color axis limits for raw dF/F

stim_fields = fieldnames(master_data);

for field_idx = 1:length(stim_fields)
    field_name = stim_fields{field_idx};
    data = master_data.(field_name);
    stim_type = data.stimulus_type;
    
    fprintf('Level 1: Generating full dFF heatmap for "%s"\n', stim_type);
    
    % Extract first complete presentation (no averaging, no concatenation)
    % This shows the full response to ONE stimulus event
    if ~isempty(data.baseline_responses)
        baseline_full = data.baseline_responses{1};  % First presentation
        n_baseline_pres = length(data.baseline_responses);
    else
        baseline_full = [];
        n_baseline_pres = 0;
    end
    
    if ~isempty(data.drug_responses)
        drug_full = data.drug_responses{1};  % First presentation
        n_drug_pres = length(data.drug_responses);
    else
        drug_full = [];
        n_drug_pres = 0;
    end
    
    if isempty(baseline_full) || isempty(drug_full)
        fprintf('  Skipping - insufficient data\n');
        continue;
    end
    
    % Create figure: baseline and drug side-by-side
    fig = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
        'Name', sprintf('Level1_FullDFF: %s', stim_type));
    
    % BASELINE - First presentation, ALL neurons, ALL timepoints
    subplot(1, 2, 1);
    imagesc(baseline_full);
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(L1_caxis_lim);
    colorbar;
    xlabel('Time (frames)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('ROI', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('Baseline - %s (Presentation 1 of %d)\n%d neurons × %d frames', ...
        stim_type, n_baseline_pres, size(baseline_full, 1), size(baseline_full, 2)), 'FontWeight', 'bold', 'FontSize', 11);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % DRUG - First presentation, ALL neurons, ALL timepoints
    subplot(1, 2, 2);
    imagesc(drug_full);
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(L1_caxis_lim);
    c = colorbar;
    ylabel(c, 'dF/F', 'FontSize', 10);
    xlabel('Time (frames)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('ROI', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('Drug - %s (Presentation 1 of %d)\n%d neurons × %d frames', ...
        stim_type, n_drug_pres, size(drug_full, 1), size(drug_full, 2)), 'FontWeight', 'bold', 'FontSize', 11);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle(sprintf('LEVEL 1 - Single Stimulus Presentation: %s', stim_type), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % Store full heatmaps
    data.full_baseline_heatmap = baseline_full;
    data.full_drug_heatmap = drug_full;
    master_data.(field_name) = data;
end

%% ====================== LEVEL 2: SELECTED ROI HEATMAPS (RAW VALUES, NO SORTING) ======================

fprintf('\n========== LEVEL 2: Selected ROI Heatmaps (raw values, unsorted) ==========\n');
fprintf('Zoom in on specific timeframes of the stimulus response for detailed inspection\n\n');

% ===== Parameters for Level 2 (ADJUST THESE TO ZOOM INTO SPECIFIC TIMEFRAMES) =====
L2_n_select = min(1000, size(baseline_full, 1));          % Number of ROIs to select
L2_frame_start = 1;                                       % Frame window START - adjust to zoom into response
L2_frame_end = size(baseline_full, 2);                    % Frame window END - adjust to zoom into response
L2_caxis_lim = [0, 5];                                    % Color axis limits

fprintf('Level 2 settings:\n');
fprintf('  ROIs: %d neurons\n', L2_n_select);
fprintf('  Frames: %d:%d (total %d frames available)\n', L2_frame_start, L2_frame_end, size(baseline_full, 2));
fprintf('  To zoom into specific response, adjust L2_frame_start and L2_frame_end above\n\n');

for field_idx = 1:length(stim_fields)
    field_name = stim_fields{field_idx};
    data = master_data.(field_name);
    stim_type = data.stimulus_type;
    
    fprintf('Level 2: Generating selected ROI heatmap for "%s"\n', stim_type);
    
    baseline_full = data.full_baseline_heatmap;
    drug_full = data.full_drug_heatmap;
    
    % Select first N ROIs (no sorting, just take first N)
    n_select = min(L2_n_select, size(baseline_full, 1));
    select_idx = 1:n_select;
    
    % Extract timeframe window
    frame_start = max(1, L2_frame_start);
    frame_end = min(size(baseline_full, 2), L2_frame_end);
    
    baseline_selected = baseline_full(select_idx, frame_start:frame_end);
    drug_selected = drug_full(select_idx, frame_start:frame_end);
    
    % Create figure
    fig = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
        'Name', sprintf('Level2_SelectedROI: %s', stim_type));
    
    % BASELINE - Selected ROIs, raw values, NOT sorted
    subplot(1, 2, 1);
    imagesc(baseline_selected);
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(L2_caxis_lim);
    colorbar;
    xlabel('Time (frames)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('ROI', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('Baseline - %s (Selected, raw)\n%d ROIs × %d timepoints (frames %d:%d)', ...
        stim_type, n_select, size(baseline_selected, 2), frame_start, frame_end), 'FontWeight', 'bold', 'FontSize', 11);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % DRUG - Same ROIs, raw values, NOT sorted
    subplot(1, 2, 2);
    imagesc(drug_selected);
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(L2_caxis_lim);
    c = colorbar;
    ylabel(c, 'dF/F', 'FontSize', 10);
    xlabel('Time (frames)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('ROI', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('Drug - %s (Selected, raw)\n%d ROIs × %d timepoints (frames %d:%d)', ...
        stim_type, n_select, size(drug_selected, 2), frame_start, frame_end), 'FontWeight', 'bold', 'FontSize', 11);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle(sprintf('LEVEL 2 - Selected ROI (raw, unsorted): %s', stim_type), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % Store for later use
    data.selected_baseline_raw = baseline_selected;
    data.selected_drug_raw = drug_selected;
    data.selected_roi_idx = select_idx;
    data.L2_frame_range = [frame_start, frame_end];
    
    master_data.(field_name) = data;
end

%% ====================== LEVEL 3: NORMALIZED HEATMAPS (SORTED BY MAX ACTIVITY) ======================

fprintf('\n========== LEVEL 3: Normalized Heatmaps (sorted by max activity) ==========\n');

% ===== Axis limits for Level 3 (adjust as needed) =====
L3_caxis_lim = [0, 15];  % Color axis limits for normalized dF/F

for field_idx = 1:length(stim_fields)
    field_name = stim_fields{field_idx};
    data = master_data.(field_name);
    stim_type = data.stimulus_type;
    
    fprintf('Level 3: Generating normalized heatmap for "%s"\n', stim_type);
    
    baseline_selected = data.selected_baseline_raw;
    drug_selected = data.selected_drug_raw;
    
    % Apply normalization: (dFF - min) / (max - min) * max per neuron (row-wise)
    % Handle edge case: if row has all same values, set entire row to that value
    baseline_min = min(baseline_selected, [], 2);
    baseline_max = max(baseline_selected, [], 2);
    baseline_range = baseline_max - baseline_min;
    baseline_range(baseline_range == 0) = 1;  % Avoid division by zero
    baseline_norm = (baseline_selected - baseline_min) ./ baseline_range .* baseline_max;
    
    drug_min = min(drug_selected, [], 2);
    drug_max = max(drug_selected, [], 2);
    drug_range = drug_max - drug_min;
    drug_range(drug_range == 0) = 1;  % Avoid division by zero
    drug_norm = (drug_selected - drug_min) ./ drug_range .* drug_max;
    
    % Sort BOTH baseline and drug by their respective max activity
    % Handle NaN values from normalization
    baseline_max_vals = max(baseline_norm, [], 2);
    baseline_max_vals(isnan(baseline_max_vals)) = 0;  % Replace NaN with 0
    [~, sort_idx_baseline] = sort(baseline_max_vals, 'descend');
    
    drug_max_vals = max(drug_norm, [], 2);
    drug_max_vals(isnan(drug_max_vals)) = 0;  % Replace NaN with 0
    [~, sort_idx_drug] = sort(drug_max_vals, 'descend');
    
    fig = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
        'Name', sprintf('Level3_Normalized: %s', stim_type));
    
    % BASELINE - Normalized, sorted by baseline max
    subplot(1, 2, 1);
    imagesc(baseline_norm(sort_idx_baseline, :));
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(L3_caxis_lim);
    colorbar;
    xlabel('Time (frames)');
    ylabel('ROI (sorted by baseline max)');
    title(sprintf('Baseline - %s (Normalized, sorted)\nFormula: (dFF-min)/(max-min)*max', stim_type), ...
        'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % DRUG - Normalized, sorted by drug max
    subplot(1, 2, 2);
    imagesc(drug_norm(sort_idx_drug, :));
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(L3_caxis_lim);
    c = colorbar;
    ylabel(c, 'dF/F (normalized)');
    xlabel('Time (frames)');
    ylabel('ROI (sorted by drug max)');
    title(sprintf('Drug - %s (Normalized, sorted)\nEach sorted independently by max', stim_type), ...
        'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle(sprintf('LEVEL 3 - Normalized (sorted independently): %s', stim_type), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % Store for later use
    data.normalized_baseline = baseline_norm;
    data.normalized_drug = drug_norm;
    data.norm_sort_idx_baseline = sort_idx_baseline;
    data.norm_sort_idx_drug = sort_idx_drug;
    
    master_data.(field_name) = data;
end

%% ====================== LEVEL 4: MATCHED ROI HEATMAPS ======================

if matched_rois_available && n_matched > 0
    fprintf('\n--- Level 4: Matched ROI Heatmaps ---\n');
    
    for field_idx = 1:length(stim_fields)
        field_name = stim_fields{field_idx};
        data = master_data.(field_name);
        stim_type = data.stimulus_type;
        
        fprintf('Generating matched ROI heatmap for "%s"\n', stim_type);
        
        % Convert response cells to 3D matrices: (n_neurons × time × n_presentations)
        % NOTE: These use LOCAL row indices (rows 1 to n_rois_in_plane)
        baseline_3d = cell_responses_to_3d(data.baseline_responses);
        drug_3d = cell_responses_to_3d(data.drug_responses);
        
        if isempty(baseline_3d) || isempty(drug_3d)
            fprintf('  (Skipping - insufficient data)\n');
            continue;
        end
        
        % Extract matched ROI responses using LOCAL indices
        % base_match_idx_local and drug_match_idx_local are row indices in base_3d/drug_3d
        if any(base_match_idx_local > size(baseline_3d, 1)) || any(drug_match_idx_local > size(drug_3d, 1))
            fprintf('ERROR: Match indices exceed data dimensions!\n');
            continue;
        end
        
        baseline_matched_3d = baseline_3d(base_match_idx_local, :, :);      % (n_matched × time × presentations)
        drug_matched_3d = drug_3d(drug_match_idx_local, :, :);             % (n_matched × time × presentations)
        
        % Average across presentations to get (n_matched × time)
        baseline_matched_avg = squeeze(mean(baseline_matched_3d, 3));
        drug_matched_avg = squeeze(mean(drug_matched_3d, 3));
        
        % Handle edge case of single presentation
        if isvector(baseline_matched_avg)
            baseline_matched_avg = baseline_matched_avg(:);
        end
        if isvector(drug_matched_avg)
            drug_matched_avg = drug_matched_avg(:);
        end
        
        fprintf('  Extracted %d matched neuron pairs\n', size(baseline_matched_avg, 1));
        fprintf('  Dimensions: %d neurons × %d timepoints\n', size(baseline_matched_avg, 1), size(baseline_matched_avg, 2));
        
        % Rank baseline by max activity (handle NaN)
        baseline_max_vals = max(baseline_matched_avg, [], 2);
        baseline_max_vals(isnan(baseline_max_vals)) = 0;
        [~, sort_idx] = sort(baseline_max_vals, 'descend');
        
        % Create figure: baseline and drug side-by-side, SAME NEURON ORDER
        fig4 = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
            'Name', sprintf('Level4_MatchedROI: %s', stim_type));
        
        % BASELINE MATCHED - Ranked by max activity
        subplot(1, 2, 1);
        imagesc(baseline_matched_avg(sort_idx, :));
        eval(['colormap(', colormap_type, ');']);
        set(gca, 'YDir', 'reverse');
        caxis(caxis_lim);
        colorbar;
        xlabel('Time (frames)');
        ylabel('Matched ROI (ranked by baseline max)');
        title(sprintf('Baseline - %s (Matched, sorted)\n%d matched neurons', stim_type, n_matched), ...
            'FontWeight', 'bold');
        set(gca, 'LineWidth', 1.5, 'FontSize', 10);
        
        % DRUG MATCHED - SAME NEURON ORDER AS BASELINE (baseline matched index i -> drug matched index i)
        subplot(1, 2, 2);
        imagesc(drug_matched_avg(sort_idx, :));
        eval(['colormap(', colormap_type, ');']);
        set(gca, 'YDir', 'reverse');
        caxis(caxis_lim);
        c = colorbar;
        ylabel(c, 'dF/F');
        xlabel('Time (frames)');
        ylabel('Matched ROI');
        title(sprintf('Drug - %s (Matched, same order)\nEach row = baseline ROI pair matched to drug', stim_type), ...
            'FontWeight', 'bold');
        set(gca, 'LineWidth', 1.5, 'FontSize', 10);
        
        sgtitle(sprintf('Level 4 - Matched ROI Activity: %s', stim_type), ...
            'FontSize', 12, 'FontWeight', 'bold');
        
        % Add annotation explaining the matching
        annotation('textbox', [0.1 0.02 0.8 0.04], 'String', ...
            sprintf('Each row shows a matched neuron pair: baseline neuron (left) and its paired drug counterpart (right). Neurons ranked by baseline max activity.'), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FitBoxToText', 'on', ...
            'EdgeColor', 'black', 'BackgroundColor', [1 1 0.8]);
        
        % Store matched neuron data for further analysis
        data.matched_baseline_avg = baseline_matched_avg;
        data.matched_drug_avg = drug_matched_avg;
        data.matched_sort_idx = sort_idx;
        data.n_matched_neurons = n_matched;
        data.base_match_idx_local = base_match_idx_local;      % Local row indices for baseline
        data.drug_match_idx_local = drug_match_idx_local;      % Local row indices for drug
        data.base_match_idx_global = base_match_idx_global;    % Global ROI IDs for reference
        data.drug_match_idx_global = drug_match_idx_global;    % Global ROI IDs for reference
        master_data.(field_name) = data;
        
        fprintf('Heatmaps created\n');
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
                data, base_match_idx_local, drug_match_idx_local, matched_rois_available);
            
 %       case 'moving_bar'
%            analyze_moving_bar_stimulus(B.Stimuli, D.Stimuli, base_dff, drug_dff, ...
 %               data, base_match_idx_local, drug_match_idx_local, matched_rois_available);
            
        case 'full_field_flash'
            analyze_full_field_flash_stimulus(B.Stimuli, D.Stimuli, base_dff, drug_dff, ...
                data, base_match_idx_local, drug_match_idx_local, matched_rois_available);

        case 'looming'
            analyze_looming_stimulus(B.Stimuli, base_dff, B.TimeCa, B.Triggers, ...
                D.Stimuli, drug_dff, D.TimeCa, D.Triggers, ...
                base_match_idx_local, drug_match_idx_local, matched_rois_available);

        case 'spontaneous'
            % Spontaneous activity: no specific trigger, but can analyze statistics
            fprintf('  (Spontaneous activity - no stimulus-triggered analysis)\n');
            
        otherwise
            fprintf('  (No specific analysis for this stimulus type)\n');
    end
end

%% ====================== HISTOGRAMS & VIOLIN PLOTS (ALL NEURONS) ======================

fprintf('\n========== GENERATING ACTIVITY DISTRIBUTIONS (ALL NEURONS) ==========\n');

hist_bins = 90;  % Number of bins

for field_idx = 1:length(stim_fields)
    field_name = stim_fields{field_idx};
    data = master_data.(field_name);
    stim_type = data.stimulus_type;
    
    fprintf('Distribution plots for "%s": all %d neurons\n', stim_type, n_rois_selected_plane);
    
    % Extract responses
    baseline_all = [];
    for i = 1:length(data.baseline_responses)
        baseline_all = [baseline_all; data.baseline_responses{i}(:)];
    end
    
    drug_all = [];
    for i = 1:length(data.drug_responses)
        drug_all = [drug_all; data.drug_responses{i}(:)];
    end
    
    % Validate data - remove NaN and Inf
    baseline_all = baseline_all(~isnan(baseline_all) & ~isinf(baseline_all));
    drug_all = drug_all(~isnan(drug_all) & ~isinf(drug_all));
    
    % Skip if insufficient data
    if length(baseline_all) < 2 || length(drug_all) < 2
        fprintf('  WARNING: Insufficient valid data for "%s" (baseline: %d, drug: %d samples)\n', ...
            stim_type, length(baseline_all), length(drug_all));
        continue;
    end
    
    % ===== Create SEPARATE FIGURES with enhanced visualization =====
    
    % ===== FIGURE 1: HISTOGRAM (Large, Full-Width) =====
    fig_hist = figure('Position', [100 150 1400 600], 'NumberTitle', 'off', ...
        'Name', sprintf('Histogram_AllNeurons_%s', stim_type));
    
    combined_data = [baseline_all; drug_all];
    p1 = prctile(combined_data, 1);
    p99 = prctile(combined_data, 99);
    
    % Ensure bin edges have width (avoid zero-width bins for constant data)
    if p1 == p99
        bin_edges = linspace(p1 - 0.5, p1 + 0.5, hist_bins + 1);
    else
        bin_edges = linspace(p1, p99, hist_bins + 1);
    end
    
    hold on;
    histogram(baseline_all, 'Normalization', 'pdf', 'FaceColor', 'r', 'FaceAlpha', 0.5, ...
        'EdgeColor', 'red', 'BinEdges', bin_edges, 'LineWidth', 2, 'DisplayName', 'Baseline');
    histogram(drug_all, 'Normalization', 'pdf', 'FaceColor', 'b', 'FaceAlpha', 0.5, ...
        'EdgeColor', 'blue', 'BinEdges', bin_edges, 'LineWidth', 2, 'DisplayName', 'Drug');
    
    % Add vertical lines for medians
    baseline_median = median(baseline_all);
    drug_median = median(drug_all);
    plot([baseline_median baseline_median], ylim, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Baseline Median');
    plot([drug_median drug_median], ylim, 'b--', 'LineWidth', 2.5, 'DisplayName', 'Drug Median');
    
    ylabel('Probability Density', 'FontSize', 13, 'FontWeight', 'bold');
    xlabel('dF/F', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Distribution of Activity - %s (n=%d neurons)', stim_type, n_rois_selected_plane), ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend('FontSize', 11, 'Location', 'best');
    grid on; set(gca, 'LineWidth', 2, 'FontSize', 11);
    hold off;
    
    % ===== FIGURE 2: VIOLIN PLOT (Large, Full-Width, Enhanced) =====
    fig_violin = figure('Position', [100 150 1400 600], 'NumberTitle', 'off', ...
        'Name', sprintf('Violin_AllNeurons_%s', stim_type));
    
    baseline_mean = mean(baseline_all);
    drug_mean = mean(drug_all);
    
    hold on;
    
    % BASELINE VIOLIN (x=0) - Solid filled with prominent outline
    if length(baseline_all) >= 2
        [f_bl, xi_bl] = ksdensity(baseline_all, 'NumPoints', 150);
        f_bl = f_bl / max(f_bl) * 0.4;  % Violin half-width (0.8 total)
        
        % Draw solid violin fill with DARK, THICK outline for clarity
        x_violin_bl = [f_bl, fliplr(-f_bl)];
        y_violin_bl = [xi_bl, fliplr(xi_bl)];
        patch(x_violin_bl, y_violin_bl, [0.2 0.6 1], 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 2.5);
        
        % Calculate mean and median
        baseline_median = median(baseline_all);
        baseline_mean = mean(baseline_all);
        
        % Add median line (thicker, darker)
        plot([-0.4 0.4], [baseline_median baseline_median], 'k-', 'LineWidth', 3, 'DisplayName', 'Median');
        
        % Add mean line (thinner, different color for distinction)
        plot([-0.4 0.4], [baseline_mean baseline_mean], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean');
        
        % Add scattered points across violin
        n_baseline = length(baseline_all);
        x_jitter_bl = randn(n_baseline, 1) * 0.12;  % Wider horizontal spread
        x_jitter_bl = max(min(x_jitter_bl, 0.4), -0.4);  % Constrain within violin
        
        scatter(x_jitter_bl, baseline_all, 30, [0.3 0.3 0.3], 'o', ...
            'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
    end
    
    % DRUG VIOLIN (x=1) - Solid filled with prominent outline
    if length(drug_all) >= 2
        [f_dr, xi_dr] = ksdensity(drug_all, 'NumPoints', 150);
        f_dr = f_dr / max(f_dr) * 0.4;  % Violin half-width
        
        % Draw solid violin fill with DARK, THICK outline for clarity
        x_violin_dr = [f_dr + 1, fliplr(-f_dr + 1)];
        y_violin_dr = [xi_dr, fliplr(xi_dr)];
        patch(x_violin_dr, y_violin_dr, [1 0.5 0.2], 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 2.5);
        
        % Calculate mean and median
        drug_median = median(drug_all);
        drug_mean = mean(drug_all);
        
        % Add median line (thicker, darker)
        plot([1 - 0.4 1 + 0.4], [drug_median drug_median], 'k-', 'LineWidth', 3, 'DisplayName', 'Median');
        
        % Add mean line (thinner, different color for distinction)
        plot([1 - 0.4 1 + 0.4], [drug_mean drug_mean], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean');
        
        % Add scattered points across violin
        n_drug = length(drug_all);
        x_jitter_dr = 1 + randn(n_drug, 1) * 0.12;  % Wider horizontal spread
        x_jitter_dr = max(min(x_jitter_dr, 1.4), 0.6);  % Constrain within violin
        
        scatter(x_jitter_dr, drug_all, 30, [0.3 0.3 0.3], 'o', ...
            'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
    end
    
    % Add legend for mean/median lines
    legend('Median', 'Mean', 'Location', 'best', 'FontSize', 11);
    
    set(gca, 'XTick', [0 1], 'XTickLabel', {'Baseline', 'Drug'}, 'FontSize', 13);
    xlabel('Condition', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('dF/F', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Violin Plot with Individual Data - %s (n=%d neurons)', stim_type, n_rois_selected_plane), ...
        'FontSize', 14, 'FontWeight', 'bold');
    xlim([-0.7 1.7]);
    grid on; set(gca, 'LineWidth', 2, 'FontSize', 11);
    hold off;
end

%% ====================== HISTOGRAMS & VIOLIN PLOTS (MATCHED NEURONS ONLY) ======================

if matched_rois_available && n_matched > 0
    fprintf('\n========== GENERATING ACTIVITY DISTRIBUTIONS (MATCHED NEURONS) ==========\n');
    
    for field_idx = 1:length(stim_fields)
        field_name = stim_fields{field_idx};
        data = master_data.(field_name);
        stim_type = data.stimulus_type;
        
        if ~isfield(data, 'matched_baseline_avg') || isempty(data.matched_baseline_avg)
            fprintf('Skipping "%s": no matched response data\n', stim_type);
            continue;
        end
        
        fprintf('Distribution plots for "%s": matched %d neuron pairs\n', stim_type, n_matched);
        
        % Extract matched neuron activity
        baseline_matched_all = [];
        for i = 1:length(data.baseline_responses)
            resp = data.baseline_responses{i};
            matched_resp = resp(data.base_match_idx_local, :);
            baseline_matched_all = [baseline_matched_all; matched_resp(:)];
        end
        
        drug_matched_all = [];
        for i = 1:length(data.drug_responses)
            resp = data.drug_responses{i};
            matched_resp = resp(data.drug_match_idx_local, :);
            drug_matched_all = [drug_matched_all; matched_resp(:)];
        end
        
        % Validate data - remove NaN and Inf
        baseline_matched_all = baseline_matched_all(~isnan(baseline_matched_all) & ~isinf(baseline_matched_all));
        drug_matched_all = drug_matched_all(~isnan(drug_matched_all) & ~isinf(drug_matched_all));
        
        % Skip if insufficient data
        if length(baseline_matched_all) < 2 || length(drug_matched_all) < 2
            fprintf('  WARNING: Insufficient valid matched data for "%s" (baseline: %d, drug: %d samples)\n', ...
                stim_type, length(baseline_matched_all), length(drug_matched_all));
            continue;
        end
        
        % ===== Create SEPARATE FIGURES with enhanced visualization =====
        
        % ===== FIGURE 1: HISTOGRAM (Large, Full-Width) =====
        fig_hist = figure('Position', [100 150 1400 600], 'NumberTitle', 'off', ...
            'Name', sprintf('Histogram_Matched_%s', stim_type));
        
        combined_matched_data = [baseline_matched_all; drug_matched_all];
        p1_m = prctile(combined_matched_data, 1);
        p99_m = prctile(combined_matched_data, 99);
        
        % Ensure bin edges have width
        if p1_m == p99_m
            bin_edges_matched = linspace(p1_m - 0.5, p1_m + 0.5, hist_bins + 1);
        else
            bin_edges_matched = linspace(p1_m, p99_m, hist_bins + 1);
        end
        
        hold on;
        histogram(baseline_matched_all, 'Normalization', 'pdf', 'FaceColor', 'r', 'FaceAlpha', 0.5, ...
            'EdgeColor', 'red', 'BinEdges', bin_edges_matched, 'LineWidth', 2, 'DisplayName', 'Baseline');
        histogram(drug_matched_all, 'Normalization', 'pdf', 'FaceColor', 'b', 'FaceAlpha', 0.5, ...
            'EdgeColor', 'blue', 'BinEdges', bin_edges_matched, 'LineWidth', 2, 'DisplayName', 'Drug');
        
        % Add vertical lines for medians
        baseline_median = median(baseline_matched_all);
        drug_median = median(drug_matched_all);
        plot([baseline_median baseline_median], ylim, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Baseline Median');
        plot([drug_median drug_median], ylim, 'b--', 'LineWidth', 2.5, 'DisplayName', 'Drug Median');
        
        ylabel('Probability Density', 'FontSize', 13, 'FontWeight', 'bold');
        xlabel('dF/F', 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('Distribution of Activity - %s (n=%d matched pairs)', stim_type, n_matched), ...
            'FontSize', 14, 'FontWeight', 'bold');
        legend('FontSize', 11, 'Location', 'best');
        grid on; set(gca, 'LineWidth', 2, 'FontSize', 11);
        hold off;
        
        % ===== FIGURE 2: VIOLIN PLOT (Large, Full-Width, Enhanced) =====
        fig_violin = figure('Position', [100 150 1400 600], 'NumberTitle', 'off', ...
            'Name', sprintf('Violin_Matched_%s', stim_type));
        
        hold on;
        
        % BASELINE VIOLIN (x=0) - Solid filled with prominent outline
        if length(baseline_matched_all) >= 2
            [f_bl, xi_bl] = ksdensity(baseline_matched_all, 'NumPoints', 150);
            f_bl = f_bl / max(f_bl) * 0.4;  % Violin half-width (0.8 total)
            
            % Draw solid violin fill with DARK, THICK outline for clarity
            x_violin_bl = [f_bl, fliplr(-f_bl)];
            y_violin_bl = [xi_bl, fliplr(xi_bl)];
            patch(x_violin_bl, y_violin_bl, [0.2 0.6 1], 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 2.5);
            
            % Calculate mean and median
            baseline_median_val = median(baseline_matched_all);
            baseline_mean_val = mean(baseline_matched_all);
            
            % Add median line (thicker, darker)
            plot([-0.4 0.4], [baseline_median_val baseline_median_val], 'k-', 'LineWidth', 3, 'DisplayName', 'Median');
            
            % Add mean line (thinner, different color for distinction)
            plot([-0.4 0.4], [baseline_mean_val baseline_mean_val], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean');
            
            % Add scattered points across violin
            n_baseline = length(baseline_matched_all);
            x_jitter_bl = randn(n_baseline, 1) * 0.12;  % Wider horizontal spread
            x_jitter_bl = max(min(x_jitter_bl, 0.4), -0.4);  % Constrain within violin
            
            scatter(x_jitter_bl, baseline_matched_all, 30, [0.3 0.3 0.3], 'o', ...
                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
        end
        
        % DRUG VIOLIN (x=1) - Solid filled with prominent outline
        if length(drug_matched_all) >= 2
            [f_dr, xi_dr] = ksdensity(drug_matched_all, 'NumPoints', 150);
            f_dr = f_dr / max(f_dr) * 0.4;  % Violin half-width
            
            % Draw solid violin fill with DARK, THICK outline for clarity
            x_violin_dr = [f_dr + 1, fliplr(-f_dr + 1)];
            y_violin_dr = [xi_dr, fliplr(xi_dr)];
            patch(x_violin_dr, y_violin_dr, [1 0.5 0.2], 'FaceAlpha', 0.8, 'EdgeColor', 'black', 'LineWidth', 2.5);
            
            % Calculate mean and median
            drug_median_val = median(drug_matched_all);
            drug_mean_val = mean(drug_matched_all);
            
            % Add median line (thicker, darker)
            plot([1 - 0.4 1 + 0.4], [drug_median_val drug_median_val], 'k-', 'LineWidth', 3, 'DisplayName', 'Median');
            
            % Add mean line (thinner, different color for distinction)
            plot([1 - 0.4 1 + 0.4], [drug_mean_val drug_mean_val], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Mean');
            
            % Add scattered points across violin
            n_drug = length(drug_matched_all);
            x_jitter_dr = 1 + randn(n_drug, 1) * 0.12;  % Wider horizontal spread
            x_jitter_dr = max(min(x_jitter_dr, 1.4), 0.6);  % Constrain within violin
            
            scatter(x_jitter_dr, drug_matched_all, 30, [0.3 0.3 0.3], 'o', ...
                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
        end
        
        % Add legend for mean/median lines
        legend('Median', 'Mean', 'Location', 'best', 'FontSize', 11);
        
        set(gca, 'XTick', [0 1], 'XTickLabel', {'Baseline', 'Drug'}, 'FontSize', 13);
        xlabel('Condition', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('dF/F', 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('Violin Plot with Individual Data - %s (n=%d matched pairs)', stim_type, n_matched), ...
            'FontSize', 14, 'FontWeight', 'bold');
        xlim([-0.7 1.7]);
        grid on; set(gca, 'LineWidth', 2, 'FontSize', 11);
        hold off;
    end
end

%% ====================== ALL STIMULI COMPARISON (MATCHED NEURONS ONLY) ======================

if matched_rois_available && n_matched > 0
    fprintf('\n========== GENERATING MULTI-STIMULUS COMPARISON (MATCHED NEURONS) ==========\n');
    
    % Create one large figure with all stimuli
    fig_all_stim = figure('Position', [100 100 1800 800], 'NumberTitle', 'off', ...
        'Name', 'AllStimuli_Matched_Comparison');
    
    n_stim = length(stim_fields);
    x_pos = 0;  % Running x-position counter for all violins
    all_y_limits = [];  % Collect all y-limits for unified scaling
    stim_labels = {};
    stim_x_positions = {};  % Store x-positions for each stimulus pair
    
    % PASS 1: Collect data and y-limits
    for field_idx = 1:length(stim_fields)
        field_name = stim_fields{field_idx};
        data = master_data.(field_name);
        stim_type = data.stimulus_type;
        
        if ~isfield(data, 'baseline_responses') || isempty(data.baseline_responses)
            continue;
        end
        
        % Extract matched neuron activity (same as before)
        baseline_matched_all = [];
        for i = 1:length(data.baseline_responses)
            resp = data.baseline_responses{i};
            matched_resp = resp(data.base_match_idx_local, :);
            baseline_matched_all = [baseline_matched_all; matched_resp(:)];
        end
        
        drug_matched_all = [];
        for i = 1:length(data.drug_responses)
            resp = data.drug_responses{i};
            matched_resp = resp(data.drug_match_idx_local, :);
            drug_matched_all = [drug_matched_all; matched_resp(:)];
        end
        
        % Validate — keep only positive values to reveal active-response spread
        baseline_matched_all = baseline_matched_all(~isnan(baseline_matched_all) & ~isinf(baseline_matched_all) & baseline_matched_all > 0);
        drug_matched_all = drug_matched_all(~isnan(drug_matched_all) & ~isinf(drug_matched_all) & drug_matched_all > 0);

        if length(baseline_matched_all) < 2 || length(drug_matched_all) < 2
            continue;
        end

        % Collect y-limits
        all_y_limits = [all_y_limits; min(baseline_matched_all); max(baseline_matched_all); ...
                        min(drug_matched_all); max(drug_matched_all)];
        
        stim_labels{end+1} = stim_type;
        stim_x_positions{end+1} = [x_pos, x_pos+1];
        x_pos = x_pos + 2.5;  % Space between stimulus pairs
    end
    
    % Unified y-limits with 5% padding
    y_min = min(all_y_limits);
    y_max = max(all_y_limits);
    y_range = y_max - y_min;
    y_min = y_min - 0.05 * y_range;
    y_max = y_max + 0.05 * y_range;
    
    % PASS 2: Plot all violins on single axes
    hold on;
    x_pos = 0;
    stim_idx_plot = 1;
    
    for field_idx = 1:length(stim_fields)
        field_name = stim_fields{field_idx};
        data = master_data.(field_name);
        stim_type = data.stimulus_type;
        
        if ~isfield(data, 'baseline_responses') || isempty(data.baseline_responses)
            continue;
        end
        
        % Extract matched neuron activity
        baseline_matched_all = [];
        for i = 1:length(data.baseline_responses)
            resp = data.baseline_responses{i};
            matched_resp = resp(data.base_match_idx_local, :);
            baseline_matched_all = [baseline_matched_all; matched_resp(:)];
        end
        
        drug_matched_all = [];
        for i = 1:length(data.drug_responses)
            resp = data.drug_responses{i};
            matched_resp = resp(data.drug_match_idx_local, :);
            drug_matched_all = [drug_matched_all; matched_resp(:)];
        end
        
        % Validate — keep only positive values to reveal active-response spread
        baseline_matched_all = baseline_matched_all(~isnan(baseline_matched_all) & ~isinf(baseline_matched_all) & baseline_matched_all > 7);
        drug_matched_all = drug_matched_all(~isnan(drug_matched_all) & ~isinf(drug_matched_all) & drug_matched_all > 7);

        if length(baseline_matched_all) < 2 || length(drug_matched_all) < 2
            continue;
        end

        % BASELINE VIOLIN at x_pos
        if length(baseline_matched_all) >= 2
            [f_bl, xi_bl] = ksdensity(baseline_matched_all, 'NumPoints', 150);
            f_bl = f_bl / max(f_bl) * 0.4;
            
            % Draw solid violin with light outline
            x_violin_bl = [f_bl + x_pos, fliplr(-f_bl + x_pos)];
            y_violin_bl = [xi_bl, fliplr(xi_bl)];
            patch(x_violin_bl, y_violin_bl, [0.2 0.6 1], 'FaceAlpha', 0.8, 'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1);
            
            % Median line
            baseline_median = median(baseline_matched_all);
            plot([x_pos - 0.4 x_pos + 0.4], [baseline_median baseline_median], 'k-', 'LineWidth', 2.5);
            
            % Scattered points
            n_baseline = length(baseline_matched_all);
            x_jitter_bl = x_pos + randn(n_baseline, 1) * 0.12;
            x_jitter_bl = max(min(x_jitter_bl, x_pos + 0.4), x_pos - 0.4);
            scatter(x_jitter_bl, baseline_matched_all, 25, [0.3 0.3 0.3], 'o', ...
                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
        end
        
        % DRUG VIOLIN at x_pos+1
        if length(drug_matched_all) >= 2
            [f_dr, xi_dr] = ksdensity(drug_matched_all, 'NumPoints', 150);
            f_dr = f_dr / max(f_dr) * 0.4;
            
            % Draw solid violin with light outline
            x_violin_dr = [f_dr + x_pos + 1, fliplr(-f_dr + x_pos + 1)];
            y_violin_dr = [xi_dr, fliplr(xi_dr)];
            patch(x_violin_dr, y_violin_dr, [1 0.5 0.2], 'FaceAlpha', 0.8, 'EdgeColor', [0.4 0.4 0.4], 'LineWidth', 1);
            
            % Median line
            drug_median = median(drug_matched_all);
            plot([x_pos + 1 - 0.4 x_pos + 1 + 0.4], [drug_median drug_median], 'k-', 'LineWidth', 2.5);
            
            % Scattered points
            n_drug = length(drug_matched_all);
            x_jitter_dr = x_pos + 1 + randn(n_drug, 1) * 0.12;
            x_jitter_dr = max(min(x_jitter_dr, x_pos + 1 + 0.4), x_pos + 1 - 0.4);
            scatter(x_jitter_dr, drug_matched_all, 25, [0.3 0.3 0.3], 'o', ...
                'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2);
        end
        
        x_pos = x_pos + 2.5;
        stim_idx_plot = stim_idx_plot + 1;
    end
    
    % Set axis properties
    set(gca, 'YLim', [-4 y_max]);
    ylabel('dF/F', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Stimulus Type', 'FontSize', 14, 'FontWeight', 'bold');
    title(sprintf('All Stimuli Comparison - Baseline vs Drug (n=%d matched pairs)', n_matched), ...
        'FontSize', 15, 'FontWeight', 'bold');
    
    % Create x-axis labels at stimulus pairs
    x_ticks = [];
    x_labels = {};
    x_pos = 0.5;  % Center between baseline and drug
    for s = 1:length(stim_labels)
        x_ticks = [x_ticks, x_pos];
        x_labels{s} = stim_labels{s};
        x_pos = x_pos + 2.5;
    end
    set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels, 'FontSize', 12);
    
    grid on; set(gca, 'LineWidth', 2, 'FontSize', 11);
    hold off;
    
    fprintf('Generated multi-stimulus comparison figure with %d stimuli\n', length(stim_labels));

    %% ==================== FIGURE: PAIRED DOT PLOT (per-neuron modulation) ====================
    fprintf('========== GENERATING PAIRED DOT PLOTS ==========\n');

    % Collect stimuli that have response data
    valid_stim_fields = {};
    for fi = 1:length(stim_fields)
        d = master_data.(stim_fields{fi});
        if isfield(d, 'baseline_responses') && ~isempty(d.baseline_responses)
            valid_stim_fields{end+1} = stim_fields{fi}; %#ok<AGROW>
        end
    end
    n_valid = length(valid_stim_fields);

    fig_paired = figure('Position', [100 50 min(420*n_valid, 1800) 540], ...
        'NumberTitle', 'off', 'Name', 'PerNeuron_PairedModulation', 'Color', 'w');

    for sp_idx = 1:n_valid
        field_name = valid_stim_fields{sp_idx};
        d = master_data.(field_name);
        stim_type = d.stimulus_type;

        % Per-neuron mean ΔF/F across ALL presentations and timepoints
        all_base_fr = [];
        for i = 1:length(d.baseline_responses)
            all_base_fr = [all_base_fr, d.baseline_responses{i}(d.base_match_idx_local, :)]; %#ok<AGROW>
        end
        base_mean = mean(all_base_fr, 2);   % n_matched × 1

        all_drug_fr = [];
        for i = 1:length(d.drug_responses)
            all_drug_fr = [all_drug_fr, d.drug_responses{i}(d.drug_match_idx_local, :)]; %#ok<AGROW>
        end
        drug_mean = mean(all_drug_fr, 2);   % n_matched × 1

        % Remove NaN/Inf
        ok = ~isnan(base_mean) & ~isinf(base_mean) & ~isnan(drug_mean) & ~isinf(drug_mean);
        base_mean = base_mean(ok);
        drug_mean = drug_mean(ok);

        n_n      = numel(base_mean);
        changes  = drug_mean - base_mean;
        pct_up   = 100 * mean(changes > 0);

        % ---- Per-neuron diverging colour (blue=suppressed, red=enhanced) ----
        max_abs = max(prctile(abs(changes), 95), eps);
        t       = max(-1, min(1, changes / max_abs));   % n×1, clamped [-1,1]
        c_blue  = [0.15 0.40 0.75];
        c_red   = [0.80 0.15 0.10];
        c_white = [1 1 1];
        neuron_colors = repmat(c_white, n_n, 1);
        neg = t < 0;  pos = ~neg;
        neuron_colors(neg,:) = c_white + abs(t(neg)) * (c_blue - c_white);  % (k×1)*(1×3)
        neuron_colors(pos,:) = c_white + t(pos)      * (c_red  - c_white);

        % ---- Reproducible jitter so each neuron's two dots align ----
        rng(sp_idx * 7);  % deterministic per stimulus
        dx = (rand(n_n, 1) - 0.5) * 0.35;

        subplot(1, n_valid, sp_idx);
        hold on;

        % Layer 1: all connecting lines in one call (grey, very transparent)
        x_l = reshape([(0+dx)'; (1+dx)'; nan(1,n_n)], [], 1);
        y_l = reshape([base_mean'; drug_mean'; nan(1,n_n)], [], 1);
        plot(x_l, y_l, '-', 'Color', [0.55 0.55 0.55 0.07], 'LineWidth', 0.5);

        % Layer 2: coloured dots at each endpoint (colour encodes direction/magnitude)
        scatter(0+dx, base_mean, 22, neuron_colors, 'o', 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0);
        scatter(1+dx, drug_mean, 22, neuron_colors, 'o', 'filled', ...
            'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0);

        % Layer 3: mean ± SEM summary — thick connecting line + capped error bars
        bm = mean(base_mean);  bs = std(base_mean) / sqrt(n_n);
        dm = mean(drug_mean);  ds = std(drug_mean)  / sqrt(n_n);
        plot([0 1], [bm dm], 'k-', 'LineWidth', 2.5);
        errorbar(0, bm, bs, 'k-', 'LineWidth', 2, 'CapSize', 7);
        errorbar(1, dm, ds, 'k-', 'LineWidth', 2, 'CapSize', 7);
        plot(0, bm, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'LineWidth', 2);
        plot(1, dm, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'w', 'LineWidth', 2);

        hold off;
        set(gca, 'XTick', [0 1], 'XTickLabel', {'Baseline','Drug'}, 'FontSize', 10);
        xlim([-0.5 1.5]);
        ylabel('Mean \DeltaF/F', 'FontSize', 10);
        title(sprintf('%s\n%.0f%% \\uparrow  |  %.0f%% \\downarrow', ...
            strrep(stim_type,'_','\_'), pct_up, 100-pct_up), ...
            'FontSize', 10, 'FontWeight', 'bold');
        grid on; box off;
        set(gca, 'LineWidth', 1.2, 'FontSize', 10);
    end

    sgtitle(sprintf('Per-Neuron Drug Modulation — Paired Mean \\DeltaF/F  (n=%d matched neurons)', n_matched), ...
        'FontSize', 13, 'FontWeight', 'bold');

    %% ==================== FIGURE: MODULATION INDEX CDF ====================
    fprintf('========== GENERATING MODULATION INDEX CDF ==========\n');

    fig_cdf = figure('Position', [100 50 min(420*n_valid, 1800) 540], ...
        'NumberTitle', 'off', 'Name', 'ModulationIndex_CDF', 'Color', 'w');

    for sp_idx = 1:n_valid
        field_name = valid_stim_fields{sp_idx};
        d = master_data.(field_name);
        stim_type = d.stimulus_type;

        % Re-use per-neuron means (same extraction as above)
        all_base_fr = [];
        for i = 1:length(d.baseline_responses)
            all_base_fr = [all_base_fr, d.baseline_responses{i}(d.base_match_idx_local, :)]; %#ok<AGROW>
        end
        base_mean = mean(all_base_fr, 2);

        all_drug_fr = [];
        for i = 1:length(d.drug_responses)
            all_drug_fr = [all_drug_fr, d.drug_responses{i}(d.drug_match_idx_local, :)]; %#ok<AGROW>
        end
        drug_mean = mean(all_drug_fr, 2);

        ok = ~isnan(base_mean) & ~isinf(base_mean) & ~isnan(drug_mean) & ~isinf(drug_mean);
        base_mean = base_mean(ok);
        drug_mean = drug_mean(ok);

        % Symmetric modulation index: (drug - base) / (|drug| + |base|)
        % Bounded [-1, 1]: positive = enhanced, negative = suppressed, 0 = no change
        denom  = abs(drug_mean) + abs(base_mean);
        ok_mi  = denom > 0;
        mi     = (drug_mean(ok_mi) - base_mean(ok_mi)) ./ denom(ok_mi);

        pct_enhanced   = 100 * mean(mi > 0);
        pct_suppressed = 100 * mean(mi < 0);

        mi_sorted = sort(mi);
        cdf_y     = (1:numel(mi_sorted))' / numel(mi_sorted);

        subplot(1, n_valid, sp_idx);
        hold on;

        % Shaded background: blue = suppressed side, red = enhanced side
        fill([-1 0 0 -1], [0 0 1 1], [0.20 0.45 0.80], 'FaceAlpha', 0.07, 'EdgeColor', 'none');
        fill([0 1 1 0],   [0 0 1 1], [0.85 0.25 0.15], 'FaceAlpha', 0.07, 'EdgeColor', 'none');

        % CDF
        plot(mi_sorted, cdf_y, 'k-', 'LineWidth', 2.5);

        % Reference lines
        xline(0,   'r--', 'LineWidth', 1.8);
        yline(0.5, '--',  'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);

        hold off;
        xlabel('Modulation Index', 'FontSize', 10);
        ylabel('Cumulative fraction', 'FontSize', 10);
        xlim([-1 1]);  ylim([0 1]);
        title(sprintf('%s\n%.0f%% \\uparrow  |  %.0f%% \\downarrow', ...
            strrep(stim_type,'_','\_'), pct_enhanced, pct_suppressed), ...
            'FontSize', 10, 'FontWeight', 'bold');
        grid on; box off;
        set(gca, 'LineWidth', 1.2, 'FontSize', 10);
    end

    sgtitle(sprintf('Modulation Index CDF — (drug \x2212 baseline) / (|drug| + |baseline|)  (n=%d neurons)', n_matched), ...
        'FontSize', 13, 'FontWeight', 'bold');

    %% ==================== PAIRED STATISTICAL TESTS ====================
    fprintf('\n========== PAIRED STATISTICS (Baseline vs Drug - Matched ROIs) ==========\n');
    fprintf('Comparing same neurons across two conditions per stimulus\n\n');
    
    for field_idx = 1:length(stim_fields)
        field_name = stim_fields{field_idx};
        data = master_data.(field_name);
        stim_type = data.stimulus_type;
        
        if ~isfield(data, 'baseline_responses') || isempty(data.baseline_responses)
            continue;
        end
        
        % Extract baseline paired data (per neuron)
        baseline_per_neuron = [];
        for i = 1:length(data.baseline_responses)
            resp = data.baseline_responses{i};
            matched_resp = resp(data.base_match_idx_local, :);
            baseline_per_neuron = [baseline_per_neuron; mean(matched_resp, 2)];
        end
        
        % Extract drug paired data (per neuron)
        drug_per_neuron = [];
        for i = 1:length(data.drug_responses)
            resp = data.drug_responses{i};
            matched_resp = resp(data.drug_match_idx_local, :);
            drug_per_neuron = [drug_per_neuron; mean(matched_resp, 2)];
        end
        
        % Remove any NaN/Inf
        valid_idx = ~(isnan(baseline_per_neuron) | isinf(baseline_per_neuron) | ...
                      isnan(drug_per_neuron) | isinf(drug_per_neuron));
        baseline_per_neuron = baseline_per_neuron(valid_idx);
        drug_per_neuron = drug_per_neuron(valid_idx);
        
        if length(baseline_per_neuron) < 2
            continue;
        end
        
        % Calculate descriptive statistics
        baseline_mean = mean(baseline_per_neuron);
        baseline_median = median(baseline_per_neuron);
        baseline_std = std(baseline_per_neuron);
        baseline_sem = baseline_std / sqrt(length(baseline_per_neuron));
        
        drug_mean = mean(drug_per_neuron);
        drug_median = median(drug_per_neuron);
        drug_std = std(drug_per_neuron);
        drug_sem = drug_std / sqrt(length(drug_per_neuron));
        
        % Calculate percent changes
        pct_change_mean = ((drug_mean - baseline_mean) / abs(baseline_mean)) * 100;
        pct_change_median = ((drug_median - baseline_median) / abs(baseline_median)) * 100;
        
        % Paired t-test (assumes normality)
        [h_ttest, p_ttest, ci_ttest, stats_ttest] = ttest(baseline_per_neuron, drug_per_neuron);
        
        % Wilcoxon signed-rank test (non-parametric, more robust)
        [p_wilcoxon, h_wilcoxon, stats_wilcoxon] = signrank(baseline_per_neuron, drug_per_neuron);
        
        % Effect size: Cohen's d for paired samples
        diff = baseline_per_neuron - drug_per_neuron;
        cohen_d = mean(diff) / std(diff);
        
        % Print results
        fprintf('-------------------------------------------\n');
        fprintf('STIMULUS: %s (n=%d neurons)\n', stim_type, length(baseline_per_neuron));
        fprintf('-------------------------------------------\n');
        
        fprintf('BASELINE:     Mean = %7.4f (SEM = %.4f), Median = %7.4f, SD = %.4f\n', ...
            baseline_mean, baseline_sem, baseline_median, baseline_std);
        fprintf('DRUG:         Mean = %7.4f (SEM = %.4f), Median = %7.4f, SD = %.4f\n', ...
            drug_mean, drug_sem, drug_median, drug_std);
        
        fprintf('\nCHANGE:       Mean = %+.4f (%+.1f%%), Median = %+.4f (%+.1f%%)\n', ...
            drug_mean - baseline_mean, pct_change_mean, ...
            drug_median - baseline_median, pct_change_median);
        
        fprintf('\nSTATISTICAL TESTS:\n');
        fprintf('  Paired t-test:           t(%d) = %7.3f, p = %.4f %s\n', ...
            stats_ttest.df, stats_ttest.tstat, p_ttest, significance_marker(p_ttest));
        fprintf('  Cohen''s d (effect size): %.4f\n', cohen_d);
        fprintf('  Wilcoxon signed-rank:    p = %.4f %s\n', ...
            p_wilcoxon, significance_marker(p_wilcoxon));
        fprintf('\n');
    end
    
    %% ==================== MODULATION PLOTS (Linear Regression) ====================
    fprintf('========== GENERATING BASELINE vs DRUG MODULATION PLOTS ==========\n');
    
    % Create one large figure with regression plots for all stimuli
    fig_modulation = figure('Position', [100 100 1400 900], 'NumberTitle', 'off', ...
        'Name', 'Baseline_vs_Drug_Modulation');
    
    n_stim_plots = 0;
    for field_idx = 1:length(stim_fields)
        field_name = stim_fields{field_idx};
        data = master_data.(field_name);
        stim_type = data.stimulus_type;
        
        if ~isfield(data, 'baseline_responses') || isempty(data.baseline_responses)
            continue;
        end
        
        % Extract baseline paired data (per neuron)
        baseline_per_neuron = [];
        for i = 1:length(data.baseline_responses)
            resp = data.baseline_responses{i};
            matched_resp = resp(data.base_match_idx_local, :);
            baseline_per_neuron = [baseline_per_neuron; mean(matched_resp, 2)];
        end
        
        % Extract drug paired data (per neuron)
        drug_per_neuron = [];
        for i = 1:length(data.drug_responses)
            resp = data.drug_responses{i};
            matched_resp = resp(data.drug_match_idx_local, :);
            drug_per_neuron = [drug_per_neuron; mean(matched_resp, 2)];
        end
        
        % Remove any NaN/Inf
        valid_idx = ~(isnan(baseline_per_neuron) | isinf(baseline_per_neuron) | ...
                      isnan(drug_per_neuron) | isinf(drug_per_neuron));
        baseline_per_neuron = baseline_per_neuron(valid_idx);
        drug_per_neuron = drug_per_neuron(valid_idx);
        
        if length(baseline_per_neuron) < 2
            continue;
        end
        
        n_stim_plots = n_stim_plots + 1;
        
        % Create subplot
        subplot(2, 2, n_stim_plots);
        
        % Scatter plot
        scatter(baseline_per_neuron, drug_per_neuron, 100, 'o', 'filled', ...
            'MarkerFaceColor', [0.3 0.5 0.8], 'MarkerEdgeColor', 'black', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.8, 'LineWidth', 1.5);
        hold on;
        
        % Linear regression line
        coeffs = polyfit(baseline_per_neuron, drug_per_neuron, 1);
        slope = coeffs(1);
        intercept = coeffs(2);
        
        % Calculate R-squared
        y_fit = polyval(coeffs, baseline_per_neuron);
        ss_res = sum((drug_per_neuron - y_fit) .^ 2);
        ss_tot = sum((drug_per_neuron - mean(drug_per_neuron)) .^ 2);
        r_squared = 1 - (ss_res / ss_tot);
        
        % Plot regression line
        x_range = [min(baseline_per_neuron) - 0.5, max(baseline_per_neuron) + 0.5];
        y_reg = slope * x_range + intercept;
        plot(x_range, y_reg, 'b-', 'LineWidth', 2.5, 'DisplayName', ...
            sprintf('y = %.3f*x + %.3f (R² = %.3f)', slope, intercept, r_squared));
        
        % Identity line (no change: drug = baseline)
        y_min = min([baseline_per_neuron; drug_per_neuron]) - 0.5;
        y_max = max([baseline_per_neuron; drug_per_neuron]) + 0.5;
        plot([y_min y_max], [y_min y_max], 'k--', 'LineWidth', 1.5, ...
            'DisplayName', 'No change (y=x)');
        
        % Labels and formatting
        xlabel('Baseline dF/F', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Drug dF/F', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('%s (n=%d neurons)', stim_type, length(baseline_per_neuron)), ...
            'FontSize', 13, 'FontWeight', 'bold');
           
        grid on; set(gca, 'LineWidth', 1.5, 'FontSize', 10);
        axis equal; 
        xlim([y_min y_max]); 
        ylim([y_min y_max]);
        hold off; 
    end
    
    % Overall title
    sgtitle(sprintf('Baseline vs Drug Modulation - Matched Neurons (n=%d pairs)', n_matched), ...
        'FontSize', 15, 'FontWeight', 'bold');
    
end



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

function marker = significance_marker(p_value)
    % Return significance marker based on p-value
    if p_value < 0.001
        marker = '***';
    elseif p_value < 0.01
        marker = '**';
    elseif p_value < 0.05
        marker = '*';
    else
        marker = 'ns';
    end
end

function [responses, frame_ranges, time_ranges, properties] = extract_full_stimulus_responses(Stimuli, dff, stim_type, TimeCa, roi_indices)
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
    %   dff: calcium data (n_neurons × n_frames, full dataset)
    %   stim_type: stimulus type to extract (string)
    %   TimeCa: time matrix from preprocessed file (2×n_frames, row 1 = time)
    %   roi_indices: (optional) row indices to subset dff to specific plane (default: use all ROIs)
    %
    % Returns:
    %   responses: cell array, each element is (n_neurons × stim_duration_frames)
    %   frame_ranges: (n_presentations × 2) [frame_start, frame_end]
    %   time_ranges: (n_presentations × 2) [time_start, time_end]
    %   properties: struct array with stimulus metadata
    
    % If roi_indices provided, subset to those ROIs only
    if nargin >= 5 && ~isempty(roi_indices)
        dff = dff(roi_indices, :);
    end
    
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

function analyze_grating_stimulus(B_Stimuli, ~, ~, ~, ~, ~, ~, ~)
    % Analyze grating stimulus: find preferred directions, average responses
    % [SCAFFOLD - Parameters available for future implementation]
    
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

function analyze_moving_bar_stimulus(~, ~, ~, ~, ...
    ~, ~, ~, ~)
    % Analyze moving bar stimulus: extract stimulus-triggered responses at onset
    % [SCAFFOLD - Parameters available for future implementation]
    
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

function analyze_full_field_flash_stimulus(~, ~, ~, ~, ...
    ~, ~, ~, ~)
    % Analyze full-field flash: plot individual trial responses + averages
    % [SCAFFOLD - Parameters available for future implementation]
    
    fprintf('    - Analyzing full-field flash responses\n');
    
    % Create stimulus-triggered average plot
    % figure('Position', [100 100 1200 700], 'NumberTitle', 'off', ...
    %     'Name', 'Full-Field Flash Analysis');
    
    % Subplot 1: Baseline individual trials + mean
    subplot(2, 2, 1);
    baseline_responses = data.baseline_responses;
    colors_base = repmat([1 0.4 0.4], length(baseline_responses), 1);
    
    % max_len = max(cellfun(@(x) size(x, 2), baseline_responses));
    
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

function analyze_looming_stimulus(base_Stimuli, base_dff, base_TimeCa, base_Triggers, ...
        drug_Stimuli, drug_dff, drug_TimeCa, drug_Triggers, ...
        base_match_idx_local, drug_match_idx_local, matched_rois_available)
    % Looming-triggered response, matched cells, baseline vs drug.
    % Looming fires exactly one red sync frame per trial, at the true
    % onset of that trial's loom (unlike the other stimulus types, which
    % have a continuous periodic red-sync rhythm). Trials are NOT evenly
    % spaced -- each has a randomized wait period beforehand -- so onsets
    % must come directly from these red frames, not from
    % stimulus_trial_t x trial_index.

    fprintf('    - Analyzing looming stimulus (matched baseline vs drug)\n');

    if ~matched_rois_available
        fprintf('      (No matched ROIs available -- skipping looming baseline/drug comparison)\n');
        return;
    end

    [Ca_base, t_com_base] = extract_looming_trials(base_Stimuli, base_dff, base_TimeCa, base_Triggers);
    [Ca_drug, t_com_drug] = extract_looming_trials(drug_Stimuli, drug_dff, drug_TimeCa, drug_Triggers);

    if isempty(Ca_base) || isempty(Ca_drug)
        fprintf('      (Looming stimulus missing in baseline or drug -- skipping)\n');
        return;
    end
    if ~isequal(size(t_com_base), size(t_com_drug)) || max(abs(t_com_base - t_com_drug)) > 1e-6
        warning('Baseline and drug looming analysis windows differ (stimulus_trial_t mismatch?) -- using baseline time axis');
    end
    t_com = t_com_base;

    % Restrict to matched cells. Use the MEDIAN across each condition's own
    % trials, not the mean -- with only a handful of trials per cell, a
    % single large-but-unreliable trial can otherwise dominate the mean
    % and flip a cell's responsive/non-responsive classification. The
    % median requires a majority of trials to be elevated instead.
    resp_base = median(Ca_base(base_match_idx_local, :, :), 3, 'omitnan');   % [n_matched x n_com]
    resp_drug = median(Ca_drug(drug_match_idx_local, :, :), 3, 'omitnan');   % [n_matched x n_com]
    n_matched = size(resp_base, 1);

    %% Plot 1: baseline vs drug, all matched cells
    plot_looming_comparison(t_com, resp_base, resp_drug, ...
        sprintf('Looming-triggered response: baseline vs drug (n = %d matched cells)', n_matched));

    %% Responsive-cell classification: mean dF/F > 0.2 from 1s post-onset to end of window
    post1s         = t_com >= 1;
    post_mean_base = mean(resp_base(:, post1s), 2);
    post_mean_drug = mean(resp_drug(:, post1s), 2);
    responsive_base = post_mean_base > 0.2;
    responsive_drug = post_mean_drug > 0.2;
    pct_base         = 100 * sum(responsive_base) / n_matched;
    pct_drug         = 100 * sum(responsive_drug) / n_matched;

    fprintf('      Responsive cells (mean dF/F > 0.2, %.1fs to %.1fs post-onset):\n', t_com(find(post1s, 1)), t_com(end));
    fprintf('        Baseline: %d / %d (%.1f%%)\n', sum(responsive_base), n_matched, pct_base);
    fprintf('        Drug:     %d / %d (%.1f%%)\n', sum(responsive_drug), n_matched, pct_drug);

    %% Plot 2: % responsive cells, baseline vs drug
    figure('Name', 'Looming: % responsive cells', 'NumberTitle', 'off', ...
        'Position', [100 100 500 500]);
    bar([pct_base, pct_drug], 'FaceColor', 'flat', 'CData', [0.2 0.6 1; 1 0.5 0.2]);
    set(gca, 'XTickLabel', {'Baseline', 'Drug'}, 'Box', 'off');
    ylabel('Responsive cells (%)', 'FontSize', 11);
    title('Looming-responsive matched cells', 'FontSize', 12);
    ylim([0, 100]);

    %% Plot 3: baseline vs drug, only cells responsive in BOTH conditions
    responsive_both = responsive_base & responsive_drug;
    n_both = sum(responsive_both);
    fprintf('      Responsive in both conditions: %d / %d (%.1f%%)\n', n_both, n_matched, 100*n_both/n_matched);

    if n_both > 0
        plot_looming_comparison(t_com, resp_base(responsive_both, :), resp_drug(responsive_both, :), ...
            sprintf('Looming-triggered response: baseline vs drug (n = %d responsive-in-both cells)', n_both));
    else
        fprintf('      (No cells responsive in both conditions -- skipping Plot 3)\n');
    end

    %% Plot 4: paired per-cell comparison -- diagnostic for mean vs. %-responsive dissociation
    cat_neither   = ~responsive_base & ~responsive_drug;
    cat_base_only =  responsive_base & ~responsive_drug;
    cat_drug_only = ~responsive_base &  responsive_drug;
    cat_both      =  responsive_base &  responsive_drug;

    figure('Name', 'Looming: per-cell baseline vs drug', 'NumberTitle', 'off', ...
        'Position', [100 100 550 550]);
    hold on;
    scatter(post_mean_base(cat_neither),   post_mean_drug(cat_neither),   25, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'Neither');
    scatter(post_mean_base(cat_base_only), post_mean_drug(cat_base_only), 25, [0.2 0.6 1],   'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Baseline only');
    scatter(post_mean_base(cat_drug_only), post_mean_drug(cat_drug_only), 25, [1 0.5 0.2],   'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Drug only');
    scatter(post_mean_base(cat_both),      post_mean_drug(cat_both),      25, [0.3 0.1 0.5], 'filled', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Both');

    ax_lim = [min([post_mean_base; post_mean_drug; 0]) - 0.05, max([post_mean_base; post_mean_drug]) + 0.05];
    plot(ax_lim, ax_lim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xline(0.2, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    yline(0.2, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    xlim(ax_lim); ylim(ax_lim); axis square;
    xlabel('Baseline: mean dF/F (1s-end post-onset)', 'FontSize', 11);
    ylabel('Drug: mean dF/F (1s-end post-onset)', 'FontSize', 11);
    title(sprintf('Per-cell post-onset response (n = %d matched cells)', n_matched), 'FontSize', 12);
    legend('Location', 'best');
    set(gca, 'Box', 'off');
end

function [Ca_trials, t_com] = extract_looming_trials(Stimuli, dff, TimeCa, Triggers)
    % Slices and aligns each looming trial in one recording to a common
    % time grid, using the per-trial red sync frame as the true onset.

    Ca_trials = [];
    t_com     = [];

    loom_idx = find(strcmp({Stimuli.type}, 'looming'), 1);
    if isempty(loom_idx)
        return;
    end

    trial_t    = Stimuli(loom_idx).stimulus_trial_t;
    time_start = Stimuli(loom_idx).TimeStimulusFrame(1);
    time_end   = Stimuli(loom_idx).TimeStimulusFrame(end);

    PT          = Triggers.TimeProjector;
    loom_onsets = PT(PT > time_start & PT < time_end)';
    n_trials    = numel(loom_onsets);

    pre_window_sec  = 2;                        % adjustable
    post_window_sec = trial_t - pre_window_sec; % rest of the trial cycle
    n_com           = 100;
    t_com           = linspace(-pre_window_sec, post_window_sec, n_com);

    n_rois    = size(dff, 1);
    Ca_trials = nan(n_rois, n_com, n_trials);
    time_vec  = TimeCa(1, :);

    for k = 1:n_trials
        idx = find(time_vec > loom_onsets(k) - pre_window_sec & ...
                   time_vec < loom_onsets(k) + post_window_sec);
        if numel(idx) < 2; continue; end
        t_rel = time_vec(idx) - loom_onsets(k);
        Ca_trials(:, :, k) = interp1(t_rel(:), dff(:, idx)', t_com(:), 'linear', 'extrap')';
    end
end

function plot_looming_comparison(t_com, resp_base, resp_drug, fig_title)
    % resp_base/resp_drug: [n_cells x n_com], already trial-averaged per cell

    mean_base = mean(resp_base, 1, 'omitnan');
    sem_base  = std(resp_base, 0, 1, 'omitnan') / sqrt(size(resp_base, 1));
    mean_drug = mean(resp_drug, 1, 'omitnan');
    sem_drug  = std(resp_drug, 0, 1, 'omitnan') / sqrt(size(resp_drug, 1));

    figure('Name', 'Looming: baseline vs drug', 'NumberTitle', 'off', ...
        'Position', [100 100 900 500]);
    hold on;
    x_sh = [t_com, fliplr(t_com)];

    y_sh_base = [(mean_base + sem_base), fliplr(mean_base - sem_base)];
    fill(x_sh, y_sh_base, [0.2 0.6 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_base, '-', 'Color', [0.2 0.6 1], 'LineWidth', 2.5, 'DisplayName', 'Baseline');

    y_sh_drug = [(mean_drug + sem_drug), fliplr(mean_drug - sem_drug)];
    fill(x_sh, y_sh_drug, [1 0.5 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t_com, mean_drug, '-', 'Color', [1 0.5 0.2], 'LineWidth', 2.5, 'DisplayName', 'Drug');

    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time from loom onset (s)', 'FontSize', 11);
    ylabel('dF/F (mean \pm SEM)', 'FontSize', 11);
    title(fig_title, 'FontSize', 12, 'Interpreter', 'none');
    legend('Location', 'best');
    set(gca, 'Box', 'off'); xlim([t_com(1), t_com(end)]);
end
