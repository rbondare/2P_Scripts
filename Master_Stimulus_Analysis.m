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
baseline_file = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1327_preprocessed.mat";
drug_file     = "Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260305_1409_preprocessed.mat";
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_V2.mat";

% Analysis parameters
ca_type = 1;                  % 1=FVDff, 2=deconvolved, 3=F
selected_plane_idx = 1;           % Plane index (1-based)
stimulus_types_to_analyze = {'spontaneous', 'moving_bar', 'sparse_local_global_flashes', 'checkers2'};

% Visualization
max_neurons_display = 500;    % For heatmap readability
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
    fprintf('  ✓ Values match - subsetting is correct!\n');
else
    fprintf('  ✗ ERROR: Values do not match - subsetting logic is broken!\n');
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
    
    % Average responses across presentations
    baseline_full = average_stimulus_presentations(data.baseline_responses);
    drug_full = average_stimulus_presentations(data.drug_responses);
    
    % Create figure: baseline and drug side-by-side
    fig = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
        'Name', sprintf('Level1_FullDFF: %s', stim_type));
    
    % BASELINE - ALL neurons, ALL timepoints
    subplot(1, 2, 1);
    imagesc(baseline_full);
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(L1_caxis_lim);
    colorbar;
    xlabel('Time (frames)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('ROI', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('Baseline - %s (Full dFF)\nAll %d neurons × %d timepoints', ...
        stim_type, size(baseline_full, 1), size(baseline_full, 2)), 'FontWeight', 'bold', 'FontSize', 11);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % DRUG - ALL neurons, ALL timepoints
    subplot(1, 2, 2);
    imagesc(drug_full);
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(L1_caxis_lim);
    c = colorbar;
    ylabel(c, 'dF/F', 'FontSize', 10);
    xlabel('Time (frames)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('ROI', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('Drug - %s (Full dFF)\nAll %d neurons × %d timepoints', ...
        stim_type, size(drug_full, 1), size(drug_full, 2)), 'FontWeight', 'bold', 'FontSize', 11);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle(sprintf('LEVEL 1 - Full dFF (No filtering): %s', stim_type), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % Store full heatmaps
    data.full_baseline_heatmap = baseline_full;
    data.full_drug_heatmap = drug_full;
    master_data.(field_name) = data;
end

%% ====================== LEVEL 2: SELECTED ROI HEATMAPS (RAW VALUES, NO SORTING) ======================

fprintf('\n========== LEVEL 2: Selected ROI Heatmaps (raw values, unsorted) ==========\n');

% ===== Parameters for Level 2 (adjust as needed) =====
L2_n_select = min(500, size(baseline_full, 1));          % Number of ROIs to select
L2_frame_start = 200;                                       % Frame window start
L2_frame_end = size(baseline_full, 2);                   % Frame window end (use all by default)
L2_caxis_lim = [0, 5];                                   % Color axis limits

fprintf('Level 2 settings: %d ROIs, frames %d:%d\n', L2_n_select, L2_frame_start, L2_frame_end);

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
            fprintf('  ✗ ERROR: Match indices exceed data dimensions!\n');
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
        
        fprintf('  ✓ Heatmaps created\n');
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
            
        case 'moving_bar'
            analyze_moving_bar_stimulus(B.Stimuli, D.Stimuli, base_dff, drug_dff, ...
                data, base_match_idx_local, drug_match_idx_local, matched_rois_available);
            
        case 'full_field_flash'
            analyze_full_field_flash_stimulus(B.Stimuli, D.Stimuli, base_dff, drug_dff, ...
                data, base_match_idx_local, drug_match_idx_local, matched_rois_available);
            
        case 'spontaneous'
            % Spontaneous activity: no specific trigger, but can analyze statistics
            fprintf('  (Spontaneous activity - no stimulus-triggered analysis)\n');
            
        otherwise
            fprintf('  (No specific analysis for this stimulus type)\n');
    end
end

%% ====================== HISTOGRAMS & VIOLIN PLOTS (ALL NEURONS) ======================

fprintf('\n========== GENERATING ACTIVITY DISTRIBUTIONS (ALL NEURONS) ==========\n');

hist_bins = 50;  % Number of bins

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
    
    % BASELINE VIOLIN (x=0) - Solid filled with scattered points
    if length(baseline_all) >= 2
        [f_bl, xi_bl] = ksdensity(baseline_all, 'NumPoints', 150);
        f_bl = f_bl / max(f_bl) * 0.4;  % Violin half-width (0.8 total)
        
        % Draw solid violin fill
        x_violin_bl = [f_bl, fliplr(-f_bl)];
        y_violin_bl = [xi_bl, fliplr(xi_bl)];
        patch(x_violin_bl, y_violin_bl, [0.2 0.6 1], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        
        % Add thick median line
        plot([-0.4 0.4], [baseline_median baseline_median], 'k-', 'LineWidth', 3);
        
        % Add scattered points across violin
        n_baseline = length(baseline_all);
        x_jitter_bl = randn(n_baseline, 1) * 0.12;  % Wider horizontal spread
        x_jitter_bl = max(min(x_jitter_bl, 0.4), -0.4);  % Constrain within violin
        
        scatter(x_jitter_bl, baseline_all, 30, [0.3 0.3 0.3], 'o', ...
            'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2, 'LineWidth', 0);
    end
    
    % DRUG VIOLIN (x=1) - Solid filled with scattered points
    if length(drug_all) >= 2
        [f_dr, xi_dr] = ksdensity(drug_all, 'NumPoints', 150);
        f_dr = f_dr / max(f_dr) * 0.4;  % Violin half-width
        
        % Draw solid violin fill
        x_violin_dr = [f_dr + 1, fliplr(-f_dr + 1)];
        y_violin_dr = [xi_dr, fliplr(xi_dr)];
        patch(x_violin_dr, y_violin_dr, [1 0.5 0.2], 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        
        % Add thick median line
        drug_median = median(drug_all);
        plot([1 - 0.4 1 + 0.4], [drug_median drug_median], 'k-', 'LineWidth', 3);
        
        % Add scattered points across violin
        n_drug = length(drug_all);
        x_jitter_dr = 1 + randn(n_drug, 1) * 0.12;  % Wider horizontal spread
        x_jitter_dr = max(min(x_jitter_dr, 1.4), 0.6);  % Constrain within violin
        
        scatter(x_jitter_dr, drug_all, 30, [0.3 0.3 0.3], 'o', ...
            'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.2, 'LineWidth', 0);
    end
    
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
        
        % BASELINE VIOLIN (x=0)
        if length(baseline_matched_all) >= 2
            [f_bl, xi_bl] = ksdensity(baseline_matched_all, 'NumPoints', 150);
            f_bl = f_bl / max(f_bl) * 0.4;  % Wider violin (0.8 total width)
            
            % Draw violin outline (no fill to see points through)
            x_violin_bl = [f_bl, fliplr(-f_bl)];
            y_violin_bl = [xi_bl, fliplr(xi_bl)];
            patch(x_violin_bl, y_violin_bl, [0.8 0.2 0.2], 'FaceAlpha', 0.3, 'EdgeColor', [0.6 0 0], 'LineWidth', 3);
            
            % Compute quartiles
            q1_bl = prctile(baseline_matched_all, 25);
            q3_bl = prctile(baseline_matched_all, 75);
            
            % Draw box plot overlay (quartile box)
            box_height_bl = q3_bl - q1_bl;
            rect_bl = rectangle('Position', [-0.15, q1_bl, 0.3, box_height_bl], ...
                'EdgeColor', [0.3 0 0], 'LineWidth', 3, 'FaceColor', 'none');
            
            % Draw whiskers
            whisker_extend = 1.5 * box_height_bl;
            whisker_lower_bl = max(min(baseline_matched_all), q1_bl - whisker_extend);
            whisker_upper_bl = min(max(baseline_matched_all), q3_bl + whisker_extend);
            plot([-0.075 0.075], [whisker_lower_bl whisker_lower_bl], 'k-', 'LineWidth', 2);
            plot([0 0], [whisker_lower_bl q1_bl], 'k-', 'LineWidth', 2);
            plot([0 0], [q3_bl whisker_upper_bl], 'k-', 'LineWidth', 2);
            plot([-0.075 0.075], [whisker_upper_bl whisker_upper_bl], 'k-', 'LineWidth', 2);
            
            % Draw thick median line
            plot([-0.15 0.15], [baseline_median baseline_median], 'k-', 'LineWidth', 4);
            
            % Add horizontal grid at quartiles
            plot(xlim, [q1_bl q1_bl], 'k--', 'LineWidth', 1, 'Alpha', 0.3);
            plot(xlim, [baseline_median baseline_median], 'k-', 'LineWidth', 1.5, 'Alpha', 0.5);
            plot(xlim, [q3_bl q3_bl], 'k--', 'LineWidth', 1, 'Alpha', 0.3);
            
            % Jittered points with density coloring
            n_baseline = length(baseline_matched_all);
            x_jitter_bl = 0 + randn(n_baseline, 1) * 0.08;
            x_jitter_bl = max(min(x_jitter_bl, 0.4), -0.4);
            
            % Compute local density for color mapping
            [density_bl, ~] = ksdensity([x_jitter_bl, baseline_matched_all], [x_jitter_bl, baseline_matched_all]);
            density_bl = density_bl / max(density_bl);
            
            scatter(x_jitter_bl, baseline_matched_all, 50, density_bl, 'o', ...
                'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.6, 'LineWidth', 0.5);
            colormap(gca, 'hot');
        end
        
        % DRUG VIOLIN (x=1)
        if length(drug_matched_all) >= 2
            [f_dr, xi_dr] = ksdensity(drug_matched_all, 'NumPoints', 150);
            f_dr = f_dr / max(f_dr) * 0.4;
            
            % Draw violin outline
            x_violin_dr = [f_dr + 1, fliplr(-f_dr + 1)];
            y_violin_dr = [xi_dr, fliplr(xi_dr)];
            patch(x_violin_dr, y_violin_dr, [0.2 0.5 0.95], 'FaceAlpha', 0.3, 'EdgeColor', [0 0.2 0.6], 'LineWidth', 3);
            
            % Compute quartiles
            q1_dr = prctile(drug_matched_all, 25);
            q3_dr = prctile(drug_matched_all, 75);
            drug_median_val = median(drug_matched_all);
            
            % Draw box plot overlay
            box_height_dr = q3_dr - q1_dr;
            rect_dr = rectangle('Position', [1 - 0.15, q1_dr, 0.3, box_height_dr], ...
                'EdgeColor', [0 0.15 0.5], 'LineWidth', 3, 'FaceColor', 'none');
            
            % Draw whiskers
            whisker_extend = 1.5 * box_height_dr;
            whisker_lower_dr = max(min(drug_matched_all), q1_dr - whisker_extend);
            whisker_upper_dr = min(max(drug_matched_all), q3_dr + whisker_extend);
            plot([1 - 0.075 1 + 0.075], [whisker_lower_dr whisker_lower_dr], 'k-', 'LineWidth', 2);
            plot([1 1], [whisker_lower_dr q1_dr], 'k-', 'LineWidth', 2);
            plot([1 1], [q3_dr whisker_upper_dr], 'k-', 'LineWidth', 2);
            plot([1 - 0.075 1 + 0.075], [whisker_upper_dr whisker_upper_dr], 'k-', 'LineWidth', 2);
            
            % Draw thick median line
            plot([1 - 0.15 1 + 0.15], [drug_median_val drug_median_val], 'k-', 'LineWidth', 4);
            
            % Add horizontal grid at quartiles
            plot(xlim, [q1_dr q1_dr], 'k--', 'LineWidth', 1, 'Alpha', 0.3);
            plot(xlim, [drug_median_val drug_median_val], 'k-', 'LineWidth', 1.5, 'Alpha', 0.5);
            plot(xlim, [q3_dr q3_dr], 'k--', 'LineWidth', 1, 'Alpha', 0.3);
            
            % Jittered points with density coloring
            n_drug = length(drug_matched_all);
            x_jitter_dr = 1 + randn(n_drug, 1) * 0.08;
            x_jitter_dr = max(min(x_jitter_dr, 1.4), 0.6);
            
            % Compute local density
            [density_dr, ~] = ksdensity([x_jitter_dr, drug_matched_all], [x_jitter_dr, drug_matched_all]);
            density_dr = density_dr / max(density_dr);
            
            scatter(x_jitter_dr, drug_matched_all, 50, density_dr, 'o', ...
                'MarkerFaceAlpha', 0.4, 'MarkerEdgeAlpha', 0.6, 'LineWidth', 0.5);
            colormap(gca, 'hot');
        end
        
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

function analyze_grating_stimulus(~, ~, ~, ~, ~, ~, ~, ~)
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
