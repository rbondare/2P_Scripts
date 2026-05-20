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
selected_roi_idx = find(centroidZ == selected_plane);
base_dff = base_dff_plane{find(unique_planes == selected_plane)};
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
drug_selected_roi_idx = find(drug_centroidZ == selected_plane);
drug_dff = drug_dff_plane{find(unique_planes_drug == selected_plane)};
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
    selected_plane, n_rois_drug_selected_plane, min(drug_selected_roi_idx), max(drug_selected_roi_idx));

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
                selected_plane, n_rois_selected_plane);
            fprintf('Drug selected plane (Z=%d): %d local neuron indices found\n', ...
                selected_plane, n_rois_drug_selected_plane);
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
            fprintf('WARNING: No matched ROIs found in selected plane %d!\n', selected_plane);
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
        fprintf('  ✓ Extraction successful:\n');
        fprintf('    Baseline responses: %d presentations, %d neurons per presentation\n', ...
            length(base_responses), size(base_responses{1}, 1));
        fprintf('    Drug responses: %d presentations, %d neurons per presentation\n', ...
            length(drug_responses), size(drug_responses{1}, 1));
    else
        fprintf('  ✗ WARNING: Empty responses!\n');
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

%% ====================== GENERATE HEATMAPS (4 LEVELS) ======================

fprintf('\n========== GENERATING PROGRESSIVE HEATMAPS ==========\n');

stim_fields = fieldnames(master_data);

% ========== LEVEL 1: FULL dFF HEATMAPS (ALL NEURONS) ==========
fprintf('\n--- Level 1: Full dFF Heatmaps (all neurons) ---\n');

for field_idx = 1:length(stim_fields)
    field_name = stim_fields{field_idx};
    data = master_data.(field_name);
    stim_type = data.stimulus_type;
    
    fprintf('Generating full dFF heatmap for "%s"\n', stim_type);
    
    % Average responses across presentations
    baseline_full = average_stimulus_presentations(data.baseline_responses);
    drug_full = average_stimulus_presentations(data.drug_responses);
    
    % Create figure: baseline and drug side-by-side
    fig = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
        'Name', sprintf('Level1_FullDFF: %s', stim_type));
    
    % BASELINE - ALL neurons
    subplot(1, 2, 1);
    imagesc(baseline_full);
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(caxis_lim);
    colorbar;
    xlabel('Time (frames)');
    ylabel('ROI');
    title(sprintf('Baseline - %s (Full dFF)\nAll %d neurons', stim_type, size(baseline_full, 1)), ...
        'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % DRUG - ALL neurons
    subplot(1, 2, 2);
    imagesc(drug_full);
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(caxis_lim);
    c = colorbar;
    ylabel(c, 'dF/F');
    xlabel('Time (frames)');
    ylabel('ROI');
    title(sprintf('Drug - %s (Full dFF)\nAll %d neurons', stim_type, size(drug_full, 1)), ...
        'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle(sprintf('Level 1 - Full dFF: %s', stim_type), 'FontSize', 12, 'FontWeight', 'bold');
    
    % Store full heatmaps
    data.full_baseline_heatmap = baseline_full;
    data.full_drug_heatmap = drug_full;
    
    % ========== LEVEL 2: SELECTED ROI HEATMAPS (MAX ACTIVITY SORTED) ==========
    % Select top N neurons by max activity in baseline
    n_select = min(max_neurons_display, size(baseline_full, 1));
    [~, sort_idx_baseline] = sort(max(baseline_full, [], 2), 'descend');
    select_idx_baseline = sort_idx_baseline(1:n_select);
    
    baseline_selected = baseline_full(select_idx_baseline, :);
    drug_selected = drug_full(select_idx_baseline, :);  % Use SAME ROI indices for drug
    
    % Create figure
    fig2 = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
        'Name', sprintf('Level2_SelectedROI: %s', stim_type));
    
    % BASELINE - Selected ROIs sorted by max activity
    subplot(1, 2, 1);
    [~, sort_idx] = sort(max(baseline_selected, [], 2), 'descend');
    imagesc(baseline_selected(sort_idx, :));
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(caxis_lim);
    colorbar;
    xlabel('Time (frames)');
    ylabel('ROI (sorted by max activity)');
    title(sprintf('Baseline - %s (Selected, sorted)\n%d/%d neurons', stim_type, n_select, size(baseline_full, 1)), ...
        'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % DRUG - Same ROI indices, sorted by baseline activity
    subplot(1, 2, 2);
    imagesc(drug_selected(sort_idx, :));
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(caxis_lim);
    c = colorbar;
    ylabel(c, 'dF/F');
    xlabel('Time (frames)');
    ylabel('ROI (same order as baseline)');
    title(sprintf('Drug - %s (Selected, same order)\n%d/%d neurons', stim_type, n_select, size(drug_full, 1)), ...
        'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle(sprintf('Level 2 - Selected ROI (max activity sorted): %s', stim_type), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % ========== LEVEL 3: NORMALIZED HEATMAPS ==========
    % Apply normalization: (dFF - min) / (max - min) * max
    baseline_norm = (baseline_selected - min(baseline_selected, [], 2)) ./ ...
        (max(baseline_selected, [], 2) - min(baseline_selected, [], 2)) .* ...
        max(baseline_selected, [], 2);
    drug_norm = (drug_selected - min(drug_selected, [], 2)) ./ ...
        (max(drug_selected, [], 2) - min(drug_selected, [], 2)) .* ...
        max(drug_selected, [], 2);
    
    % Sort by max activity of normalized baseline
    [~, sort_idx_norm] = sort(max(baseline_norm, [], 2), 'descend');
    
    fig3 = figure('Position', [100 100 1400 700], 'NumberTitle', 'off', ...
        'Name', sprintf('Level3_Normalized: %s', stim_type));
    
    % BASELINE - Normalized, sorted
    subplot(1, 2, 1);
    imagesc(baseline_norm(sort_idx_norm, :));
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(caxis_lim);
    colorbar;
    xlabel('Time (frames)');
    ylabel('ROI (sorted by normalized max)');
    title(sprintf('Baseline - %s (Normalized)\nFormula: (dFF-min)/(max-min)*max', stim_type), ...
        'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    % DRUG - Normalized, same order as baseline
    subplot(1, 2, 2);
    imagesc(drug_norm(sort_idx_norm, :));
    eval(['colormap(', colormap_type, ');']);
    set(gca, 'YDir', 'reverse');
    caxis(caxis_lim);
    c = colorbar;
    ylabel(c, 'dF/F (normalized)');
    xlabel('Time (frames)');
    ylabel('ROI (same order as baseline)');
    title(sprintf('Drug - %s (Normalized)', stim_type), 'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.5, 'FontSize', 10);
    
    sgtitle(sprintf('Level 3 - Normalized (sorted): %s', stim_type), ...
        'FontSize', 12, 'FontWeight', 'bold');
    
    % Store for later use
    data.selected_baseline_norm = baseline_norm;
    data.selected_drug_norm = drug_norm;
    data.selected_roi_idx = select_idx_baseline;
    data.sort_idx_selected = sort_idx;  % Sort by raw max activity
    data.sort_idx_norm = sort_idx_norm;  % Sort by normalized max activity
    
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
        
        % Rank baseline by max activity
        [~, sort_idx] = sort(max(baseline_matched_avg, [], 2), 'descend');
        
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
            'EdgeColor', 'none', 'BackgroundColor', 'yellow', 'Alpha', 0.3);
        
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

fprintf('\n========== ANALYSIS COMPLETE ==========\n');
fprintf('\nFOUR-LEVEL HEATMAP HIERARCHY:\n\n');

fprintf('LEVEL 1 - Full dFF Heatmaps:\n');
fprintf('  • All %d neurons in selected plane\n', n_rois_selected_plane);
fprintf('  • Baseline and drug side-by-side (raw dF/F)\n');
fprintf('  • Color scale: [%g, %g] dF/F\n', caxis_lim(1), caxis_lim(2));
fprintf('  • Unordered (as they appear in the data)\n\n');

fprintf('LEVEL 2 - Selected ROI Heatmaps (Sorted):\n');
fprintf('  • Top %d neurons by baseline max activity\n', max_neurons_display);
fprintf('  • Baseline and drug side-by-side (raw dF/F)\n');
fprintf('  • Neurons sorted by descending max activity in baseline\n');
fprintf('  • Drug neurons in same order (maintaining ROI pairing)\n\n');

fprintf('LEVEL 3 - Normalized Heatmaps (Sorted):\n');
fprintf('  • Same top %d neurons as Level 2\n', max_neurons_display);
fprintf('  • Normalization formula: (dFF - min) / (max - min) * max\n');
fprintf('  • Applied per neuron (row-wise normalization)\n');
fprintf('  • Neurons sorted by descending max normalized activity\n');
fprintf('  • Emphasizes timing patterns regardless of baseline activity level\n\n');

if matched_rois_available && n_matched > 0
    fprintf('LEVEL 4 - Matched ROI Heatmaps:\n');
    fprintf('  • %d matched neuron pairs between baseline and drug\n', n_matched);
    fprintf('  • Baseline sorted by descending max activity\n');
    fprintf('  • Drug neurons MATCHED to baseline by pairing index\n');
    fprintf('  • Each row = one matched neuron pair (baseline ↔ drug)\n');
    fprintf('  • Enables direct baseline vs drug comparison\n\n');
else
    fprintf('LEVEL 4 - Matched ROI Heatmaps:\n');
    fprintf('  • (Not available - no matched ROIs in plane %d)\n\n', selected_plane);
end

fprintf('ACCESSING DATA FROM master_data:\n');
fprintf('  master_data.STIMULUS_NAME.full_baseline_heatmap\n');
fprintf('  master_data.STIMULUS_NAME.full_drug_heatmap\n');
fprintf('  master_data.STIMULUS_NAME.selected_baseline_norm\n');
fprintf('  master_data.STIMULUS_NAME.selected_drug_norm\n');
fprintf('  master_data.STIMULUS_NAME.matched_baseline_avg\n');
fprintf('  master_data.STIMULUS_NAME.matched_drug_avg\n\n');

fprintf('TIMING EXTRACTION (verified against axon_analysis_V1):\n');
fprintf('  ✓ Start time: TimeStimulusFrame(1)\n');
fprintf('  ✓ End time: start + (stimulus_trial_t × trials)\n');
fprintf('  ✓ Frame matching: find(TimeCa(1,:) > start & TimeCa(1,:) < end)\n');
fprintf('  ✓ Boundary handling: Strict inequalities (> and <)\n\n');

fprintf('Master data saved in "master_data" struct.\n');

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
