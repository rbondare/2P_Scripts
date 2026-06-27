%% Preprocessing: load baseline/drug recordings, match ROIs, set shared variables
%
% Script 1 of 3 (Preprocess_Stimulus_Data -> analysis_functions/ -> Plot_Stimulus_Results).
% This is the known/verified loading pipeline -- ported as-is. Run standalone
% to sanity-check loading/matching alone, or let Plot_Stimulus_Results.m call
% it directly (scripts share the caller's workspace, so all variables below
% land directly in the calling script).

clear; clc; close all;

%% ====================== CONFIGURATION ======================

% Data files
baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1444_preprocessed.mat";
drug_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260320_1609_preprocessed.mat";
roi_match_file = "C:\Users\rbondarenko\projects\2P_Scripts\roi_matching\roiMatch_plane0_1444_V1.mat";

% Analysis parameters
ca_type = 1;                  % 1=FVDff, 2=deconvolved, 3=F
selected_plane_idx = 1;        % Plane index (1-based)
stimulus_types_to_analyze = {'spontaneous','sparse_local_global_flashes', 'checkers2', 'grating', 'moving_bar', 'looming'};

% Output folder for this recording pair (Plot_Stimulus_Results.m appends
% per-stimulus and cross_stimulus subfolders under this).
[~, base_name, ~] = fileparts(baseline_file); base_name = regexprep(base_name, '_preprocessed$', '');
[~, drug_name, ~] = fileparts(drug_file);     drug_name = regexprep(drug_name, '_preprocessed$', '');
pair_outdir = fullfile(fileparts(mfilename('fullpath')), 'analysis_figures', sprintf('%s_vs_%s', base_name, drug_name));
if ~exist(pair_outdir, 'dir'); mkdir(pair_outdir); end

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
fprintf('\nPreprocessing complete. Figures for this pair will be saved under:\n  %s\n', pair_outdir);
