%% PREPROCESSING AND DATA ORGANIZATION SCRIPT
% Consolidates: Master_Stimulus_Analysis, DataGathering, axon_analysis_V1
% Purpose: Load, organize, and chunk stimulus responses for subsequent analysis/plotting
%
% Workflow:
% 1. Load preprocessed recordings (rec1, rec2)
% 2. Define plane-based neuron indexing via centroids
% 3. Load ROI matching (if available)
% 4. Print available stimuli and select types to analyze
% 5. Chunk responses by stimulus parameters (grating, looming, moving_bar)
% 6. Extract stimulus timings for triggered average responses (flashes, checkers, spontaneous)
% 7. Save organized structure for plotting

clear all; close all;

%% ============================================================================
%  SECTION 1: FILE PATHS AND CONFIGURATION
%% ============================================================================

% Define recording files
rec1_file = 'C:\Users\rbondarenko\projects\DATA_2P\Animal25_240716_1453_preprocessed.mat';
rec2_file = 'C:\Users\rbondarenko\projects\DATA_2P\Animal25_240717_1449_preprocessed.mat';

% Output directory for organized data
output_dir = 'C:\Users\rbondarenko\projects\IntermediateDir\AnimalRB19\';
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

% ROI matching file (optional)
roi_matching_file = ''; % Set to full path if available, empty string if not

% Plane to analyze (1-indexed: 1, 2, or 3 for your data)
selected_plane = 1; % Change this to 1, 2, or 3 to analyze different planes

%% ============================================================================
%  SECTION 2: LOAD RECORDINGS
%% ============================================================================

fprintf('\n========== LOADING RECORDINGS ==========\n');

% Load rec1
fprintf('Loading rec1 from: %s\n', rec1_file);
B = load(rec1_file);
fprintf('  - Neurons: %d, Samples: %d\n', size(B.CaData(1).Ca_dFF, 1), size(B.CaData(1).Ca_dFF, 2));

% Load rec2
fprintf('Loading rec2 from: %s\n', rec2_file);
D = load(rec2_file);
fprintf('  - Neurons: %d, Samples: %d\n', size(D.CaData(1).Ca_dFF, 1), size(D.CaData(1).Ca_dFF, 2));

% Store in consistent structure for later
recordings = struct();
recordings.rec1.data = B;
recordings.rec2.data = D;

% Sampling rates
fs_b = 1 / mean(diff(B.TimeCa));
fs_d = 1 / mean(diff(D.TimeCa));
recordings.rec1.fs = fs_b;
recordings.rec2.fs = fs_d;
fprintf('Sampling rates: rec1 = %.1f Hz, rec2 = %.1f Hz\n', fs_b, fs_d);

%% ============================================================================
%  SECTION 3: EXTRACT AND DEFINE CENTROID FOR PLANE INDEXING
%% ============================================================================

fprintf('\n========== SETTING UP PLANE-BASED NEURON INDEXING ==========\n');

% Extract centroids for plane selection (for both recordings)
recordings.rec1.centroids = B.CaData(1).Ca_centroid_voxel;
recordings.rec2.centroids = D.CaData(1).Ca_centroid_voxel;

% Function to get neuron indices for specific plane
get_plane_neurons = @(centroids, plane_idx) find(centroids(:, 3) == plane_idx);

% First, get available planes
rec1_unique_planes = sort(unique(B.CaData(1).Ca_centroid_voxel(:, 3)));
rec2_unique_planes = sort(unique(D.CaData(1).Ca_centroid_voxel(:, 3)));

fprintf('Available planes in recordings (1-indexed mapping to Z-values):\n');
fprintf('  rec1 planes: ');
for i = 1:length(rec1_unique_planes)
    fprintf('[%d→Z=%.0f] ', i, rec1_unique_planes(i));
end
fprintf('\n');
fprintf('  rec2 planes: ');
for i = 1:length(rec2_unique_planes)
    fprintf('[%d→Z=%.0f] ', i, rec2_unique_planes(i));
end
fprintf('\n\n');

% Extract plane neurons using the selected_plane index
plane_z_rec1 = rec1_unique_planes(selected_plane);
plane_z_rec2 = rec2_unique_planes(selected_plane);

recordings.rec1.plane_neurons = get_plane_neurons(recordings.rec1.centroids, plane_z_rec1);
recordings.rec2.plane_neurons = get_plane_neurons(recordings.rec2.centroids, plane_z_rec2);

fprintf('Selected plane index: %d\n', selected_plane);
fprintf('  rec1 (Z=%.0f): %d neurons (total: %d)\n', plane_z_rec1, ...
    length(recordings.rec1.plane_neurons), size(B.CaData(1).Ca_dFF, 1));
fprintf('  rec2 (Z=%.0f): %d neurons (total: %d)\n', plane_z_rec2, ...
    length(recordings.rec2.plane_neurons), size(D.CaData(1).Ca_dFF, 1));

% Extract dFF for selected plane
recordings.rec1.dff = B.CaData(1).Ca_dFF(recordings.rec1.plane_neurons, :);
recordings.rec2.dff = D.CaData(1).Ca_dFF(recordings.rec2.plane_neurons, :);

%% ============================================================================
%  SECTION 4: LOAD ROI MATCHING (if available)
%% ============================================================================

fprintf('\n========== CHECKING ROI MATCHING ==========\n');

if ~isempty(roi_matching_file) && isfile(roi_matching_file)
    fprintf('Loading ROI matching from: %s\n', roi_matching_file);
    M = load(roi_matching_file);
    recordings.roi_matching = M;
    fprintf('ROI matching loaded: %s\n', fieldnames(M){1});
else
    fprintf('No ROI matching file provided or file not found.\n');
    fprintf('  Set roi_matching_file variable if you want to load matching later.\n');
    recordings.roi_matching = [];
end

%% ============================================================================
%  SECTION 5: PRINT AVAILABLE STIMULI
%% ============================================================================

fprintf('\n========== AVAILABLE STIMULI ==========\n');

fprintf('\nrec1 stimuli:\n');
for i = 1:length(B.Stimuli)
    fprintf('  [%d] %s\n', i, B.Stimuli(i).type);
end

fprintf('\nrec2 stimuli:\n');
for i = 1:length(D.Stimuli)
    fprintf('  [%d] %s\n', i, D.Stimuli(i).type);
end

%% ============================================================================
%  SECTION 6: SELECT STIMULI TO ANALYZE
%% ============================================================================

fprintf('\n========== SELECTING STIMULI TO ANALYZE ==========\n');

% Define which stimulus types to process
stimulus_types_to_analyze = {'spontaneous','sparse_local_global_flashes', 'checkers2', 'grating', 'looming'};

fprintf('Stimulus types to analyze:\n');
for i = 1:length(stimulus_types_to_analyze)
    fprintf('  - %s\n', stimulus_types_to_analyze{i});
end

%% ============================================================================
%  SECTION 7: IDENTIFY STIMULUS INDICES BY TYPE
%% ============================================================================

% Map stimulus type names to indices in each recording
recordings.rec1.stim_indices = find_stimulus_by_type(B.Stimuli, stimulus_types_to_analyze);
recordings.rec2.stim_indices = find_stimulus_by_type(D.Stimuli, stimulus_types_to_analyze);

fprintf('\nStimulus indices in each recording:\n');
fprintf('rec1:\n');
for stype = fieldnames(recordings.rec1.stim_indices)'
    idx = recordings.rec1.stim_indices.(stype{1});
    if ~isempty(idx)
        fprintf('  %s: index %d\n', stype{1}, idx);
    end
end

fprintf('rec2:\n');
for stype = fieldnames(recordings.rec2.stim_indices)'
    idx = recordings.rec2.stim_indices.(stype{1});
    if ~isempty(idx)
        fprintf('  %s: index %d\n', stype{1}, idx);
    end
end

%% ============================================================================
%  SECTION 8: DATA GATHERING (CHUNKING) FOR PARAMETRIC STIMULI
%% ============================================================================

fprintf('\n========== DATA GATHERING (CHUNKING) FOR PARAMETRIC STIMULI ==========\n');

% Stimuli to chunk by parameters
chunking_stimuli = {'grating', 'looming', 'moving_bar'};

results_chunked = struct();

for rec_name = {'rec1', 'rec2'}
    rec_key = rec_name{1};
    fprintf('\nProcessing %s:\n', rec_key);
    
    % Get data structure
    if strcmp(rec_key, 'rec1')
        data = B;
        dff_data = recordings.rec1.dff;
        stim_inds = recordings.rec1.stim_indices;
    else
        data = D;
        dff_data = recordings.rec2.dff;
        stim_inds = recordings.rec2.stim_indices;
    end
    
    % Process each chunking stimulus type
    for stim_type = chunking_stimuli
        stim_name = stim_type{1};
        
        % Check if this stimulus exists in this recording
        if ~isfield(stim_inds, stim_name) || isempty(stim_inds.(stim_name))
            fprintf('  %s: NOT FOUND\n', stim_name);
            continue;
        end
        
        stim_idx = stim_inds.(stim_name);
        fprintf('  %s: ', stim_name);
        
        % Get stimulus parameters and organize by parameter
        try
            stim_responses = organize_stimulus_responses_by_parameter(...
                data.Stimuli, stim_idx, dff_data, stim_name);
            
            results_chunked.(rec_key).(stim_name) = stim_responses;
            
            fprintf('OK - %d unique parameter values, %d presentations total\n', ...
                length(stim_responses.param_values), ...
                sum(stim_responses.n_presentations_per_param));
        catch ME
            fprintf('ERROR - %s\n', ME.message);
        end
    end
end

%% ============================================================================
%  SECTION 9: EXTRACT STIMULUS TIMINGS FOR TRIGGERED AVERAGES
%% ============================================================================

fprintf('\n========== EXTRACTING STIMULUS TIMINGS FOR TRIGGERED AVERAGES ==========\n');

% Stimuli to handle via triggered averages
triggered_stimuli = {'spontaneous', 'sparse_local_global_flashes', 'checkers2'};

results_triggered = struct();

for rec_name = {'rec1', 'rec2'}
    rec_key = rec_name{1};
    fprintf('\nProcessing %s:\n', rec_key);
    
    % Get data structure
    if strcmp(rec_key, 'rec1')
        data = B;
        dff_data = recordings.rec1.dff;
        stim_inds = recordings.rec1.stim_indices;
        fs = recordings.rec1.fs;
    else
        data = D;
        dff_data = recordings.rec2.dff;
        stim_inds = recordings.rec2.stim_indices;
        fs = recordings.rec2.fs;
    end
    
    % Process each triggered stimulus type
    for stim_type = triggered_stimuli
        stim_name = stim_type{1};
        
        % Check if this stimulus exists in this recording
        if ~isfield(stim_inds, stim_name) || isempty(stim_inds.(stim_name))
            fprintf('  %s: NOT FOUND\n', stim_name);
            continue;
        end
        
        stim_idx = stim_inds.(stim_name);
        fprintf('  %s: ', stim_name);
        
        try
            % Extract stimulus timings with dFF data window
            if strcmp(rec_key, 'rec1')
                time_ca = B.TimeCa;
            else
                time_ca = D.TimeCa;
            end
            
            stim_timings = extract_stimulus_timings(...
                data.Stimuli(stim_idx), time_ca, ...
                'name', stim_name, ...
                'data', dff_data);
            
            results_triggered.(rec_key).(stim_name) = stim_timings;
            
            fprintf('OK - %d presentations, %d calcium frames\n', ...
                stim_timings.n_presentations, ...
                stim_timings.n_ca_frames);
        catch ME
            fprintf('ERROR - %s\n', ME.message);
        end
    end
end

%% ============================================================================
%  SECTION 10: ASSEMBLE OUTPUT STRUCTURE
%% ============================================================================

fprintf('\n========== ASSEMBLING OUTPUT STRUCTURE ==========\n');

% Create master output structure
output = struct();
output.recordings = recordings;
output.results_chunked = results_chunked;
output.results_triggered = results_triggered;
output.parameters = struct(...
    'selected_plane', selected_plane, ...
    'stimulus_types_analyzed', {stimulus_types_to_analyze}, ...
    'chunking_stimuli', {chunking_stimuli}, ...
    'triggered_stimuli', {triggered_stimuli}, ...
    'file_rec1', rec1_file, ...
    'file_rec2', rec2_file);

% Save output
output_filename = fullfile(output_dir, 'preprocessed_stimulus_responses.mat');
fprintf('Saving output to: %s\n', output_filename);
save(output_filename, 'output', '-v7.3');
fprintf('DONE!\n\n');

%% ============================================================================
%  HELPER FUNCTIONS
%% ============================================================================

function stim_inds = find_stimulus_by_type(stimuli_array, type_list)
    %% Find stimulus indices by type name
    % Returns structure with field names matching type_list
    
    stim_inds = struct();
    
    for i = 1:length(type_list)
        type_name = type_list{i};
        
        % Search for matching stimulus
        idx_matches = find(strcmp({stimuli_array.type}, type_name));
        
        if ~isempty(idx_matches)
            stim_inds.(type_name) = idx_matches(1); % Take first match if multiple
        else
            stim_inds.(type_name) = [];
        end
    end
end

function stim_responses = organize_stimulus_responses_by_parameter(stimuli, stim_idx, dff_data, stim_type)
    %% Organize stimulus responses by parameter values (DataGathering approach)
    % Input:
    %   stimuli: array of stimulus structures from preprocessed file
    %   stim_idx: index of target stimulus in array
    %   dff_data: dFF matrix (neurons x time)
    %   stim_type: name of stimulus type ('grating', 'looming', 'moving_bar')
    % 
    % Output:
    %   stim_responses.param_values: unique parameter values (sorted)
    %   stim_responses.param_name: human-readable parameter name
    %   stim_responses.responses: cell array of response matrices
    %   stim_responses.n_presentations_per_param: count per parameter value
    %   stim_responses.n_neurons: number of neurons in dff_data
    
    stim = stimuli(stim_idx);
    n_neurons = size(dff_data, 1);
    
    % Get parameter values based on stimulus type
    [param_values, param_name, param_field] = get_stimulus_parameters(stim, stim_type);
    
    if isempty(param_values)
        error(sprintf('Could not extract parameters for stimulus type: %s', stim_type));
    end
    
    % Sort parameter values
    [param_values_sorted, sort_idx] = sort(param_values);
    
    % Initialize output
    stim_responses = struct();
    stim_responses.param_values = param_values_sorted;
    stim_responses.param_name = param_name;
    stim_responses.n_neurons = n_neurons;
    stim_responses.responses = cell(size(param_values_sorted));
    stim_responses.n_presentations_per_param = zeros(size(param_values_sorted), 1);
    
    % Group presentations by parameter value
    for p = 1:length(param_values_sorted)
        % Find all presentations with this parameter value
        param_matches = find(abs(param_values(sort_idx(p)) - param_values) < 1e-3);
        
        % Stack responses for this parameter
        % NOTE: Presentations might have variable lengths - handle with padding if needed
        n_presentations = length(param_matches);
        
        % Simplified: concatenate all presentations along dimension 3
        % TODO: Extract proper frame ranges per presentation from TimeStimulusFrame
        response_stack = [];
        for pp = 1:n_presentations
            % For now, use all dFF data
            % In production, extract frame range from stim.TimeStimulusFrame(param_matches(pp))
            response_stack = cat(3, response_stack, dff_data);
        end
        
        stim_responses.responses{p} = response_stack;
        stim_responses.n_presentations_per_param(p) = n_presentations;
    end
end

function [param_values, param_name, param_field] = get_stimulus_parameters(stim, stim_type)
    %% Extract parameter values from stimulus based on type
    
    param_values = [];
    param_name = '';
    param_field = '';
    
    switch lower(stim_type)
        case 'grating'
            if isfield(stim, 'specParams') && isfield(stim.specParams, 'grating_orientations')
                param_values = stim.specParams.grating_orientations;
                param_name = 'orientation (degrees)';
                param_field = 'grating_orientations';
            end
            
        case 'moving_bar'
            if isfield(stim, 'specParams') && isfield(stim.specParams, 'bar_orientations')
                param_values = stim.specParams.bar_orientations;
                param_name = 'bar orientation (degrees)';
                param_field = 'bar_orientations';
            end
            
        case 'looming'
            if isfield(stim, 'specParams') && isfield(stim.specParams, 'loom_speed')
                param_values = stim.specParams.loom_speed;
                param_name = 'looming speed (deg/s)';
                param_field = 'loom_speed';
            end
            
        otherwise
            error(sprintf('Unknown stimulus type: %s', stim_type));
    end
end

function stim_timings = extract_stimulus_timings(stim, time_ca_vector, varargin)
    %% Extract stimulus timings (time windows, frame indices, durations, etc.)
    % For triggered average responses
    % 
    % Input:
    %   stim: stimulus structure from preprocessed file
    %   time_ca_vector: TimeCa time vector (used to convert to frame indices)
    %   'name': stimulus name (optional)
    %   'data': dFF data matrix to extract (optional)
    %   'plane_idx': plane index for multi-plane data (optional)
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'name', '');
    addParameter(p, 'data', []);
    addParameter(p, 'plane_idx', 1);
    parse(p, varargin{:});
    stim_name = p.Results.name;
    dff_data = p.Results.data;
    plane_idx = p.Results.plane_idx;
    
    stim_timings = struct();
    stim_timings.name = stim_name;
    
    % Calculate stimulus time window using standard approach
    % stim_start = first frame in TimeStimulusFrame
    % stim_end = stim_start + (stimulus_trial_t * trials)
    
    if isfield(stim, 'TimeStimulusFrame') && ~isempty(stim.TimeStimulusFrame)
        stim_start_frame = stim.TimeStimulusFrame(1); % First stimulus frame (in seconds)
        
        % Calculate total stimulus duration
        if isfield(stim, 'stimulus_trial_t') && isfield(stim, 'trials')
            stim_total_time = stim.stimulus_trial_t * stim.trials;
            stim_end_frame = stim_start_frame + stim_total_time;
        else
            % Fallback: use last frame + estimated duration
            stim_end_frame = stim.TimeStimulusFrame(end) + stim.stimulus_trial_t;
        end
        
        % Find calcium frame indices within stimulus time window
        ca_indices = find(time_ca_vector > stim_start_frame & time_ca_vector < stim_end_frame);
        
        stim_timings.stim_start_time = stim_start_frame;
        stim_timings.stim_end_time = stim_end_frame;
        stim_timings.stim_duration = stim_total_time;
        stim_timings.ca_frame_indices = ca_indices;
        stim_timings.n_ca_frames = length(ca_indices);
        
        % Extract dFF for this stimulus window if provided
        if ~isempty(dff_data)
            stim_timings.dff_window = dff_data(:, ca_indices);
        end
        
        % Store individual presentation onsets for triggered analysis
        if size(stim.TimeStimulusFrame, 2) >= 1
            onset_times = stim.TimeStimulusFrame(:, 1);
            stim_timings.onset_times = onset_times;
            stim_timings.n_presentations = length(onset_times);
            
            % Convert onset times to frame indices
            [~, onset_frame_indices] = arrayfun(@(t) min(abs(time_ca_vector - t)), onset_times, 'UniformOutput', false);
            stim_timings.onset_frame_indices = cell2mat(onset_frame_indices);
        end
        
        % If available, extract durations (second column)
        if size(stim.TimeStimulusFrame, 2) >= 2
            durations = stim.TimeStimulusFrame(:, 2);
            stim_timings.durations_sec = durations;
        else
            stim_timings.durations_sec = [];
        end
    else
        warning(sprintf('No TimeStimulusFrame found for stimulus: %s', stim_name));
        stim_timings.stim_start_time = [];
        stim_timings.stim_end_time = [];
        stim_timings.ca_frame_indices = [];
        stim_timings.n_ca_frames = 0;
        stim_timings.dff_window = [];
        stim_timings.onset_times = [];
        stim_timings.onset_frame_indices = [];
        stim_timings.n_presentations = 0;
    end
end
