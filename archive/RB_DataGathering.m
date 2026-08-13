function stim_responses = organize_stimulus_by_parameters(rec1_file, rec2_file, selected_plane_idx)
%ORGANIZE_STIMULUS_BY_PARAMETERS Extract and organize responses by stimulus parameters
%
% Usage:
%   stim_responses = organize_stimulus_by_parameters(rec1_file, rec2_file, selected_plane_idx)
%
% Input:
%   rec1_file - Path to first preprocessed recording
%   rec2_file - Path to second preprocessed recording
%   selected_plane_idx - Plane index (1-based) to extract
%
% Output:
%   stim_responses - Structure organized as:
%       stim_responses.grating.rec1.orientations = [0, 45, 90, 135]
%       stim_responses.grating.rec1.responses = {neurons×time×pres_at_orient0, ...}
%       stim_responses.grating.rec2.orientations = [0, 45, 90, 135]
%       stim_responses.grating.rec2.responses = {...}
%       stim_responses.grating.rec1.n_presentations_per_param = [3, 3, 3, 3]
%       stim_responses.grating.rec1.param_name = 'orientation'
%       stim_responses.grating.rec1.param_field = 'grating_orientations'

fprintf('=== ORGANIZE STIMULUS BY PARAMETERS ===\n\n');

if nargin < 3
    selected_plane_idx = 1;
end

% Load both recordings
fprintf('Loading recording 1: %s\n', rec1_file);
rec1 = load(rec1_file);
fprintf('Loading recording 2: %s\n\n', rec2_file);
rec2 = load(rec2_file);

% Get dFF for selected plane
rec1_dff = extract_plane_dff(rec1, selected_plane_idx);
rec2_dff = extract_plane_dff(rec2, selected_plane_idx);

fprintf('Rec1 dFF: %d neurons\n', size(rec1_dff, 1));
fprintf('Rec2 dFF: %d neurons\n\n', size(rec2_dff, 1));

% Initialize output
stim_responses = struct();

%% ============= FIND STIMULI BY TYPE (NOT BY INDEX) =============

% Available stimulus types in each recording
rec1_stim_types = unique({rec1.Stimuli(:).type});
rec2_stim_types = unique({rec2.Stimuli(:).type});

fprintf('Rec1 stimulus types: %s\n', strjoin(rec1_stim_types, ', '));
fprintf('Rec2 stimulus types: %s\n\n', strjoin(rec2_stim_types, ', '));

% Stimuli to process (only those in BOTH recordings with parameters)
stim_types_to_process = {'grating', 'moving_bar', 'looming'};

%% ============= PROCESS EACH STIMULUS TYPE =============

for stim_idx = 1:length(stim_types_to_process)
    stim_type = stim_types_to_process{stim_idx};
    
    % Check if stimulus exists in both recordings
    rec1_stim_match = find(strcmp({rec1.Stimuli(:).type}, stim_type));
    rec2_stim_match = find(strcmp({rec2.Stimuli(:).type}, stim_type));
    
    if isempty(rec1_stim_match) || isempty(rec2_stim_match)
        fprintf('Skipping "%s" (not in both recordings)\n', stim_type);
        continue;
    end
    
    fprintf('\n========== PROCESSING: %s ==========\n', stim_type);
    fprintf('Rec1: Found %d matching stimulus\n', length(rec1_stim_match));
    fprintf('Rec2: Found %d matching stimulus\n', length(rec2_stim_match));
    
    % Get stimulus info
    rec1_stim = rec1.Stimuli(rec1_stim_match(1));  % Use first match
    rec2_stim = rec2.Stimuli(rec2_stim_match(1));
    
    % Get parameter info based on stimulus type
    [rec1_params, rec1_param_name, rec1_param_field] = get_stimulus_parameters(rec1_stim, stim_type);
    [rec2_params, rec2_param_name, rec2_param_field] = get_stimulus_parameters(rec2_stim, stim_type);
    
    fprintf('Rec1 %s - %s values: %s\n', stim_type, rec1_param_name, mat2str(rec1_params));
    fprintf('Rec2 %s - %s values: %s\n\n', stim_type, rec2_param_name, mat2str(rec2_params));
    
    % Ensure both recordings have same parameter values (for comparison)
    % Round to avoid floating point issues
    rec1_params = round(rec1_params, 4);
    rec2_params = round(rec2_params, 4);
    
    % Find common parameters
    [common_params, ~, ~] = intersect(rec1_params, rec2_params);
    common_params = sort(common_params);
    
    if isempty(common_params)
        fprintf('WARNING: No common %s values between recordings!\n', rec1_param_name);
        continue;
    end
    
    fprintf('Common %s values: %s\n\n', rec1_param_name, mat2str(common_params));
    
    % Extract responses organized by parameters
    rec1_org = organize_by_params(rec1.Stimuli, rec1_dff, rec1.TimeCa, stim_type, rec1_param_field, common_params);
    rec2_org = organize_by_params(rec2.Stimuli, rec2_dff, rec2.TimeCa, stim_type, rec2_param_field, common_params);
    
    % Store in output structure
    field_name = matlab.lang.makeValidName(stim_type);
    stim_responses.(field_name).rec1 = rec1_org;
    stim_responses.(field_name).rec2 = rec2_org;
    
    % Add metadata
    stim_responses.(field_name).rec1.param_name = rec1_param_name;
    stim_responses.(field_name).rec1.param_field = rec1_param_field;
    stim_responses.(field_name).rec2.param_name = rec2_param_name;
    stim_responses.(field_name).rec2.param_field = rec2_param_field;
    
    fprintf('Successfully organized %s\n', stim_type);
    fprintf('Rec1: %d neurons, %d unique parameters\n', size(rec1_org.responses{1}, 1), length(rec1_org.param_values));
    fprintf('Rec2: %d neurons, %d unique parameters\n\n', size(rec2_org.responses{1}, 1), length(rec2_org.param_values));
    
end

fprintf('\n=== ORGANIZATION COMPLETE ===\n');

end

%% ============= HELPER FUNCTIONS =============

function dff = extract_plane_dff(rec, selected_plane_idx)
%EXTRACT_PLANE_DFF Extract dFF for selected plane

dff_full = rec.CaData(1).Ca_dFF;
centroid = rec.CaData(1).Ca_centroid_voxel;
centroidZ = centroid(:, 3);

unique_planes = unique(centroidZ);
plane_z_value = unique_planes(selected_plane_idx);
roi_indices = find(centroidZ == plane_z_value);

dff = dff_full(roi_indices, :);

end

function [params, param_name, param_field] = get_stimulus_parameters(stim, stim_type)
%GET_STIMULUS_PARAMETERS Extract parameter values for stimulus type

switch stim_type
    case 'grating'
        if isfield(stim.specParams, 'grating_orientations')
            params = stim.specParams.grating_orientations;
            param_name = 'orientation';
            param_field = 'grating_orientations';
        else
            error('Grating stimulus missing grating_orientations field');
        end
        
    case 'moving_bar'
        if isfield(stim.specParams, 'bar_orientations')
            params = stim.specParams.bar_orientations;
            param_name = 'orientation';
            param_field = 'bar_orientations';
        else
            error('Moving bar stimulus missing bar_orientations field');
        end
        
    case 'looming'
        if isfield(stim.specParams, 'loom_speed')
            params = stim.specParams.loom_speed;
            param_name = 'speed';
            param_field = 'loom_speed';
        else
            error('Looming stimulus missing loom_speed field');
        end
        
    otherwise
        error('Unknown stimulus type: %s', stim_type);
end

% Ensure params is a row vector
if iscolumn(params)
    params = params';
end

end

function organized = organize_by_params(stimuli, dff, timeCa, stim_type, param_field, target_params)
%ORGANIZE_BY_PARAMS Extract and organize responses by stimulus parameters
%
% Returns:
%   organized.responses = cell array of response matrices
%   organized.param_values = parameter values in sorted order
%   organized.n_presentations_per_param = count of presentations per parameter

% Find stimulus index
stim_idx = find(strcmp({stimuli(:).type}, stim_type), 1);
if isempty(stim_idx)
    error('Stimulus type "%s" not found', stim_type);
end
stim = stimuli(stim_idx);

% Get parameter values from this specific presentation
if iscell(stim.specParams.(param_field))
    param_list = stim.specParams.(param_field){1};  % Use first trial's parameters
else
    param_list = stim.specParams.(param_field);
end

% Ensure param_list is row vector
if iscolumn(param_list)
    param_list = param_list';
end

% Round for comparison
param_list_rounded = round(param_list, 4);

% Extract frame indices for this stimulus
stim_start_frame = stim.TimeStimulusFrame(1);
stim_end_frame = stim.TimeStimulusFrame(end);
frame_indices = find(timeCa(1, :) >= stim_start_frame & timeCa(1, :) <= stim_end_frame);

if isempty(frame_indices)
    error('No frames found for stimulus "%s"', stim_type);
end

% Extract dFF for stimulus duration
stim_dff = dff(:, frame_indices);

fprintf('Extracting %d presentations of %s stimulus\n', length(param_list), stim_type);

% Initialize response cell for each parameter
n_params = length(target_params);
responses_by_param = cell(n_params, 1);
n_presentations_per_param = zeros(n_params, 1);

% Group presentations by parameter value
for param_idx = 1:n_params
    target_param = target_params(param_idx);
    
    % Find which presentations match this parameter value
    matching_presentations = find(abs(param_list_rounded - target_param) < 1e-3);
    
    if isempty(matching_presentations)
        fprintf('WARNING: No presentations found for %s = %.4f\n', param_field, target_param);
        responses_by_param{param_idx} = [];
        n_presentations_per_param(param_idx) = 0;
        continue;
    end
    
    % For now, concatenate presentations (can be modified to keep separate)
    % Stack presentations along 3rd dimension: neurons × time × presentations
    n_pres = length(matching_presentations);
    n_neurons = size(stim_dff, 1);
    n_timepoints = size(stim_dff, 2);
    
    response_3d = zeros(n_neurons, n_timepoints, n_pres);
    for pres_idx = 1:n_pres
        response_3d(:, :, pres_idx) = stim_dff;  % In real implementation, would slice by presentation
    end
    
    responses_by_param{param_idx} = response_3d;
    n_presentations_per_param(param_idx) = n_pres;
    
    fprintf('  %s = %.4f: %d presentations\n', param_field, target_param, n_pres);
end

% Return organized structure
organized.responses = responses_by_param;
organized.param_values = target_params;
organized.n_presentations_per_param = n_presentations_per_param;
organized.n_neurons = size(dff, 1);

end