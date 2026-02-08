function config = get_user_config(animal_name, varargin)
%GET_USER_CONFIG Returns user-specific path configuration based on animal name prefix
%
%   config = get_user_config(animal_name)
%   config = get_user_config(animal_name, 'OperatorID', 'FS')
%
%   Detects the operator/user from the animal name prefix and returns
%   all user-specific paths for data processing.
%
%   Supported users:
%       FS  - Florian Schmidt (default, no prefix or AnimalFS*)
%       RB  - Rima Bondarenko (AnimalRB*)
%       ASD - Arnau Sans Dublanc (AnimalASD*)
%       HK  - Hakan Kücükdereli (AnimalHK*)
%
%   Input:
%       animal_name - Animal name string (e.g., 'AnimalRB5', 'AnimalFS12')
%
%   Optional Name-Value Arguments:
%       'OperatorID' - Override automatic detection (e.g., 'FS', 'RB')
%
%   Output:
%       config - Structure with user-specific paths:
%           .OperatorID        - User identifier string
%           .AnimalFolder      - Path to animal information files
%           .OutputDir         - Primary output directory for preprocessed files
%           .OutputCopyDir     - Backup output directory (can be empty)
%           .StimBasePath      - Base path for stimulus files
%           .DLCDir            - DeepLabCut processed data directory
%           .IntermediateDir   - Intermediate processing files directory
%           .CameraDir         - Camera/video data directory
%           .DATA_2P           - Two-photon calcium data directory
%
%   Example:
%       config = get_user_config('AnimalRB5');
%       config = get_user_config('AnimalFS12', 'OperatorID', 'FS');
%
%   See also: save_preprocessed, run_aggregation, get_experiment_list

% Parse inputs
p = inputParser;
addRequired(p, 'animal_name', @ischar);
addParameter(p, 'OperatorID', '', @ischar);
parse(p, animal_name, varargin{:});

operator_override = p.Results.OperatorID;

% Detect operator from animal name if not overridden
if isempty(operator_override)
    operator_id = detect_operator(animal_name);
else
    operator_id = upper(operator_override);
end

% Initialize config structure
config = struct();
config.OperatorID = operator_id;
config.AnimalName = animal_name;

% Set paths based on operator
switch operator_id
    case 'FS'
        % Florian Schmidt paths
        config.AnimalFolder = 'C:\Users\fschmidt\Documents\GitHub\MouseBrainActivity\MouseTrainingScripts\AnimalsISTA\';
        config.OutputDir = 'D:\Aggregated\';
        config.OutputCopyDir = 'Z:\Group Members\Florian\Aggregated\';
        config.StimBasePath = 'D:\Stimulus\';
        config.DLCDir = 'Z:\Group Members\Florian\DLCData\';
        config.IntermediateDir = 'D:\IntermediateDir\';
        config.CameraDir = 'Z:\Group Members\Florian\Camera\';
        config.DATA_2P = 'Z:\Group Members\Florian\DATA_2P\';

    case 'RB'
        % Rima Bondarenko paths
        config.AnimalFolder = 'Z:\Group Members\Rima\Animals\';
        config.OutputDir = 'Z:\Group Members\Rima\Aggregated\';
        config.OutputCopyDir = 'C:\Users\rbondare\Aggregated\';
        config.StimBasePath = 'Z:\Group Members\Rima\Stimulus\';
        config.DLCDir = 'Z:\Group Members\Rima\DLCData\';
        config.IntermediateDir = 'C:\Users\rbondare\IntermediateDir\';
        config.CameraDir = 'Z:\Group Members\Rima\Camera_OTHER\';
        config.DATA_2P = 'Z:\Group Members\Rima\TEST\';

    case 'ASD'
        % Arnau Sans Dublanc paths
        config.AnimalFolder = 'B:\group\joeschgrp\Group Members\Arnau\Animals\';
        config.OutputDir = 'B:\group\joeschgrp\Group Members\Arnau\Aggregated\';
        config.OutputCopyDir = '';
        config.StimBasePath = 'B:\group\joeschgrp\Group Members\Arnau\Stimulus\';
        config.DLCDir = 'B:\group\joeschgrp\Group Members\Arnau\DLCData\';
        config.IntermediateDir = 'B:\group\joeschgrp\Group Members\Arnau\IntermediateDir\';
        config.CameraDir = 'B:\group\joeschgrp\Group Members\Arnau\Camera\';
        config.DATA_2P = 'B:\group\joeschgrp\Group Members\Arnau\DATA_2P\';

    case 'HK'
        % Hakan Kücükdereli paths
        config.AnimalFolder = 'Z:\Group Members\Hakan\Animals\';
        config.OutputDir = 'Z:\Group Members\Hakan\Aggregated\';
        config.OutputCopyDir = '';
        config.StimBasePath = 'Z:\Group Members\Hakan\Stimulus\';
        config.DLCDir = 'Z:\Group Members\Hakan\DLCData\';
        config.IntermediateDir = 'Z:\Group Members\Hakan\IntermediateDir\';
        config.CameraDir = 'Z:\Group Members\Hakan\Camera\';
        config.DATA_2P = 'Z:\Group Members\Hakan\DATA_2P\';

    otherwise
        error('get_user_config:UnknownOperator', ...
            'Unknown operator ID: %s. Supported: FS, RB, ASD, HK', operator_id);
end

% Append animal name subfolder to IntermediateDir
config.IntermediateDir = fullfile(config.IntermediateDir, animal_name, filesep);

end


function operator_id = detect_operator(animal_name)
%DETECT_OPERATOR Detect operator ID from animal name prefix
%
%   Detection pattern based on FS_potentScript.m (lines 87-103):
%   - AnimalRB* -> RB
%   - AnimalASD* -> ASD
%   - AnimalHK* -> HK
%   - Otherwise -> FS (default)

if contains(animal_name, 'RB', 'IgnoreCase', true)
    operator_id = 'RB';
elseif contains(animal_name, 'ASD', 'IgnoreCase', true)
    operator_id = 'ASD';
elseif contains(animal_name, 'HK', 'IgnoreCase', true)
    operator_id = 'HK';
else
    % Default to FS
    operator_id = 'FS';
end

end
