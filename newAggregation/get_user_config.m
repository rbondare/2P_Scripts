function config = get_user_config(animal_name, varargin)
%GET_USER_CONFIG Returns RB-specific path configuration for data processing
%
%   config = get_user_config(animal_name)
%
%   Returns all RB-specific paths for data processing.
%   Simplified version that only supports RB (Rima Bondarenko) pathways.
%
%   Input:
%       animal_name - Animal name string (e.g., 'AnimalRB5')
%
%   Output:
%       config - Structure with RB-specific paths:
%           .OperatorID        - User identifier string ('RB')
%           .AnimalFolder      - Path to animal information files
%           .OutputDir         - Primary output directory for preprocessed files
%           .OutputCopyDir     - Backup output directory (can be empty)
%           .StimBasePath      - Base path for stimulus files
%           .DLCDir            - DeepLabCut processed data directory
%           .IntermediateDir   - Intermediate processing files directory
%           .CameraDir         - Folder with the DLC-processed camera sessions (<session>\success.log,
%                                eye_results, body_results, whisker_results); import_video_behavior archives
%                                them to the sibling folder preCamera afterwards
%           .DATA_2P           - Two-photon calcium data directory
%
%   Example:
%       config = get_user_config('AnimalRB5');
%
%   See also: save_preprocessed, run_aggregation, get_experiment_list

% Parse inputs
p = inputParser;
addRequired(p, 'animal_name', @ischar);
parse(p, animal_name, varargin{:});

% Initialize config structure with RB settings
config = struct();
config.OperatorID = 'RB';
config.AnimalName = animal_name;

% Rima's folder on the group share. Z: has been mapped both to the group root (Z:\joeschgrp\Group Members\...)
% and to joeschgrp itself (Z:\Group Members\...); use whichever exists.
persistent rima_root;
if isempty(rima_root)
    if exist('Z:\joeschgrp\Group Members\Rima', 'dir') == 7
        rima_root = 'Z:\joeschgrp\Group Members\Rima\';
    else
        rima_root = 'Z:\Group Members\Rima\';
    end
end

% Set all RB-specific paths
config.AnimalFolder = [rima_root 'Animals\'];
config.OutputDir = [rima_root 'Aggregated\'];
config.OutputCopyDir = '';  % Backup directory (empty for now)
%config.StimBasePath = [rima_root 'Stimulus\'];          % RB19 and earlier
config.StimBasePath = [rima_root 'Stimulus_Axons\'];     % RB22-24 (axon experiments)
config.DLCDir = [rima_root 'DLCData\'];
config.IntermediateDir = 'C:\Users\rbondarenko\projects\IntermediateDir\';
%config.CameraDir = [rima_root 'preCamera\'];            % archive of imported sessions, not the input
config.CameraDir = [rima_root 'Camera\'];                % DLC + whisker output (run_dlc.bat) lands here
%config.DATA_2P = [rima_root 'DATA_2P_OTHER\AnimalRB6\'];
config.DATA_2P = [rima_root 'DATA_2P\'];

% Set all DG-specific paths

% config.AnimalFolder = 'Z:\joeschgrp\Group Members\Rima\Animals\';
% config.OutputDir = 'Z:\joeschgrp\Group Members\Rima\Aggregated\';
% config.OutputCopyDir = '';  % Backup directory (empty for now)
% %config.StimBasePath = 'Z:\Group Members\Rima\Stimulus';
% config.StimBasePath = 'Z:\joeschgrp\Group Members\Rima\Stimulus_Axons\';
% config.DLCDir = 'Z:\joeschgrp\Group Members\Rima\DLCData\';
% config.IntermediateDir = 'C:\Users\rbondarenko\projects\IntermediateDir\';
% config.CameraDir = 'Z:\joeschgrp\Group Members\Rima\preCamera\';
% config.DATA_2P = 'Z:\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\';
% %config.DATA_2P = 'C:\Users\rbondarenko\projects\DATA_2P\TEST\';

% Append animal name subfolder to IntermediateDir
config.IntermediateDir = fullfile(config.IntermediateDir, animal_name, filesep);

% Display warning if paths don't exist (helps with debugging)
persistent warned_paths;
if isempty(warned_paths)
    warned_paths = containers.Map;
end

key = sprintf('%s_%s', animal_name, config.StimBasePath);
if ~warned_paths.isKey(key)
    if ~exist(config.StimBasePath, 'dir')
        warning('get_user_config:PathNotFound', ...
            'StimBasePath not found: %s (looking for %s data)', ...
            config.StimBasePath, animal_name);
    end
    warned_paths(key) = true;
end

end
