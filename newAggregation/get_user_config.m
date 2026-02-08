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
%           .CameraDir         - Camera/video data directory
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

% Set all RB-specific paths
config.AnimalFolder = 'Z:\Group Members\Rima\Animals\';
config.OutputDir = 'Z:\Group Members\Rima\Aggregated\';
config.OutputCopyDir = '';  % Backup directory (empty for now)
config.StimBasePath = 'C:\Users\rbondare\Stimulus_test';
config.DLCDir = 'Z:\Group Members\Rima\DLCData\';
config.IntermediateDir = 'C:\Users\rbondare\IntermediateDir\';
config.CameraDir = 'Z:\Group Members\Rima\Camera_OTHER\';
config.DATA_2P = 'Z:\Group Members\Rima\TEST\';

% Append animal name subfolder to IntermediateDir
config.IntermediateDir = fullfile(config.IntermediateDir, animal_name, filesep);

end
