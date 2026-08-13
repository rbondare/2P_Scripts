%RUN_AGGREGATION Universal script for aggregating 2P preprocessing data
%
%   This script replaces RB_call_2SP_extraction.m with a more flexible approach:
%   - User auto-detection from animal names
%   - Manual experiment specification option
%   - Stimulus folder as default source (not DATA_2P)
%   - Support for behavior-only processing (require_calcium = false by default)
%
%   Usage:
%       1. Edit the configuration section below
%       2. Run the script
%
%   See also: save_preprocessed, get_experiment_list, get_user_config, import_video_behavior

%% ===== CONFIGURATION ===== 
% Edit these settings to customize processing

% --- Experiment Discovery ---
% 'source' options: 'stimulus', 'calcium', 'camera', 'manual'
config_source = 'manual';          % Where to look for experiments

% --- Manual Recording Specification ---
% Only used when config_source = 'manual'
config_recordings = {};              % e.g., {'AnimalRB5_20231201_001', 'AnimalRB14_20231215_002'}

% --- Animal Filtering ---
% Leave empty {} to process all animals found
config_animals = {'AnimalRB5'};                 % e.g., {'AnimalRB5', 'AnimalRB14'}

% --- Processing Options ---
config_require_calcium = false;      % false = allow behavior-only processing
config_require_stimulus = true;      % true = require stimulus files
config_check_processed = true;       % true = skip already processed recordings

% --- Strict Source Mode ---
% When true: ONLY discover recordings from the source folder (e.g., stimulus folder)
%            Does NOT look for additional recordings in DATA_2P or Camera folders
%            Calcium data is still loaded from DATA_2P, but only for recordings found in stimulus folder
%
% When false: Discover recordings from all available locations
%             (stimulus folder, DATA_2P, Camera folders)
config_strict_source = true;         % true = ONLY discover from source folder
                                     % false = discover from all locations

config_ca_type = 1;                  % Calcium data type: 1=Raw, 2=di_masks, 3=di
config_ca_format = 'auto';           % 'auto', 'suite2p', 'caiman', 'none'

% --- Import DLC After Aggregation ---
config_import_dlc = true;            % Run import_video_behavior after aggregation

%% ===== SETUP =====
% Save user configuration before clearing workspace
user_config_source = config_source;
user_config_recordings = config_recordings;
user_config_animals = config_animals;
user_config_require_calcium = config_require_calcium;
user_config_require_stimulus = config_require_stimulus;
user_config_check_processed = config_check_processed;
user_config_strict_source = config_strict_source;
user_config_ca_type = config_ca_type;
user_config_ca_format = config_ca_format;
user_config_import_dlc = config_import_dlc;

% Clear workspace except user configuration variables
clearvars -except user_config_*

% Add required paths
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

% Restore user configuration
config_source = user_config_source;
config_recordings = user_config_recordings;
config_animals = user_config_animals;
config_require_calcium = user_config_require_calcium;
config_require_stimulus = user_config_require_stimulus;
config_check_processed = user_config_check_processed;
config_strict_source = user_config_strict_source;
config_ca_type = user_config_ca_type;
config_ca_format = user_config_ca_format;
config_import_dlc = user_config_import_dlc;

%% ===== GET EXPERIMENT LIST =====
fprintf('=== Universal Aggregation Script ===\n');
fprintf('Source: %s\n', config_source);
fprintf('Require calcium: %d\n', config_require_calcium);
fprintf('Require stimulus: %d\n', config_require_stimulus);
fprintf('Strict source mode: %d', config_strict_source);
if config_strict_source
    fprintf(' (ONLY data in source folder will be used)\n');
else
    fprintf(' (will check multiple locations for data)\n');
end
if ~isempty(config_animals)
    fprintf('Animal filter: %s (ONLY these animals will be processed)\n', strjoin(config_animals, ', '));
    % Show config paths for first animal to verify settings
    if ~strcmp(config_source, 'manual') && numel(config_animals) > 0
        sample_config = get_user_config(config_animals{1});
        fprintf('  Config paths for %s:\n', config_animals{1});
        fprintf('    StimBasePath: %s\n', sample_config.StimBasePath);
        fprintf('    DATA_2P: %s\n', sample_config.DATA_2P);
        fprintf('    DLCDir: %s\n', sample_config.DLCDir);
        fprintf('    CameraDir: %s\n', sample_config.CameraDir);
    end
else
    fprintf('Animal filter: NONE (ALL animals will be included)\n');
    fprintf('WARNING: No animal filter set! All animals in the source folder will be processed.\n');
end
fprintf('\n');

% Get experiments to process
if strcmp(config_source, 'manual')
    if isempty(config_recordings)
        error('run_aggregation:NoRecordings', ...
            'source=''manual'' requires config_recordings to be set');
    end
    experiments = get_experiment_list(...
        'source', 'manual', ...
        'recordings', config_recordings, ...
        'require_calcium', config_require_calcium, ...
        'require_stimulus', config_require_stimulus, ...
        'check_processed', config_check_processed);
else
    experiments = get_experiment_list(...
        'source', config_source, ...
        'animals', config_animals, ...
        'require_calcium', config_require_calcium, ...
        'require_stimulus', config_require_stimulus, ...
        'check_processed', config_check_processed, ...
        'strict_source', config_strict_source);
end

if isempty(experiments)
    fprintf('No experiments found to process.\n');
    return;
end

%% ===== DISPLAY EXPERIMENTS =====
fprintf('\n=== Experiments to Process ===\n');
fprintf('%-40s %-12s %-10s %-10s %-10s %-10s\n', 'Recording', 'Animal', 'Calcium', 'Stimulus', 'DLC', 'Processed');
fprintf('%s\n', repmat('-', 1, 92));

% Validate all experiments match animal filter
if ~isempty(config_animals)
    for i = 1:numel(experiments)
        if ~ismember(experiments(i).animal_name, config_animals)
            error('run_aggregation:FilterError', ...
                'Found experiment from %s but only %s should be processed!', ...
                experiments(i).animal_name, strjoin(config_animals, ', '));
        end
    end
end

for i = 1:numel(experiments)
    exp = experiments(i);
    fprintf('%-40s %-12s %-10s %-10s %-10s %-10s\n', ...
        exp.name, ...
        exp.animal_name, ...
        yesno(exp.has_calcium), ...
        yesno(exp.has_stimulus), ...
        yesno(exp.has_dlc), ...
        yesno(exp.is_processed));
end
fprintf('\n');

%% ===== CONFIRM PROCESSING =====
reply = input(sprintf('Process %d recording(s)? [Y/n]: ', numel(experiments)), 's');
if ~isempty(reply) && ~strcmpi(reply, 'y')
    fprintf('Cancelled.\n');
    return;
end

%% ===== PROCESS RECORDINGS =====
fprintf('\n=== Starting Processing ===\n');

success_count = 0;
error_count = 0;
error_recordings = {};

for i = 1:numel(experiments)
    exp = experiments(i);

    fprintf('\n[%d/%d] Processing: %s\n', i, numel(experiments), exp.name);
    fprintf('  Animal: %s, Operator: %s\n', exp.animal_name, exp.config.OperatorID);
    
    % Safety check: Validate animal is in filter (if filter is set)
    if ~isempty(config_animals) && ~ismember(exp.animal_name, config_animals)
        fprintf('  SKIPPING: Animal %s not in filter [%s]\n', ...
            exp.animal_name, strjoin(config_animals, ', '));
        error_count = error_count + 1;
        error_recordings{end+1} = sprintf('%s (wrong animal)', exp.name); %#ok<SAGROW>
        continue;
    end

    % Determine calcium format based on data availability
    if ~exp.has_calcium
        ca_format = 'none';
        fprintf('  Mode: Behavior-only (no calcium data)\n');
    else
        ca_format = config_ca_format;
        fprintf('  Mode: Full processing (calcium format: %s)\n', ca_format);
    end


    % Determine recording path for processing
    if exp.has_calcium
        % Use DATA_2P folder (contains headers/TIF/suite2p)
        recording_path = fullfile(exp.config.DATA_2P, exp.name);
        fprintf('  Using DATA_2P path (has calcium data)\n');
        
    elseif exp.has_stimulus
        % Use stimulus folder for behavior-only processing
        recording_path = fullfile(exp.config.StimBasePath, exp.animal_name, exp.name);
        fprintf('  Using stimulus path (behavior-only mode)\n');
        
    else
        error('run_aggregation:NoData', 'No data found for recording: %s', exp.name);
    end

    % Process recording
    try
        OutputFilename = save_preprocessed(recording_path, ...
            'ca_type', config_ca_type, ...
            'ca_format', ca_format, ...
            'config', exp.config, ...
            'overwrite_intermediate', false);

        fprintf('  Output: %s\n', OutputFilename);
        success_count = success_count + 1;

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('    at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
        end
        error_count = error_count + 1;
        error_recordings{end+1} = exp.name; %#ok<SAGROW>
    end
end

%% ===== IMPORT DLC DATA =====
if config_import_dlc
    fprintf('\n=== Importing DLC Data ===\n');

    try
        % Pass animal filter to respect the same filtering as aggregation
        if ~isempty(config_animals)
            fprintf('Importing DLC only for: %s\n', strjoin(config_animals, ', '));
            import_video_behavior('animals', config_animals);
        else
            fprintf('Importing DLC for all animals\n');
            import_video_behavior();
        end
        fprintf('DLC import completed.\n');
    catch ME
        fprintf('DLC import error: %s\n', ME.message);
    end
end

%% ===== SUMMARY =====
fprintf('\n=== Processing Summary ===\n');
fprintf('Total recordings: %d\n', numel(experiments));
fprintf('Successful: %d\n', success_count);
fprintf('Errors: %d\n', error_count);

if ~isempty(error_recordings)
    fprintf('\nFailed recordings:\n');
    for i = 1:numel(error_recordings)
        fprintf('  - %s\n', error_recordings{i});
    end
end

fprintf('\n=== DONE ===\n');

%% ===== HELPER FUNCTIONS =====
function str = yesno(val)
%YESNO Convert boolean to Yes/No string
if val
    str = 'Yes';
else
    str = 'No';
end
end
