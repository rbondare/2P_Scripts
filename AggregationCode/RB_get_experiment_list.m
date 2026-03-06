function experiments = RB_get_experiment_list(varargin)
%GET_EXPERIMENT_LIST Flexible experiment discovery for aggregation pipeline
%
%   experiments = get_experiment_list()
%   experiments = get_experiment_list('source', 'stimulus')
%   experiments = get_experiment_list('animals', {'AnimalRB5', 'AnimalRB14'})
%   experiments = get_experiment_list('recordings', {'AnimalRB5_20231201_001'})
%
%   Discovers experiments/recordings from various sources for processing.
%   By default, uses stimulus folder as the source (not DATA_2P), allowing
%   behavior-only processing without requiring calcium data.
%
%   Optional Name-Value Arguments:
%       'source' - Where to look for experiments:
%           'stimulus'  - Scan stimulus folder (default)
%           'calcium'   - Scan DATA_2P folder (requires calcium data)
%           'camera'    - Scan camera folder
%           'manual'    - Use 'recordings' parameter only
%
%       'config' - User config structure from get_user_config
%           If not provided, will use default FS paths
%
%       'animals' - Cell array of animal names to filter by
%           e.g., {'AnimalRB5', 'AnimalRB14'}
%
%       'recordings' - Cell array of specific recording names
%           e.g., {'AnimalRB5_20231201_001'}
%           When source='manual', this is required
%
%       'require_calcium' - Whether to require calcium data exists (default: false)
%           false: Include recordings even without suite2p folder
%           true:  Only include recordings with suite2p folder
%
%       'require_stimulus' - Whether to require stimulus data exists (default: true)
%           true: Only include recordings with stimulus files
%           false: Include recordings even without stimulus data
%
%       'check_processed' - Whether to check for existing preprocessed files (default: true)
%           true: Skip already processed recordings
%           false: Include all recordings regardless of processing status
%
%       'strict_source' - Whether to only look in the source folder for data (default: false)
%           true: ONLY process recordings found in source folder, don't check other locations
%                 (e.g., if source='stimulus', don't check DATA_2P folder for calcium data)
%           false: Check all configured locations for data availability
%
%   Output:
%       experiments - Structure array with fields:
%           .name           - Recording name (e.g., 'AnimalRB5_20231201_001')
%           .animal_name    - Animal name (e.g., 'AnimalRB5')
%           .path           - Full path to recording folder
%           .has_calcium    - Boolean: suite2p folder exists
%           .has_stimulus   - Boolean: stimulus files exist
%           .has_dlc        - Boolean: DLC data exists
%           .is_processed   - Boolean: preprocessed file exists
%           .config         - User configuration for this animal
%
%   Example:
%       % Get all unprocessed RB experiments from stimulus folder
%       exps = get_experiment_list('source', 'stimulus', ...
%           'animals', {'AnimalRB5', 'AnimalRB14'});
%
%       % Get specific recordings for manual processing
%       exps = get_experiment_list('source', 'manual', ...
%           'recordings', {'AnimalRB5_20231201_001', 'AnimalRB14_20231215_002'});
%
%       % Get all experiments including those without calcium data
%       exps = get_experiment_list('require_calcium', false);
%
%   See also: get_user_config, save_preprocessed, run_aggregation

% Parse inputs
p = inputParser;
addParameter(p, 'source', 'stimulus', @(x) ismember(x, {'stimulus', 'calcium', 'camera', 'manual'}));
addParameter(p, 'config', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'animals', {}, @iscell);
addParameter(p, 'recordings', {}, @iscell);
addParameter(p, 'require_calcium', false, @islogical);
addParameter(p, 'require_stimulus', true, @islogical);
addParameter(p, 'check_processed', true, @islogical);
addParameter(p, 'strict_source', false, @islogical);
parse(p, varargin{:});

source = p.Results.source;
config = p.Results.config;
animals_filter = p.Results.animals;
recordings_filter = p.Results.recordings;
require_calcium = p.Results.require_calcium;
require_stimulus = p.Results.require_stimulus;
check_processed = p.Results.check_processed;
strict_source = p.Results.strict_source;

% Initialize empty experiments array
experiments = struct('name', {}, 'animal_name', {}, 'path', {}, ...
    'has_calcium', {}, 'has_stimulus', {}, 'has_dlc', {}, ...
    'is_processed', {}, 'config', {});

% Handle manual source
if strcmp(source, 'manual')
    if isempty(recordings_filter)
        warning('get_experiment_list:NoRecordings', ...
            'source=''manual'' requires ''recordings'' parameter');
        return;
    end

    for i = 1:numel(recordings_filter)
        rec_name = recordings_filter{i};
        exp = create_experiment_entry(rec_name, config);
        if ~isempty(exp)
            experiments(end+1) = exp; %#ok<AGROW>
        end
    end
    return;
end

% Get base directory based on source
if isempty(config)
    % Use default FS config to get base paths
    config = get_user_config('AnimalRB0');
end

switch source
    case 'stimulus'
        base_dir = config.StimBasePath;
    case 'calcium'
        base_dir = config.DATA_2P;
    case 'camera'
        base_dir = config.CameraDir;
end

if ~exist(base_dir, 'dir')
    warning('get_experiment_list:DirNotFound', ...
        'Base directory not found: %s', base_dir);
    return;
end

% Scan for experiments
if strcmp(source, 'stimulus')
    % Stimulus folder structure: StimBasePath/AnimalName/RecordingName/
    experiments = scan_stimulus_folder(base_dir, animals_filter, recordings_filter);
else
    % DATA_2P and Camera structure: BaseDir/RecordingName/
    experiments = scan_flat_folder(base_dir, animals_filter, recordings_filter);
end

% Enrich experiments with additional info and filter
valid_idx = true(numel(experiments), 1);

if strict_source
    fprintf('Checking data availability for %d experiment(s)... (strict source mode)\n', numel(experiments));
else
    fprintf('Checking data availability for %d experiment(s)...\n', numel(experiments));
end

for i = 1:numel(experiments)
    exp = experiments(i);
    
    % Double-check animal filter (safety check)
    if ~isempty(animals_filter) && ~ismember(exp.animal_name, animals_filter)
        fprintf('  WARNING: Filtering out %s (animal %s not in filter)\n', exp.name, exp.animal_name);
        valid_idx(i) = false;
        continue;
    end

    % Get user config for this animal
    exp.config = get_user_config(exp.animal_name);

    % Check for calcium data (always in DATA_2P folder)
    calcium_path = fullfile(exp.config.DATA_2P, exp.name, 'suite2p');
    exp.has_calcium = exist(calcium_path, 'dir') == 7;

    % Check for stimulus data (always in StimBasePath folder)
    stim_path = fullfile(exp.config.StimBasePath, exp.animal_name, exp.name);
    stim_files = dir(fullfile(stim_path, 'stim*.mat'));
    exp.has_stimulus = ~isempty(stim_files);

    % Check for DLC data (always in DLCDir folder)
    dlc_file = fullfile(exp.config.DLCDir, [exp.name '_DLC_filtered.mat']);
    exp.has_dlc = exist(dlc_file, 'file') == 2;

    % Check if already processed
    preprocessed_file = fullfile(exp.config.OutputDir, [exp.name '_preprocessed.mat']);
    exp.is_processed = exist(preprocessed_file, 'file') == 2;

    experiments(i) = exp;

    % Apply filters
    if require_calcium && ~exp.has_calcium
        valid_idx(i) = false;
        continue;
    end

    if require_stimulus && ~exp.has_stimulus
        valid_idx(i) = false;
        continue;
    end

    if check_processed && exp.is_processed
        valid_idx(i) = false;
        continue;
    end
end

% Apply filter
experiments = experiments(valid_idx);

% Get unique animals in final list
unique_animals = unique({experiments.animal_name});

% Print summary
fprintf('Found %d experiment(s) to process\n', numel(experiments));
if ~isempty(experiments)
    fprintf('  Animals: %s\n', strjoin(unique_animals, ', '));
    fprintf('  With calcium: %d, Without calcium: %d\n', ...
        sum([experiments.has_calcium]), sum(~[experiments.has_calcium]));
end

% Final validation: ensure all experiments match animal filter
if ~isempty(animals_filter)
    for i = 1:numel(experiments)
        if ~ismember(experiments(i).animal_name, animals_filter)
            warning('get_experiment_list:FilterViolation', ...
                'Experiment %s from animal %s passed filter but should not have!', ...
                experiments(i).name, experiments(i).animal_name);
        end
    end
end

end


function experiments = scan_stimulus_folder(base_dir, animals_filter, recordings_filter)
%SCAN_STIMULUS_FOLDER Scan stimulus folder structure (BaseDir/Animal/Recording/)

experiments = struct('name', {}, 'animal_name', {}, 'path', {}, ...
    'has_calcium', {}, 'has_stimulus', {}, 'has_dlc', {}, ...
    'is_processed', {}, 'config', {});

% List animal folders
animal_dirs = dir(base_dir);
animal_dirs = animal_dirs([animal_dirs.isdir]);
animal_dirs = animal_dirs(~ismember({animal_dirs.name}, {'.', '..'}));

fprintf('Scanning stimulus folder: %s\n', base_dir);
fprintf('Found %d animal folder(s)\n', numel(animal_dirs));

for a = 1:numel(animal_dirs)
    animal_name = animal_dirs(a).name;

    % Apply animal filter
    if ~isempty(animals_filter) && ~ismember(animal_name, animals_filter)
        fprintf('  Skipping %s (not in filter)\n', animal_name);
        continue;
    end
    
    fprintf('  Processing %s...\n', animal_name);

    % List recording folders
    animal_path = fullfile(base_dir, animal_name);
    rec_dirs = dir(animal_path);
    rec_dirs = rec_dirs([rec_dirs.isdir]);
    rec_dirs = rec_dirs(~ismember({rec_dirs.name}, {'.', '..'}));
    
    fprintf('    Found %d recording(s)\n', numel(rec_dirs));

    for r = 1:numel(rec_dirs)
        rec_name = rec_dirs(r).name;

        % Apply recording filter
        if ~isempty(recordings_filter) && ~ismember(rec_name, recordings_filter)
            continue;
        end

        % Create experiment entry
        exp = struct();
        exp.name = rec_name;
        exp.animal_name = animal_name;
        exp.path = fullfile(animal_path, rec_name);
        exp.has_calcium = false;
        exp.has_stimulus = true;  % We found it in stimulus folder
        exp.has_dlc = false;
        exp.is_processed = false;
        exp.config = struct();

        experiments(end+1) = exp; %#ok<AGROW>
    end
end

fprintf('\n');

end


function experiments = scan_flat_folder(base_dir, animals_filter, recordings_filter)
%SCAN_FLAT_FOLDER Scan flat folder structure (BaseDir/Recording/)

experiments = struct('name', {}, 'animal_name', {}, 'path', {}, ...
    'has_calcium', {}, 'has_stimulus', {}, 'has_dlc', {}, ...
    'is_processed', {}, 'config', {});

% List recording folders
fprintf('Scanning flat folder: %s\n', base_dir);
rec_dirs = dir(base_dir);
rec_dirs = rec_dirs([rec_dirs.isdir]);
rec_dirs = rec_dirs(~ismember({rec_dirs.name}, {'.', '..'}));
fprintf('Found %d recording folder(s)\n', numel(rec_dirs));

for r = 1:numel(rec_dirs)
    rec_name = rec_dirs(r).name;

    % Apply recording filter
    if ~isempty(recordings_filter) && ~ismember(rec_name, recordings_filter)
        continue;
    end

    % Extract animal name from recording name (first part before underscore)
    name_parts = strsplit(rec_name, '_');
    animal_name = name_parts{1};

    % Apply animal filter
    if ~isempty(animals_filter) && ~ismember(animal_name, animals_filter)
        fprintf('  Skipping %s (animal %s not in filter)\n', rec_name, animal_name);
        continue;
    end
    
    fprintf('  Found %s (animal: %s)\n', rec_name, animal_name);

    % Create experiment entry
    exp = struct();
    exp.name = rec_name;
    exp.animal_name = animal_name;
    exp.path = fullfile(base_dir, rec_name);
    exp.has_calcium = false;
    exp.has_stimulus = false;
    exp.has_dlc = false;
    exp.is_processed = false;
    exp.config = struct();

    experiments(end+1) = exp; %#ok<AGROW>
end

fprintf('\n');

end


function exp = create_experiment_entry(rec_name, config)
%CREATE_EXPERIMENT_ENTRY Create experiment entry for a specific recording name

% Extract animal name from recording name
name_parts = strsplit(rec_name, '_');
animal_name = name_parts{1};

% Get config for this animal
if isempty(config)
    exp_config = get_user_config(animal_name);
else
    exp_config = config;
end

% Create experiment entry
exp = struct();
exp.name = rec_name;
exp.animal_name = animal_name;
exp.path = '';  % Will be determined based on what data exists
exp.has_calcium = false;
exp.has_stimulus = false;
exp.has_dlc = false;
exp.is_processed = false;
exp.config = exp_config;

% Check for calcium data path
calcium_path = fullfile(exp_config.DATA_2P, rec_name);
if exist(calcium_path, 'dir')
    exp.path = calcium_path;
    exp.has_calcium = exist(fullfile(calcium_path, 'suite2p'), 'dir') == 7;
end

% Check for stimulus data
stim_path = fullfile(exp_config.StimBasePath, animal_name, rec_name);
stim_files = dir(fullfile(stim_path, 'stim*.mat'));
exp.has_stimulus = ~isempty(stim_files);
if isempty(exp.path) && exp.has_stimulus
    exp.path = stim_path;
end

% Check for DLC data
dlc_file = fullfile(exp_config.DLCDir, [rec_name '_DLC_filtered.mat']);
exp.has_dlc = exist(dlc_file, 'file') == 2;

% Check if already processed
preprocessed_file = fullfile(exp_config.OutputDir, [rec_name '_preprocessed.mat']);
exp.is_processed = exist(preprocessed_file, 'file') == 2;

end
