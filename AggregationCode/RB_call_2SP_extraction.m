% Setup paths
addpath('Z:\Group Members\Rima\Stimulus'); 
addpath('C:\Users\rbondare\2P_scripts\MouseBrainActivity-main\2p Pre-processing');

maindir = 'Z:\Group Members\Rima\TEST\'; %DATA_2P
ca_type = 1;

% Get all experiment folders
explist = dir(maindir);
explist = explist(cat(1,explist.isdir));
explist = explist(~ismember({explist.name},{'. ','..'}));  % FIXED directory filtering

if isempty(explist)
    error('No experiments found in %s', maindir);
end

fprintf('Found %d recording(s) to process\n', numel(explist));

AnimalFolder = 'Z:\Group Members\Rima\Animals\';
parallel. gpu.enableCUDAForwardCompatibility(true);
params = FS_load_default_settings();

params. Neurons.use_deconvolved = true;
params.Neurons. activity_type = 'both';

%% Process ALL recordings
for e = 1:numel(explist)  
    
    fprintf('\n=== [%d/%d] Processing:  %s ===\n', e, numel(explist), explist(e).name);
    
    recording_path = [explist(e).folder '\' explist(e).name '\'];
    
    % Check if suite2p folder exists
    if ~exist([recording_path 'suite2p'], 'dir')
        fprintf('No suite2p folder, skipping\n');
        continue;
    end
    
    % Check if already processed
    preprocessed_file = ['Z:\Group Members\Rima\Aggregated\' explist(e).name '_preprocessed.mat'];
    if exist(preprocessed_file, 'file')
        fprintf('Already processed, skipping\n');
        continue;
    end
    
    % Process this recording
    try
        OutputFilename = RB_save_S2Pout_and_headers(recording_path, ca_type, 'overwrite_intermediate', false);
        fprintf('Created:  %s\n', OutputFilename);
    catch ME
        fprintf('Error:  %s\n', ME.message);
        fprintf('at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
    end
end

% import DLC data for all processed recordings
fprintf('Starting DLC import');
RB_import_VideoBehav();

fprintf('ALL DONE');