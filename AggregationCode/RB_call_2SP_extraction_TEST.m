clear;

% Setup paths
addpath('Z:\Group Members\Rima\Stimulus'); 
addpath('C:\Users\rbondare\2P_scripts\MouseBrainActivity-main\2p Pre-processing');

maindir = 'Z:\Group Members\Rima\TEST\';
ca_type = 1;

explist = dir(maindir);
explist = explist(cat(1,explist.isdir));
explist = explist(~ismember({explist.name},{'.','..'}));  % FIXED

AnimalFolder = 'Z:\Group Members\Rima\Animals\';
e = 1;
parallel. gpu. enableCUDAForwardCompatibility(true);
params = FS_load_default_settings();

params. Neurons. use_deconvolved = true;
params.Neurons. activity_type = 'both';

%% STEP 1: Process Suite2P (creates preprocessed file)
fprintf('\n=== STEP 1: Processing Suite2P ===\n');
fprintf('Processing %s\n', explist(e).name);
OutputFilename = RB_save_S2Pout_and_headers([explist(e).folder '\' explist(e).name '\'], ca_type, 'overwrite_intermediate', false);
fprintf('✓ Created:  %s\n', OutputFilename);

%% STEP 2: Import DLC (adds to preprocessed file)
fprintf('\n=== STEP 2: Importing DLC ===\n');
RB_import_VideoBehav_TEST();  % Use test version

fprintf('\n=== ALL DONE ===\n');