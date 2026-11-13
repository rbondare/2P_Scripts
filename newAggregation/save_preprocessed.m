function OutputFilename = save_preprocessed(recording_path, varargin)
%SAVE_PREPROCESSED Universal aggregation function for 2P preprocessing
%
%   OutputFilename = save_preprocessed(recording_path)
%   OutputFilename = save_preprocessed(recording_path, 'ca_type', 1)
%   OutputFilename = save_preprocessed(recording_path, 'ca_format', 'none')
%
%   Universal aggregation function that combines calcium imaging data,
%   behavioral data, and stimulus information into a single preprocessed file.
%
%   Key improvements over RB_save_S2Pout_and_headers.m:
%   - Optional calcium processing (behavior-only mode with ca_format='none')
%   - Config-based paths (no hardcoded user paths)
%   - Multiple calcium data formats (suite2p, caiman, auto-detect)
%   - Bug fixes from RB version (size(FN) -> size(FN,1), datetime format)
%
%   Input:
%       recording_path - Path to recording folder (containing suite2p/, tif files, etc.)
%                        Can also be just the recording name if config is provided
%
%   Optional Name-Value Arguments:
%       'ca_type' - Calcium data type index (default: 1)
%           1 = Raw, 2 = di_masks, 3 = di
%
%       'ca_format' - Calcium data format (default: 'auto')
%           'auto'    - Auto-detect from files present
%           'suite2p' - Load from suite2p output
%           'caiman'  - Load from Caiman HDF5 files
%           'none'    - Skip calcium processing (behavior-only)
%
%       'process_calcium' - Whether to process calcium data (default: true)
%           Set to false for behavior-only processing
%
%       'process_stimulus' - Whether to process stimulus data (default: true)
%
%       'process_behavior' - Whether to process DLC behavioral data (default: true)
%
%       'config' - User configuration from get_user_config (default: auto-detect)
%
%       'overwrite_intermediate' - Overwrite intermediate files (default: false)
%
%   Output:
%       OutputFilename - Path to the created preprocessed.mat file
%
%   Stimuli(k).TimeStimulusFrame is the time (s, session clock) of every stimulus frame of block k, taken from
%   the projector pulses of the red sync frames; Stimuli(k).RedFrameSynchronized says whether that worked
%   (false = stimulus PC clock, which is 0-35 ms off in RB22-24) and Stimuli(k).RedSyncInfo says how
%   (pulses found/lost, PC-clock offset, ...; fields described in make_sync_info).
%   A training stimulus holds many trials in one Stimuli entry: its TimeStimulusFrame is then frames x trials like
%   Frameinfo (NaN where there is no real frame), RedFrameSynchronized is true only if every trial was
%   synchronized, and RedFrameSynchronizedTrial, FrameValid, RedSyncInfo(trial) and RedSyncSession give the
%   detail (synchronize_red_frames_training). Its ball camera log is time stamped by training_ball_data
%   (see LocomotionCal.SyncInfo).
%
%   Example:
%       % Full processing with auto-detection
%       OutputFilename = save_preprocessed('D:\DATA_2P\AnimalRB5_20231201_001\');
%
%       % Behavior-only processing (no calcium data required)
%       OutputFilename = save_preprocessed('D:\DATA_2P\AnimalRB14_20231215_002\', ...
%           'ca_format', 'none');
%
%       % Process with specific config
%       config = get_user_config('AnimalRB5');
%       OutputFilename = save_preprocessed('AnimalRB5_20231201_001', 'config', config);
%
%   See also: get_user_config, get_experiment_list, load_calcium_data, run_aggregation

%% Parse inputs
p = inputParser;
addRequired(p, 'recording_path', @ischar);
addParameter(p, 'ca_type', 1, @(x) isnumeric(x) && x >= 0 && x <= 3);
addParameter(p, 'ca_format', 'auto', @(x) ismember(x, {'auto', 'suite2p', 'caiman', 'none'}));
addParameter(p, 'process_calcium', true, @islogical);
addParameter(p, 'process_stimulus', true, @islogical);
addParameter(p, 'process_behavior', true, @islogical);
addParameter(p, 'config', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'overwrite_intermediate', false, @islogical);
addParameter(p, 'stim_file', '', @ischar);
addParameter(p, 'cal_idx', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
parse(p, recording_path, varargin{:});

ca_type = p.Results.ca_type;
ca_format = p.Results.ca_format;
process_calcium = p.Results.process_calcium;
process_stimulus = p.Results.process_stimulus;
process_behavior = p.Results.process_behavior;
config = p.Results.config;
overwrite_intermediate = p.Results.overwrite_intermediate;
stim_file = p.Results.stim_file;
cal_idx = p.Results.cal_idx;

% If ca_format is 'none', disable calcium processing
if strcmp(ca_format, 'none')
    process_calcium = false;
end

%% Setup paths and directories
% Ensure recording_path ends with filesep
basedir = recording_path;
if ~endsWith(basedir, filesep)
    basedir = [basedir filesep];
end

% Extract animal and recording names from path
pathsplit = strsplit(basedir, filesep);
pathsplit = pathsplit(~cellfun(@isempty, pathsplit));  % Remove empty entries
RecordingName = pathsplit{end};
if isempty(RecordingName)
    RecordingName = pathsplit{end-1};
end

namesplit = strsplit(RecordingName, '_');
AnimalName = namesplit{1};

% Get user config if not provided
if isempty(config)
    config = get_user_config(AnimalName);
end

% Store original config for looking up AnimalInfo in original location
original_config = get_user_config(AnimalName);

fprintf('=== Processing: %s ===\n', RecordingName);
fprintf('Animal: %s, Operator: %s\n', AnimalName, config.OperatorID);
fprintf('Processing: Calcium=%d, Stimulus=%d, Behavior=%d\n', ...
    process_calcium, process_stimulus, process_behavior);

%% Setup output paths
OutputDir = config.OutputDir;
OutputCopyDir = config.OutputCopyDir;
OutputFilename = fullfile(OutputDir, [RecordingName '_preprocessed.mat']);
AnimalFolder = fullfile(config.AnimalFolder, AnimalName, filesep);
StimPath = fullfile(config.StimBasePath, AnimalName, RecordingName, filesep);
DLCdir = config.DLCDir;
IntermediateDir = config.IntermediateDir;
IntermediateFilename = fullfile(IntermediateDir, [RecordingName '_S2Presult.mat']);

% Create directories if needed
if ~exist(OutputDir, 'dir'); mkdir(OutputDir); end
if ~exist(IntermediateDir, 'dir'); mkdir(IntermediateDir); end

%% Initialize output file
Data = matfile(OutputFilename, 'Writable', true);

%% Process headers and calcium data
if process_calcium
    [Headers, SI_Info, Fall, S2Pdirs, multicolor] = process_calcium_headers(...
        basedir, ca_format, IntermediateFilename, overwrite_intermediate);

    % Save SI_Info to output
    Data.SI_Info = SI_Info;

    % Handle animal info (pass original_config to check original location)
    [AnimalInfo, Data] = handle_animal_info(AnimalFolder, AnimalName, Headers, Data, original_config);

else
    % Behavior-only mode: create minimal structures
    fprintf('Behavior-only mode: skipping calcium/header processing\n');

    % Try to load stimulus files to get timing info
    Headers = struct();
    SI_Info = struct();
    SI_Info.Clockstart = datetime('now');  % Placeholder, will be updated from stimulus
    Fall = [];
    S2Pdirs = [];
    multicolor = false;

    % Create placeholder AnimalInfo
    AnimalInfo = struct();
    AnimalInfo.AnimalName = AnimalName;
    AnimalInfo.Dates.Imaging = datetime('now');
end

%% Process frame timestamps
if process_calcium && ~isempty(fieldnames(Headers))
    [FT, multi, dtA, dt] = process_frame_timestamps(Headers, SI_Info);
    Data.TimeCa = FT;
else
    FT = [];
    multi = false;
    dtA = [];
    dt = [];
end

%% Process calcium data
if process_calcium && ~isempty(Fall)
    [CaData, Img_max_proj] = process_calcium_signals(Fall, S2Pdirs, SI_Info, FT, multicolor, ca_type);
    Data.CaData = CaData;
    Data.Img_max = Img_max_proj;
else
    % Create empty CaData structure
    CaData = create_empty_cadata();
    Data.CaData = CaData;
    Data.dFFsettings = create_default_dff_settings();
end

%% Compute acquisition timestamps (needed for locomotion, before trigger processing strips Headers)
if process_calcium && ~isempty(fieldnames(Headers))
    acq_timestamps = nan(numel(Headers), 2);
    for h = 1:numel(Headers)
        acq_timestamps(h,:) = Headers(h).frameTimestamps_sec([1 end]);
        if multi; acq_timestamps(h,:) = acq_timestamps(h,:) + dtA(h); end
    end
else
    acq_timestamps = [];
end

%% Process triggers
if process_calcium && ~isempty(fieldnames(Headers))
    [Triggers, Headers, TimingProblems] = process_triggers(Headers, multi, dt, dtA);
    Data.Triggers = Triggers;
    Data.Headers = Headers;
else
    Triggers = struct();
    Triggers.TimeProjector = [];
    Triggers.TimeCamera = [];
    Triggers.TimeBall = {};
    Triggers.TimingProblems = false;
    Triggers.Stimulus_acq_timestamps = [];
    Data.Triggers = Triggers;
end

%% Process stimulus data
if process_stimulus
    [Stimulusoptions, BallCamData, LocomotionCal] = process_stimulus_data(...
        StimPath, Headers, SI_Info, Triggers, multi, dtA, acq_timestamps, stim_file, cal_idx);

    if ~isempty(Stimulusoptions)
        Data.Stimuli = Stimulusoptions;
    end
    if ~isempty(LocomotionCal)
        Data.LocomotionNonCal = BallCamData;
        Data.LocomotionCal = LocomotionCal;
    end
end

%% Process DLC behavioral data
if process_behavior
    dlc_file = fullfile(DLCdir, [RecordingName '_DLC_filtered.mat']);
    if exist(dlc_file, 'file')
        fprintf('Loading DLC behavioral data\n');
        load(dlc_file, 'DLCSupporting', 'Behav', 'VideoParams', 'VideoFrameTimes');

        % Dropped camera frames: if the video and the 2P camera trigger log differ by more than one frame, pairing
        % frame i with trigger i (below) times every frame after a drop too early and cuts off the trigger tail.
        % The camera's hardware time stamps (frames.csv) then tell which trigger each saved frame belongs to.
        ticks_used = false;
        if isfield(Triggers, 'TimeCamera') && ~isempty(Triggers.TimeCamera) && ...
                abs(numel(Triggers.TimeCamera) - numel(Behav.FrameNum)) > 1
            [ticks_used, Triggers, selected_frameIdx] = match_video_frames_by_ticks(...
                Triggers, Behav, config, AnimalName, RecordingName);
            if ticks_used
                Data.Triggers = Triggers;
                VideoParams.selected_frameIdx = selected_frameIdx;
            end
        end

        % Correct video frame timing
        if ticks_used
            % selected_frameIdx and Triggers.TimeCamera were set from the camera time stamps above
        elseif isfield(VideoParams, 'selected_frameIdx')
            selected_frameIdx = VideoParams.selected_frameIdx;
        else
            if exist('Triggers', 'var') && isfield(Triggers, 'TimeCamera')
                [selected_frameIdx, Triggers] = correct_videoFrames_OS(Triggers, VideoFrameTimes, Behav);
                Data.Triggers = Triggers;
                VideoParams.selected_frameIdx = selected_frameIdx;
                save(dlc_file, 'VideoParams', '-append');
            else
                selected_frameIdx = 1:size(Behav.VideoNum, 1);
            end
        end

        % Apply frame selection to DLC data
        DLCSupporting = apply_frame_selection(DLCSupporting, selected_frameIdx);
        Behav = apply_frame_selection_behav(Behav, selected_frameIdx);

        Data.DLCSupporting = DLCSupporting;
        Data.Behav = Behav;
        Data.VideoParams = VideoParams;
    else
        fprintf('DLC file not found: %s\n', dlc_file);
    end
end

%% Copy to backup location
if ~isempty(OutputCopyDir)
    if ~exist(OutputCopyDir, 'dir'); mkdir(OutputCopyDir); end
    copyfile(OutputFilename, fullfile(OutputCopyDir, [RecordingName '_preprocessed.mat']));
    fprintf('Copied to backup: %s\n', OutputCopyDir);
end

fprintf('=== Completed: %s ===\n', OutputFilename);

end

%% ===== HELPER FUNCTIONS =====

function [Headers, SI_Info, Fall, S2Pdirs, multicolor] = process_calcium_headers(...
    basedir, ca_format, IntermediateFilename, overwrite_intermediate)
%PROCESS_CALCIUM_HEADERS Load and process calcium data and headers

suite2Pdir = fullfile(basedir, 'suite2p', filesep);
S2Pdirs = dir([suite2Pdir 'plane*']);

% Sort plane directories
if ~isempty(S2Pdirs)
    temp = struct2cell(S2Pdirs);
    temp = temp(1,:)';
    temp = cellfun(@(x) str2double(x(6:end)), temp);
    [~, sort_idx] = sort(temp);
    S2Pdirs = S2Pdirs(sort_idx);
end

% Find and process header files
files = dir(fullfile(basedir, '*header.mat'));
files_tif = dir(fullfile(basedir, '*.tif'));

if ~isempty(files_tif)
    files = dir(fullfile(basedir, '*.tif'));
end

if isempty(files)
    error('save_preprocessed:NoHeaders', 'No header or TIF files found in: %s', basedir);
end

% Load or create headers
if ~exist(IntermediateFilename, 'file') || overwrite_intermediate
    Headers = load_headers_from_files(basedir, files, files_tif);

    % Get SI_Info
    k1 = 1;
    if isempty(Headers(k1).SI); k1 = 2; end
    SI_Info = Headers(k1).SI;
    SI_Info.Clockstart = datetime(Headers(k1).epoch{1});

    % Determine multicolor (channels saved to the TIFFs, not channelsActive - see is_multicolor.m)
    multicolor = is_multicolor(SI_Info);

    % Load calcium data
    [Fall, S2Pdirs, ~] = load_calcium_data(basedir, 'format', ca_format, 'SI_Info', SI_Info);

    % Save intermediate
    fprintf('Saving intermediate file\n');
    save(IntermediateFilename, 'Fall', 'Headers', '-v7.3');
else
    fprintf('Loading from intermediate file\n');
    load(IntermediateFilename, 'Fall', 'Headers');

    k1 = 1;
    if isempty(fieldnames(Headers(k1).ImgInfo)); k1 = 2; end
    SI_Info = Headers(k1).SI;
    SI_Info.Clockstart = datetime(Headers(k1).epoch{1});

    multicolor = is_multicolor(SI_Info);
end

% Fix numFramesPerVolume if missing
if ~isfield(SI_Info.hFastZ, 'numFramesPerVolume')
    if isfield(SI_Info.hStackManager, 'numFramesPerVolumeWithFlyback')
        SI_Info.hFastZ.numFramesPerVolume = SI_Info.hStackManager.numFramesPerVolumeWithFlyback;
    else
        error('save_preprocessed:PlaneNumNotFound', 'Plane number not found in SI_Info');
    end
end

end


function Headers = load_headers_from_files(basedir, files, files_tif)
%LOAD_HEADERS_FROM_FILES Load header information from TIF or header files

numFiles = length(files);
Headers = struct();

w = waitbar(0, 'Loading headers');
for File = 1:numFiles
    fullname = files(File).name;
    [~, file_name, ext] = fileparts(fullname);

    if strcmp(ext, '.tif')
        header_file = fullfile(basedir, [file_name '_header.mat']);
        if ~exist(header_file, 'file')
            [header, ~, ImgInfo] = fast_get_scanimage_header_notMROI(fullfile(basedir, fullname), 'slice', 1, 'volume', 1);
            save(header_file, 'header', 'ImgInfo', '-v7.3');
        else
            load(header_file, 'header', 'ImgInfo');
        end
    else
        load(fullfile(basedir, [file_name '.mat']), 'header', 'ImgInfo');
    end

    Headers(File).ImgInfo = ImgInfo;
    if ~isempty(header)
        FN = fieldnames(header);
        % BUG FIX: Use size(FN,1) instead of size(FN) which returns [rows,cols]
        for n = 1:size(FN, 1)
            Headers(File).(FN{n}) = header.(FN{n});
        end
    end
    waitbar(File/numFiles, w);
end
close(w);

% A TIFF whose ScanImage frame description was cut off (first frame of a session with many pre-start aux triggers)
% has no parsed header; rebuild it, otherwise its frame is missing and calcium / video / ball are shifted
tif_paths = cell(1, numFiles);
for File = 1:numFiles
    if strcmp(files(File).name(max(end-3, 1):end), '.tif')
        tif_paths{File} = fullfile(basedir, files(File).name);
    end
end
Headers = recover_truncated_headers(Headers, tif_paths);

% Clean up headers for last file if needed
k1 = 1;
if isempty(Headers(k1).SI); k1 = 2; end
SI_Info = Headers(k1).SI;

if ~isfield(Headers(end).ImgInfo, 'numVolumes') && Headers(end).SI.hStackManager.numSlices > 1
    numS = Headers(end-1).ImgInfo.numSlices + SI_Info.hFastZ.numDiscardFlybackFrames;
    FN = fieldnames(Headers);
    last = numel(Headers(end).frameNumbers) - rem(numel(Headers(end).frameNumbers), numS);
    for n = 3:numel(FN)
        Headers(end).(FN{n}) = Headers(end).(FN{n})(1:last);
    end
end

end


function [AnimalInfo, Data] = handle_animal_info(AnimalFolder, AnimalName, Headers, Data, original_config)
%HANDLE_ANIMAL_INFO Load or create animal information
%   Also checks original user folder if local folder doesn't have the file

k1 = 1;
if isempty(Headers(k1).SI); k1 = 2; end

animal_file = fullfile(AnimalFolder, [AnimalName '.mat']);
animal_info_found = false;

% First check local/test folder
if exist(animal_file, 'file')
    load(animal_file, 'AnimalInfo');
    animal_info_found = true;
    fprintf('Loaded AnimalInfo from: %s\n', animal_file);
end

% If not found and we have original config, check original location
if ~animal_info_found && nargin >= 5 && ~isempty(original_config)
    original_animal_folder = fullfile(original_config.AnimalFolder, AnimalName, filesep);
    original_animal_file = fullfile(original_animal_folder, [AnimalName '.mat']);
    if exist(original_animal_file, 'file')
        load(original_animal_file, 'AnimalInfo');
        animal_info_found = true;
        fprintf('Loaded AnimalInfo from original location: %s\n', original_animal_file);
    end
end

% If still not found, create default (without dialog for batch mode compatibility)
if ~animal_info_found
    % Check if we're in batch/non-interactive mode
    is_batch_mode = batchStartupOptionUsed || ~usejava('desktop');

    if is_batch_mode
        % Create default AnimalInfo without dialog
        fprintf('Creating default AnimalInfo (batch mode)\n');
        AnimalInfo = struct();
        AnimalInfo.AnimalName = AnimalName;
        AnimalInfo.Genotype = 'Unknown';
        AnimalInfo.PCF_AnimalNumber = 'xx/xxx';
        AnimalInfo.Sex = 'unknown';
        AnimalInfo.Notes = 'Auto-created in batch mode';
        AnimalInfo.Dates.Birth = datetime('now') - years(1);
        AnimalInfo.Dates.Injection = datetime('now') - days(30);
        AnimalInfo.Dates.Implantation = datetime('now') - days(14);
        AnimalInfo.Indicator.Type = 'GCamP6f';
        AnimalInfo.Indicator.Expression = '';
        AnimalInfo.Dates.Imaging = [];
    else
        % Interactive mode - show dialog
        N = {'AnimalName'; 'Genotype'; 'PCF_AnimalNumber'; 'Sex'; 'DateOfBirth'; ...
            'DateOfInjection'; 'DateOfImplantation'; 'Indicator'; 'Expression'; 'Notes'};
        def = {AnimalName; 'Ai148/Rorb-IRES2-cre ki'; 'xx/xxx'; 'm'; '20190101'; ...
            '20190101'; '20190101'; 'GCamP6f'; ''; ''};
        answerE = inputdlg(N, 'Add animal', 1, def);

        infmt = 'yyyyMMdd';
        AnimalInfo = struct();
        AnimalInfo.AnimalName = answerE{1};
        AnimalInfo.Genotype = answerE{2};
        AnimalInfo.PCF_AnimalNumber = answerE{3};
        AnimalInfo.Sex = answerE{4};
        AnimalInfo.Notes = answerE{10};
        AnimalInfo.Dates.Birth = datetime(answerE{5}, 'InputFormat', infmt);
        AnimalInfo.Dates.Injection = datetime(answerE{6}, 'InputFormat', infmt);
        AnimalInfo.Dates.Implantation = datetime(answerE{7}, 'InputFormat', infmt);
        AnimalInfo.Indicator.Type = answerE{8};
        AnimalInfo.Indicator.Expression = answerE{9};
        AnimalInfo.Dates.Imaging = [];
    end
end

% Add imaging date (not again if this session is already listed, e.g. when it is re-aggregated)
imaging_date = datetime(Headers(k1).epoch{1});
if isempty(AnimalInfo.Dates.Imaging) || ~any(abs(AnimalInfo.Dates.Imaging - imaging_date) < seconds(1))
    AnimalInfo.Dates.Imaging = [AnimalInfo.Dates.Imaging; imaging_date];
end

% Save animal info to local folder
if ~exist(AnimalFolder, 'dir'); mkdir(AnimalFolder); end
save(animal_file, 'AnimalInfo');
Data.AnimalInfo = AnimalInfo;

end


function [FT, multi, dtA, dt] = process_frame_timestamps(Headers, SI_Info)
%PROCESS_FRAME_TIMESTAMPS Process and validate frame timestamps

k1 = 1;
if isempty(Headers(k1).SI); k1 = 2; end

FT = cat(2, Headers(:).frameTimestamps_sec);

if any(diff(FT) < 0)
    % Handle non-monotonic timestamps (multi-acquisition)
    dt = cellfun(@(x) datetime(x), cat(2, Headers(:).epoch), 'UniformOutput', false);
    dt = cat(2, dt{:});
    dt = seconds(dt - dt(1));

    dtA = zeros(size(Headers));
    stE = datetime(Headers(k1).epoch{1});
    for h = 1:numel(Headers)
        dtA(h) = seconds(datetime(Headers(h).epoch{1}) - stE);
    end

    tempf = figure('Name', 'Timestamp Correction');
    plot(FT, 'k'); hold on;
    FT = FT + dt;
    plot(FT, 'r'); hold off;
    legend('timestamps raw', 'timestamps corrected');

    multi = any(input('Frame Timestamps not monotonically increasing, multi-acquisition? 1 for yes/0 for cancel: '));
    close(tempf);

    if ~multi
        error('save_preprocessed:TimestampCancel', 'Cancelling due to timestamp issues');
    end
else
    multi = false;
    dt = [];
    dtA = [];
end

FT = unique(FT, 'stable');

if isfield(Headers(k1).ImgInfo, 'numSlices')
    FT = reshape(FT, Headers(k1).ImgInfo.numSlices + SI_Info.hFastZ.numDiscardFlybackFrames, []);
    if SI_Info.hFastZ.numDiscardFlybackFrames > 0
        FT = FT(1:end-SI_Info.hFastZ.numDiscardFlybackFrames, :);
    end
end

end


function [CaData, Img_max_proj] = process_calcium_signals(Fall, S2Pdirs, SI_Info, FT, multicolor, ca_type)
%PROCESS_CALCIUM_SIGNALS Process calcium fluorescence signals

proj = exist([S2Pdirs(1).folder filesep 'proj'], 'dir') == 7;

if proj
    FT = FT(floor(size(FT, 1)/2), :);
end

if ~proj
    num_planes = SI_Info.hFastZ.numFramesPerVolume - SI_Info.hFastZ.numDiscardFlybackFrames;
    Img_max_proj = cell(num_planes, 1);

    for f = 1:min(numel(S2Pdirs), num_planes)
        Ncell = find(Fall(f).iscell(:, 1));
        S2Presult(f).deconvolved = Fall(f).spks(logical(Fall(f).iscell(:, 1)), 1:size(FT, 2));
        S2Presult(f).F = Fall(f).F(logical(Fall(f).iscell(:, 1)), 1:size(FT, 2));
        S2Presult(f).Fneu = Fall(f).Fneu(logical(Fall(f).iscell(:, 1)), 1:size(FT, 2));
        S2Presult(f).center = nan(numel(Ncell), 2);
        S2Presult(f).cellstats = reshape(Fall(f).stat(Ncell), [], 1);

        for n = 1:numel(Ncell)
            if iscell(Fall(f).stat)
                S2Presult(f).center(n, :) = Fall(f).stat{Ncell(n)}.med;
            else
                S2Presult(f).center(n, :) = Fall(f).stat(Ncell(n)).med;
            end
        end

        if multicolor
            S2Presult(f).redcell = Fall(f).redcell(Ncell, :);
        end

        Img_max_proj{f, 1} = Fall(f).ops.max_proj;
        if multicolor
            Img_max_proj{f, 2} = single(Fall(f).ops.meanImg_chan2);
        end
    end
else
    Img_max_proj = {Fall.ops.max_proj};
    Ncell = find(Fall.iscell(:, 1));
    S2Presult.deconvolved = Fall.spks(logical(Fall.iscell(:, 1)), 1:size(FT, 2));
    S2Presult.F = Fall.F(logical(Fall.iscell(:, 1)), 1:size(FT, 2));
    S2Presult.Fneu = Fall.Fneu(logical(Fall.iscell(:, 1)), 1:size(FT, 2));
    S2Presult.center = nan(numel(Ncell), 2);
    S2Presult.cellstats = reshape(Fall.stat(Ncell), [], 1);

    if multicolor
        S2Presult.redcell = Fall.redcell(Ncell, :);
    end

    for n = 1:numel(Ncell)
        S2Presult.center(n, :) = Fall.stat{Ncell(n)}.med;
    end
end

% Add plane index to center coordinates
for n = 1:numel(S2Presult)
    S2Presult(n).center(:, end+1) = n;
end

% Calculate dF/F
fprintf('Calculating dF/F\n');
dFFsettings = create_default_dff_settings();

Ca_F_neuropil = cat(1, S2Presult(:).Fneu);
Ca_F = cat(1, S2Presult(:).F);
Ca_deconvolved = cat(1, S2Presult(:).deconvolved);
Ca_centroid_voxel = cat(1, S2Presult(:).center);
Ca_cellstats = cat(1, S2Presult(:).cellstats);

% Convert cellstats from cell to struct array
% First union all field names so structs from different planes are compatible
all_stat_fields = fieldnames(Ca_cellstats{1});
for n = 2:numel(Ca_cellstats)
    all_stat_fields = union(all_stat_fields, fieldnames(Ca_cellstats{n}));
end
for n = 1:numel(Ca_cellstats)
    for fi = 1:numel(all_stat_fields)
        fn = all_stat_fields{fi};
        if ~isfield(Ca_cellstats{n}, fn)
            Ca_cellstats{n}.(fn) = [];
        end
    end
end
Ca_cellstats2 = Ca_cellstats{1};
for n = 1:numel(Ca_cellstats)
    Ca_cellstats2(n) = Ca_cellstats{n};
end
Ca_cellstats = Ca_cellstats2;

if multicolor
    Red_cell_prob = cat(1, S2Presult(:).redcell);
    for n = 1:numel(Ca_cellstats)
        Ca_cellstats(n).redprob = Red_cell_prob(n, end);
    end
end

[Ca_dFF, Ca_baseline] = estimate_continuous_dff(FT, Ca_F, Ca_F_neuropil, ...
    dFFsettings.neuropil_factor, dFFsettings.window, ...
    dFFsettings.medianmode, dFFsettings.baselinemode, dFFsettings.prclevel);

% Create CaData structure
CaData = create_empty_cadata();

if ca_type > 0
    CaData(ca_type).Ca_dFF = Ca_dFF;
    CaData(ca_type).Ca_deconvolved = Ca_deconvolved;
    CaData(ca_type).Ca_centroid_voxel = Ca_centroid_voxel;
    CaData(ca_type).Ca_F = Ca_F;
    CaData(ca_type).Ca_F_neuropil = Ca_F_neuropil;
    CaData(ca_type).Ca_baseline = Ca_baseline;
    CaData(ca_type).Ca_cellstats = Ca_cellstats;
end

end


function CaData = create_empty_cadata()
%CREATE_EMPTY_CADATA Create empty CaData structure

for fe = 1:3
    CaData(fe).Ca_dFF = [];
    CaData(fe).Ca_deconvolved = [];
    CaData(fe).Ca_centroid_voxel = [];
    CaData(fe).Ca_F = [];
    CaData(fe).Ca_F_neuropil = [];
    CaData(fe).Ca_baseline = [];
    CaData(fe).Ca_cellstats = [];
end

end


function dFFsettings = create_default_dff_settings()
%CREATE_DEFAULT_DFF_SETTINGS Create default dF/F calculation settings

dFFsettings.window = 15;
dFFsettings.neuropil_factor = 0.5;
dFFsettings.medianmode = 'windowed';
dFFsettings.baselinemode = 'prctile';
dFFsettings.prclevel = 8;

end


function [Triggers, Headers, TimingProblems] = process_triggers(Headers, multi, dt, dtA)
%PROCESS_TRIGGERS Extract and process timing triggers from headers

TimingProblems = false;
k1 = 1;
if isempty(Headers(k1).SI); k1 = 2; end

% Remove unnecessary fields from Headers
Headers = rmfield(Headers, {'SI'; 'acqTriggerTimestamps_sec'; 'acquisitionNumbers'; ...
    'frameNumberAcquisition'; 'nextFileMarkerTimestamps_sec'; 'endOfAcquisition'; ...
    'endOfAcquisitionMode'; 'dcOverVoltage'; 'I2CData'; 'epoch'});

fprintf('Processing timing triggers\n');

FrameNum = cat(2, Headers(:).frameNumbers);

if multi
    newAcqIdx = [find(diff(FrameNum) < 0) + 1, numel(FrameNum) + 1];
    for m = 1:numel(newAcqIdx) - 1
        FrameNum(newAcqIdx(m):newAcqIdx(m+1)-1) = ...
            FrameNum(newAcqIdx(m):newAcqIdx(m+1)-1) + FrameNum(newAcqIdx(m)-1);
    end
end

[~, ia, ~] = unique(FrameNum, 'stable');

% Projector triggers
PT = cat(2, Headers(:).auxTrigger0);
if multi
    for m = 1:numel(PT)
        PT{m} = PT{m} + dt(m);
    end
end
PT = cat(2, PT{ia});
dPT = [inf diff(PT)];
if any(dPT < 0)
    warning('Projector Timestamps not monotonically increasing');
    TimingProblems = true;
end
PT(dPT < 1/60) = [];
Triggers.TimeProjector = PT;

% Camera triggers
CT = cat(2, Headers(:).auxTrigger1);
if multi
    for m = 1:numel(CT)
        CT{m} = CT{m} + dt(m);
    end
end
CT = cat(2, CT{ia(ia < numel(CT))});
dCT = [inf diff(CT)];
if any(dCT < 0)
    warning('Camera Timestamps not monotonically increasing');
    TimingProblems = true;
end
Triggers.TimeCamera = CT;

% Ball triggers
BT = cell(numel(Headers), 1);
for b = 1:numel(BT)
    if ~isempty(Headers(b).auxTrigger2)
        [~, temp_fn, ~] = unique(Headers(b).frameNumbers);
        BT{b} = cat(2, Headers(b).auxTrigger2{temp_fn});
        if multi
            BT{b} = BT{b} + dtA(b);
        end
    end
end
Triggers.TimeBall = BT;
Triggers.TimingProblems = TimingProblems;

% Extract stimulus acquisition timestamps before removing frameTimestamps_sec
Stimulus_acq_timestamps = nan(numel(Headers), 2);
for n = 1:numel(Headers)
    if isfield(Headers(n), 'frameTimestamps_sec') && ~isempty(Headers(n).frameTimestamps_sec)
        Stimulus_acq_timestamps(n,:) = Headers(n).frameTimestamps_sec([1 end]);
        if multi && ~isempty(dtA) && n <= numel(dtA)
            Stimulus_acq_timestamps(n,:) = Stimulus_acq_timestamps(n,:) + dtA(n);
        end
    end
end
Triggers.Stimulus_acq_timestamps = Stimulus_acq_timestamps;

% Remove trigger fields if no problems
if ~TimingProblems
    Headers = rmfield(Headers, {'auxTrigger0'; 'auxTrigger1'; 'auxTrigger2'; 'auxTrigger3'});
end

Headers = rmfield(Headers, {'frameNumbers'; 'frameTimestamps_sec'});

end


function [used, Triggers, selected_frameIdx] = match_video_frames_by_ticks(Triggers, Behav, config, AnimalName, RecordingName)
%MATCH_VIDEO_FRAMES_BY_TICKS Time the saved video frames with the camera hardware time stamps (frames.csv)
%   used = false (Triggers unchanged, a warning says why) when there is no usable frames.csv or the match fails
%   its checks; then the caller falls back to pairing frame i with trigger i.
%   used = true: Triggers.TimeCamera holds the trigger time of every video frame (frame k <-> element k) and
%   Triggers.CameraMatch documents the match (dropped frames, residuals); see match_camera_ticks_to_triggers.

used = false;
selected_frameIdx = [];
nB = numel(Behav.FrameNum);
nT = numel(Triggers.TimeCamera);
csvf = find_camera_frames_csv(config, AnimalName, RecordingName);
if isempty(csvf)
    warning('save_preprocessed:CameraFrameCount', ...
        ['The video has %d frames but the 2P logged %d camera triggers, and no frames.csv with the camera time stamps ' ...
        'was found for %s (Camera or preCamera folder): frame i is paired with trigger i, so frames after a dropped ' ...
        'frame are timed too early.'], nB, nT, RecordingName);
    return
end
Tb = readtable(csvf);
if height(Tb) ~= nB || ~ismember('camera_timestamp_ticks', Tb.Properties.VariableNames)
    warning('save_preprocessed:CameraFrameCount', ...
        '%s has %d rows (video: %d frames) or no camera_timestamp_ticks column: not used for the frame timing.', ...
        csvf, height(Tb), nB);
    return
end
ticks = double(Tb.camera_timestamp_ticks - Tb.camera_timestamp_ticks(1)) * 1e-9;
trig = double(Triggers.TimeCamera(:));
[t, info] = match_camera_ticks_to_triggers(trig, ticks);
if ~info.ok
    warning('save_preprocessed:CameraTickMatch', ...
        ['Matching the video frames to the camera triggers with the time stamps of %s failed its checks (re-anchored %d x, ' ...
        '%d frames more than a quarter period from their trigger, %d frames after the last trigger): ' ...
        'falling back to pairing frame i with trigger i.'], csvf, info.n_reanchor, info.n_bad, info.n_untimed_tail);
    return
end

CM = struct();
CM.method = 'camera hardware time stamps (frames.csv) matched to the 2P camera triggers';
[~, nm, ext] = fileparts(csvf);
CM.frames_csv = [nm ext];
CM.n_frames = nB;
CM.n_triggers_raw = nT;
CM.n_dropped = info.n_dropped;
CM.dropped_time_s = trig(info.dropped_idx)';
CM.n_before_first = info.n_before_first;
CM.n_after_last = info.n_after_last;
CM.n_untimed_tail = info.n_untimed_tail;
CM.residual_sd_ms = info.residual_sd_ms;
CM.max_abs_residual_ms = info.max_abs_residual_ms;
CM.drift_ppm = info.drift_ppm;
CM.period_s = info.period_s;
Triggers.CameraMatch = CM;
Triggers.TimeCamera = t(:)';
selected_frameIdx = info.frame_idx(:)';
used = true;
fprintf(['  Camera frame drops: %d of %d camera triggers (%.2f %%) have no video frame; every video frame is matched to its ' ...
    'own trigger with the camera time stamps (residual SD %.3f ms, clock drift %.1f ppm)\n'], ...
    info.n_dropped, nT, 100 * info.n_dropped / nT, info.residual_sd_ms, info.drift_ppm);
if info.n_untimed_tail > 0 || info.n_before_first > 0 || info.n_after_last > 0
    fprintf('  Camera frames: %d video frames after the last trigger left out; %d triggers before the first / %d after the last video frame\n', ...
        info.n_untimed_tail, info.n_before_first, info.n_after_last);
end
end


function csvf = find_camera_frames_csv(config, AnimalName, RecordingName)
%FIND_CAMERA_FRAMES_CSV <recording>'s *_frames.csv (camera time stamps) in the Camera folder or its preCamera archive
csvf = '';
if ~isfield(config, 'CameraDir') || isempty(config.CameraDir)
    return
end
cam = regexprep(config.CameraDir, '[\\/]+$', '');
parent = fileparts(cam);
cand = {fullfile(cam, RecordingName), ...
        fullfile(parent, 'preCamera', AnimalName, RecordingName), ...
        fullfile(parent, 'preCamera', RecordingName)};
for k = 1:numel(cand)
    L = dir(fullfile(cand{k}, '*_frames.csv'));
    if numel(L) == 1
        csvf = fullfile(L(1).folder, L(1).name);
        return
    elseif numel(L) > 1
        warning('save_preprocessed:CameraFrameCount', '%d *_frames.csv files in %s: not used', numel(L), cand{k});
        return
    end
end
end


function [Stimulusoptions, BallCamData, LocomotionCal] = process_stimulus_data(...
    StimPath, Headers, SI_Info, Triggers, multi, dtA, acq_timestamps, stim_file, cal_idx)
%PROCESS_STIMULUS_DATA Load and process stimulus files

Stimulusoptions = [];
BallCamData = [];
LocomotionCal = [];

% Find stimulus files
if nargin >= 8 && ~isempty(stim_file)
    % Use specific stim file (e.g. individual Training session)
    StimFiles = dir(stim_file);
else
    StimFiles = dir(fullfile(StimPath, 'stim*.mat'));
    if isempty(StimFiles)
        % Try alternative path pattern
        [parent_path, rec_name, ~] = fileparts(StimPath(1:end-1));
        StimFiles = dir(fullfile(parent_path, [rec_name '*'], 'stim*.mat'));
    end
end

if isempty(StimFiles)
    fprintf('No stimulus files found in: %s\n', StimPath);
    return;
end

fprintf('Loading %d stimulus file(s)\n', size(StimFiles, 1));

% Load stimulus files
Stimulusoptions = [];
Stimuli = [];
for F = 1:size(StimFiles, 1)
    stim_data = load(fullfile(StimFiles(F).folder, StimFiles(F).name));
    if isfield(stim_data, 'options')
        if isempty(Stimulusoptions)
            Stimulusoptions = stim_data.options;
        else
            Stimulusoptions(F) = stim_data.options;
        end
    end
    if isempty(Stimuli)
        Stimuli = stim_data;
    else
        Stimuli(F) = stim_data;
    end
end

if isempty(Stimulusoptions)
    fprintf('No stimulus options found in files\n');
    return;
end

% Extract ball camera data
BallCamData_temp = cat(1, {Stimuli(:).BallCamData})';
Stimuli = rmfield(Stimuli, {'BallCamData', 'options'});

% Transfer fields from Stimuli to Stimulusoptions
stim_fields = fieldnames(Stimuli);
for i = 1:length(stim_fields)
    for j = 1:size(Stimuli, 2)
        Stimulusoptions(j).(stim_fields{i}) = Stimuli(j).(stim_fields{i});
    end
end

% Synchronize stimulus frames using projector triggers (red sync frames)
% Check stimulus type to determine synchronization method
if ~isempty(Stimulusoptions) && isfield(Stimulusoptions(1), 'type')
    PT = Triggers.TimeProjector;
    Stimulus_acquisition_start_stop_timestamp = Triggers.Stimulus_acq_timestamps;
    
    if ~strcmp(Stimulusoptions(1).type, 'training')
        % Non-training stimulus: synchronize using projector triggers
        fprintf('Synchronizing stimulus frames for non-training paradigm\n');
        % Map each stimulus block to its correct acquisition file via nearest
        % end-time match (same helper already used below for training stimuli).
        % Without this, Stimnumber was used directly as a row index into
        % Stimulus_acquisition_start_stop_timestamp, which has one row per TIF
        % acquisition file -- not one per stimulus block -- so any inter-block
        % dead time recorded as its own acquisition broke the mapping for every block.
        StimAcqNumber = match_stim_to_acq(Stimulusoptions, Stimulus_acquisition_start_stop_timestamp, SI_Info);
        for Stimnumber = 1:numel(Stimulusoptions)
            stim_sync = synchronize_red_frames(...
                Stimulusoptions(Stimnumber), PT, Stimulus_acquisition_start_stop_timestamp, ...
                StimAcqNumber(Stimnumber), SI_Info);
            
            % Copy synchronized fields back to stimulus struct
            Stimulusoptions(Stimnumber).TimeStimulusFrame = stim_sync.TimeStimulusFrame;
            Stimulusoptions(Stimnumber).RedFrameSynchronized = stim_sync.RedFrameSynchronized;
            Stimulusoptions(Stimnumber).TimeStimulusFrameVBL = stim_sync.TimeStimulusFrameVBL;
            Stimulusoptions(Stimnumber).RedSyncInfo = stim_sync.RedSyncInfo;
        end
        % Report every block and compare the blocks with each other (falls back to PC clock
        % times for a block whose pulses were numbered inconsistently with the others)
        Stimulusoptions = check_red_sync_consistency(Stimulusoptions, PT, ...
            Stimulus_acquisition_start_stop_timestamp, StimAcqNumber, SI_Info);
    else
        % Training stimulus: every trial (column of Frameinfo) is synchronized on its own and finds its
        % acquisition by overlap in time (match_stim_to_acq maps by end time, which fails for the trials
        % that share an acquisition with a break or neighbour a stub acquisition)
        fprintf('Synchronizing stimulus frames for training paradigm\n');
        for Stimnumber = 1:numel(Stimulusoptions)
            stim_sync = synchronize_red_frames_training(...
                Stimulusoptions(Stimnumber), PT, Stimulus_acquisition_start_stop_timestamp, SI_Info);

            % Copy synchronized fields back to stimulus struct
            Stimulusoptions(Stimnumber).TimeStimulusFrame = stim_sync.TimeStimulusFrame;
            Stimulusoptions(Stimnumber).RedFrameSynchronized = stim_sync.RedFrameSynchronized;
            Stimulusoptions(Stimnumber).TimeStimulusFrameVBL = stim_sync.TimeStimulusFrameVBL;
            Stimulusoptions(Stimnumber).RedFrameSynchronizedTrial = stim_sync.RedFrameSynchronizedTrial;
            Stimulusoptions(Stimnumber).FrameValid = stim_sync.FrameValid;
            Stimulusoptions(Stimnumber).RedSyncInfo = stim_sync.RedSyncInfo;
            Stimulusoptions(Stimnumber).RedSyncSession = stim_sync.RedSyncSession;
        end
    end
end

% --- Process ball camera data into locomotion ---
BT = Triggers.TimeBall;
if isempty(BT) || all(cellfun(@isempty, BT))
    % No ball triggers available
    return;
end

% Determine number of stimuli and mapping mode
if numel(BT) == numel(BallCamData_temp)
    N = numel(BallCamData_temp);
    StimAcqNumber = [];  % direct 1:1 mapping
else
    N = numel(Stimulusoptions);
    % Compute StimAcqNumber: map each stimulus to its acquisition block
    if ~isempty(acq_timestamps)
        ClockAcqEnd = seconds(acq_timestamps(:,2)) + SI_Info.Clockstart;
        ClockStimEnd = cat(1, Stimulusoptions(:).timeend);
        StimAcqNumber = zeros(size(ClockStimEnd));
        for s = 1:numel(ClockStimEnd)
            [~, I] = sort(abs(seconds(ClockAcqEnd - ClockStimEnd(s))));
            StimAcqNumber(s) = I(1);
        end
    else
        N = min(numel(BT), numel(BallCamData_temp));
        StimAcqNumber = [];
    end
end

% A training stimulus file holds ONE ball camera log for the whole session (all its trials) while
% Triggers.TimeBall has one cell per acquisition: training_ball_data time stamps it with all of them
train_ball = isfield(Stimulusoptions(1), 'type') && strcmp(Stimulusoptions(1).type, 'training') && ...
    numel(BallCamData_temp) == 1 && numel(BT) > 1;
if train_ball
    N = 0;
end

BallCamData = struct;
for n = 1:N
    % Find bad readings (cells that don't have exactly 5 values)
    Ca_f = find(cellfun(@numel, BallCamData_temp{n}) ~= 5);

    if numel(Ca_f) <= 0
        % Normal case: all readings valid
        BallCamData(n).Values = cell2mat(BallCamData_temp{n}');

        if isempty(StimAcqNumber)
            % Direct mapping
            if size(BallCamData(n).Values,1) < numel(BT{n}); BT{n}(Ca_f) = []; end
            BallCamData(n).Triggers_sec = BT{n}';
        else
            % Multi-acquisition mapping
            if StimAcqNumber(n) > 1
                temp_BT = [BT{StimAcqNumber(n)-1}'; BT{StimAcqNumber(n)}'];
                if size(BallCamData(n).Values,1) < numel(temp_BT); temp_BT(Ca_f) = []; end
                BallCamData(n).Triggers_sec = temp_BT;
            else
                if size(BallCamData(n).Values,1) < numel(BT{StimAcqNumber(n)})
                    BT{StimAcqNumber(n)}(Ca_f) = [];
                end
                BallCamData(n).Triggers_sec = BT{StimAcqNumber(n)}';
            end
        end
    else
        % Bad readings present: fill with NaN, splice good values
        warning('Ball camera data contains %d bad readings for stimulus %d', numel(Ca_f), n);
        BallCamData_temp{n}(Ca_f) = [];

        if isempty(StimAcqNumber)
            if numel(BallCamData_temp{n}) < numel(BT{n}); BT{n}(Ca_f) = []; end
            BallCamData(n).Triggers_sec = BT{n}';
        else
            if StimAcqNumber(n) > 1
                temp_BT = [BT{StimAcqNumber(n)-1}'; BT{StimAcqNumber(n)}'];
                if numel(BallCamData_temp{n}) < numel(temp_BT); temp_BT(Ca_f) = []; end
                BallCamData(n).Triggers_sec = temp_BT;
            else
                if numel(BallCamData_temp{n}) < numel(BT{StimAcqNumber(n)})
                    BT{StimAcqNumber(n)}(Ca_f) = [];
                end
                BallCamData(n).Triggers_sec = BT{StimAcqNumber(n)}';
            end
        end

        BallCamData(n).Values = nan(numel(BallCamData(n).Triggers_sec), 5);
        good_data = cell2mat(BallCamData_temp{n}');
        BallCamData(n).Values(1:size(good_data,1), :) = good_data;
    end
end

if train_ball
    % a session that was stopped early (cockpit_aborted_training*.mat) keeps the ball readings it had not yet
    % written into the stim file in <stim folder>\*_BallCamTail.mat (RB24_261104_1651: 30,403 readings in the stim
    % file + 26,936 in the tail = 57,339 for 57,340 ball triggers); without it the second half has no locomotion
    tail_readings = load_ball_tail(StimFiles(1), numel(BallCamData_temp{1}), sum(cellfun(@numel, BT)));
    main_readings = reshape(BallCamData_temp{1}, [], 1);
    use_tail = ~isempty(tail_readings);
    if use_tail && ~training_ball_data([main_readings; tail_readings], BT, true).SyncInfo.ok && ...
            training_ball_data(main_readings, BT, true).SyncInfo.ok
        % The pairing of readings and triggers verifies only WITHOUT the tail: the last reading(s) were taken after
        % the last trigger the headers hold (RB24_261104_1540: 1 reading, counter 241, after the last 2P frame;
        % both trigger trains end just before it), so they have no time stamp and are left out
        use_tail = false;
        fprintf('  Ball camera tail: %d reading(s) have no trigger in the headers, not used\n', numel(tail_readings));
    end
    if use_tail
        BallCamData = training_ball_data([main_readings; tail_readings], BT);
        BallCamData.SyncInfo.n_tail_readings = numel(tail_readings);
    else
        BallCamData = training_ball_data(main_readings, BT);
        BallCamData.SyncInfo.n_tail_readings = 0;
    end
    BallCamData.SyncInfo.n_tail_unused = numel(tail_readings) * ~use_tail;
end

% Calibrate locomotion — pre-load model to avoid uigetfile in batch mode
cal_paths = {
   % 'A:\Florian\BallCalibration\CalibrationModel.mat'
   % 'S:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat'
   % 'K:\Toni\BallCalibration\CalibrationModel.mat'
   % 'K:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat'
   % 'B:\fs3-joeschgrp\Toni\BallCalibration\CalibrationModel.mat'
   'C:\Users\rbondarenko\projects\2P_Scripts\newAggregation\\CalibrationModel.mat'
    fullfile(fileparts(mfilename('fullpath')), 'CalibrationModel.mat')
};
ModelParams = [];
for cp = 1:numel(cal_paths)
    if exist(cal_paths{cp}, 'file')
        tmp = load(cal_paths{cp}, 'ModelParams');
        ModelParams = tmp.ModelParams;
        break;
    end
end
if ~isempty(cal_idx) && ~isempty(ModelParams)
    if cal_idx > numel(ModelParams)
        warning('save_preprocessed:calIdxOOB', ...
            'cal_idx %d exceeds available models (1-%d); using auto-selection.', ...
            cal_idx, numel(ModelParams));
    else
        fprintf('  Calibration override: forcing model #%d (%s)\n', ...
            cal_idx, char(datetime(ModelParams(cal_idx).Date, 'Format', 'yyyy-MMM-dd')));
        ModelParams = ModelParams(cal_idx);
    end
end
try
    if ~isempty(ModelParams)
        LocomotionCal = transform_locomotion(BallCamData, SI_Info, ModelParams);
    else
        LocomotionCal = transform_locomotion(BallCamData, SI_Info);
    end
catch ME
    warning('Locomotion calibration failed: %s', ME.message);
    LocomotionCal = [];
end

end


function Stimopts = synchronize_red_frames(Stimopts, PT, timestamps, Stimnumber, SI_Info, pc_clock_reason, lag_center_s)
%SYNCHRONIZE_RED_FRAMES Synchronize projected frames using red sync frames
%   Aligns stimulus frame times from stimulus PC with camera acquisition times
%   using projector trigger times as reference
%
%   Every red frame should give one projector pulse. When the pulse count matches the number of
%   red frames they are paired in order. When pulses are missing (periodic red frames only),
%   match_red_pulses pairs them by time instead, and its result is used only if all of its
%   consistency checks pass. Otherwise the stimulus PC clock is used (RedFrameSynchronized = false).
%   pc_clock_reason (optional text) forces the PC clock, lag_center_s (optional, s) centres the window that
%   anchors the first pulse on a known PC-clock offset instead of the default one; both are used by
%   check_red_sync_consistency.
%   Stimopts.RedSyncInfo records what was done (see make_sync_info for the fields).

if nargin < 6; pc_clock_reason = ''; end
if nargin < 7; lag_center_s = NaN; end

Timestamp_Projected_Sync_Frames = PT(PT >= timestamps(Stimnumber,1) & PT <= timestamps(Stimnumber,2)+1);
PT_window = Timestamp_Projected_Sync_Frames;
valid_frame = ~isnan(Stimopts.Frameinfo.frame_count);
num_frames = sum(valid_frame);
is_flag_based = Stimopts.red_sync_nth_frame == 0;
% stimulus PC clock time of every (valid) frame on the session clock
pc_time = reshape(seconds(Stimopts.Frameinfo.frametime(valid_frame) - SI_Info.Clockstart), 1, []);
if ~is_flag_based
    Red_frame_number = find(mod(0:num_frames-1, Stimopts.red_sync_nth_frame) == 0);
else
    % red_sync_nth_frame==0 means no periodic sync (e.g. looming): red frames are
    % flagged once per trial instead, at subtrial_frame_count==1 (the trial's first
    % frame, per present_stimulus_addnew.m's redframe_flag, which only fires for this case)
    sfc = Stimopts.Frameinfo.subtrial_frame_count(~isnan(Stimopts.Frameinfo.frame_count));
    Red_frame_number = find(sfc == 1);
end
n_red_expected = numel(Red_frame_number);

% Handle edge cases where red frame count doesn't match projector sync count
if numel(Red_frame_number) == numel(Timestamp_Projected_Sync_Frames) + 1
    Red_frame_number(end) = [];
end

if numel(Red_frame_number) ~= numel(Timestamp_Projected_Sync_Frames)
    % Try filtering projector times by larger gaps (likely sync point)
    dt = diff([-inf Timestamp_Projected_Sync_Frames]);
    tc = Timestamp_Projected_Sync_Frames(dt > .04);
    if numel(Red_frame_number) == numel(tc) + 1
        Red_frame_number(end) = [];
    end
    if numel(Red_frame_number) == numel(tc)
        Timestamp_Projected_Sync_Frames = tc;
    end
end

% Pulses are missing and the simple recoveries above did not fix it: pair pulses and red frames by time
sync_info = make_sync_info('pulse count', true, '');
sync_info.n_pulses = numel(PT_window);
sync_info.n_expected = n_red_expected;
sync_info.n_lost = n_red_expected - numel(PT_window);
if isempty(pc_clock_reason) && ~is_flag_based && ...
        numel(Red_frame_number) ~= numel(Timestamp_Projected_Sync_Frames)
    [rf, pt, match] = match_red_pulses(PT_window, num_frames, Stimopts.red_sync_nth_frame, ...
        Stimopts.ifi, pc_time, lag_center_s);
    if match.ok
        Red_frame_number = rf;
        Timestamp_Projected_Sync_Frames = pt;
        sync_info = match;
    else
        sync_info.reason = match.reason;
    end
end

% If still mismatched or insufficient points, fall back to using PC clock times
if ~isempty(pc_clock_reason) || numel(Red_frame_number) < 2 || numel(Timestamp_Projected_Sync_Frames) < 2 || numel(Red_frame_number) ~= numel(Timestamp_Projected_Sync_Frames)
    why = pc_clock_reason;
    if isempty(why); why = sync_info.reason; end
    if isempty(why); why = sprintf('%d pulses for %d red frames', numel(PT_window), n_red_expected); end
    warning('Projector sync frames not matching! Using PC Clock Times instead (%s)', why);
    Stimopts.TimeStimulusFrame = seconds(Stimopts.Frameinfo.frametime - SI_Info.Clockstart);
    Stimopts.RedFrameSynchronized = false;
    sync_info = make_sync_info('PC clock', false, why);
    sync_info.n_pulses = numel(PT_window);
    sync_info.n_expected = n_red_expected;
    sync_info.n_lost = n_red_expected - numel(PT_window);
elseif num_frames < 1
    warning('Stimulus has no valid frames! Setting TimeStimulusFrame to NaN');
    Stimopts.TimeStimulusFrame = nan(1, numel(Stimopts.Frameinfo.frametime));
    Stimopts.RedFrameSynchronized = false;
    sync_info = make_sync_info('none', false, 'no valid frames');
elseif is_flag_based
    % Trial-anchored red frames are sparse and separated by an inter-trial
    % WaitSecs() that advances real time but not frame_count, so interpolating
    % or extrapolating across anchors (like the periodic case below) would
    % smear that wait into the per-frame mapping. Instead anchor each frame to
    % its own trial's onset and propagate using the actually-recorded
    % frame-to-frame time (Frameinfo.frametime) -- never bridging across a
    % trial boundary with an assumed constant rate.
    ft = Stimopts.Frameinfo.frametime(~isnan(Stimopts.Frameinfo.frame_count));
    Stimopts.TimeStimulusFrame = nan(1, num_frames);
    for a = 1:numel(Red_frame_number)
        if a < numel(Red_frame_number)
            seg = Red_frame_number(a):(Red_frame_number(a+1)-1);
        else
            seg = Red_frame_number(a):num_frames;
        end
        Stimopts.TimeStimulusFrame(seg) = Timestamp_Projected_Sync_Frames(a) + ...
            seconds(ft(seg) - ft(Red_frame_number(a)));
    end
    Stimopts.RedFrameSynchronized = true;
else
    % Interpolate to map stimulus frame numbers to camera-synchronized times
    Stimopts.TimeStimulusFrame = interp1(Red_frame_number, Timestamp_Projected_Sync_Frames, ...
        1:num_frames, 'linear', 'extrap');
    Stimopts.RedFrameSynchronized = true;
    if ~isempty(sync_info.InStall) && any(sync_info.InStall)
        % across a stall of the stimulus PC the pulses are further apart than a straight line through them
        % suggests; the PC clock saw the stall, so keep its steps and interpolate only the offset
        off_red = reshape(Timestamp_Projected_Sync_Frames, 1, []) - pc_time(Red_frame_number);
        off_all = interp1(Red_frame_number, off_red, 1:num_frames, 'linear', 'extrap');
        Stimopts.TimeStimulusFrame(sync_info.InStall) = pc_time(sync_info.InStall) + off_all(sync_info.InStall);
    end
end

if Stimopts.RedFrameSynchronized && ~strcmp(sync_info.method, 'gap-tolerant')
    % how far the pulses are from the PC clock times of the red frames they were paired with
    off = reshape(Timestamp_Projected_Sync_Frames, 1, []) - pc_time(Red_frame_number);
    sync_info.offset_median_ms = 1e3 * median(off);
    sync_info.offset_p05_ms = 1e3 * pctl(off, 5);
    sync_info.offset_p95_ms = 1e3 * pctl(off, 95);
end
Stimopts.RedSyncInfo = sync_info;

% Compute VBL-relative timing (time relative to first stimulus frame)
Stimopts.TimeStimulusFrameVBL = seconds(Stimopts.Frameinfo.frametime - Stimopts.Frameinfo.frametime(1)) + ...
    Stimopts.TimeStimulusFrame(1);
Stimopts.TimeStimulusFrameVBL = Stimopts.TimeStimulusFrameVBL(~isnan(Stimopts.TimeStimulusFrameVBL))';

end


function [Red_frame_number, Sync_times, info] = match_red_pulses(pulses, num_frames, nth, ifi, pc_time, lag_center_s, opts)
%MATCH_RED_PULSES Pair projector pulses with red-frame numbers when some pulses are missing
%   For periodic red frames (every nth stimulus frame). Pulses can be lost in bursts (seen in
%   RB24: about 2% of them, in bursts every ~73 s), so pulses and red frames cannot be paired by
%   position. Instead each pulse is numbered by counting the pulse periods since the previous pulse
%   (a lost pulse leaves a hole in the numbering instead of shifting everything after it), and the
%   first pulse is anchored with the stimulus PC clock. The pairs then go to the usual interpolation.
%
%   Accepted only if every check passes, otherwise info.ok is false and info.reason says which one
%   failed (the caller then uses the PC clock):
%     - at least min_pulses pulses, and a measured pulse period within period_tol of nth*ifi
%     - every interval between neighbouring pulses within max_int_dev periods of a whole number
%       of periods (otherwise the count is ambiguous); no two pulses less than one period apart
%     - at most max_edge red frames before the first or after the last pulse
%     - anchoring: with n red frames before/after the pulse train the first pulse can be red frame
%       1..n+1; the PC-clock offset (pulse time - PC clock time) must fall in a window exactly one
%       period wide, [lag_lo_s, lag_lo_s + period), which only one candidate can satisfy (lag_center_s,
%       if given, centres the window on a known offset instead). With no red frame left over (n = 0)
%       the anchoring is unambiguous and the window is not needed
%     - pulse time vs red-frame number is a straight line to within max_fit_s
%     - no frame ends up more than max_pc_s away from its PC clock time (catches pulses of another
%       block, which are tens of seconds off; real clock offsets are tens of ms)
%   info.InGap marks frames whose time is interpolated across lost pulses or extrapolated before
%   the first / after the last pulse (the others lie between two neighbouring pulses).
%
%   opts (optional struct, used for training trials) changes limits that are fixed above by default:
%     max_int_dev       replaces the limit of the interval check (up to < 0.5, where the rounding of an
%                       interval to a whole number of periods is still unambiguous)
%     outlier_s         > 0: a pulse further than this from the median PC-clock offset of its 11
%                       neighbouring pulses is dropped as an outlier (an isolated late pulse; the frames
%                       around it are then interpolated from its neighbours). Level steps that last, like
%                       the one-frame steps of RB24, are not outliers
%     max_outlier_frac  more than this share of outliers rejects the block
%     max_off_dev_s     every remaining pulse must be within this of the block's median PC-clock offset;
%                       a pulse numbered one period wrong would be 83 ms off
%     stall_fix         true (default without opts) corrects the period count where the stimulus PC stalled, see
%                       the comment at the rounding of the intervals; info.n_stall_fixes counts the places
%   info.n_outliers counts the dropped pulses.

% --- limits of the checks (named so a rig that needs different ones can change them here)
min_pulses = 20;
period_tol = 0.02;
max_int_dev = 0.35;
max_edge = 8;
lag_lo_s = -0.020;
max_fit_s = 2 * ifi;
max_pc_s = 1;
outlier_s = 0;
max_outlier_frac = 0.02;
max_off_dev_s = inf;
stall_fix = nargin < 7 || isempty(opts);    % off for training trials (opts given) unless opts.stall_fix says otherwise
jump_tol_s = 0.2 * nth * ifi;               % how far an offset jump may be from a whole number of periods (one frame)
if nargin >= 7 && ~isempty(opts)
    if isfield(opts, 'stall_fix'); stall_fix = opts.stall_fix; end
    if isfield(opts, 'max_int_dev'); max_int_dev = opts.max_int_dev; end
    if isfield(opts, 'outlier_s'); outlier_s = opts.outlier_s; end
    if isfield(opts, 'max_outlier_frac'); max_outlier_frac = opts.max_outlier_frac; end
    if isfield(opts, 'max_off_dev_s'); max_off_dev_s = opts.max_off_dev_s; end
end

if nargin < 6 || ~isfinite(lag_center_s)
    lag_lo = lag_lo_s;
else
    lag_lo = lag_center_s - nth * ifi / 2;
end
Red_frame_number = [];
Sync_times = [];
info = make_sync_info('gap-tolerant', false, '');
red_nums = find(mod(0:num_frames-1, nth) == 0);
M = numel(red_nums);
w = reshape(pulses, 1, []);
info.n_pulses = numel(w);
info.n_expected = M;
info.n_lost = M - numel(w);

if isempty(ifi) || ~isfinite(ifi) || ifi <= 0
    info.reason = 'stimulus frame interval (ifi) not available';
    return;
end
period_nom = nth * ifi;
if numel(w) < min_pulses
    info.reason = sprintf('only %d pulses, need %d to number them reliably', numel(w), min_pulses);
    return;
end

d = diff(w);
P = median(d);
if abs(P / period_nom - 1) > period_tol
    info.reason = sprintf('pulse period %.1f ms, expected %.1f ms', 1e3 * P, 1e3 * period_nom);
    return;
end
q = d / P;
step = round(q);
if any(step < 1)
    info.reason = 'two pulses less than one period apart';
    return;
end
idx0 = 1 + [0 cumsum(step)];            % red-frame number of every pulse if the first pulse is red frame 1

% A stall of the stimulus PC (frames in its own log that take several frame intervals) does not move the
% projector's rhythm, but the interval between the two pulses around it is no whole number of periods and the
% rounding above can count one period too many or few (RB24_261023_1538 moving_bar: 183 ms across a 107 ms stall
% counted as 2 periods instead of 1). From there on the pulse - PC clock offset sits a whole number of periods off.
% Such a jump of the offset between neighbouring pulses is undone, so the offset stays continuous.
fixed = false(size(d));
if stall_fix
    pc_at = @(rows) interp1(1:num_frames, pc_time, rows, 'linear', 'extrap');
    doff = diff(w - pc_at(1 + nth * (idx0 - 1)));
    m = round(doff / period_nom);
    fixed = m ~= 0 & abs(doff - m * period_nom) < jump_tol_s & step + m >= 1;
    if any(fixed)
        step(fixed) = step(fixed) + m(fixed);
        idx0 = 1 + [0 cumsum(step)];
        info.n_stall_fixes = nnz(fixed);
    end
end
info.interval_dev = max(abs(q(~fixed) - step(~fixed)));
if isempty(info.interval_dev); info.interval_dev = 0; end
if info.interval_dev > max_int_dev
    info.reason = sprintf('pulse intervals are not whole numbers of periods (off by up to %.2f periods)', info.interval_dev);
    return;
end

n_edge = M - idx0(end);                 % red frames left over before the first / after the last pulse
if n_edge < 0
    info.reason = sprintf('the pulses span %d periods but there are only %d red frames', idx0(end), M);
    return;
elseif n_edge > max_edge
    info.reason = sprintf('%d red frames before the first or after the last pulse (limit %d)', n_edge, max_edge);
    return;
end

if n_edge == 0
    shift = 0;
else
    cand = 0:n_edge;
    med = arrayfun(@(c) median(w - pc_time(red_nums(idx0 + c))), cand);
    in_window = med >= lag_lo & med < lag_lo + period_nom;
    if nnz(in_window) ~= 1
        info.reason = sprintf(['cannot tell which red frame the first pulse is: PC-clock offset would be ' ...
            '%s ms for first pulse = red frame %s (allowed %.0f to %.0f ms)'], ...
            mat2str(round(1e3 * med)), mat2str(1 + cand), 1e3 * lag_lo, 1e3 * (lag_lo + period_nom));
        return;
    end
    shift = cand(in_window);
    info.anchor_borderline = med(in_window) < lag_lo + 0.010 || med(in_window) >= lag_lo + period_nom - 0.010;
end

rf = red_nums(idx0 + shift);
in_stall = false(1, num_frames);        % frames between two pulses around a stall of the stimulus PC
for g = find(fixed)
    in_stall(rf(g)+1:rf(g+1)-1) = true;
end
if outlier_s > 0
    off0 = w - pc_time(rf);
    out = abs(off0 - movmedian(off0, 11)) > outlier_s;
    if nnz(out) > max_outlier_frac * numel(w)
        info.reason = sprintf('%d of %d pulses disagree with the PC-clock offset of their neighbours by more than %.0f ms (limit %.0f%%)', ...
            nnz(out), numel(w), 1e3 * outlier_s, 100 * max_outlier_frac);
        return;
    end
    info.n_outliers = nnz(out);
    w = w(~out);
    rf = rf(~out);
end
if isfinite(max_off_dev_s)
    off0 = w - pc_time(rf);
    if max(abs(off0 - median(off0))) > max_off_dev_s
        info.reason = sprintf('PC-clock offset of a pulse is %.0f ms from the median of the block (limit %.0f ms): pulse numbering is not reliable', ...
            1e3 * max(abs(off0 - median(off0))), 1e3 * max_off_dev_s);
        return;
    end
end
if info.n_stall_fixes > 0
    % the stall is a real step of the pulse train, so the straight line is fitted to the offset to the PC clock
    off1 = w - pc_time(rf);
    c = polyfit(rf, off1, 1);
    res = off1 - polyval(c, rf);
else
    c = polyfit(rf, w, 1);
    res = w - polyval(c, rf);
end
info.fit_res_max_ms = 1e3 * max(abs(res));
if max(abs(res)) > max_fit_s
    info.reason = sprintf('pulse times deviate up to %.0f ms from a straight line (limit %.0f ms)', ...
        1e3 * max(abs(res)), 1e3 * max_fit_s);
    return;
end
mapped = interp1(rf, w, 1:num_frames, 'linear', 'extrap');
dev = max(abs(mapped - pc_time));
if dev > max_pc_s
    info.reason = sprintf('mapped frame times are up to %.0f ms from the PC clock (limit %.0f ms): pulses of another block?', 1e3 * dev, 1e3 * max_pc_s);
    return;
end

d = diff(w);                 % intervals between the pulses that are used
gap = diff(rf) > nth;        % at least one red frame without a pulse in between
in_gap = false(1, num_frames);
in_gap(1:rf(1)-1) = true;
in_gap(rf(end)+1:end) = true;
for g = find(gap)
    in_gap(rf(g)+1:rf(g+1)-1) = true;
end
off = w - pc_time(rf);
Red_frame_number = rf;
Sync_times = w;
info.ok = true;
if info.n_stall_fixes > 0
    info.note = sprintf('period count corrected at %d place(s) where the stimulus PC log shows a stall', info.n_stall_fixes);
end
info.InStall = in_stall;
info.shift = shift;
info.n_gaps = nnz(gap);
info.longest_gap_s = max([0 d(gap)]);
info.offset_median_ms = 1e3 * median(off);
info.offset_p05_ms = 1e3 * pctl(off, 5);
info.offset_p95_ms = 1e3 * pctl(off, 95);
info.InGap = in_gap;

end


function v = pctl(x, p)
%PCTL Percentile p (0-100) of the vector x without the Statistics Toolbox
x = sort(x(~isnan(x)));
v = x(min(numel(x), max(1, round(p / 100 * numel(x)))));
end


function info = make_sync_info(method, ok, reason)
%MAKE_SYNC_INFO Stimulus.RedSyncInfo: how the stimulus frames of one block were synchronized
%   method            'pulse count' (pulses and red frames paired in order), 'gap-tolerant' (some
%                     pulses missing, paired by time with match_red_pulses), 'PC clock' (fallback,
%                     times are the stimulus PC clock, not the projector) or 'none'
%   ok                true if the frame times come from the projector pulses
%   reason            why the PC clock was used (empty when ok)
%   note              anything worth knowing about an accepted result (e.g. the first pulse was
%                     re-anchored with the PC-clock offset of the other blocks)
%   n_pulses          projector pulses found in the block's acquisition window
%   n_expected        red frames shown in the block
%   n_lost            n_expected - n_pulses
%   n_outliers        pulses dropped as isolated outliers (training trials only, see match_red_pulses)
%   n_stall_fixes     places where the period count was corrected for a stall of the stimulus PC
%   n_gaps           number of places where pulses are missing in the middle of the block
%   longest_gap_s     longest stretch without a pulse (s)
%   shift             'gap-tolerant' only: the first pulse is red frame 1 + shift
%   offset_median_ms  median of (pulse time - stimulus PC clock time of the same red frame), ms; the
%   offset_p05_ms     blocks of a session should agree on it within a few ms; two levels one stimulus
%   offset_p95_ms     frame apart are normal
%   fit_res_max_ms    'gap-tolerant' only: largest distance of a pulse from the straight line, ms
%   interval_dev      'gap-tolerant' only: largest distance of a pulse interval from a whole number
%                     of periods (0 = perfect, 0.5 = ambiguous)
%   anchor_borderline 'gap-tolerant' only: the PC-clock offset is within 10 ms of the edge of the
%                     window used to anchor the first pulse
%   InGap             'gap-tolerant' only: 1 x frames logical, true where the frame time is
%                     interpolated across lost pulses or extrapolated before/after the pulse train
%   InStall           'gap-tolerant' only: 1 x frames logical, true for the frames between the two pulses
%                     around a stall of the stimulus PC (n_stall_fixes); their time is the PC clock plus
%                     the pulse - PC clock offset interpolated between those pulses, not a straight line
info = struct('method', method, 'ok', ok, 'reason', reason, 'note', '', 'n_pulses', NaN, 'n_expected', NaN, ...
    'n_lost', NaN, 'n_outliers', 0, 'n_stall_fixes', 0, 'n_gaps', NaN, 'longest_gap_s', NaN, 'shift', NaN, 'offset_median_ms', NaN, ...
    'offset_p05_ms', NaN, 'offset_p95_ms', NaN, 'fit_res_max_ms', NaN, 'interval_dev', NaN, ...
    'anchor_borderline', false, 'InGap', [], 'InStall', []);
end


function Stimulusoptions = check_red_sync_consistency(Stimulusoptions, PT, timestamps, StimAcqNumber, SI_Info)
%CHECK_RED_SYNC_CONSISTENCY Report the red-frame sync of every block and compare the blocks
%   The offset between the stimulus PC clock and the projector pulses changes by only a few ms
%   over a session, so every synchronized block should show about the same median
%   (pulse time - PC clock time). A block that differs by more than tol_ms was most likely
%   numbered one red frame (83 ms) off. A block synchronized by match_red_pulses that disagrees
%   with the median of the blocks (needs >= 3 synchronized blocks to tell which one is off) is
%   matched again with its first pulse anchored on that median offset; if that does not give an
%   agreeing result either it goes back to PC clock times. Any other disagreement is only reported.

tol_ms = 33;     % 0.4 of a red-frame period (every 5th frame at 60 Hz); real blocks differ by < 15 ms
n = numel(Stimulusoptions);
if ~all(arrayfun(@(s) isfield(s, 'RedSyncInfo') && ~isempty(s.RedSyncInfo), Stimulusoptions))
    return;
end
info = [Stimulusoptions.RedSyncInfo];
fprintf('  Red-frame sync per stimulus block:\n');
fprintf('    %3s %-28s %-13s %11s %5s %5s %6s %10s\n', '#', 'type', 'method', 'pulses/red', 'lost', 'gaps', 'anchor', 'PC offset');
for k = 1:n
    anchor = '-';
    if strcmp(info(k).method, 'gap-tolerant'); anchor = sprintf('+%d', info(k).shift); end
    offs = '-';
    if isfinite(info(k).offset_median_ms); offs = sprintf('%+.0f ms', info(k).offset_median_ms); end
    fprintf('    %3d %-28s %-13s %5d/%-5d %5d %5s %6s %10s\n', k, char(string(Stimulusoptions(k).type)), ...
        info(k).method, info(k).n_pulses, info(k).n_expected, info(k).n_lost, num2str(info(k).n_gaps), anchor, offs);
end

off = [info.offset_median_ms];
synced = [Stimulusoptions.RedFrameSynchronized] & isfinite(off);
if nnz(synced) < 2
    return;
end
ref = median(off(synced));
far = synced & abs(off - ref) > tol_ms;
if ~any(far)
    fprintf('  Red-frame sync: the %d synchronized blocks agree on the PC-clock offset (median %+.0f ms, spread %.0f ms)\n', ...
        nnz(synced), ref, max(off(synced)) - min(off(synced)));
    return;
end
for k = find(far)
    if nnz(synced) >= 3 && strcmp(info(k).method, 'gap-tolerant')
        % anchor the first pulse on the median offset of the blocks that agree with each other
        ref_k = median(off(synced & ~far));
        stim_sync = synchronize_red_frames(Stimulusoptions(k), PT, timestamps, StimAcqNumber(k), SI_Info, '', ref_k / 1e3);
        new_off = stim_sync.RedSyncInfo.offset_median_ms;
        if stim_sync.RedFrameSynchronized && abs(new_off - ref_k) <= tol_ms
            stim_sync.RedSyncInfo.note = sprintf('first pulse re-anchored: PC-clock offset %+.0f ms (was %+.0f ms) like the other blocks (%+.0f ms)', new_off, off(k), ref_k);
            fprintf('  Stimulus %d: %s\n', k, stim_sync.RedSyncInfo.note);
        else
            why = sprintf('PC-clock offset %+.0f ms but %+.0f ms in the other blocks, and anchoring it like them did not work either', off(k), ref_k);
            stim_sync = synchronize_red_frames(Stimulusoptions(k), PT, timestamps, StimAcqNumber(k), SI_Info, why);
        end
        Stimulusoptions(k).TimeStimulusFrame = stim_sync.TimeStimulusFrame;
        Stimulusoptions(k).RedFrameSynchronized = stim_sync.RedFrameSynchronized;
        Stimulusoptions(k).TimeStimulusFrameVBL = stim_sync.TimeStimulusFrameVBL;
        Stimulusoptions(k).RedSyncInfo = stim_sync.RedSyncInfo;
    else
        warning('check_red_sync_consistency:offset', ...
            ['Stimulus %d: PC-clock offset %+.0f ms differs from the median of the blocks (%+.0f ms) by more than %d ms; ' ...
            'check its red-frame sync'], k, off(k), ref, tol_ms);
    end
end

end


function Stimopts = synchronize_red_frames_training(Stimopts, PT, acq_timestamps, SI_Info)
%SYNCHRONIZE_RED_FRAMES_TRAINING Synchronize every trial of a training stimulus with the projector pulses
%   A training stimulus file holds all the trials of a session. Frameinfo.* are frames x trials matrices
%   (padded to the longest trial, and in some trials with placeholder rows), every trial is imaged in its
%   own acquisition, and the red frame (every red_sync_nth_frame-th frame, counted over the whole trial)
%   gives one projector pulse. Each trial is synchronized on its own:
%     1. its real frames are picked out: rows with a frame number >= 1 and a time inside the trial
%        (padding and placeholder rows are ignored and stay NaN in the output). The frame numbers must be
%        1..n in order and the times must increase, otherwise the trial keeps PC clock times
%     2. at least min_cover of the trial must lie inside imaged acquisitions
%     3. the projector pulses from the trial's first frame to its last red frame are paired with its red
%        frames by match_red_pulses, which tolerates lost pulses (see there for its consistency checks)
%     4. the PC-clock offset (pulse time - PC clock time) of a trial must agree with that of the trials around
%        it (the two clocks drift apart slowly and smoothly, but not necessarily linearly over a session): it is
%        compared with a line through its n_neighbours nearest synchronized trials in time. A trial that
%        disagrees, or one that could not be matched, is matched again with its first pulse anchored on that
%        expectation; if that does not agree either, it keeps PC clock times
%     5. the frame times are interpolated between the pulses and must increase
%   A trial that cannot be synchronized safely keeps PC clock times and says why in its RedSyncInfo.
%
%   Output, added to Stimopts (frames x trials like Frameinfo, NaN where there is no real frame):
%     TimeStimulusFrame          time (s, session clock = TimeCa) of every real frame
%     TimeStimulusFrameVBL       time of the trial's first frame + stimulus PC clock intervals
%     FrameValid                 logical, true for the real frames (Frameinfo rows to use)
%     RedFrameSynchronizedTrial  1 x trials logical, true if the trial's times come from the pulses
%     RedFrameSynchronized       true only if every trial was synchronized
%     RedSyncInfo(trial)         how the trial was done: make_sync_info fields + trial, acquisition (the
%                                acquisition it overlaps most), acq_cover (share of the trial inside imaged
%                                acquisitions), n_frames (real frames), n_dropped_rows (rows with a time that
%                                are not real frames, e.g. placeholders)
%     RedSyncSession             session summary: n_trials, n_synced, trials_pc_clock, offset_min_ms,
%                                offset_max_ms, drift_ms_per_min and offset_start_ms (overall straight-line
%                                trend of the trials' PC-clock offsets, information only), max_dev_ms (largest
%                                distance of a trial from the offset expected from its neighbouring trials)

% --- limits of the checks (named so a rig that needs different ones can change them here)
lag_tol_s = 0.020;      % a trial's PC-clock offset may differ this much from that of its neighbouring trials
n_neighbours = 8;       % neighbouring trials (nearest in time) that give the expected offset of a trial
min_trials_fit = 5;     % synchronized trials needed for that comparison
min_cover = 0.9;        % share of a trial that must lie inside imaged acquisitions
win_margin_s = 0.005;   % pulse window of a trial: see match_training_trial
% limits of match_red_pulses for training trials: the projector chain shows isolated pulses 2 frames (33 ms) late and
% level steps of 1-2 frames in RB24 (see the opts of match_red_pulses), so intervals may be off by up to 0.45 period
opts = struct('max_int_dev', 0.45, 'outlier_s', 0.025, 'max_outlier_frac', 0.02, 'max_off_dev_s', 0.0375);

FI = Stimopts.Frameinfo;
nth = Stimopts.red_sync_nth_frame;
periodic = isscalar(nth) && nth >= 1;
[nrow, ntr] = size(FI.frametime);
PT = reshape(PT, 1, []);
pc_all = seconds(FI.frametime - SI_Info.Clockstart);        % stimulus PC clock on the session clock

% --- 1. real frames of every trial
t_lo = -inf(1, ntr);
t_hi = inf(1, ntr);
if isfield(Stimopts, 'timestart') && isfield(Stimopts, 'timeend') && isdatetime(Stimopts.timestart) && ...
        isdatetime(Stimopts.timeend) && numel(Stimopts.timestart) == ntr && numel(Stimopts.timeend) == ntr
    t_lo = reshape(seconds(Stimopts.timestart - SI_Info.Clockstart), 1, []) - 1;
    t_hi = reshape(seconds(Stimopts.timeend - SI_Info.Clockstart), 1, []) + 1;
end
valid = false(nrow, ntr);
pc = cell(1, ntr);
why = repmat({''}, 1, ntr);
n_drop = zeros(1, ntr);
for q = 1:ntr
    fc = FI.frame_count(:, q);
    t = pc_all(:, q);
    v = ~isnan(fc) & fc >= 1 & ~isnan(t) & t >= t_lo(q) & t <= t_hi(q);
    valid(:, q) = v;
    pc{q} = reshape(t(v), 1, []);
    n_drop(q) = nnz(~isnan(t)) - nnz(v);
    if ~any(v)
        why{q} = 'no real frames';
    elseif ~isequal(reshape(fc(v), 1, []), 1:nnz(v))
        why{q} = 'frame numbers are not 1..n in order';
    elseif any(diff(pc{q}) <= 0)
        why{q} = 'frame times do not increase';
    end
end
if isfield(Stimopts, 'ifi') && ~isempty(Stimopts.ifi) && Stimopts.ifi(1) > 0
    ifi = Stimopts.ifi(1);
else
    dts = cellfun(@diff, pc(cellfun(@isempty, why)), 'UniformOutput', false);
    ifi = median([dts{:}]);       % frame interval of the stimulus PC (NaN if there are no frames)
end

% --- 2./3. every trial on its own
acq = nan(1, ntr);
cover = nan(1, ntr);
rf = cell(1, ntr);
pt = cell(1, ntr);
info = cell(1, ntr);
tried = false(1, ntr);
for q = 1:ntr
    if ~isempty(pc{q})
        ov = max(0, min(acq_timestamps(:, 2), pc{q}(end)) - max(acq_timestamps(:, 1), pc{q}(1)));
        [~, acq(q)] = max(ov);
        cover(q) = sum(ov) / max(pc{q}(end) - pc{q}(1), eps);
    end
    if ~periodic
        fail = 'no periodic red frames (red_sync_nth_frame < 1)';
    elseif ~isempty(why{q})
        fail = why{q};
    elseif cover(q) < min_cover
        fail = sprintf('only %.0f%% of the trial lies inside imaged acquisitions', 100 * cover(q));
    else
        fail = '';
    end
    if isempty(fail)
        tried(q) = true;
        [rf{q}, pt{q}, info{q}] = match_training_trial(PT, pc{q}, nth, ifi, NaN, win_margin_s, opts);
        if ~info{q}.ok
            info{q} = pc_clock_info(info{q}.reason, info{q});
        end
    else
        info{q} = make_sync_info('PC clock', false, fail);
    end
end

% --- 4. the PC-clock offset of a trial must agree with the trials around it
% A trial whose first pulse was numbered one pulse period (or one frame) off has an offset that jumps away from its
% neighbours. The stimulus PC clock does not drift linearly against the 2P clock over a whole session
% (RB24_261104_1540: +0.4 ms per trial for 65 trials, then -1.5 ms per trial down to +4 ms), so a trial is compared
% with a robust line through its n_neighbours nearest synchronized trials in time (itself left out), not with one
% line for the session. Trials that disagree are rejected from the reference and matched again anchored on it.
% Limit: this finds trials that disagree with their neighbours. If the true offset drifted across the edge of the
% anchoring window of match_red_pulses for ALL later trials, they would agree with each other (anchor_borderline,
% reported below, warns when a trial is near that edge).
ok = cellfun(@(s) s.ok, info);
off = cellfun(@(s) s.offset_median_ms, info);
tmid = nan(1, ntr);
for q = 1:ntr
    if ~isempty(pc{q}); tmid(q) = mean(pc{q}([1 end])); end
end
sel = ok & isfinite(off) & isfinite(tmid);
pred = nan(1, ntr);                     % offset expected from the neighbouring trials
if nnz(sel) >= min_trials_fit
    good = sel;
    for it = 1:6
        pred = neighbour_offset_reference(tmid, off, good, n_neighbours);
        newgood = sel & abs(off - pred) <= 1e3 * lag_tol_s;
        if nnz(newgood) < min_trials_fit; pred(:) = NaN; break; end
        if isequal(newgood, good); break; end
        good = newgood;
    end
end
if any(isfinite(pred))
    redo = find(tried & isfinite(pred) & (~ok | (sel & abs(off - pred) > 1e3 * lag_tol_s)));
    for q = redo
        [rf2, pt2, m2] = match_training_trial(PT, pc{q}, nth, ifi, pred(q) / 1e3, win_margin_s, opts);
        if m2.ok && abs(m2.offset_median_ms - pred(q)) <= 1e3 * lag_tol_s
            m2.note = sprintf('first pulse anchored on the neighbouring trials: PC-clock offset %+.0f ms (neighbours %+.0f ms)', ...
                m2.offset_median_ms, pred(q));
            rf{q} = rf2;
            pt{q} = pt2;
            info{q} = m2;
            ok(q) = true;
        elseif ok(q)
            info{q} = pc_clock_info(sprintf(['PC-clock offset %+.0f ms but the neighbouring trials give %+.0f ms, ' ...
                'and anchoring on them did not work either'], off(q), pred(q)), info{q});
            ok(q) = false;
        end
    end
end

% --- 5. frame times of every trial, frames x trials like Frameinfo
TSF = nan(nrow, ntr);
VBL = nan(nrow, ntr);
for q = 1:ntr
    v = valid(:, q);
    n = nnz(v);
    if n == 0; continue; end
    if ok(q)
        % between the pulses interpolate; before the first / after the last one continue with the slope of the
        % whole trial (the last two pulses alone would pass their timing jitter or a level step on)
        k = 1:n;
        t = interp1(rf{q}, pt{q}, k, 'linear');
        slope = polyfit(rf{q}, pt{q}, 1);
        t(k < rf{q}(1)) = pt{q}(1) + slope(1) * (k(k < rf{q}(1)) - rf{q}(1));
        t(k > rf{q}(end)) = pt{q}(end) + slope(1) * (k(k > rf{q}(end)) - rf{q}(end));
        if any(diff(t) <= 0)
            info{q} = pc_clock_info('frame times from the pulses do not increase', info{q});
            ok(q) = false;
            t = pc{q};
        end
    else
        t = pc{q};
    end
    TSF(v, q) = t(:);
    ft = FI.frametime(v, q);
    VBL(v, q) = t(1) + seconds(ft - ft(1));
end
for q = 1:ntr
    info{q}.trial = q;
    info{q}.acquisition = acq(q);
    info{q}.acq_cover = cover(q);
    info{q}.n_frames = numel(pc{q});
    info{q}.n_dropped_rows = n_drop(q);
end
info = [info{:}];

sess = struct('n_trials', ntr, 'n_synced', nnz(ok), 'trials_pc_clock', find(~ok), 'offset_min_ms', NaN, ...
    'offset_max_ms', NaN, 'drift_ms_per_min', NaN, 'offset_start_ms', NaN, 'max_dev_ms', NaN, ...
    'tolerance_ms', 1e3 * lag_tol_s);
offF = [info.offset_median_ms];
if any(ok)
    sess.offset_min_ms = min(offF(ok));
    sess.offset_max_ms = max(offF(ok));
end
if any(isfinite(pred(ok)))
    % overall trend of the offset (information only; the decision above is made against the neighbouring trials)
    tm = tmid(ok & isfinite(offF)); mu = mean(tm);
    p = polyfit(tm - mu, offF(ok & isfinite(offF)), 1);
    sess.drift_ms_per_min = 60 * p(1);
    sess.offset_start_ms = polyval(p, -mu);
    sess.max_dev_ms = max(abs(offF(ok) - pred(ok)));    % largest distance of a trial from its neighbouring trials
end

Stimopts.TimeStimulusFrame = TSF;
Stimopts.TimeStimulusFrameVBL = VBL;
Stimopts.FrameValid = valid;
Stimopts.RedFrameSynchronizedTrial = ok;
Stimopts.RedFrameSynchronized = all(ok);
Stimopts.RedSyncInfo = info;
Stimopts.RedSyncSession = sess;

% --- report
fprintf('  Red-frame sync of the %d training trials: %d synchronized from the projector pulses\n', ntr, nnz(ok));
if any(ok)
    n_exp = sum([info(ok).n_expected]);
    n_lost = sum([info(ok).n_lost]);
    fprintf('    pulses lost: %d of %d red frames (%.1f%%), trials with a gap: %d; isolated late pulses dropped: %d\n', n_lost, n_exp, ...
        100 * n_lost / n_exp, nnz([info(ok).n_gaps] > 0), sum([info(ok).n_outliers]));
    fprintf('    PC-clock offset (pulse - PC clock): %+.0f to %+.0f ms', sess.offset_min_ms, sess.offset_max_ms);
    if isfinite(sess.max_dev_ms)
        fprintf('; overall trend %+.2f ms/min, largest deviation of a trial from its neighbouring trials %.1f ms (limit %.0f ms)', ...
            sess.drift_ms_per_min, sess.max_dev_ms, sess.tolerance_ms);
    end
    fprintf('\n');
    re = find(ok & ~cellfun(@isempty, {info.note}));
    if ~isempty(re)
        fprintf('    first pulse anchored on the neighbouring trials for trial(s) %s\n', mat2str(re));
    end
    bl = find(ok & [info.anchor_borderline]);
    if ~isempty(bl)
        fprintf('    first pulse anchoring is within 10 ms of the edge of its window for trial(s) %s\n', mat2str(bl));
    end
end
bad = find(~ok);
if ~isempty(bad)
    fprintf('    PC clock times used for %d trial(s):\n', numel(bad));
    for q = bad(1:min(numel(bad), 10))
        fprintf('      trial %d: %s\n', q, info(q).reason);
    end
    warning('Projector sync frames not matching for training! Using PC Clock Times instead for trial(s) %s (see RedSyncInfo)', ...
        mat2str(bad));
end

end


function ref = neighbour_offset_reference(tmid, off, good, k)
%NEIGHBOUR_OFFSET_REFERENCE PC-clock offset expected for every trial from the trials around it
%   tmid  trial times (s), off  PC-clock offsets (ms), good  logical, trials that may serve as reference,
%   k     number of reference trials. For every trial (also those not in good) the k nearest trials in time
%   that are in good, the trial itself left out, are taken and a straight line through their offsets is
%   evaluated at the trial's time. The line is a Theil-Sen line (median of the slopes between all pairs, median
%   intercept), so that a few wrong neighbours do not pull it: with a least-squares line, 3 wrong trials at the end
%   of a session vouch for each other and the correct ones next to them are rejected instead (tested on the
%   offsets of RB24_261104_1540 with 83 ms errors injected: a robust line over 8 neighbours finds up to 4 wrong
%   trials at either end or in the middle). With fewer than 4 reference trials their median is used, with fewer
%   than 3 the result is NaN.
ref = nan(size(off));
idx = find(good & isfinite(tmid) & isfinite(off));
for q = find(isfinite(tmid))
    cand = idx(idx ~= q);
    if numel(cand) < 3; continue; end
    [~, o] = sort(abs(tmid(cand) - tmid(q)));
    nb = cand(o(1:min(k, numel(cand))));
    x = reshape(tmid(nb) - tmid(q), [], 1);
    y = reshape(off(nb), [], 1);
    if numel(nb) >= 4 && range(x) > 0
        [i1, i2] = find(triu(true(numel(nb), numel(nb)), 1));
        dx = x(i2) - x(i1);
        use = dx ~= 0;
        slope = median((y(i2(use)) - y(i1(use))) ./ dx(use));
        ref(q) = median(y - slope * x);   % value of the line at the trial's own time (x = 0)
    else
        ref(q) = median(y);
    end
end
end


function [rf, pt, info] = match_training_trial(PT, pc, nth, ifi, lag_center_s, win_margin_s, opts)
%MATCH_TRAINING_TRIAL Pair the projector pulses of one training trial with its red frames
%   pc are the PC-clock times of the trial's frames. The pulse of a red frame comes lag_center_s after its
%   PC-clock time (default: the middle of the anchoring window of match_red_pulses), and pulses are taken
%   from half a period (+ win_margin_s) before that for the first frame to the same after it for the last
%   red frame. That keeps out the pulses of the neighbouring trials, and also the pulse of a red frame that
%   was projected right after the last logged frame, which comes one period after the last one.
%   See match_red_pulses.

n = numel(pc);
last_red = pc(1 + nth * floor((n - 1) / nth));
lag_in = lag_center_s;
if ~isfinite(lag_center_s)
    lag_center_s = -0.020 + nth * ifi / 2;          % centre of the default anchoring window
end
half = nth * ifi / 2 + win_margin_s;
w = PT(PT >= pc(1) + lag_center_s - half & PT <= last_red + lag_center_s + half);
[rf, pt, info] = match_red_pulses(w, n, nth, ifi, pc, lag_in, opts);

end


function info = pc_clock_info(why, prev)
%PC_CLOCK_INFO RedSyncInfo of a trial that keeps PC clock times, with the pulse counts of its failed match

info = make_sync_info('PC clock', false, why);
info.n_pulses = prev.n_pulses;
info.n_expected = prev.n_expected;
info.n_lost = prev.n_lost;

end


function tail = load_ball_tail(StimFile, n_main, n_trig)
%LOAD_BALL_TAIL Ball camera readings saved next to the stim file of a session that was stopped early
%   Returns the readings of every non-empty *_BallCamTail.mat in the stim file's folder (cell, one 1x5 row per
%   reading, as in the stim file), or {} when there is none or when appending it would NOT bring the number of
%   readings closer to the number of ball triggers (then it is reported and ignored).
tail = {};
tf = dir(fullfile(StimFile.folder, '*_BallCamTail.mat'));
tf = tf([tf.bytes] > 0);
for k = 1:numel(tf)
    S = load(fullfile(tf(k).folder, tf(k).name));
    if ~isfield(S, 'BallCamData') || isempty(S.BallCamData); continue; end
    t = reshape(cellfun(@(r) double(r(:))', S.BallCamData, 'UniformOutput', false), [], 1);
    if abs(n_main + numel(tail) + numel(t) - n_trig) < abs(n_main + numel(tail) - n_trig)
        tail = [tail; t]; %#ok<AGROW>
        fprintf('  Ball camera tail %s: %d readings appended (%d in the stim file, %d ball triggers)\n', tf(k).name, numel(t), n_main, n_trig);
    else
        warning('Ball camera tail %s (%d readings) does not fit %d triggers for %d readings: ignored', tf(k).name, numel(t), n_trig, n_main + numel(tail));
    end
end
end


function Ball = training_ball_data(readings, BT, quiet)
%TRAINING_BALL_DATA Ball camera readings of a training session, each with a time stamp
%   A training stimulus file holds ONE ball camera log for the whole session while Triggers.TimeBall has one
%   cell per acquisition (BT), so the time stamps of the log are all the triggers in acquisition order -
%   except for triggers that have no reading. The camera counter (column 5 of a reading, 1..L and then
%   wrapping) shows where readings are missing, and in RB24 it jumps exactly where the trigger list has
%   a gap (ScanImage starting a new acquisition), and the trigger list has one trigger more than the log
%   at the start of the session and after each gap. With those triggers removed, readings and
%   triggers pair one to one. This is accepted only if
%     - the counter jumps and the trigger gaps (interval > 5 x the median) are equally many,
%     - the trigger index minus the reading index at the h-th jump / gap is h, and
%     - there are 1 + (number of gaps) triggers more than readings,
%   (or if there are exactly as many triggers as readings and no gap / jump); otherwise Triggers_sec is
%   the first numel(readings) triggers and SyncInfo.ok is false: the pairing is then not verified.
%   Readings that are not 5 numbers are NaN.
%   Ball.SyncInfo: ok, method, reason, n_readings, n_triggers, n_dropped_triggers, counter_jumps (reading index
%   before each jump), readings_missing (per jump), trigger_gaps (trigger index before each gap), gap_s.
%   quiet (optional, default false): print / warn nothing (to try a variant of the readings).

if nargin < 3; quiet = false; end
nr = numel(readings);
good = reshape(cellfun(@numel, readings) == 5, [], 1);
Values = nan(nr, 5);
Values(good, :) = cell2mat(reshape(readings(good), [], 1));
if any(~good) && ~quiet
    warning('Ball camera data contains %d bad readings (set to NaN)', nnz(~good));
end
tb = cell2mat(cellfun(@(x) reshape(x, [], 1), reshape(BT, [], 1), 'UniformOutput', false));

% holes in the camera counter
L = max(Values(:, 5));
idx = find(good & ~isnan(Values(:, 5)));
miss = mod(diff(Values(idx, 5)) - diff(idx), L);          % 0 where no reading is missing
hole_after = reshape(idx(miss ~= 0), [], 1);              % reading before each jump
hole_missing = reshape(miss(miss ~= 0), [], 1);           % readings missing at each jump
% gaps in the trigger list
iv = diff(tb);
gap_after = reshape(find(iv > 5 * median(iv)), [], 1);
nh = numel(hole_after);

info = struct('ok', false, 'method', '', 'reason', '', 'n_readings', nr, 'n_triggers', numel(tb), ...
    'n_dropped_triggers', 0, 'counter_jumps', hole_after, 'readings_missing', hole_missing, ...
    'trigger_gaps', gap_after, 'gap_s', iv(gap_after));
if numel(tb) == nr && nh == 0 && isempty(gap_after)
    info.ok = true;
    info.method = 'one to one';
elseif L >= 10 && nh == numel(gap_after) && all(gap_after - hole_after == (1:nh)') && ...
        numel(tb) - nr == 1 + nh
    drop = [1; gap_after + 1];
    tb(drop) = [];
    info.ok = true;
    info.method = 'trigger at the session start and after each gap removed';
    info.n_dropped_triggers = numel(drop);
else
    info.reason = sprintf(['%d triggers for %d readings, %d counter jump(s) and %d trigger gap(s) do not fit ' ...
        'the expected pattern, the pairing is not verified'], numel(tb), nr, nh, numel(gap_after));
    info.method = 'first triggers';
    tb = tb(1:min(numel(tb), nr));
    tb = [tb; nan(nr - numel(tb), 1)];
    if ~quiet; warning('Ball camera time stamps: %s', info.reason); end
end
if ~quiet
    fprintf('  Ball camera log: %d readings, %d triggers (%s)\n', nr, info.n_triggers, info.method);
end

Ball = struct('Values', Values, 'Triggers_sec', tb, 'SyncInfo', info);

end


function StimAcqNumber = match_stim_to_acq(Stimulusoptions, acq_timestamps, SI_Info)
%MATCH_STIM_TO_ACQ Map each training stimulus to its acquisition block
%   Uses stimulus end time to find closest acquisition block end time

ClockAcqEnd = seconds(acq_timestamps(:,2)) + SI_Info.Clockstart;
ClockStimEnd = cat(1, Stimulusoptions(:).timeend);
StimAcqNumber = zeros(size(ClockStimEnd));
for s = 1:numel(ClockStimEnd)
    [~, I] = sort(abs(seconds(ClockAcqEnd - ClockStimEnd(s))));
    StimAcqNumber(s) = I(1);
end

end


function DLCSupporting = apply_frame_selection(DLCSupporting, selected_frameIdx)
%APPLY_FRAME_SELECTION Apply frame selection to DLC supporting data

if isfield(DLCSupporting, 'Probability')
    FN = fieldnames(DLCSupporting.Probability);
    for fn = 1:numel(FN)
        DLCSupporting.Probability.(FN{fn}) = DLCSupporting.Probability.(FN{fn})(selected_frameIdx, :, :);
    end
end

if isfield(DLCSupporting, 'Pupil')
    FN = fieldnames(DLCSupporting.Pupil);
    for fn = 1:numel(FN)
        DLCSupporting.Pupil.(FN{fn}) = DLCSupporting.Pupil.(FN{fn})(selected_frameIdx, :, :);
    end
end

if isfield(DLCSupporting, 'Ball')
    FN = fieldnames(DLCSupporting.Ball);
    for fn = 1:numel(FN)
        DLCSupporting.Ball.(FN{fn}) = DLCSupporting.Ball.(FN{fn})(selected_frameIdx, :, :);
    end
end

end


function Behav = apply_frame_selection_behav(Behav, selected_frameIdx)
%APPLY_FRAME_SELECTION_BEHAV Apply frame selection to behavioral data

FN = fieldnames(Behav);
for fn = 1:numel(FN)
    Behav.(FN{fn}) = Behav.(FN{fn})(selected_frameIdx, :, :);
end

end
