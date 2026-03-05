function import_video_behavior(varargin)
%IMPORT_VIDEO_BEHAVIOR Universal DLC data importer with config-based paths
%
%   import_video_behavior()
%   import_video_behavior('config', config)
%   import_video_behavior('animals', {'AnimalRB5'})
%
%   Config-based version of RB_import_VideoBehav.m that uses get_user_config
%   for paths. Includes bug fixes from the RB version.
%
%   Bug fixes compared to RB_import_VideoBehav.m:
%   - Fixed datetime format 'yyyy_MM_dd HH:mm:ss: SSS' -> 'yyyy_MM_dd HH:mm:ss:SSS'
%     (removed extra space before SSS that caused parsing failures)
%   - Fixed extra spaces in property access (e.g., 'animallist. isdir')
%
%   Optional Name-Value Arguments:
%       'config' - User configuration from get_user_config
%           If not provided, will process all users based on animals found
%
%       'animals' - Cell array of animal names to filter
%           e.g., {'AnimalRB5', 'AnimalRB14'}
%
%       'overwrite' - Overwrite existing DLC files (default: false)
%
%       'delete_camera_tifs' - Delete camera TIFs after processing (default: false)
%
%   Parameters (can be customized):
%       movmed_filt   - Median filter window (default: 3)
%       pupil_min_p   - Minimum pupil probability threshold (default: 0.95)
%       body_min_p    - Minimum body probability threshold (default: 0.8)
%       eye_r_mm      - Eye radius in mm (default: 1.5)
%       eye_px_per_mm - Pixels per mm for eye tracking (default: 23.4)
%       body_vid      - Body video width (default: 480)
%
%   Example:
%       % Import all available DLC data
%       import_video_behavior();
%
%       % Import for specific animals
%       import_video_behavior('animals', {'AnimalRB5', 'AnimalRB14'});
%
%       % Force overwrite existing files
%       import_video_behavior('overwrite', true);
%
%   See also: get_user_config, save_preprocessed, run_aggregation

%% Parse inputs
p = inputParser;
addParameter(p, 'config', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'animals', {}, @iscell);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'delete_camera_tifs', false, @islogical);
addParameter(p, 'movmed_filt', 3, @isnumeric);
addParameter(p, 'pupil_min_p', 0.95, @isnumeric);
addParameter(p, 'body_min_p', 0.8, @isnumeric);
addParameter(p, 'eye_r_mm', 1.5, @isnumeric);
addParameter(p, 'eye_px_per_mm', 23.4, @isnumeric);
addParameter(p, 'body_vid', 480, @isnumeric);
parse(p, varargin{:});

config = p.Results.config;
animals_filter = p.Results.animals;
overwrite = p.Results.overwrite;
delete_camera_tifs = p.Results.delete_camera_tifs;

% Create params structure
params.movmed_filt = p.Results.movmed_filt;
params.pupil_min_p = p.Results.pupil_min_p;
params.body_min_p = p.Results.body_min_p;
params.eye_r_mm = p.Results.eye_r_mm;
params.eye_px_per_mm = p.Results.eye_px_per_mm;
params.body_vid = p.Results.body_vid;
params.overwrite = overwrite;
params.delete_camera_tifs = delete_camera_tifs;

%% Get directories based on config or detect from animals
if isempty(config)
    if ~isempty(animals_filter)
        config = get_user_config(animals_filter{1});
    else
        config = get_user_config('AnimalFS0');
    end
end

parentdir = config.CameraDir;
archivedir = fullfile(fileparts(config.CameraDir(1:end-1)), 'preCamera', filesep);
DLCdir = config.DLCDir;
preprocesseddir = config.OutputDir;

% Ensure directories exist
if ~exist(DLCdir, 'dir'); mkdir(DLCdir); end
if ~exist(archivedir, 'dir'); mkdir(archivedir); end

%% Start processing
fprintf('=== Starting DLC import scan ===\n');
fprintf('Camera directory: %s\n', parentdir);
fprintf('DLC output: %s\n', DLCdir);
fprintf('Preprocessed: %s\n', preprocesseddir);
fprintf('\n');

if ~exist(parentdir, 'dir')
    fprintf('Camera directory not found: %s\n', parentdir);
    return;
end

% Get folder list - FIXED: Properly filter directories
animallist = dir(parentdir);
animallist = animallist([animallist.isdir]);  % Fixed: removed space
animallist = animallist(~ismember({animallist.name}, {'.', '..'}));

if isempty(animallist)
    fprintf('No folders found in %s\n', parentdir);
    return;
end

% Detect folder structure: animal subfolders vs flat recording folders
folder_names = {animallist.name};
if ~isempty(animals_filter)
    has_animal_dirs = any(ismember(folder_names, animals_filter));
else
    % Folders without date suffixes (e.g. _YYMMDD_HHMM) are animal dirs
    has_animal_dirs = any(cellfun(@isempty, regexp(folder_names, '_\d{6}_\d{4}$')));
end
flat_mode = ~has_animal_dirs;

% Set up animal list depending on folder structure
if flat_mode
    fprintf('Flat folder structure detected (no animal subfolders).\n');
    tokens = regexp(folder_names, '^(.+?)_\d{6}_\d{4}$', 'tokens');
    valid = ~cellfun(@isempty, tokens);
    rec_animal_names = repmat({''}, size(folder_names));
    rec_animal_names(valid) = cellfun(@(t) t{1}{1}, tokens(valid), 'UniformOutput', false);
    unique_animals = unique(rec_animal_names(valid), 'stable');

    if ~isempty(animals_filter)
        unique_animals = unique_animals(ismember(unique_animals, animals_filter));
    end

    all_rec_folders = animallist;
    animallist = [];
    for i = 1:numel(unique_animals)
        animallist(i).name = unique_animals{i};
        animallist(i).folder = parentdir;
    end
else
    all_rec_folders = [];
    rec_animal_names = {};
    if ~isempty(animals_filter)
        animallist = animallist(ismember(folder_names, animals_filter));
    end
end

fprintf('Found %d animal folder(s)\n', numel(animallist));

processed_count = 0;
skipped_count = 0;

%% Process each animal
for a = 1:numel(animallist)
    animal_name = animallist(a).name;
    fprintf('\n--- Processing animal: %s ---\n', animal_name);

    % Get config for this animal
    animal_config = get_user_config(animal_name);

    % Update paths for this animal
    DLCdir = animal_config.DLCDir;
    preprocesseddir = animal_config.OutputDir;

    % Get experiment directories - mode-aware
    if flat_mode
        % Recording folders sit directly in parentdir
        mask = strcmp(rec_animal_names, animal_name);
        explist = all_rec_folders(mask);
    else
        explist = dir(fullfile(animallist(a).folder, animallist(a).name));
        explist = explist([explist.isdir]);  % Fixed: removed space
        explist = explist(~ismember({explist.name}, {'.', '..'}));
    end

    if isempty(explist)
        if ~flat_mode
            fprintf('  No experiments found, removing empty folder\n');
            rmdir(fullfile(animallist(a).folder, animallist(a).name));
        end
        continue;
    end

    fprintf('  Found %d experiment(s)\n', numel(explist));

    %% Process each experiment
    for e = 1:numel(explist)
        exp_name = explist(e).name;
        exp_path = fullfile(explist(e).folder, exp_name);

        fprintf('    Checking: %s\n', exp_name);

        % Check for success.log (indicates DLC processing complete)
        if ~exist(fullfile(exp_path, 'success.log'), 'file')
            fprintf('      No success.log, skipping\n');
            skipped_count = skipped_count + 1;
            continue;
        end

        fprintf('      success.log found\n');

        % Check preprocessed file
        finalFile = fullfile(preprocesseddir, [exp_name '_preprocessed.mat']);
        preprocessed_exists = exist(finalFile, 'file') == 2;
        fprintf('      Preprocessed file exists: %d\n', preprocessed_exists);

        if preprocessed_exists
            DLCnotimported = isempty(who('-file', finalFile, 'VideoParams'));
        else
            DLCnotimported = true;
        end

        % Decide whether to process
        should_process = ~preprocessed_exists || ...
            ((preprocessed_exists && DLCnotimported) || params.overwrite);

        if ~should_process
            fprintf('      Already imported, skipping\n');
            skipped_count = skipped_count + 1;
            continue;
        end

        % Check for or create DLC filtered file
        dlc_filtered_file = fullfile(DLCdir, [exp_name '_DLC_filtered.mat']);
        if ~exist(dlc_filtered_file, 'file') || params.overwrite
            [success, DLCSupporting, Behav, VideoParams, VideoFrameTimes] = ...
                process_dlc_data(exp_path, exp_name, DLCdir, params);

            if ~success
                skipped_count = skipped_count + 1;
                continue;
            end
        else
            fprintf('      Loading existing DLC data\n');
            load(dlc_filtered_file, 'DLCSupporting', 'Behav', 'VideoParams', 'VideoFrameTimes');
        end

        % Add to preprocessed file if it exists
        if preprocessed_exists
            success = add_dlc_to_preprocessed(finalFile, dlc_filtered_file, ...
                DLCSupporting, Behav, VideoParams, VideoFrameTimes);

            if success
                fprintf('      Successfully aggregated!\n');
                processed_count = processed_count + 1;
            else
                skipped_count = skipped_count + 1;
            end
        else
            fprintf('      No preprocessed file found, skipping aggregation\n');
            skipped_count = skipped_count + 1;
        end

        % Delete camera TIFs if requested
        if params.delete_camera_tifs && exist(fullfile(exp_path, 'camera'), 'dir')
            rmdir(fullfile(exp_path, 'camera'), 's');
        end

        % Archive processed files
        archive_animal_dir = fullfile(archivedir, animal_name);
        archive_exp_dir = fullfile(archive_animal_dir, exp_name);

        if ~exist(archive_animal_dir, 'dir'); mkdir(archive_animal_dir); end
        if ~exist(archive_exp_dir, 'dir'); mkdir(archive_exp_dir); end

        stat = movefile(fullfile(exp_path, '*'), archive_exp_dir);
        if stat
            rmdir(exp_path);
            fprintf('      Archived to: %s\n', archive_exp_dir);
        end
    end
end

%% Summary
fprintf('\n=== Scan Complete ===\n');
fprintf('Processed: %d\n', processed_count);
fprintf('Skipped: %d\n', skipped_count);
fprintf('Exiting.\n');

end


%% ===== LOCAL FUNCTIONS =====

function [success, DLCSupporting, Behav, VideoParams, VideoFrameTimes] = ...
    process_dlc_data(exp_path, exp_name, DLCdir, params)
%PROCESS_DLC_DATA Process DLC CSV files and create filtered output

success = false;
DLCSupporting = struct();
Behav = struct();
VideoParams = struct();
VideoFrameTimes = [];

% Find eye CSV file
eyename = fullfile(exp_path, 'eye_results', [exp_name '_eye.csv']);
if ~exist(eyename, 'file')
    tempf = dir(fullfile(exp_path, 'eye_results', '*_eye.csv'));
    if isempty(tempf)
        fprintf('      Eye CSV not found, skipping\n');
        return;
    end
    eyename = fullfile(tempf.folder, tempf.name);
end

% Find body CSV file
bodyname = fullfile(exp_path, 'body_results', [exp_name '_body.csv']);
if ~exist(bodyname, 'file')
    tempf = dir(fullfile(exp_path, 'body_results', '*_body.csv'));
    if isempty(tempf)
        fprintf('      Body CSV not found, skipping\n');
        return;
    end
    bodyname = fullfile(tempf.folder, tempf.name);
end

fprintf('      Importing DLC data\n');

try
    [DLC, VideoFrameTimes] = import_DLC_eye(eyename, params);
    [DLCbody, DLCraw] = import_DLC_body(bodyname, params);
catch ME
    fprintf('      Error importing DLC: %s\n', ME.message);
    return;
end

% Populate structures - FIXED: removed extra spaces in property access
DLCraw.PupilPoints = DLC.Pupil.PointsRaw;

VideoParams.Eye.Radius_mm = DLC.eye_r_mm;
VideoParams.Eye.Center = DLC.Eye.center;
VideoParams.Eye.EllipseParams = DLC.Eye.A;
VideoParams.Eye.Semiaxes = DLC.Eye.semiaxes;
VideoParams.Eye.Rotation = DLC.Eye.phi;
VideoParams.Eye.medfilt = DLC.Eye.denoise.movmed_points;
VideoParams.Eye.CornealReflection = DLC.cornealReflection;
VideoParams.Eye.px_per_mm = DLC.eye_px_per_mm;
VideoParams.Pupil.p_thresh = DLC.Pupil.denoise.p_thresh;
VideoParams.Body.HeadplatRBC = DLCbody.headplateRBC;
VideoParams.Body.BallLoc_RAC = DLCbody.BallRAC;
VideoParams.Body.EyeLocFull = DLCbody.EyeLocFull;
VideoParams.Body.p_trhesh = DLCbody.prob_cutoff;
VideoParams.Body.medfilt = DLCbody.medfilt;

DLCSupporting.Probability.EyeLid = DLC.Eye.p;
DLCSupporting.Probability.Pupil = DLC.Pupil.p;
DLCSupporting.Probability.EllipseRMSE = DLC.Pupil.RMSE;
DLCSupporting.Probability.BallApex = DLCbody.ballApex(:, 3);
DLCSupporting.Probability.BallRostral = DLCbody.ballRostral(:, 3);
DLCSupporting.Probability.BallCaudal = DLCbody.ballCaudal(:, 3);

DLCSupporting.Pupil.Eccentricity = DLC.Pupil.eccentricity;
DLCSupporting.Pupil.EllipsePhi = DLC.Pupil.Rotation;
DLCSupporting.Pupil.EllipseSemiAxes = DLC.Pupil.semiaxes;
DLCSupporting.Pupil.EllipseCenter = DLC.Pupil.center;
DLCSupporting.Pupil.Points = DLC.Pupil.PointsRaw;
DLCSupporting.CropParams = DLC.CropParams;
DLCSupporting.CropParams.Body_vid_width = params.body_vid;
DLCSupporting.CropParams.Body_scale = DLCbody.bodyvid_scale;
DLCSupporting.Ball.Apex = DLCbody.ballApex(:, 1:2);
DLCSupporting.Ball.Rostral = DLCbody.ballRostral(:, 1:2);
DLCSupporting.Ball.Caudal = DLCbody.ballCaudal(:, 1:2);

Behav.VideoNum = DLC.Vid(:, 1);
Behav.FrameNum = DLC.Vid(:, 2);
Behav.PupilAzEl = DLC.Pupil.az_el;
Behav.PupilArea = DLC.Pupil.area;

FN = fieldnames(DLCbody);
% Use actual field count (robust to different DLC versions)
num_fields = min(19, numel(FN));
for fn = 1:num_fields
    Behav.(FN{fn}) = DLCbody.(FN{fn})(:, 1:2);
    DLCSupporting.Probability.(FN{fn}) = DLCbody.(FN{fn})(:, 3);
end

Behav.EyeLidArea = DLC.Eye.area;
Behav.EyeLidCenter = DLC.Eye.position;

% Save DLC files
fprintf('      Saving DLC files\n');
save(fullfile(DLCdir, [exp_name '_DLC_raw.mat']), 'DLCraw');
save(fullfile(DLCdir, [exp_name '_DLC_filtered.mat']), ...
    'DLCSupporting', 'Behav', 'VideoParams', 'VideoFrameTimes');

success = true;

end


function success = add_dlc_to_preprocessed(finalFile, dlc_file, ...
    DLCSupporting, Behav, VideoParams, VideoFrameTimes)
%ADD_DLC_TO_PREPROCESSED Add DLC data to existing preprocessed file

success = false;

try
    fprintf('      Adding DLC to preprocessed file\n');

    Data = matfile(finalFile, 'Writable', true);
    Triggers = Data.Triggers;

    if isfield(VideoParams, 'selected_frameIdx')
        selected_frameIdx = VideoParams.selected_frameIdx;
    else
        [selected_frameIdx, Triggers] = correct_videoFrames_OS(Triggers, VideoFrameTimes, Behav);
        Data.Triggers = Triggers;
        VideoParams.selected_frameIdx = selected_frameIdx;
        save(dlc_file, 'VideoParams', '-append');
    end

    % Apply frame selection
    FN = fieldnames(DLCSupporting.Probability);
    for fn = 1:numel(FN)
        DLCSupporting.Probability.(FN{fn}) = DLCSupporting.Probability.(FN{fn})(selected_frameIdx, :, :);
    end

    FN = fieldnames(DLCSupporting.Pupil);
    for fn = 1:numel(FN)
        DLCSupporting.Pupil.(FN{fn}) = DLCSupporting.Pupil.(FN{fn})(selected_frameIdx, :, :);
    end

    FN = fieldnames(DLCSupporting.Ball);
    for fn = 1:numel(FN)
        DLCSupporting.Ball.(FN{fn}) = DLCSupporting.Ball.(FN{fn})(selected_frameIdx, :, :);
    end

    FN = fieldnames(Behav);
    for fn = 1:numel(FN)
        Behav.(FN{fn}) = Behav.(FN{fn})(selected_frameIdx, :, :);
    end

    Data.DLCSupporting = DLCSupporting;
    Data.Behav = Behav;
    Data.VideoParams = VideoParams;

    success = true;

catch ME
    fprintf('      Error aggregating: %s\n', ME.message);
end

end


function [DLC, VideoFrameTimes] = import_DLC_eye(filename, params)
%IMPORT_DLC_EYE Import eye tracking data from DLC CSV

delimiter = ',';

fileID = fopen(filename, 'r');
if fileID == -1
    error('Cannot open eye CSV file: %s', filename);
end

% Read general info (row 2)
startRow = 2;
endRow = 2;
formatSpec = '%C%f%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, ...
    'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, ...
    'ReturnOnError', false, 'EndOfLine', '\r\n');

output_general = table2struct(table(dataArray{1:end-1}, 'VariableNames', ...
    {'Experiment_Name', 'eye_center_x', 'eye_center_y', 'eye_width', 'eye_height', ...
    'eye_theta', 'pupil_jackknife_hold', 'eye_R', 'eye_area', 'crop_size', ...
    'original_width', 'original_height', 'original_eye_center_x', 'original_eye_center_y'}));
fclose(fileID);

% Read frame data (rows 4+)
fileID = fopen(filename, 'r');
startRow = 4;
endRow = inf;

formatSpec = '%s%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%C%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, ...
    'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, ...
    'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

output_frame = table2struct(table(dataArray{1:end-1}, 'VariableNames', ...
    {'videoref', 'frame', 'tiff_name', 'time_stamp', 'eye3D_x', 'eye3D_y', 'eye3D_z', ...
    'eye_width_2', 'eye_heigth', 'eye_theta', 'eye_area', 'pupil3d_x', 'pupil3d_y', ...
    'pupil3d_z', 'pupil_width', 'pupil_height', 'pupil_theta', 'pupil_el', 'pupil_az', ...
    'pupil_std_x', 'pupil_std_y', 'jackknife_inds', 'pupilN', 'pupilN1', 'pupilN2', ...
    'pupilNE', 'pupilNE1', 'pupilNE2', 'pupilE', 'pupilE1', 'pupilE2', 'pupilSE', ...
    'pupilSE1', 'pupilSE2', 'pupilS', 'pupilS1', 'pupilS2', 'pupilSW', 'pupilSW1', ...
    'pupilSW2', 'pupilW', 'pupilW1', 'pupilW2', 'pupilNW', 'pupilNW1', 'pupilNW2', ...
    'cornealReflection', 'cornealReflection1', 'cornealReflection2', 'eyeN', 'eyeN1', ...
    'eyeN2', 'eyeNE', 'eyeNE1', 'eyeNE2', 'eyeE', 'eyeE1', 'eyeE2', 'eyeSE', 'eyeSE1', ...
    'eyeSE2', 'eyeS', 'eyeS1', 'eyeS2', 'eyeSW', 'eyeSW1', 'eyeSW2', 'eyeW', 'eyeW1', ...
    'eyeW2', 'eyeNW', 'eyeNW1', 'eyeNW2'}));

% Process data
temp = struct2cell(output_frame);
frame_temp = temp(1, :)';
DLC.Vid = cat(2, str2double(cellfun(@(x) x{1}(end-4:end), frame_temp, 'UniformOutput', false)), ...
    cat(1, output_frame(:).frame));

DLC.Eye.denoise.movmed_points = params.movmed_filt;
DLC.Eye.center = [output_general.eye_center_x, output_general.eye_center_y];
DLC.Eye.semiaxes = [output_general.eye_width, output_general.eye_height] / 2;
DLC.Eye.phi = output_general.eye_theta;
DLC.Eye.area = cat(1, output_frame(:).eye_area);
DLC.Eye.position = cat(2, cat(1, output_frame(:).eye3D_x), cat(1, output_frame(:).eye3D_y));
DLC.Eye.p = mean(cell2mat(temp(end-21:3:end, :)'), 2);
DLC.cornealReflection = median(cell2mat(temp(end-26:end-25, :)'), 1);

DLC.Pupil.center = cat(2, cat(1, output_frame(:).pupil3d_x), cat(1, output_frame(:).pupil3d_y));
DLC.Pupil.semiaxes = cat(2, cat(1, output_frame(:).pupil_width), cat(1, output_frame(:).pupil_height)) / 2;
DLC.Pupil.Rotation = cat(1, output_frame(:).pupil_theta);
DLC.Pupil.RMSE = mean(cat(2, cat(1, output_frame(:).pupil_std_x), cat(1, output_frame(:).pupil_std_y)), 2);
DLC.Pupil.p = cell2mat(temp(25:3:46, :)');
DLC.Pupil.denoise.p_thresh = params.pupil_min_p;
DLC.Pupil.denoise.movmed_points = params.movmed_filt;
DLC.Pupil.PointsRaw = reshape(cell2mat(temp(23:46, :)'), size(temp, 2), 3, 8);
DLC.Pupil.eccentricity = sqrt(1 - (DLC.Pupil.semiaxes(:, 1).^2) ./ (DLC.Pupil.semiaxes(:, 2).^2));
DLC.Pupil.area = prod(DLC.Pupil.semiaxes, 2) * pi;

DLC.eye_r_mm = params.eye_r_mm;
DLC.eye_px_per_mm = params.eye_px_per_mm;

DLC.Pupil.az_el = cat(2, cat(1, output_frame(:).pupil_az), cat(1, output_frame(:).pupil_el));

DLC.CropParams.Size = output_general.crop_size;
DLC.CropParams.original_width = output_general.original_width;
DLC.CropParams.original_height = output_general.original_height;
DLC.CropParams.original_eye_center_x = output_general.original_eye_center_x;
DLC.CropParams.original_eye_center_y = output_general.original_eye_center_y;
DLC.Eye.A = nan(1, 8);

% Parse video frame times
% BUG FIX: Removed extra space before SSS in datetime format
% Original RB: 'yyyy_MM_dd HH:mm:ss: SSS' (fails)
% Fixed:       'yyyy_MM_dd HH:mm:ss:SSS'  (correct)
videoframetimestemp = cellfun(@(x) x{1}, temp(4, :), 'UniformOutput', false)';
if str2double(videoframetimestemp{1}(1:4)) == 0
    VideoFrameTimes = nan;
else
    VideoFrameTimes = datetime(videoframetimestemp, 'InputFormat', 'yyyy_MM_dd HH:mm:ss:SSS');
end

end


function [DLC, DLCraw] = import_DLC_body(filename, params)
%IMPORT_DLC_BODY Import body tracking data from DLC CSV

delimiter = ',';
startRow = 4;
endRow = inf;

formatSpec = '%C%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename, 'r');
if fileID == -1
    error('Cannot open body CSV file: %s', filename);
end

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, ...
    'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, ...
    'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

dataArray = dataArray(1:end-1);
NameOut = import_body_part_name(filename)';
[NameUn, ~, ib] = unique(NameOut, 'stable');

for n = 1:numel(NameUn)
    DLC.(NameUn{n}) = cat(2, dataArray{ib == n});
end

DLC.headplateRBC = [median(DLC.headplateRostralCorner(:, 1:2), 1); ...
    median(DLC.headplateBendCorner(:, 1:2), 1); nan(1, 2)];
DLC.BallRAC = [median(DLC.ballRostral(:, 1:2), 1); ...
    median(DLC.ballApex(:, 1:2), 1); ...
    median(DLC.ballCaudal(:, 1:2), 1)];
DLC.EyeLocFull = median(DLC.eye(:, 1:2));

FN = fieldnames(DLC);
DLC = rmfield(DLC, FN([1:4, 8:9]));
FN = fieldnames(DLC);
DLCraw = DLC;
DLC.prob_cutoff = params.body_min_p;
DLC.medfilt = params.movmed_filt;
DLC.bodyvid_scale = import_DLC_eyescale(filename);

fprintf('      Filtering body data\n');

% Use actual field count (robust to different DLC versions)
num_fields = min(19, numel(FN));
for fn = 1:num_fields
    DLC.(FN{fn})(DLC.(FN{fn})(:, end) < DLC.prob_cutoff, 1:end-1) = nan;
    DLC.(FN{fn})(:, 1:end-1) = fillmissing(DLC.(FN{fn})(:, 1:end-1), 'linear');
    DLC.(FN{fn})(:, 1:end-1) = movmedian(DLC.(FN{fn})(:, 1:end-1), DLC.medfilt);
    DLC.(FN{fn})(DLC.(FN{fn})(:, end) < DLC.prob_cutoff, 1:end-1) = nan;
end

end


function NameOut = import_body_part_name(filename)
%IMPORT_BODY_PART_NAME Import body part names from DLC CSV header

delimiter = ',';
startRow = 3;
endRow = 3;
formatSpec = repmat('%q', 1, 100);
formatSpec = [formatSpec, '%[^\n\r]'];

fileID = fopen(filename, 'r');
if fileID == -1
    error('Cannot open file for body part names: %s', filename);
end

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, ...
    'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, ...
    'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

nonemptystrings = cellfun(@(x) strlength(x), dataArray) > 0;
NameOut = [dataArray{nonemptystrings}];

end


function scale = import_DLC_eyescale(filename)
%IMPORT_DLC_EYESCALE Import eye scale factor from DLC body CSV

startRow = 2;
endRow = 2;

delimiter = ',';
formatSpec = '%*s%*s%*s%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

fileID = fopen(filename, 'r');
if fileID == -1
    error('Cannot open file for eye scale: %s', filename);
end

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, ...
    'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, ...
    'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

scale = dataArray{1};

end
