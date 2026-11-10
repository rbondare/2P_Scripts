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
%   Whisking (whisker_results\<recording>_whisker.csv + .json from DLC\face_model\face_motion.py, written by
%   run_dlc.bat next to eye_results/body_results) is added when present; recordings without whisker_results
%   are imported as before. The pupil mask (Behav.PupilBad) is computed here as well. Everything is read from
%   the session folder, which is archived right after the import, so this is the only place to add it.
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
%       'sessions' - Cell array of recording names (e.g. {'AnimalRB24_261030_1744'}); all other
%           session folders are left untouched (default: all sessions of the animals)
%
%       'overwrite' - Overwrite existing DLC files (default: false)
%
%       'delete_camera_tifs' - Delete camera TIFs after processing (default: false)
%
%   Parameters (can be customized):
%       movmed_filt      - Median filter window in frames for the body traces (default: 3)
%       pupil_min_p      - Pupil point likelihood below which the pupil is masked (default: 0.8)
%       lid_open_frac    - Eyelid opening below this fraction of the recording's median masks the pupil
%                          (default: 0.75)
%       pupil_pad_frames - Masked pupil stretches are widened by this many frames on both sides
%                          (default: 5 = 0.1 s)
%       body_min_p       - Minimum body probability threshold (default: 0.8)
%       eye_r_mm         - Eye radius in mm (default: 1.5)
%       eye_px_per_mm    - Pixels per mm for eye tracking (default: 23.4)
%       body_vid         - Body video width (default: 480)
%
%   ====================================================================================================
%   OUTPUT VARIABLES
%   ====================================================================================================
%   Saved to <DLCDir>\<recording>_DLC_filtered.mat as Behav, DLCSupporting, VideoParams and VideoFrameTimes
%   (unfiltered body points and pupil points: <recording>_DLC_raw.mat as DLCraw). save_preprocessed then
%   copies Behav, DLCSupporting and VideoParams into the aggregated file.
%
%   N = number of video frames (50 fps). Row k of every Behav and DLCSupporting trace is the same video
%   frame. In the aggregated file only the frames matched to the 2P camera triggers are kept
%   (VideoParams.selected_frameIdx), so there N = numel(Triggers.TimeCamera).
%   Positions are image coordinates in pixels: x to the right, y downward. NaN = no valid value on that frame.
%   Eye variables are in pixels of the eye crop (DLCSupporting.CropParams.Size, 150 px); body variables are in
%   pixels of the body video (480 px wide = CropParams.Body_scale x the 1856 x 1020 camera frame).
%
%   ---- Behav : one row per video frame --------------------------------------------------------------
%   Frame bookkeeping
%     VideoNum         Nx1  Number of the video file the frame belongs to (0 when the session is one mp4).
%     FrameNum         Nx1  Frame number inside that video, counting from 0.
%
%   Pupil and eye (eye camera)
%     PupilAzEl        Nx2  [azimuth elevation] of the pupil in degrees; [0 0] = pupil centred on the eye
%                           centre. Azimuth > 0: pupil displaced towards +x of the eye image; elevation > 0:
%                           towards -y (up in the image). NaN on PupilBad frames.
%     PupilArea        Nx1  Area of the fitted pupil ellipse in pixels^2 (divide by
%                           VideoParams.Eye.px_per_mm^2 for mm^2). NaN on PupilBad frames.
%     PupilBad         Nx1  logical. true = the pupil fit cannot be trusted on this frame (rule below).
%     PupilAzElRaw     Nx2  PupilAzEl before masking (equal to PupilAzEl wherever PupilBad is false). It is
%                           NaN only on the few frames where the eye csv itself has no value (PupilAzEl too).
%     PupilAreaRaw     Nx1  PupilArea before masking.
%     EyeLidOpening    Nx1  Distance between the lower and the upper eyelid point (eyeS_y - eyeN_y), pixels.
%                           Larger = eye more open, near 0 = shut, negative = the model swapped the two points
%                           (a bad frame). Never masked.
%     EyeLidArea       Nx1  Area of the fitted eye-opening ellipse in pixels^2; shrinks when the eye closes.
%     EyeLidCenter     Nx2  [x y] centre of the fitted eye (eye3D_x and eye3D_y of the eye csv), pixels.
%
%   PupilBad is true on a frame when
%     (a) any of its 8 pupil points has a likelihood < pupil_min_p (DLCSupporting.Probability.Pupil), or
%     (b) EyeLidOpening < lid_open_frac x the median EyeLidOpening of the whole recording (blink, eye nearly
%         shut; the median is VideoParams.Pupil.lid_ref_px), or
%     (c) a pupil likelihood or the eyelid opening is missing (NaN),
%   and every flagged stretch is then widened by pupil_pad_frames frames on both sides. Long stretches with the
%   eye shut are flagged as well: they are an eye state rather than tracking noise, so look at
%   VideoParams.Pupil.frac_bad before using a session. Recompute (a) and (b) separately from the variables
%   above to tell the causes apart.
%
%   Body (body camera): [x y] in body-video pixels, one Nx2 field per tracked body part. These are the first 19
%   body parts of the body csv in csv order; the two headplate corners are not included (their medians are in
%   VideoParams.Body.HeadplatRBC).
%     ballRostral, ballApex, ballCaudal    rostral (head side), top and caudal (tail side) point of the ball
%     noseTip, whiskerFollicleA1, eye      nose tip, whisker follicle A1, eye
%     PinnaRightBase, PinnaRightTip, PinnaRightBottom                right ear
%     shoulderApex, hipApex                shoulder and hip
%     foreLimbRightTip, foreLimbRightWrist, foreLimbRightElbow       right forelimb
%     hindLimbRightTip, hindLimbRightAnkle, hindLimbRightKnee        right hindlimb
%     foreLimbLeftTip, foreLimbLeftWrist                             left forelimb
%   NaN where the point's likelihood is < body_min_p (DLCSupporting.Probability.<same name>); all other frames
%   are median filtered over movmed_filt frames.
%
%   Whisking and face motion (only if whisker_results\<recording>_whisker.csv exists; made by face_motion.py
%   from the raw face video, the method is described at the top of that file)
%     WhiskIndex       Nx1  Whisker motion relative to rest, unitless: motion energy of the UPPER whisker region
%                           (front_top, above the lick spout, so licking does not drive it) divided by its
%                           resting level. 1 = whiskers at rest, larger = faster whisker movement (it follows the
%                           speed of the whiskers, not their angle, and cannot tell whisker movement from head
%                           movement). The resting level is set per video codec (MPEG-4 up to 2026-10-29, H.264
%                           from 2026-10-30; the H.264 files keep much more sensor noise), see
%                           VideoParams.Whisker.rest_level / .codec. Always a number, but not meaningful on
%                           WhiskArtifact frames. Before 2026-11-10 this was the mean of two boxes (one of them
%                           the jaw): lick-driven and on a different scale from 2026-10-30 on; VideoParams.Whisker
%                           has the field whisker_box ('front_top') only in the new version.
%     WhiskAmp         Nx1  WhiskIndex - 1, floored at 0: 0 = at rest, graded amount of whisking.
%                           NaN on WhiskArtifact frames.
%     WhiskOn          Nx1  Whisking state: 1 = whisking, 0 = not. Hysteresis on WhiskIndex, levels per codec: MPEG-4
%                           on above 1.76 / off below 1.44 (1.6 x rest +-10%), H.264 on above 1.45 / off below 1.40
%                           (VideoParams.Whisker.whisk_on_above_index / whisk_off_below_index; since 2026-11-10 the
%                           H.264 levels were set from 51 blind human labels, with 1.76 the state missed sustained
%                           moderate whisking; before that H.264 sessions also used 1.76 / 1.44). Pauses < 0.3 s inside
%                           a bout are closed, bouts < 0.3 s are dropped. NaN on WhiskArtifact frames. A binary state
%                           is a convention: WhiskIndex / WhiskAmp carry the graded information.
%     WhiskArtifact    Nx1  1 = motion that is not whisking (grooming, a forepaw in front of the face):
%                           WhiskIndex, or the whole-face motion (the former two-box index of the upper whisker
%                           region and the jaw region), above 2 x its own 95th percentile in the recording,
%                           widened by 0.5 s on both sides. 0 otherwise. Grooming shows mostly in the jaw / pad
%                           boxes, hence the second criterion.
%     JawME            Nx1  Motion energy of the jaw / lick region (chin, lower jaw, tongue and the shafts hanging
%                           below the pad; the former second whisker region): mean |frame(t) - frame(t-1)| inside
%                           the region, gray levels per frame, not normalized. Licking and chewing dominate it
%                           (45% of the variance of the old index was lick-driven, 87% of that through this
%                           region); use it as the movement covariate for lick-related analyses. Added 2026-11-10.
%     FacePadME        Nx1  Motion energy of the whisker-pad skin: mean |frame(t) - frame(t-1)| inside the
%                           pad region, in gray levels per frame. Also moves with chewing and breathing.
%                           Not normalized to rest.
%
%   ---- VideoParams : numbers that describe the whole recording ----------------------------------------
%     Eye.Radius_mm          1.5 = eyeball radius of the 3D eye model, mm (parameter eye_r_mm)
%     Eye.px_per_mm          23.4 = eye-crop pixels per mm (parameter eye_px_per_mm)
%     Eye.Center             [x y] centre of the eye ellipse of the recording (header row of the eye csv), px
%     Eye.Semiaxes           [a b] half width and half height of that ellipse, px
%     Eye.Rotation           rotation of that ellipse, rad
%     Eye.CornealReflection  [x y] median position of the corneal reflection, px
%     Eye.medfilt            movmed_filt
%     Eye.EllipseParams      NaN(1,8) placeholder, unused
%     Pupil.p_thresh         pupil_min_p that was used for PupilBad
%     Pupil.lid_frac         lid_open_frac that was used for PupilBad
%     Pupil.lid_ref_px       median EyeLidOpening of the whole recording, px (= 100% of lid_frac)
%     Pupil.pad_frames       pupil_pad_frames that was used for PupilBad
%     Pupil.frac_bad         fraction of the N video frames with PupilBad = true (before frame selection)
%     Body.HeadplatRBC       3x2 [x y] median headplate rostral corner, bend corner, then a NaN row
%     Body.BallLoc_RAC       3x2 [x y] median ball rostral, apex and caudal point
%     Body.EyeLocFull        [x y] median position of the 'eye' body point
%     Body.p_trhesh          body_min_p (the name is misspelled, kept for compatibility)
%     Body.medfilt           movmed_filt
%     Whisker                contents of whisker_results\<recording>_whisker.json (absent without whisker data):
%         session, frames, fps               recording name, number of video frames, frame rate
%         fraction_whisking                  fraction of the non-artifact frames with WhiskOn = 1
%         fraction_artifact                  fraction of all frames with WhiskArtifact = 1
%         n_bouts, median_bout_s             number and median duration (s) of whisking bouts
%         index_p5, index_median, index_p95  percentiles of WhiskIndex over the non-artifact frames
%         codec, whisker_box                 video codec ('mpeg4' / 'h264') and the region the index is made from
%                                            ('front_top'; the field is absent in the version before 2026-11-10)
%         rest_level                         resting motion energy of that region for this codec, a struct with
%                                            one field named like the region (= WhiskIndex 1)
%         whisk_off_below_index, whisk_on_above_index   the WhiskOn hysteresis levels on WhiskIndex (field names of
%                                            the 2026-11-10 H.264 version; results saved before have threshold_x_rest
%                                            = 1.6 and params.hysteresis = 0.1, i.e. 1.44 / 1.76, which is what
%                                            every MPEG-4 session still uses)
%         ctrl_me_median                     median motion energy of the background region (noise level)
%         anchor_median_xy_px, anchor_range_px   median nose position the regions are placed on, and how far
%                                            it moved over the recording, full-frame px
%         boxes_rel_to_nose_px               the five regions [dx0 dy0 dx1 dy1] relative to the nose, px
%         face_track_frames, params          frames of the face-model track used; smoothing/artifact settings
%     selected_frameIdx      added later by save_preprocessed: the video frames kept to match the triggers
%
%   ---- DLCSupporting : tracking quality and the rest of the pupil fit (not masked by PupilBad) ---------
%     Probability.EyeLid        Nx1   mean likelihood of the 8 eyelid points
%     Probability.Pupil         Nx8   likelihood of each pupil point, order N NE E SE S SW W NW (1 = certain)
%     Probability.EllipseRMSE   Nx1   spread of the pupil ellipse fit (mean of pupil_std_x and pupil_std_y of
%                                     the csv), pixels. Not a probability; lower = better
%     Probability.BallApex, .BallRostral, .BallCaudal, and one field per body part (same name as in Behav)
%                               Nx1   likelihood of that point
%     Pupil.EllipseSemiAxes     Nx2   [width/2 height/2] of the fitted pupil ellipse, pixels
%     Pupil.EllipseCenter       Nx2   [x y] centre of the fitted pupil ellipse, pixels
%     Pupil.EllipsePhi          Nx1   rotation of the pupil ellipse, rad
%     Pupil.Eccentricity        Nx1   sqrt(1 - (width/height)^2). COMPLEX whenever width > height, which is
%                                     most frames (96% in RB24_261023_1453): do not use it
%     Pupil.Points              Nx3x8 [x y likelihood] of the 8 pupil points (same order), pixels
%     Ball.Apex, .Rostral, .Caudal  Nx2   [x y] of the ball points, body-video pixels (as in Behav)
%     CropParams                Size (eye crop, px); original_width, original_height (camera frame);
%                               original_eye_center_x, _y (centre of the eye crop in the camera frame);
%                               Body_vid_width (480); Body_scale (body video / camera frame)
%
%   ---- VideoFrameTimes -----------------------------------------------------------------------------
%     NaN for mp4 recordings (no per-frame clock; frames are matched to the camera triggers by frame count),
%     otherwise an Nx1 datetime with the clock time of every video frame.
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
addParameter(p, 'sessions', {}, @iscell);
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'delete_camera_tifs', false, @islogical);
addParameter(p, 'movmed_filt', 3, @isnumeric);
addParameter(p, 'pupil_min_p', 0.8, @isnumeric);
addParameter(p, 'lid_open_frac', 0.75, @isnumeric);
addParameter(p, 'pupil_pad_frames', 5, @isnumeric);
addParameter(p, 'body_min_p', 0.8, @isnumeric);
addParameter(p, 'eye_r_mm', 1.5, @isnumeric);
addParameter(p, 'eye_px_per_mm', 23.4, @isnumeric);
addParameter(p, 'body_vid', 480, @isnumeric);
parse(p, varargin{:});

config = p.Results.config;
animals_filter = p.Results.animals;
sessions_filter = p.Results.sessions;
overwrite = p.Results.overwrite;
delete_camera_tifs = p.Results.delete_camera_tifs;

% Create params structure
params.movmed_filt = p.Results.movmed_filt;
params.pupil_min_p = p.Results.pupil_min_p;
params.lid_open_frac = p.Results.lid_open_frac;
params.pupil_pad_frames = p.Results.pupil_pad_frames;
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

        if ~isempty(sessions_filter) && ~ismember(exp_name, sessions_filter)
            continue;
        end

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
[pupil_bad, lid_ref_px] = pupil_bad_frames(DLC.Pupil.p, DLC.Eye.opening, params);
VideoParams.Pupil.p_thresh = DLC.Pupil.denoise.p_thresh;
VideoParams.Pupil.lid_frac = params.lid_open_frac;
VideoParams.Pupil.lid_ref_px = lid_ref_px;
VideoParams.Pupil.pad_frames = params.pupil_pad_frames;
VideoParams.Pupil.frac_bad = mean(pupil_bad);
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
Behav.PupilBad = pupil_bad;
Behav.PupilAzElRaw = Behav.PupilAzEl;
Behav.PupilAreaRaw = Behav.PupilArea;
Behav.PupilAzEl(pupil_bad, :) = NaN;
Behav.PupilArea(pupil_bad) = NaN;

FN = fieldnames(DLCbody);
% Use actual field count (robust to different DLC versions)
num_fields = min(19, numel(FN));
for fn = 1:num_fields
    Behav.(FN{fn}) = DLCbody.(FN{fn})(:, 1:2);
    DLCSupporting.Probability.(FN{fn}) = DLCbody.(FN{fn})(:, 3);
end

Behav.EyeLidArea = DLC.Eye.area;
Behav.EyeLidCenter = DLC.Eye.position;
Behav.EyeLidOpening = DLC.Eye.opening;

[Behav, VideoParams] = add_whisker_data(exp_path, exp_name, Behav, VideoParams);

% Save DLC files
fprintf('      Saving DLC files\n');
save(fullfile(DLCdir, [exp_name '_DLC_raw.mat']), 'DLCraw');
save(fullfile(DLCdir, [exp_name '_DLC_filtered.mat']), ...
    'DLCSupporting', 'Behav', 'VideoParams', 'VideoFrameTimes');

success = true;

end


function [bad, lid_ref_px] = pupil_bad_frames(pupil_p, lid_opening, params)
%PUPIL_BAD_FRAMES Frames on which the pupil fit cannot be trusted (rule: PupilBad in the file header)
%   pupil_p      N-by-8 likelihoods of the pupil points
%   lid_opening  N-by-1 lower minus upper eyelid point (y), pixels
%   bad          N-by-1 logical
%   lid_ref_px   median lid opening of the whole recording, pixels

lid_ref_px = median(lid_opening, 'omitnan');
bad = ~all(pupil_p >= params.pupil_min_p, 2);   % low likelihood, or NaN
if lid_ref_px > 0
    bad = bad | ~(lid_opening >= params.lid_open_frac * lid_ref_px);   % eye nearly shut, or NaN
else
    warning('import_video_behavior:lidOpening', ...
        'Median eyelid opening is %g px; the eyelid criterion of PupilBad is not applied', lid_ref_px);
end
bad = movmax(double(bad), 2 * params.pupil_pad_frames + 1) > 0;

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
DLC.Eye.opening = cat(1, output_frame(:).eyeS1) - cat(1, output_frame(:).eyeN1);  % lower minus upper eyelid point (y)
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
