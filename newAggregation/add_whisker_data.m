function [Behav, VideoParams] = add_whisker_data(session_dir, exp_name, Behav, VideoParams)
%ADD_WHISKER_DATA Add whisking on/off, amplitude and jaw / face motion from <session_dir>\whisker_results\
%
%   [Behav, VideoParams] = add_whisker_data(session_dir, exp_name, Behav, VideoParams)
%
%   Reads whisker_results\<exp_name>_whisker.csv (+ .json) written by DLC\face_model\face_motion.py and adds them to
%   Behav / VideoParams.Whisker. Every trace is an N-by-1 column (N = video frames) in Behav, so the generic row
%   selection by selected_frameIdx in save_preprocessed / add_dlc_to_preprocessed keeps it aligned with the eye and
%   body traces. (A field in DLCSupporting would NOT be selected: only Probability/Pupil/Ball are.)
%   Nothing is changed (and a warning is given where the frame counts differ) when the csv is missing.
%
%   Behav fields added (definitions: import_video_behavior.m, section "Whisking and face motion"):
%     WhiskIndex, WhiskAmp, WhiskOn, WhiskArtifact   whisking, from the upper whisker region (front_top, above the spout)
%     JawME                                          motion energy of the jaw / lick region (the former fan_down box)
%     FacePadME                                      motion energy of the whisker-pad skin
%   VideoParams.Whisker = contents of the json (codec, whisker_box, rest_level, ...).
%
%   Used by import_video_behavior (new sessions) and refresh_whisker_behavior (sessions that are already aggregated).

wfile = fullfile(session_dir, 'whisker_results', [exp_name '_whisker.csv']);
if ~exist(wfile, 'file')
    fprintf('      No whisker_results found, skipping whisker import\n');
    return;
end

T = readtable(wfile);
if height(T) ~= numel(Behav.FrameNum)
    warning('add_whisker_data:whiskerFrames', ...
        'Whisker csv has %d frames but the DLC csv has %d; whisker data NOT imported', ...
        height(T), numel(Behav.FrameNum));
    return;
end

Behav.WhiskIndex    = T.whisk_index;   % 1 = whiskers at rest (rest per video codec), upper whisker region only
Behav.WhiskAmp      = T.whisk_amp;     % index - 1; NaN on artifact frames
Behav.WhiskOn       = T.whisk_on;      % 0/1; NaN on artifact frames
Behav.WhiskArtifact = T.artifact;      % 1 = grooming etc.
Behav.FacePadME     = T.me_pad;        % whisker-pad (face) motion energy
if ismember('me_fan_down', T.Properties.VariableNames)
    Behav.JawME = T.me_fan_down;       % jaw / lick motion energy (chin, lower jaw, tongue)
end

jfile = fullfile(session_dir, 'whisker_results', [exp_name '_whisker.json']);
if exist(jfile, 'file')
    VideoParams.Whisker = jsondecode(fileread(jfile));
end
fprintf('      Imported whisker data (%d frames)\n', height(T));

end
