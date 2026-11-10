function refresh_whisker_behavior(sessions, varargin)
%REFRESH_WHISKER_BEHAVIOR Replace the whisking / jaw-motion traces of sessions that are already aggregated
%
%   refresh_whisker_behavior({'AnimalRB24_261030_1744', ...})
%   refresh_whisker_behavior(sessions, 'whisker_root', 'C:\...\DLC\work', 'dry_run', true)
%
%   For every session: reads the (new) whisker_results\<session>_whisker.csv/.json with add_whisker_data, writes the
%   whisking fields (WhiskIndex, WhiskAmp, WhiskOn, WhiskArtifact, FacePadME, JawME) into
%       <DLCDir>\<session>_DLC_filtered.mat   (Behav, VideoParams; all video frames)
%       <OutputDir>\<session>_preprocessed.mat (Behav rows selected by VideoParams.selected_frameIdx, VideoParams.Whisker)
%   and leaves every other variable untouched. The old whisking fields are saved first to
%   <backup_dir>\<session>_whisk_v1.mat (struct 'old': aggregated Behav fields + VideoParams.Whisker, and the
%   full-length DLC_filtered Behav fields), so the previous version can be restored.
%
%   Name-Value arguments:
%       'whisker_root' - folder that holds <session>\whisker_results (default: the archived Camera folder of the
%                        session, <rima_root>\preCamera\<animal>\<session>, then <CameraDir>\<session>)
%       'backup_dir'   - default C:\Users\rbondarenko\2P_cache\backup_whisk_v1
%       'dry_run'      - only check and print what would be changed (default false)
%
%   A session is skipped with a message (the others go on) when the aggregated file, the DLC_filtered file or the new
%   whisker csv is missing, the frame counts differ, or the new csv is not the front_top version.
%
%   See also: add_whisker_data, import_video_behavior

p = inputParser;
addParameter(p, 'whisker_root', '', @ischar);
addParameter(p, 'backup_dir', 'C:\Users\rbondarenko\2P_cache\backup_whisk_v1', @ischar);
addParameter(p, 'dry_run', false, @islogical);
parse(p, varargin{:});
whisker_root = p.Results.whisker_root;
backup_dir = p.Results.backup_dir;
dry_run = p.Results.dry_run;
if ischar(sessions) || isstring(sessions); sessions = cellstr(sessions); end
if ~exist(backup_dir, 'dir') && ~dry_run; mkdir(backup_dir); end

fields = {'WhiskIndex', 'WhiskAmp', 'WhiskOn', 'WhiskArtifact', 'FacePadME', 'JawME'};
nok = 0;
failed = {};
for s = 1:numel(sessions)
    session = sessions{s};
    fprintf('--- %s\n', session);
    try
        refresh_one(session, whisker_root, backup_dir, dry_run, fields);
        nok = nok + 1;
    catch ME
        fprintf('    SKIPPED: %s\n', ME.message);
        failed{end + 1} = session; %#ok<AGROW>
    end
end
fprintf('\nrefreshed %d of %d sessions%s\n', nok, numel(sessions), ternary(dry_run, ' (dry run, nothing written)', ''));
if ~isempty(failed)
    fprintf('not refreshed: %s\n', strjoin(failed, ', '));
end

end


function refresh_one(session, whisker_root, backup_dir, dry_run, fields)

tok = regexp(session, '^(.+?)_\d{6}_\d{4}$', 'tokens', 'once');
if isempty(tok); error('cannot read the animal name from %s', session); end
animal = tok{1};
config = get_user_config(animal);
finalFile = fullfile(config.OutputDir, [session '_preprocessed.mat']);
dlcFile = fullfile(config.DLCDir, [session '_DLC_filtered.mat']);
if ~exist(finalFile, 'file'); error('no aggregated file %s', finalFile); end
if ~exist(dlcFile, 'file'); error('no DLC_filtered file %s', dlcFile); end

% where the new whisker_results are
cands = {};
if ~isempty(whisker_root); cands{end + 1} = fullfile(whisker_root, session); end
cands{end + 1} = fullfile(fileparts(config.CameraDir(1:end-1)), 'preCamera', animal, session);
cands{end + 1} = fullfile(config.CameraDir, session);
sdir = '';
for c = 1:numel(cands)
    if exist(fullfile(cands{c}, 'whisker_results', [session '_whisker.csv']), 'file')
        sdir = cands{c};
        break;
    end
end
if isempty(sdir); error('no whisker_results\%s_whisker.csv in %s', session, strjoin(cands, ' | ')); end
fprintf('    new whisker data from %s\n', sdir);

D = load(dlcFile, 'Behav', 'VideoParams');
if ~isfield(D.VideoParams, 'selected_frameIdx')
    error('DLC_filtered has no VideoParams.selected_frameIdx (session not aggregated with a DLC import)');
end
sel = D.VideoParams.selected_frameIdx;

[B2, VP2] = add_whisker_data(sdir, session, D.Behav, D.VideoParams);
if ~isfield(VP2, 'Whisker') || ~isfield(VP2.Whisker, 'whisker_box')
    error('the whisker csv/json in %s is not the front_top version (VideoParams.Whisker.whisker_box missing)', sdir);
end
for f = 1:numel(fields)
    if ~isfield(B2, fields{f}); error('new Behav has no %s', fields{f}); end
    if numel(B2.(fields{f})) ~= numel(D.Behav.FrameNum); error('%s has the wrong length', fields{f}); end
end

Data = matfile(finalFile, 'Writable', ~dry_run);
A = Data.Behav;
VPa = Data.VideoParams;
if numel(A.FrameNum) ~= numel(sel)
    error('aggregated Behav has %d rows but selected_frameIdx has %d', numel(A.FrameNum), numel(sel));
end
if ~isequal(A.FrameNum, D.Behav.FrameNum(sel))
    error('aggregated Behav.FrameNum differs from DLC_filtered Behav.FrameNum(selected_frameIdx)');
end

old = struct();
old.date = char(datetime('now'));
old.aggregated_Behav = struct();
old.dlc_Behav = struct();
for f = 1:numel(fields)
    if isfield(A, fields{f}); old.aggregated_Behav.(fields{f}) = A.(fields{f}); end
    if isfield(D.Behav, fields{f}); old.dlc_Behav.(fields{f}) = D.Behav.(fields{f}); end
end
if isfield(VPa, 'Whisker'); old.aggregated_VideoParams_Whisker = VPa.Whisker; end
if isfield(D.VideoParams, 'Whisker'); old.dlc_VideoParams_Whisker = D.VideoParams.Whisker; end

fprintf('    codec %s, REST %.3f, whisking %.0f%%, artifact %.1f%%, index p5/50/95 = %.2f / %.2f / %.2f   (old p5/50/95 of the aggregated trace: %.2f / %.2f / %.2f)\n', ...
    VP2.Whisker.codec, VP2.Whisker.rest_level.front_top, 100 * VP2.Whisker.fraction_whisking, 100 * VP2.Whisker.fraction_artifact, ...
    VP2.Whisker.index_p5, VP2.Whisker.index_median, VP2.Whisker.index_p95, ...
    prctile(A.WhiskIndex, 5), prctile(A.WhiskIndex, 50), prctile(A.WhiskIndex, 95));
if dry_run
    return;
end

bfile = fullfile(backup_dir, [session '_whisk_v1.mat']);
if exist(bfile, 'file')
    fprintf('    backup %s exists already (first backup kept: it is the original v1)\n', bfile);
else
    save(bfile, 'old', '-v7.3');
    fprintf('    backed up the old whisking fields to %s\n', bfile);
end

% DLC_filtered: full-length traces + params
Behav = B2;
VideoParams = VP2;
save(dlcFile, 'Behav', 'VideoParams', '-append');

% aggregated file: only the whisking fields and VideoParams.Whisker
for f = 1:numel(fields)
    A.(fields{f}) = B2.(fields{f})(sel, :);
end
VPa.Whisker = VP2.Whisker;
Data.Behav = A;
Data.VideoParams = VPa;
fprintf('    updated %s\n', finalFile);

end


function out = ternary(c, a, b)
if c; out = a; else; out = b; end
end
