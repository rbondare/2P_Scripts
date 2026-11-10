function export_session_meta(out_dir, varargin)
% EXPORT_SESSION_META  stage A of the recordings database: what only MATLAB can read cheaply.
%
%   export_session_meta                      -> all RB22-24 sessions found in DATA_2P, written to <Recordings_DB>\_build
%   export_session_meta(out_dir, 'sessions', {'AnimalRB24_261030_1744'}, 'force', true)
%
% Writes two CSV files (read by build_recdb.py, which assembles the final tables):
%   m_sessions.csv  one row per DATA_2P session: stimulus folder, aggregated file status, #planes, #aggregated ROIs, duration,
%                   timing / sync / ball / camera QC taken from the aggregated file, task type + reward settings
%   m_blocks.csv    one row per stimulus block (stim*.mat): type, key, loom direction, trials, wall-clock start, duration,
%                   interrupted, 2P-clock span + red-sync status (aggregated sessions only), compact parameter string
% Blocks are read from the RAW stim*.mat of every session (so un-aggregated sessions are covered too) and enriched from the
% aggregated Stimuli struct when it exists. Block key / loom direction / trial-onset logic is the one of
% axon analysis/cross_session/export_xs.m (same onsets as the analyses, so the trial counts agree with session_qc.csv).
% Incremental: a session is re-read only when its aggregated file or a stim file changed since the last run ('force' overrides).
% Read-only on all data; the only thing written is <out_dir>.

p = inputParser;  addParameter(p, 'sessions', {});  addParameter(p, 'force', false);  parse(p, varargin{:});
here = fileparts(mfilename('fullpath'));  repo = fileparts(fileparts(here));
addpath(repo, fullfile(repo, 'newAggregation'), fullfile(repo, 'analysis_functions'));
cfg = get_user_config('AnimalRB24');                       % only for the root paths (Z: share, DATA_2P, Stimulus_Axons, Aggregated)
rima_root = fileparts(fileparts(cfg.OutputDir));  rima_root = [rima_root filesep];
if nargin < 1 || isempty(out_dir); out_dir = fullfile(rima_root, 'Recordings_DB', '_build'); end
if ~isfolder(out_dir); mkdir(out_dir); end

% DATA_2P folder name -> stimulus folder name when they differ (typo in the 2P folder name: recorded 2026-10-28, see fov_matching/fov_common.py)
stim_override = containers.Map({'AnimalRB24_261021_1826'}, {'AnimalRB24_261028_1826'});

d = dir(cfg.DATA_2P);  names = {d([d.isdir]).name};
names = names(~cellfun(@isempty, regexp(names, '^AnimalRB2[234]_\d{6}_\d{4}$', 'once')));
if ~isempty(p.Results.sessions); names = names(ismember(names, p.Results.sessions)); end
names = sort(names);

fS = fullfile(out_dir, 'm_sessions.csv');  fB = fullfile(out_dir, 'm_blocks.csv');
oldS = table();  oldB = table();
if ~p.Results.force && isfile(fS) && isfile(fB)
    try; oldS = readtable(fS, 'TextType', 'char', 'DatetimeType', 'text', 'VariableNamingRule', 'preserve');  oldB = readtable(fB, 'TextType', 'char', 'DatetimeType', 'text', 'VariableNamingRule', 'preserve'); catch; oldS = table(); oldB = table(); end
end

rowsS = {};  rowsB = {};
for i = 1:numel(names)
    nm = names{i};  t0 = tic;
    sf = nm;  if isKey(stim_override, nm); sf = stim_override(nm); end
    us = strfind(nm, '_');  animal = nm(1:us(1) - 1);
    stim_dir = fullfile(cfg.StimBasePath, animal, sf);
    agg_file = fullfile(cfg.OutputDir, [nm '_preprocessed.mat']);
    sig = signature_(stim_dir, agg_file);
    if ~isempty(oldS) && any(strcmp(oldS.session, nm)) && strcmp(oldS.signature{strcmp(oldS.session, nm)}, sig)
        rowsS{end+1} = table2struct(oldS(strcmp(oldS.session, nm), :)); %#ok<AGROW>
        ob = oldB(strcmp(oldB.session, nm), :);
        for k = 1:height(ob); rowsB{end+1} = table2struct(ob(k, :)); end %#ok<AGROW>
        fprintf('%s: unchanged, reused (%d/%d)\n', nm, i, numel(names));  continue;
    end
    try
        [S, B] = one_session_(nm, sf, stim_dir, agg_file, sig);
    catch ME
        warning('%s FAILED: %s', nm, ME.message);
        S = empty_session_(nm, sf, sig);  S.meta_error = strtrim(ME.message);  B = {};
    end
    rowsS{end+1} = S; %#ok<AGROW>
    for k = 1:numel(B); rowsB{end+1} = B{k}; end %#ok<AGROW>
    fprintf('%s: %d blocks, aggregated=%d (%.1f s) (%d/%d)\n', nm, numel(B), S.has_aggregated, toc(t0), i, numel(names));
end

writetable(rows_to_table_(rowsS), fS);  writetable(rows_to_table_(rowsB), fB);
fprintf('EXPORT_SESSION_META done: %d sessions, %d blocks -> %s\n', numel(rowsS), numel(rowsB), out_dir);
end


% ------------------------------------------------------------------------------------------------------------------------------
function [S, B] = one_session_(nm, sf, stim_dir, agg_file, sig)
S = empty_session_(nm, sf, sig);
has_agg = isfile(agg_file);
S.has_aggregated = double(has_agg);
files = dir(fullfile(stim_dir, 'stim*.mat'));  [~, o] = sort({files.name});  files = files(o);
S.stim_folder_found = double(isfolder(stim_dir));  S.n_stim_files = numel(files);

D = struct();
if has_agg
    fi = dir(agg_file);  S.agg_mtime = datestr(fi.datenum, 'yyyy-mm-dd HH:MM:SS');  S.agg_size_mb = round(fi.bytes / 1e6, 1);
    m = matfile(agg_file);  vars = who(m);
    for v = {'Stimuli', 'Triggers', 'VideoParams', 'LocomotionCal', 'TimeCa', 'SI_Info'}
        if ismember(v{1}, vars); D.(v{1}) = m.(v{1}); end
    end
    if ismember('CaData', vars)
        try; c = m.CaData(1, 1);  S.n_agg_rois = size(c.Ca_dFF, 1); catch; end    % row 1 of CaData = the populated entry (accepted ROIs)
    end
    if isfield(D, 'TimeCa'); S.n_planes = size(D.TimeCa, 1);  S.duration_min_2p = round((D.TimeCa(1, end) - D.TimeCa(1, 1)) / 60, 2);  S.n_frames_2p = size(D.TimeCa, 2); end
    if isfield(D, 'SI_Info') && isfield(D.SI_Info, 'Clockstart'); S.clock_start = char(D.SI_Info.Clockstart, 'yyyy-MM-dd HH:mm:ss'); end
    if isfield(D, 'Triggers') && isfield(D.Triggers, 'TimingProblems'); S.timing_problems = double(any(D.Triggers.TimingProblems(:))); end
    if isfield(D, 'Triggers') && isfield(D.Triggers, 'TimeCamera'); S.n_cam_triggers = numel(D.Triggers.TimeCamera); end
    if isfield(D, 'VideoParams')
        V = D.VideoParams;
        if isfield(V, 'Pupil') && isfield(V.Pupil, 'frac_bad'); S.pupil_frac_bad = V.Pupil.frac_bad; end
        if isfield(V, 'Whisker') && isstruct(V.Whisker)
            W = V.Whisker;
            if isfield(W, 'fraction_whisking'); S.whisk_fraction = W.fraction_whisking; end
            if isfield(W, 'fraction_artifact'); S.whisk_artifact_fraction = W.fraction_artifact; end
            if isfield(W, 'codec'); S.video_codec = W.codec; end
            if isfield(W, 'fps'); S.video_fps = W.fps; end
        end
    end
    if isfield(D, 'LocomotionCal')
        L = D.LocomotionCal;  nrd = 0;  ntr = 0;  oks = [];
        for k = 1:numel(L)
            if isfield(L, 'Forward_mm') && ~isempty(L(k).Forward_mm); nrd = nrd + numel(L(k).Forward_mm); end
            if isfield(L, 'Triggers_sec') && ~isempty(L(k).Triggers_sec); ntr = ntr + numel(L(k).Triggers_sec); end
            if isfield(L, 'SyncInfo') && isstruct(L(k).SyncInfo) && isfield(L(k).SyncInfo, 'ok'); oks(end+1) = double(L(k).SyncInfo.ok); end %#ok<AGROW>
        end
        S.ball_n_entries = numel(L);  S.ball_readings = nrd;  S.ball_triggers = ntr;
        if ~isempty(oks); S.ball_sync_ok = double(all(oks)); end    % empty = open-loop sessions (no SyncInfo)
    end
end

% ---- blocks: raw stim files (all sessions), enriched from the aggregated Stimuli struct
B = {};
agg_blocks = [];  if isfield(D, 'Stimuli'); agg_blocks = D.Stimuli; end
if ~isempty(agg_blocks) && numel(agg_blocks) ~= numel(files)
    S.meta_error = sprintf('aggregated Stimuli has %d blocks, folder has %d stim files', numel(agg_blocks), numel(files));
end
for k = 1:numel(files)
    R = load(fullfile(stim_dir, files(k).name));
    opt = R.options;  bl = empty_block_(nm, k, files(k).name);
    bl.type = opt.type;  if isfield(opt, 'addnote'); bl.addnote = char(string(opt.addnote)); end
    sp = struct();  if isfield(opt, 'specParams') && isstruct(opt.specParams); sp = opt.specParams; end
    [bl.key, bl.loom_dir] = block_key_(bl.type, sp, bl.addnote);
    if isfield(opt, 'trials') && isnumeric(opt.trials) && ~isempty(opt.trials); bl.trials_requested = double(opt.trials(1)); end
    % timestart / timeend are one datetime per block, or one per trial in the training blocks -> block span = earliest start .. latest end
    if isfield(R, 'timestart') && isdatetime(R.timestart) && ~isempty(R.timestart); t_s = min(R.timestart(:));  bl.start_wallclock = char(t_s, 'yyyy-MM-dd HH:mm:ss'); end
    if isfield(R, 'timestart') && isfield(R, 'timeend') && isdatetime(R.timestart) && isdatetime(R.timeend) && ~isempty(R.timeend)
        bl.duration_s = round(seconds(max(R.timeend(:)) - min(R.timestart(:))), 1);
    end
    if isfield(R, 'interrupted') && ~isempty(R.interrupted); bl.interrupted = double(any(R.interrupted(:))); end
    if isfield(opt, 'projectortype'); bl.projector = char(string(opt.projectortype)); end
    if isfield(opt, 'SetUp'); bl.setup = double(opt.SetUp(1)); end
    bl.params = params_str_(bl.type, sp);
    Fi = [];  if isfield(R, 'Frameinfo'); Fi = R.Frameinfo; end
    s = [];
    if ~isempty(agg_blocks) && numel(agg_blocks) == numel(files); s = agg_blocks(k); end
    try
        if ~isempty(s)
            [bl.n_trials, bl.t_start_2p, bl.t_end_2p] = agg_trials_(s, bl.key);
            if isfield(s, 'RedFrameSynchronized') && ~isempty(s.RedFrameSynchronized); bl.red_sync_ok = double(all(s.RedFrameSynchronized(:))); end
            if isfield(s, 'RedSyncInfo') && isstruct(s.RedSyncInfo) && isfield(s.RedSyncInfo, 'method') && numel(s.RedSyncInfo) == 1; bl.red_sync_method = char(string(s.RedSyncInfo.method)); end
        elseif ~isempty(Fi) && any(strcmp(bl.key, {'loom_expand', 'loom_recede', 'bar', 'fff'})) && isfield(Fi, 'subtrial_frame_count') && isfield(Fi, 'frame_count')
            sfc = Fi.subtrial_frame_count(~isnan(Fi.frame_count));  bl.n_trials = sum(sfc == 1);
        end
    catch ME
        bl.note = ['trial count failed: ' strtrim(ME.message)];
    end
    if strcmp(bl.key, 'square') && isfield(R, 'RewardMode'); S = task_fields_(S, sp, R); end
    B{end+1} = bl; %#ok<AGROW>
end
S.n_blocks = numel(B);
if ~isempty(B)
    S.n_blocks_interrupted = nansum_(cellfun(@(b) b.interrupted, B));
    keys = cellfun(@(b) b.key, B, 'UniformOutput', false);  nt = cellfun(@(b) b.n_trials, B);
    u = unique(keys, 'stable');  parts = cell(1, numel(u));
    for j = 1:numel(u); parts{j} = sprintf('%s x%d', u{j}, nansum_(nt(strcmp(keys, u{j})))); if all(isnan(nt(strcmp(keys, u{j})))); parts{j} = u{j}; end; end
    S.stimuli = strjoin(parts, ', ');
    S.block_sequence = strjoin(keys, ' > ');
    rs = cellfun(@(b) b.red_sync_ok, B);  S.n_blocks_red_sync_ok = nansum_(rs);  S.n_blocks_red_sync_checked = sum(~isnan(rs));
    S.task_type = task_type_(keys, S);
end
end


function S = task_fields_(S, sp, R)
% lick task vs passive square: RewardMode 1 = spout + reward, 0 = square without reward (axon_analysis_task_analysis_V1.m:270)
if isfield(sp, 'RewardMode'); S.reward_mode = double(sp.RewardMode); elseif isfield(R, 'RewardMode'); S.reward_mode = double(R.RewardMode(1)); end
if isfield(sp, 'level'); S.reward_level = double(sp.level); end
if isfield(sp, 'RewardDelay'); S.reward_delay_s = double(sp.RewardDelay); end
if isfield(sp, 'StimulusContrast'); S.square_contrast = double(sp.StimulusContrast); end
end


function t = task_type_(keys, S)
if any(strcmp(keys, 'square'))
    if isnan(S.reward_mode); t = 'square';
    elseif S.reward_mode == 1; t = 'lick task';
    else; t = 'passive square'; end
elseif any(strcmp(keys, 'alleyCL'))
    if any(ismember(keys, {'flash', 'fff', 'bar', 'loom_expand', 'loom_recede'})); t = 'open-loop visual + corridor'; else; t = 'corridor (visuomotor)'; end
elseif any(ismember(keys, {'flash', 'fff', 'bar', 'loom_expand', 'loom_recede'})); t = 'open-loop visual';
elseif all(strcmp(keys, 'spont')); t = 'spontaneous';
else; t = 'other'; end
end


function [n, t0, t1] = agg_trials_(s, key)
% trial onsets exactly as in cross_session/export_xs.m (same counts as the analyses)
n = NaN;  t0 = NaN;  t1 = NaN;  o = [];
TF = s.TimeStimulusFrame;
if isempty(TF); return; end
f = isfinite(TF);  if any(f(:)); t0 = round(min(TF(f)), 2);  t1 = round(max(TF(f)), 2); end
switch key
    case 'flash'; o = get_flash_onsets(s);
    case {'fff', 'bar', 'loom_expand', 'loom_recede'}
        sfc = s.Frameinfo.subtrial_frame_count(~isnan(s.Frameinfo.frame_count));  o = TF(find(sfc == 1));
    case 'square'
        Fi = s.Frameinfo;  n_w1 = size(Fi.frame_count_wait1, 1);  n_st = size(Fi.frame_count_stim, 1);
        FV = s.FrameValid & isfinite(TF);  o = nan(size(TF, 2), 1);
        for t = 1:size(TF, 2); j = find(FV(n_w1 + (1:n_st), t), 1); if ~isempty(j); o(t) = TF(n_w1 + j, t); end; end
    otherwise; return;
end
n = sum(isfinite(o(:)));
end


function [key, ld] = block_key_(type, sp, addnote)
ld = NaN;
switch type
    case 'sparse_local_global_flashes'; key = 'flash';
    case 'full_field_flash';            key = 'fff';
    case 'moving_bar';                  key = 'bar';
    case 'looming'
        key = 'loom_expand';  ld = 1;
        if isfield(sp, 'loom_maxexpand') && ~isempty(sp.loom_maxexpand); ld = sign(sp.loom_maxexpand(1));
        elseif contains(lower(addnote), 'receding'); ld = -1; end
        if ld < 0; key = 'loom_recede'; end
    case 'training'
        key = 'square';  if isfield(sp, 'stimulusType') && ~strcmp(char(string(sp.stimulusType)), 'square'); key = 'training'; end
    case 'spontaneous'; key = 'spont';
    case 'alleyCL';     key = 'alleyCL';
    otherwise;          key = ['other_' type];
end
end


function s = params_str_(type, sp)
% compact "name=value;" string of the scalar / short-vector stimulus parameters (training blocks: whitelist only, the struct has ~50 fields)
if isempty(fieldnames(sp)); s = ''; return; end
fn = fieldnames(sp);
if strcmp(type, 'training')
    fn = intersect({'level', 'RewardMode', 'RewardDelay', 'StimulusContrast', 'StimSize', 'locationMode', 'minWait', 'maxWait', 'Tend', 'lickWindowMin', 'lickWindowMax', 'subsessionwait', 'repetitions'}, fn, 'stable');
end
parts = {};
for i = 1:numel(fn)
    v = sp.(fn{i});
    if isnumeric(v) || islogical(v)
        if numel(v) <= 4 && ~isempty(v); parts{end+1} = sprintf('%s=%s', fn{i}, strtrim(sprintf('%g ', v))); end %#ok<AGROW>
    elseif ischar(v) && numel(v) <= 24 && ~isempty(v)
        parts{end+1} = sprintf('%s=%s', fn{i}, v); %#ok<AGROW>
    elseif isstring(v) && isscalar(v) && strlength(v) <= 24
        parts{end+1} = sprintf('%s=%s', fn{i}, char(v)); %#ok<AGROW>
    end
end
s = strjoin(strrep(parts, ',', ' '), ';');
end


function sig = signature_(stim_dir, agg_file)
% change detector for the incremental mode: aggregated file (mtime, size) + stim files (names, mtimes)
parts = {SCHEMA_VERSION_};     % bump when the extracted columns / logic change: forces a re-read of every session in the incremental mode
if isfile(agg_file); fi = dir(agg_file); parts{end+1} = sprintf('agg:%.5f:%d', fi.datenum, fi.bytes); else; parts{end+1} = 'agg:none'; end
if isfolder(stim_dir); f = dir(fullfile(stim_dir, 'stim*.mat')); for k = 1:numel(f); parts{end+1} = sprintf('%s:%.5f', f(k).name, f(k).datenum); end; end %#ok<AGROW>
sig = sprintf('s%08x', mod(sum(double(strjoin(parts, '|')) .* (1:numel(strjoin(parts, '|')))), 2^32));
end


function v = SCHEMA_VERSION_(); v = 'schema3'; end


function n = nansum_(x); x = x(~isnan(x)); if isempty(x); n = 0; else; n = sum(x); end; end


function S = empty_session_(nm, sf, sig)
S = struct('session', nm, 'stim_folder', sf, 'stim_folder_found', NaN, 'n_stim_files', NaN, 'has_aggregated', 0, 'agg_mtime', '', 'agg_size_mb', NaN, ...
    'n_planes', NaN, 'n_agg_rois', NaN, 'n_frames_2p', NaN, 'duration_min_2p', NaN, 'clock_start', '', 'n_blocks', NaN, 'n_blocks_interrupted', NaN, ...
    'stimuli', '', 'block_sequence', '', 'task_type', '', 'reward_mode', NaN, 'reward_level', NaN, 'reward_delay_s', NaN, 'square_contrast', NaN, ...
    'timing_problems', NaN, 'n_blocks_red_sync_ok', NaN, 'n_blocks_red_sync_checked', NaN, 'n_cam_triggers', NaN, 'video_codec', '', 'video_fps', NaN, ...
    'pupil_frac_bad', NaN, 'whisk_fraction', NaN, 'whisk_artifact_fraction', NaN, 'ball_n_entries', NaN, 'ball_readings', NaN, 'ball_triggers', NaN, 'ball_sync_ok', NaN, ...
    'meta_error', '', 'signature', sig);
end


function bl = empty_block_(nm, k, fname)
bl = struct('session', nm, 'block', k, 'file', fname, 'type', '', 'key', '', 'loom_dir', NaN, 'addnote', '', 'trials_requested', NaN, 'n_trials', NaN, ...
    'start_wallclock', '', 'duration_s', NaN, 'interrupted', NaN, 't_start_2p', NaN, 't_end_2p', NaN, 'red_sync_ok', NaN, 'red_sync_method', '', ...
    'projector', '', 'setup', NaN, 'params', '', 'note', '');
end


function T = rows_to_table_(rows)
if isempty(rows); T = table(); return; end
fn = fieldnames(rows{1});
for i = 2:numel(rows); fn = union(fn, fieldnames(rows{i}), 'stable'); end
for i = 1:numel(rows)
    for f = 1:numel(fn)
        if ~isfield(rows{i}, fn{f}); rows{i}.(fn{f}) = []; end
    end
    rows{i} = orderfields(rows{i}, fn);
end
T = struct2table([rows{:}], 'AsArray', true);
for f = 1:numel(fn)
    v = T.(fn{f});
    if iscell(v)
        T.(fn{f}) = cellfun(@(x) cell_to_char_(x), v, 'UniformOutput', false);
    end
end
end


function s = cell_to_char_(x)
if isempty(x); s = ''; elseif isnumeric(x) || islogical(x); s = sprintf('%g', x); elseif isstring(x); s = char(x); else; s = x; end
end
