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
%                   interrupted, 2P-clock span + red-sync status (aggregated sessions only), compact parameter string,
%                   a readable description (stim_desc / stim_short) and, for the square (training) blocks, the task facts:
%                   stimulus shape / contrast / size / location, GUI level (2.5 = adaptive) + the level of every trial (levelHistory),
%                   reward mode / delay / lick window, and SPOUT yes/no = whether any lick was logged in lickData at all
%
% Also (option 'xs', default true): the ROI classes of the FOV pictures (signal-like / noise-like / neuropil-limited) are computed from the xs export
% of cross_session/export_xs.m, so every aggregated session without an xs file gets one made here. Existing xs files are never overwritten.
% Blocks are read from the RAW stim*.mat of every session (so un-aggregated sessions are covered too) and enriched from the
% aggregated Stimuli struct when it exists. Block key / loom direction / trial-onset logic is the one of
% axon analysis/cross_session/export_xs.m (same onsets as the analyses, so the trial counts agree with session_qc.csv).
% Incremental: a session is re-read only when its aggregated file or a stim file changed since the last run ('force' overrides).
% Read-only on all data; the only thing written is <out_dir>.

p = inputParser;  addParameter(p, 'sessions', {});  addParameter(p, 'force', false);  addParameter(p, 'xs', true);  parse(p, varargin{:});
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

% ---- xs exports for the ROI classes (aggregated sessions that have none yet; never overwrites)
if p.Results.xs
    xs_dir = 'C:\Users\rbondarenko\2P_cache\xs';
    addpath(fullfile(repo, 'axon analysis', 'cross_session'));
    for i = 1:numel(names)
        nm = names{i};
        if ~isfile(fullfile(cfg.OutputDir, [nm '_preprocessed.mat'])) || isfile(fullfile(xs_dir, [nm '_xs.mat'])); continue; end
        try
            t0 = tic;  export_xs(nm, xs_dir);  fprintf('  xs export made for %s (%.0f s)\n', nm, toc(t0));
        catch ME
            warning('xs export of %s FAILED (its FOV picture stays unclassified): %s', nm, strtrim(ME.message));
        end
    end
end
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
    if strcmp(bl.type, 'training'); bl = training_fields_(bl, R, sp); end
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
    bl.stim_desc = describe_block_(bl, sp, false);  bl.stim_short = describe_block_(bl, sp, true);
    B{end+1} = bl; %#ok<AGROW>
end
S.n_blocks = numel(B);
if ~isempty(B)
    S.n_blocks_interrupted = nansum_(cellfun(@(b) b.interrupted, B));
    keys = cellfun(@(b) b.key, B, 'UniformOutput', false);  nt = cellfun(@(b) b.n_trials, B);
    % session summary of the stimuli: one entry per distinct short description (blocks that look alike are pooled), with the completed trials
    shorts = cellfun(@(b) b.stim_short, B, 'UniformOutput', false);
    nt_s = nt;  for q = 1:numel(B); if isnan(nt_s(q)) && strcmp(B{q}.key, 'square'); nt_s(q) = B{q}.trials_logged; end; end    % raw-only session: trials in the stim file
    u = unique(shorts, 'stable');  parts = cell(1, numel(u));
    for j = 1:numel(u)
        sel = strcmp(shorts, u{j});  parts{j} = sprintf('%s x%d', u{j}, nansum_(nt_s(sel)));
        if all(isnan(nt_s(sel))); parts{j} = u{j}; end
    end
    S.stimuli = strjoin(parts, '; ');
    S.block_sequence = strjoin(keys, ' > ');
    S = square_session_fields_(S, B);
    rs = cellfun(@(b) b.red_sync_ok, B);  S.n_blocks_red_sync_ok = nansum_(rs);  S.n_blocks_red_sync_checked = sum(~isnan(rs));
    S.task_type = task_type_(keys, S);
end
end


function S = task_fields_(S, sp, R)
% lick task vs passive square: RewardMode 1 = spout + reward, 0 = square without reward (axon_analysis_task_analysis_V1.m:270)
if isfield(sp, 'RewardMode'); S.reward_mode = double(sp.RewardMode); elseif isfield(R, 'RewardMode'); S.reward_mode = double(R.RewardMode(1)); end
if isfield(sp, 'RewardDelay'); S.reward_delay_s = double(sp.RewardDelay); end
if isfield(sp, 'StimulusContrast'); S.square_contrast = double(sp.StimulusContrast); end
% the level is NOT taken from specParams.level here: in the adaptive mode (GUI level 2.5) it is only the level of the LAST trial (see square_session_fields_)
end


function bl = training_fields_(bl, R, sp)
% task facts of a training (square) block, read from the raw stim file (so un-aggregated sessions are covered too)
bl.stim_shape = char(string(pick_(sp, R, 'stimulusType', '')));
bl.stim_contrast_pct = num_(pick_(sp, R, 'StimulusContrast', NaN));
bl.stim_size = num_(pick_(sp, R, 'StimSize', NaN));
bl.stim_location = char(string(pick_(sp, R, 'locationMode', '')));
bl.reward_mode = num_(pick_(sp, R, 'RewardMode', NaN));
bl.reward_delay_s = num_(pick_(sp, R, 'RewardDelay', NaN));
lw = [num_(pick_(sp, R, 'lickWindowMin', NaN)) num_(pick_(sp, R, 'lickWindowMax', NaN))];  if all(isfinite(lw)); bl.lick_window_s = sprintf('%g-%g', lw); end
mw = [num_(pick_(sp, R, 'minWait', NaN)) num_(pick_(sp, R, 'maxWait', NaN))];  if all(isfinite(mw)); bl.wait_range_s = sprintf('%g-%g', mw); end

% level: the GUI level is the top-level variable `level` (1, 2, 3, or 2.5 = ADAPTIVE: every trial runs level 2 or 3, drawn from the agent state; the level of each trial is in levelHistory).
% specParams.level is only the level of the last trial, so it must not be used as the session level (axon_analysis_task_analysis_V1.m, reward section).
gui = num_(pick_(struct(), R, 'level', NaN));  if ~isfinite(gui); gui = num_(pick_(sp, struct(), 'level', NaN)); end
bl.gui_level = gui;
lh = [];  if isfield(R, 'levelHistory') && isnumeric(R.levelHistory); lh = double(R.levelHistory(:)); lh = lh(isfinite(lh)); end
nts = {};
if ~isempty(lh)
    bl.level_mode = 'adaptive';  bl.trials_logged = numel(lh);
    bl.n_L1 = sum(lh == 1);  bl.n_L2 = sum(lh == 2);  bl.n_L3 = sum(lh == 3);
    if gui ~= 2.5; nts{end+1} = sprintf('levelHistory present but GUI level is %g', gui); end %#ok<AGROW>
else
    bl.level_mode = 'fixed';  nt = NaN;
    if isfield(R, 'timestart') && isdatetime(R.timestart); nt = numel(R.timestart); end
    bl.trials_logged = nt;
    if ismember(gui, [1 2 3]) && isfinite(nt); bl.(sprintf('n_L%d', gui)) = nt; end
    if gui == 2.5; nts{end+1} = 'GUI level 2.5 (adaptive) but no levelHistory: per-trial levels unknown'; end %#ok<AGROW>
end

% spout = a lick was detected at all: lickData (cell of Arduino ms strings, one column per trial) holds at least one finite entry. No spout / reward mode 0: lickData is empty (0 x 0).
n_l = 0;  n_lt = 0;  ncol = 0;
if isfield(R, 'lickData') && ~isempty(R.lickData)
    L = R.lickData;  ncol = size(L, 2);
    if isnumeric(L) || islogical(L); v = double(L); else; v = str2double(L); end
    fin = isfinite(v);  n_l = sum(fin(:));  n_lt = sum(any(fin, 1));
end
bl.n_lick_entries = n_l;  bl.n_trials_with_licks = n_lt;  bl.lick_log_trials = ncol;
if n_l > 0; bl.spout = 'yes'; else; bl.spout = 'no'; end
if bl.reward_mode == 1 && n_l == 0; nts{end+1} = 'RewardMode 1 but no lick in the log'; end %#ok<AGROW>
if bl.reward_mode == 0 && n_l > 0; nts{end+1} = 'RewardMode 0 but licks in the log'; end %#ok<AGROW>
if n_l > 0 && isfinite(bl.trials_logged) && ncol ~= bl.trials_logged; nts{end+1} = sprintf('lick log has %d trial columns for %d trials', ncol, bl.trials_logged); end %#ok<AGROW>
bl.spout_note = strjoin(nts, '; ');
end


function v = pick_(sp, R, name, dflt)
% parameter `name` from specParams, else from the top level of the stim file, else the default
if isfield(sp, name) && ~isempty(sp.(name)); v = sp.(name); elseif isfield(R, name) && ~isempty(R.(name)); v = R.(name); else; v = dflt; end
end


function x = num_(v)
if (isnumeric(v) || islogical(v)) && ~isempty(v); x = double(v(1)); else; x = NaN; end
end


function s = describe_block_(bl, sp, short)
% readable one-line description of a block: short = for the session summary (`stimuli`), else the full one (blocks.csv stim_desc). Only facts that are in the stim file.
SEP = ' | ';
switch bl.key
    case 'square'
        shape = bl.stim_shape;  if isempty(shape); shape = 'square'; end
        lvl = '';
        if isfinite(bl.gui_level)
            lvl = sprintf('level %g', bl.gui_level);
            if strcmp(bl.level_mode, 'adaptive'); lvl = sprintf('%s adaptive (%s)', lvl, level_counts_(bl)); end
        end
        if strcmp(bl.spout, 'yes'); sp_s = 'spout'; elseif strcmp(bl.spout, 'no'); sp_s = 'no spout'; else; sp_s = ''; end
        if short
            head = shape;  if isfinite(bl.stim_contrast_pct); head = sprintf('%s %g%%', shape, bl.stim_contrast_pct); end
            parts = {head, bl.stim_location, lvl, sp_s};
            s = strjoin(parts(~cellfun(@isempty, parts)), SEP);
        else
            parts = {shape};
            if isfinite(bl.stim_contrast_pct); parts{end+1} = sprintf('%g%% contrast', bl.stim_contrast_pct); end %#ok<AGROW>
            if ~isempty(bl.stim_location); parts{end+1} = bl.stim_location; end %#ok<AGROW>
            if isfinite(bl.stim_size); parts{end+1} = sprintf('size %g', bl.stim_size); end %#ok<AGROW>
            parts{end+1} = lvl;  %#ok<AGROW>
            if bl.reward_mode == 0; parts{end+1} = 'reward mode 0 (square shown without reward)';  %#ok<AGROW>
            elseif bl.reward_mode == 1; parts{end+1} = sprintf('reward mode 1 (reward delay %g s)', bl.reward_delay_s); end %#ok<AGROW>
            if ~isempty(bl.lick_window_s); parts{end+1} = ['lick window ' bl.lick_window_s ' s']; end %#ok<AGROW>
            if ~isempty(bl.wait_range_s); parts{end+1} = ['wait ' bl.wait_range_s ' s before the square']; end %#ok<AGROW>
            if strcmp(bl.spout, 'yes'); parts{end+1} = sprintf('spout (%d lick entries in %d trials)', bl.n_lick_entries, bl.n_trials_with_licks); %#ok<AGROW>
            elseif strcmp(bl.spout, 'no'); parts{end+1} = 'no spout (no lick in the log)'; end %#ok<AGROW>
            s = strjoin(parts(~cellfun(@isempty, parts)), SEP);
        end
    case {'loom_expand', 'loom_recede'}
        if strcmp(bl.key, 'loom_expand'); nm = 'expanding loom'; else; nm = 'receding loom'; end
        has_grey = isfield(sp, 'pre_grey_t') || isfield(sp, 'post_grey_t');
        if short
            s = bl.key;  if has_grey; s = [s ' (+grey)']; end
        else
            parts = {nm};
            v = fld_(sp, 'loom_start_radius');  if ~isempty(v); parts{end+1} = sprintf('start radius %g', v(1)); end %#ok<AGROW>
            v = fld_(sp, 'loom_maxexpand');  if ~isempty(v); parts{end+1} = sprintf('max expand %g', v(1)); end %#ok<AGROW>
            v = fld_(sp, 'loom_speed');  if ~isempty(v); parts{end+1} = sprintf('speed %g', v(1)); end %#ok<AGROW>
            v = fld_(sp, 'loom_origin');  if ~isempty(v); parts{end+1} = ['origin ' strtrim(sprintf('%g ', v))]; end %#ok<AGROW>
            if has_grey; parts{end+1} = sprintf('grey %g s before / %g s after (outside the stimulus frames)', num_(pick_(sp, struct(), 'pre_grey_t', NaN)), num_(pick_(sp, struct(), 'post_grey_t', NaN))); end %#ok<AGROW>
            s = strjoin(parts, SEP);
        end
    case 'flash'
        s = bl.key;
        if ~short
            parts = {'flash (sparse local/global)'};
            v = fld_(sp, 'duration');  if ~isempty(v); parts{end+1} = sprintf('%g s duration', v(1)); end %#ok<AGROW>
            v = fld_(sp, 'interval');  if ~isempty(v); parts{end+1} = sprintf('every %s s', strjoin(arrayfun(@(x) sprintf('%g', x), v, 'UniformOutput', false), '-')); end %#ok<AGROW>
            s = strjoin(parts, SEP);
        end
    case 'fff'
        s = bl.key;
        if ~short
            parts = {'full-field flash (luminance modulation)'};
            v = fld_(sp, 'fff_sineamp');  if ~isempty(v); parts{end+1} = sprintf('sine amplitude %g', v(1)); end %#ok<AGROW>
            v = fld_(sp, 'fff_freq_factor');  if ~isempty(v); parts{end+1} = sprintf('frequency factor %g', v(1)); end %#ok<AGROW>
            s = strjoin(parts, SEP);
        end
    case 'bar'
        s = bl.key;
        if ~short
            parts = {'moving bar'};
            v = fld_(sp, 'bar_width');  if ~isempty(v); parts{end+1} = sprintf('width %g', v(1)); end %#ok<AGROW>
            v = fld_(sp, 'bar_moving_duration');  if ~isempty(v); parts{end+1} = sprintf('%g s per sweep', v(1)); end %#ok<AGROW>
            if ~isempty(bl.addnote); parts{end+1} = ['note: ' bl.addnote]; end %#ok<AGROW>
            s = strjoin(parts, SEP);
        end
    case 'spont'
        s = bl.key;  if ~short; s = 'spontaneous (grey screen)'; end
    case 'alleyCL'
        s = bl.key;  if ~short && ~isempty(bl.addnote); s = ['alleyCL | ' bl.addnote]; end
    otherwise
        s = bl.key;
end
s = strrep(s, ',', ' ');
end


function v = fld_(sp, name)
if isfield(sp, name) && isnumeric(sp.(name)) && ~isempty(sp.(name)); v = double(sp.(name)); else; v = []; end
end


function t = level_counts_(bl)
parts = {};
for l = 1:3
    n = bl.(sprintf('n_L%d', l));
    if isfinite(n) && n > 0; parts{end+1} = sprintf('L%d %d', l, n); end %#ok<AGROW>
end
t = strjoin(parts, ' / ');
end


function S = square_session_fields_(S, B)
% session-level summary of the square (training) blocks: level, spout, lick entries
sq = B(cellfun(@(b) strcmp(b.key, 'square'), B));
if isempty(sq); return; end
sps = cellfun(@(b) b.spout, sq, 'UniformOutput', false);
if all(strcmp(sps, 'yes')); S.spout = 'yes'; elseif all(strcmp(sps, 'no')); S.spout = 'no'; else; S.spout = 'mixed'; end
S.n_lick_entries = nansum_(cellfun(@(b) b.n_lick_entries, sq));
lv = cellfun(@(b) b.gui_level, sq);  md = cellfun(@(b) b.level_mode, sq, 'UniformOutput', false);
if all(lv == lv(1)); S.task_level = lv(1); else; S.task_level = NaN; end
if all(strcmp(md, md{1})); S.level_mode = md{1}; else; S.level_mode = 'mixed'; end
n = zeros(1, 3);  for l = 1:3; n(l) = nansum_(cellfun(@(b) b.(sprintf('n_L%d', l)), sq)); end
parts = arrayfun(@(l) sprintf('L%d %d', l, n(l)), find(n > 0), 'UniformOutput', false);
S.level_counts = strjoin(parts, ' / ');
nts = cellfun(@(b) b.spout_note, sq, 'UniformOutput', false);  nts = nts(~cellfun(@isempty, nts));
S.spout_note = strjoin(unique(nts, 'stable'), '; ');
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


function v = SCHEMA_VERSION_(); v = 'schema4'; end     % schema4: block descriptions + task level / spout fields (2026-11-11)


function n = nansum_(x); x = x(~isnan(x)); if isempty(x); n = 0; else; n = sum(x); end; end


function S = empty_session_(nm, sf, sig)
S = struct('session', nm, 'stim_folder', sf, 'stim_folder_found', NaN, 'n_stim_files', NaN, 'has_aggregated', 0, 'agg_mtime', '', 'agg_size_mb', NaN, ...
    'n_planes', NaN, 'n_agg_rois', NaN, 'n_frames_2p', NaN, 'duration_min_2p', NaN, 'clock_start', '', 'n_blocks', NaN, 'n_blocks_interrupted', NaN, ...
    'stimuli', '', 'block_sequence', '', 'task_type', '', 'reward_mode', NaN, 'task_level', NaN, 'level_mode', '', 'level_counts', '', 'reward_delay_s', NaN, 'square_contrast', NaN, ...
    'spout', '', 'n_lick_entries', NaN, 'spout_note', '', ...
    'timing_problems', NaN, 'n_blocks_red_sync_ok', NaN, 'n_blocks_red_sync_checked', NaN, 'n_cam_triggers', NaN, 'video_codec', '', 'video_fps', NaN, ...
    'pupil_frac_bad', NaN, 'whisk_fraction', NaN, 'whisk_artifact_fraction', NaN, 'ball_n_entries', NaN, 'ball_readings', NaN, 'ball_triggers', NaN, 'ball_sync_ok', NaN, ...
    'meta_error', '', 'signature', sig);
end


function bl = empty_block_(nm, k, fname)
bl = struct('session', nm, 'block', k, 'file', fname, 'type', '', 'key', '', 'loom_dir', NaN, 'addnote', '', 'trials_requested', NaN, 'n_trials', NaN, ...
    'start_wallclock', '', 'duration_s', NaN, 'interrupted', NaN, 't_start_2p', NaN, 't_end_2p', NaN, 'red_sync_ok', NaN, 'red_sync_method', '', ...
    'projector', '', 'setup', NaN, 'params', '', 'note', '', 'stim_desc', '', 'stim_short', '', ...
    'stim_shape', '', 'stim_contrast_pct', NaN, 'stim_size', NaN, 'stim_location', '', 'gui_level', NaN, 'level_mode', '', 'n_L1', NaN, 'n_L2', NaN, 'n_L3', NaN, 'trials_logged', NaN, ...
    'reward_mode', NaN, 'reward_delay_s', NaN, 'lick_window_s', '', 'wait_range_s', '', 'spout', '', 'n_lick_entries', NaN, 'n_trials_with_licks', NaN, 'lick_log_trials', NaN, 'spout_note', '');
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
