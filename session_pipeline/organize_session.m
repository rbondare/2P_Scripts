function S = organize_session(preprocessed_file, cfg)
% ORGANIZE_SESSION  Chunk, align and average one aggregated recording into
% a canonical per-session struct (see session_pipeline/PIPELINE.md).
%
% This is the working replacement for RB_DataGathering.m: it does the real
% per-trial slicing and trial averaging RB_DataGathering only stubbed, and
% reuses the already-validated timing/RNG-reconstruction extractors in
% analysis_functions/ verbatim (no timing logic is re-derived here).
%
%   S = organize_session(preprocessed_file)
%   S = organize_session(preprocessed_file, cfg)
%
% cfg fields (all optional):
%   .ca_type            1=Ca_dFF (default), 2=Ca_deconvolved, 3=Ca_F
%   .selected_plane_idx imaging plane, 1-based (default 1)
%   .fs_range           plausible sampling-rate bounds, Hz (default [2 30])
%
% Returns S with S.meta and S.stimuli.<type> (see PIPELINE.md for the full
% field contract). QC is NOT run here -- call qc_session(S, B, cfg)
% separately (run_organize.m does both).

    if nargin < 2 || isempty(cfg); cfg = struct(); end
    if ~isfield(cfg, 'ca_type');            cfg.ca_type = 1;            end
    if ~isfield(cfg, 'selected_plane_idx'); cfg.selected_plane_idx = 1; end
    if ~isfield(cfg, 'fs_range');           cfg.fs_range = [2 30];      end

    % Make the reused extractors + reconstructFlashes visible on the path.
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);                       % 2P_Scripts/
    addpath(fullfile(root, 'analysis_functions'));
    addpath(root);

    assert(isfile(preprocessed_file), 'File not found: %s', preprocessed_file);
    fprintf('\n=== organize_session: %s ===\n', preprocessed_file);
    B = load(preprocessed_file);

    % ---- resolve calcium field + plane ------------------------------------
    ca_fields = {'Ca_dFF', 'Ca_deconvolved', 'Ca_F'};
    ca_field = ca_fields{cfg.ca_type};
    assert(isfield(B.CaData(1), ca_field), ...
        'CaData(1).%s not present in %s', ca_field, preprocessed_file);

    dff_full  = B.CaData(1).(ca_field);
    centroidZ = B.CaData(1).Ca_centroid_voxel(:, 3);
    unique_planes = unique(centroidZ);
    assert(cfg.selected_plane_idx <= numel(unique_planes), ...
        'Requested plane %d but only %d planes present', ...
        cfg.selected_plane_idx, numel(unique_planes));
    plane_z = unique_planes(cfg.selected_plane_idx);
    roi_idx = find(centroidZ == plane_z);
    dff = dff_full(roi_idx, :);

    % ---- resolve this plane's time vector ---------------------------------
    % TimeCa is [n_planes x n_volumes]; use this plane's own row when it
    % exists (sub-volume-accurate), else fall back to row 1 / vector.
    if isvector(B.TimeCa)
        time_vec = B.TimeCa(:)'; time_row = 1;
    elseif size(B.TimeCa, 1) >= cfg.selected_plane_idx
        time_row = cfg.selected_plane_idx;
        time_vec = B.TimeCa(time_row, :);
    else
        time_row = 1;
        time_vec = B.TimeCa(1, :);
        warning('TimeCa has %d rows < plane %d; using row 1', size(B.TimeCa,1), cfg.selected_plane_idx);
    end

    % ---- meta -------------------------------------------------------------
    [~, sid] = fileparts(preprocessed_file);
    sid = regexprep(sid, '_preprocessed$', '');

    S = struct();
    S.meta = struct();
    S.meta.session_id         = sid;
    S.meta.source_file        = preprocessed_file;
    S.meta.created            = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    S.meta.git_commit         = local_git_commit(root);
    S.meta.ca_type            = cfg.ca_type;
    S.meta.ca_field           = ca_field;
    S.meta.selected_plane_idx = cfg.selected_plane_idx;
    S.meta.plane_z            = plane_z;
    S.meta.time_row           = time_row;
    S.meta.n_rois             = size(dff, 1);
    % plane-local ROI bookkeeping: row r of dff is global ROI roi_idx(r) with
    % centroid roi_centroid(r,:). Kept so downstream layers can verify an ROI
    % match (matched-cell centroid distance vs shuffled) without re-loading
    % the ~GB preprocessed file.
    S.meta.roi_index_global   = roi_idx(:);
    S.meta.roi_centroid       = B.CaData(1).Ca_centroid_voxel(roi_idx, :);
    S.meta.n_frames           = size(dff, 2);
    S.meta.fs_hz              = 1 / median(diff(time_vec), 'omitnan');
    S.meta.time_span          = [time_vec(1), time_vec(end)];
    S.meta.qc_passed          = false;     % set by qc_session
    S.meta.qc                 = struct();   % filled by qc_session

    fprintf('  plane %d (Z=%g): %d ROIs, %d frames, fs=%.2f Hz\n', ...
        cfg.selected_plane_idx, plane_z, S.meta.n_rois, S.meta.n_frames, S.meta.fs_hz);

    % ---- per-stimulus organization ----------------------------------------
    Stimuli = B.Stimuli;
    types   = unique({Stimuli.type}, 'stable');
    S.meta.stim_types = types;
    S.stimuli = struct();

    for ti = 1:numel(types)
        raw_type = types{ti};
        fld = matlab.lang.makeValidName(raw_type);
        plan = stim_plan(raw_type);
        fprintf('  - %s [%s]\n', raw_type, plan.paradigm);

        switch plan.paradigm
            case 'trial_onset'
                E = organize_trial_onset(Stimuli, dff, time_vec, B.Triggers, raw_type, plan);
            case 'direction'
                E = organize_direction(Stimuli, dff, time_vec, raw_type, plan);
            case 'block'
                E = organize_block(Stimuli, dff, time_vec, raw_type, plan);
            otherwise
                E = struct('type', raw_type, 'paradigm', 'unknown');
        end

        if isempty(E)
            fprintf('      (no usable data -- skipped)\n');
            continue;
        end
        E.stim_index = find(strcmp({Stimuli.type}, raw_type), 1);
        s0 = Stimuli(E.stim_index);
        if isfield(s0, 'RedFrameSynchronized')
            E.RedFrameSynchronized = logical(s0.RedFrameSynchronized);
        else
            E.RedFrameSynchronized = NaN;
        end
        E.n_rois = S.meta.n_rois;
        S.stimuli.(fld) = E;
    end
end

% ======================================================================
function E = organize_trial_onset(Stimuli, dff, time_vec, Triggers, raw_type, plan)
    E = [];
    idx = find(strcmp({Stimuli.type}, raw_type), 1);
    if isempty(idx); return; end
    s = Stimuli(idx);

    switch raw_type
        case 'looming'
            [Ca, t_com] = extract_looming_trials(Stimuli, dff, time_vec, Triggers);
            PT = Triggers.TimeProjector;
            onsets = PT(PT >= s.TimeStimulusFrame(1) & PT <= s.TimeStimulusFrame(end));
        case 'sparse_local_global_flashes'
            onsets = get_flash_onsets(Stimuli);
            Ca = align_trials_to_onsets(onsets, dff, time_vec, ...
                plan.pre, plan.post, plan.n_com);
            t_com = linspace(-plan.pre, plan.post, plan.n_com);
        otherwise
            return;
    end
    if isempty(Ca); return; end

    E = struct();
    E.type        = raw_type;
    E.paradigm    = 'trial_onset';
    E.t_com       = t_com;
    E.onsets      = onsets(:);
    E.trials      = Ca;                                  % rois x n_com x n_trials
    E.n_trials    = size(Ca, 3);
    E.avg         = median(Ca, 3, 'omitnan');            % rois x n_com (median convention)
    E.sem         = std(Ca, 0, 3, 'omitnan') ./ sqrt(max(E.n_trials, 1));
    E.peak_window = plan.peak_window;
end

% ======================================================================
function E = organize_direction(Stimuli, dff, time_vec, raw_type, plan)
    E = [];
    idx = find(strcmp({Stimuli.type}, raw_type), 1);
    if isempty(idx); return; end
    s = Stimuli(idx);

    [Ca, t_com, direction_of] = extract_direction_trials(Stimuli, dff, time_vec, ...
        raw_type, plan.orientation_field, plan.pre, plan.post, plan.n_com);
    if isempty(Ca); return; end

    % subtrial onsets (same anchor extract_direction_trials uses internally)
    sfc = s.Frameinfo.subtrial_frame_count(~isnan(s.Frameinfo.frame_count));
    subtrial_starts = find(sfc == 1);
    onsets = s.TimeStimulusFrame(subtrial_starts);

    [Ca_by_dir, dir_values] = chunk_by_direction(Ca, direction_of);
    peak_idx = t_com >= plan.peak_window(1) & t_com <= plan.peak_window(2);
    amp_by_dir = squeeze(mean(Ca_by_dir(:, peak_idx, :), 2, 'omitnan'));   % rois x n_dir
    if isrow(amp_by_dir); amp_by_dir = amp_by_dir'; end                    % guard 1-cell case
    [~, pref_dir_idx] = max(amp_by_dir, [], 2);

    n_rois = size(Ca_by_dir, 1);
    pref_resp = nan(n_rois, numel(t_com));
    for c = 1:n_rois
        pref_resp(c, :) = Ca_by_dir(c, :, pref_dir_idx(c));
    end

    % independent QC of the RNG direction reconstruction on this session
    [r_same, r_shifted] = verify_direction_reconstruction(Stimuli, dff, time_vec, ...
        raw_type, plan.orientation_field);

    E = struct();
    E.type          = raw_type;
    E.paradigm      = 'direction';
    E.t_com         = t_com;
    E.subtrials     = Ca;                     % rois x n_com x n_subtrials
    E.direction_of  = direction_of(:)';
    E.onsets        = onsets(:);
    E.dir_values    = dir_values(:)';
    E.by_dir        = Ca_by_dir;              % rois x n_com x n_dir (median of repeats)
    E.amp_by_dir    = amp_by_dir;             % rois x n_dir
    E.pref_dir_idx  = pref_dir_idx;
    E.pref_resp     = pref_resp;              % rois x n_com
    E.param_name    = plan.param_name;
    E.param_field   = plan.orientation_field;
    E.peak_window   = plan.peak_window;
    % 'applicable' is false when the stimulus has <2 repeat cycles (e.g.
    % grating shows each orientation exactly once, trials=1) -- the
    % cross-repeat correlation is undefined then, so verification is N/A,
    % NOT a failure. The RNG-replay mechanism itself is still validated on
    % stimuli that DO have repeats (moving_bar, 3 reps).
    applicable = ~isnan(r_same) && ~isnan(r_shifted);
    E.verify        = struct('r_same', r_same, 'r_shifted', r_shifted, ...
                             'applicable', applicable, ...
                             'pass', ~applicable || (r_same > r_shifted + 0.05));
end

% ======================================================================
function E = organize_block(Stimuli, dff, time_vec, raw_type, ~)
    E = [];
    resp = extract_full_stimulus_responses(Stimuli, dff, raw_type, time_vec);
    if isempty(resp); return; end
    block = resp{1};                          % rois x n_frames_block

    idx = find(strcmp({Stimuli.type}, raw_type), 1);
    s = Stimuli(idx);
    fidx = find(time_vec >= s.TimeStimulusFrame(1) & time_vec <= max(s.TimeStimulusFrame(:)));
    % block width is defined by extract_full_stimulus_responses; take the
    % matching real-seconds axis for it.
    if numel(fidx) >= size(block, 2)
        block_t = time_vec(fidx(1:size(block,2))) - time_vec(fidx(1));
    else
        dt = median(diff(time_vec), 'omitnan');
        block_t = (0:size(block,2)-1) * dt;
    end

    E = struct();
    E.type         = raw_type;
    E.paradigm     = 'block';
    E.block        = block;
    E.block_t      = block_t;
    E.duration_sec = block_t(end);
end

% ======================================================================
function plan = stim_plan(raw_type)
    % Central map: stimulus type -> paradigm + extraction parameters. One
    % place to add a new stimulus type.
    plan = struct('paradigm', 'block', 'orientation_field', '', ...
        'param_name', '', 'pre', 1, 'post', 6, 'n_com', 100, ...
        'peak_window', [0.5 4]);
    switch raw_type
        case 'looming'
            plan.paradigm = 'trial_onset';
            plan.peak_window = [0.5 1.5];
        case 'sparse_local_global_flashes'
            plan.paradigm = 'trial_onset';
            plan.pre = 1; plan.post = 4; plan.peak_window = [0 0.8];
        case 'grating'
            plan.paradigm = 'direction';
            plan.orientation_field = 'grating_orientations';
            plan.param_name = 'orientation';
            plan.pre = 1; plan.post = 6; plan.peak_window = [0.5 4];
        case 'moving_bar'
            plan.paradigm = 'direction';
            plan.orientation_field = 'bar_orientations';
            plan.param_name = 'orientation';
            plan.pre = 1; plan.post = 6; plan.peak_window = [0.5 4];
        case {'spontaneous', 'checkers2'}
            plan.paradigm = 'block';
        otherwise
            plan.paradigm = 'block';   % safe default for unknown types
    end
end

% ======================================================================
function h = local_git_commit(root)
    h = 'unknown';
    try
        old = pwd; cleaner = onCleanup(@() cd(old)); %#ok<NASGU>
        cd(root);
        [st, out] = system('git rev-parse --short HEAD');
        if st == 0; h = strtrim(out); end
    catch
    end
end
