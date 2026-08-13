function R = compare_sessions(src1, src2, roi_match_file, opts)
% COMPARE_SESSIONS  Thin pairwise layer over two cached *_organized.mat
% structs (see PIPELINE.md). Re-derives NOTHING about timing/chunking --
% it only maps matched cells and computes matched-cell raw/z metrics +
% modulation index per common stimulus, reusing the same metric functions
% (compute_baseline_zscore[_from_resp], compute_whole_block_zscore,
% compute_modulation_index) the validated pipeline uses.
%
%   R = compare_sessions(org1, org2, roi_match_file)
%   R = compare_sessions(org1, org2, roi_match_file, opts)
%
% src1/src2 : path to *_organized.mat OR an in-memory organized struct.
%             src1 must correspond to column 1 of the ROI-match mapping,
%             src2 to column 2 (same ordering roiMatchPub used to build it).
% opts.verbose (default true)
%
% Returns R.<stim>.{raw_base,raw_drug,z_base,z_drug,mi_raw,mi_z,p_signrank,n_matched}
% and R.meta.

    if nargin < 4 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'verbose'); opts.verbose = true; end

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(fileparts(here), 'analysis_functions'));

    S1 = load_org(src1);
    S2 = load_org(src2);

    if ~S1.meta.qc_passed || ~S2.meta.qc_passed
        warning('compare_sessions: one/both sessions did not pass QC (%s: %d, %s: %d)', ...
            S1.meta.session_id, S1.meta.qc_passed, S2.meta.session_id, S2.meta.qc_passed);
    end

    % ---- matched ROI indices (per-plane-local; same convention as Preprocess_Stimulus_Data) ----
    M = load(roi_match_file);
    assert(isfield(M, 'roiMatchData') && isfield(M.roiMatchData, 'allSessionMapping'), ...
        'roi_match_file missing roiMatchData.allSessionMapping');
    idx1 = M.roiMatchData.allSessionMapping(:, 1);
    idx2 = M.roiMatchData.allSessionMapping(:, 2);
    valid = idx1 >= 1 & idx1 <= S1.meta.n_rois & idx2 >= 1 & idx2 <= S2.meta.n_rois;
    idx1 = idx1(valid); idx2 = idx2(valid);
    n_matched = numel(idx1);
    assert(n_matched > 0, 'No matched ROIs fall within the selected plane of both sessions');

    R = struct();
    R.meta = struct('session1', S1.meta.session_id, 'session2', S2.meta.session_id, ...
        'roi_match_file', roi_match_file, 'n_matched', n_matched, ...
        'qc1', S1.meta.qc_passed, 'qc2', S2.meta.qc_passed, ...
        'idx1', idx1(:), 'idx2', idx2(:));   % plane-local matched ROI rows,
        % kept so a figure can go back to the organized structs (e.g. per-repeat
        % subtrials) for the SAME cells without redoing the matching.

    if opts.verbose
        fprintf('\ncompare_sessions: %s  vs  %s  (%d matched cells)\n', ...
            S1.meta.session_id, S2.meta.session_id, n_matched);
    end

    common = intersect(fieldnames(S1.stimuli), fieldnames(S2.stimuli));
    for i = 1:numel(common)
        nm = common{i};
        E1 = S1.stimuli.(nm); E2 = S2.stimuli.(nm);
        if ~strcmp(E1.paradigm, E2.paradigm); continue; end

        switch E1.paradigm
            case 'trial_onset'
                res = cmp_trial_onset(E1, E2, idx1, idx2);
            case 'direction'
                res = cmp_direction(E1, E2, idx1, idx2);
            case 'block'
                res = cmp_block(E1, E2, idx1, idx2);
            otherwise
                continue;
        end
        res.n_matched = n_matched;
        res.stim_label = E1.type;
        res.paradigm   = E1.paradigm;
        R.(nm) = res;
        if opts.verbose
            fprintf('  %-28s  MI(z) med=%+.3f  signrank p=%.3g\n', E1.type, ...
                median(res.mi_z, 'omitnan'), res.p_signrank);
        end
    end
end

% ----------------------------------------------------------------------
function res = cmp_trial_onset(E1, E2, i1, i2)
    t = E1.t_com; pw = E1.peak_window;
    pk = t >= pw(1) & t <= pw(2);
    raw_base = mean(E1.avg(i1, pk), 2, 'omitnan');
    raw_drug = mean(E2.avg(i2, pk), 2, 'omitnan');
    z_base = compute_baseline_zscore(E1.trials, i1, t, pw);
    z_drug = compute_baseline_zscore(E2.trials, i2, t, pw);
    res = pack(raw_base, raw_drug, z_base, z_drug);
    res = add_traces(res, t, pw, E1.avg(i1, :), E2.avg(i2, :));
end

% ----------------------------------------------------------------------
function res = cmp_direction(E1, E2, i1, i2)
    % Mirror analyze_direction_stimulus's matched logic: preferred direction
    % from session1 (baseline), carried to session2, matched by angle.
    t = E1.t_com; pw = E1.peak_window;
    pk = t >= pw(1) & t <= pw(2);

    by1 = E1.by_dir(i1, :, :); d1 = E1.dir_values;
    by2 = E2.by_dir(i2, :, :); d2 = E2.dir_values;
    common = intersect(d1, d2);
    [~, ia] = ismember(common, d1); [~, ib] = ismember(common, d2);
    by1 = by1(:, :, ia); by2 = by2(:, :, ib);

    amp1 = squeeze(mean(by1(:, pk, :), 2, 'omitnan'));   % cells x dirs
    if isrow(amp1); amp1 = amp1'; end
    [~, pref] = max(amp1, [], 2);

    n = size(by1, 1);
    resp1 = nan(n, numel(t)); resp2 = nan(n, numel(t));
    for c = 1:n
        resp1(c, :) = by1(c, :, pref(c));
        resp2(c, :) = by2(c, :, pref(c));
    end
    raw_base = mean(resp1(:, pk), 2, 'omitnan');
    raw_drug = mean(resp2(:, pk), 2, 'omitnan');
    z_base = compute_baseline_zscore_from_resp(resp1, t, pw);
    z_drug = compute_baseline_zscore_from_resp(resp2, t, pw);
    res = pack(raw_base, raw_drug, z_base, z_drug);
    res = add_traces(res, t, pw, resp1, resp2);
    res.dir_values_common = common(:)';
    res.pref_dir_idx      = pref(:);
end

% ----------------------------------------------------------------------
function res = cmp_block(E1, E2, i1, i2)
    raw_base = mean(E1.block(i1, :), 2, 'omitnan');
    raw_drug = mean(E2.block(i2, :), 2, 'omitnan');
    z_base = compute_whole_block_zscore(E1.block, i1);
    z_drug = compute_whole_block_zscore(E2.block, i2);
    res = pack(raw_base, raw_drug, z_base, z_drug);
end

% ----------------------------------------------------------------------
function res = pack(raw_base, raw_drug, z_base, z_drug)
    res = struct();
    res.raw_base = raw_base; res.raw_drug = raw_drug;
    res.z_base   = z_base;   res.z_drug   = z_drug;
    res.mi_raw = compute_modulation_index(raw_base, raw_drug);
    res.mi_z   = compute_modulation_index(z_base, z_drug);
    valid = ~isnan(z_base) & ~isnan(z_drug);
    if nnz(valid) >= 2
        res.p_signrank = signrank(z_base(valid), z_drug(valid));
    else
        res.p_signrank = NaN;
    end
end

% ----------------------------------------------------------------------
function res = add_traces(res, t, pw, resp_base, resp_drug)
% Attach the matched-cell response traces the metrics above were computed
% from, so a group layer can average time courses across recordings without
% re-deriving anything (same cells, same preferred-direction choice).
    res.t_com       = t(:)';
    res.peak_window = pw;
    res.resp_base   = resp_base;      % n_matched x n_com
    res.resp_drug   = resp_drug;
end

% ----------------------------------------------------------------------
function S = load_org(src)
    if ischar(src) || isstring(src)
        S = load(char(src));
    else
        S = src;
    end
    assert(isfield(S, 'stimuli') && isfield(S, 'meta'), 'Not an organized-session struct');
end
