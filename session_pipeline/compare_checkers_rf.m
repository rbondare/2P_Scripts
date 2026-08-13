function R = compare_checkers_rf(rf1, rf2, roi_match_file, opts)
% COMPARE_CHECKERS_RF  Baseline vs drug checkers2 receptive fields, matched cells.
%
%   R = compare_checkers_rf(rf1, rf2, roi_match_file)
%   R = compare_checkers_rf(rf1, rf2, roi_match_file, opts)
%
%   rf1/rf2        outputs of checkers_rf_session, or paths to saved ones.
%                  rf1 must be column 1 of the ROI-match mapping (baseline).
%   roi_match_file roiMatch_plane<P>_<date>_<t1>_<t2>_V1.mat
%
% opts:
%   .snr_min     dB above the circular-shift null for a cell to count as having
%                an RF (default 3). Applied to the SUM of the two conditions'
%                evidence, never to one side -- selecting on baseline alone and
%                then testing the drug manufactures a difference.
%   .require_gate  refuse to report a comparison when either session failed its
%                own null gate (default true). Turn this off only to inspect.
%
% This is deliberately the same shape as compare_sessions.m: it computes, it does
% not plot (checkers_rf_plots.m does that), so a figure cannot disagree with the
% numbers.
%
% THE GATE EXISTS BECAUSE THE OBVIOUS FAILURE MODE IS SILENT. Every cell has a
% peak somewhere in its RF volume whether or not it has a receptive field, and
% those peaks have positions, widths and latencies that will happily differ
% between two sessions and produce small p-values. If checkers_rf_session's null
% gate says a session has no recoverable RFs, the per-cell numbers below are
% descriptions of noise and the comparison means nothing.

    if nargin < 4 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'snr_min');      opts.snr_min = 3;         end
    if ~isfield(opts, 'require_gate'); opts.require_gate = true; end

    A = load_rf(rf1);
    Bs = load_rf(rf2);

    % ---- ROI matching (same convention as compare_sessions) ----------------
    M = load(roi_match_file);
    assert(isfield(M, 'roiMatchData') && isfield(M.roiMatchData, 'allSessionMapping'), ...
        'roi_match_file missing roiMatchData.allSessionMapping');
    idx1 = M.roiMatchData.allSessionMapping(:, 1);
    idx2 = M.roiMatchData.allSessionMapping(:, 2);
    n1 = numel(A.props.snr_db); n2 = numel(Bs.props.snr_db);
    valid = idx1 >= 1 & idx1 <= n1 & idx2 >= 1 & idx2 <= n2;
    idx1 = idx1(valid); idx2 = idx2(valid);
    n_matched = numel(idx1);
    assert(n_matched > 0, 'no valid matched ROIs -- wrong plane for this roi_match_file?');

    R = struct();
    R.meta = struct('session1', A.meta.session_id, 'session2', Bs.meta.session_id, ...
        'plane1', A.meta.selected_plane_idx, 'plane2', Bs.meta.selected_plane_idx, ...
        'roi_match_file', string(roi_match_file), 'n_matched', n_matched, ...
        'idx1', idx1(:), 'idx2', idx2(:), 'snr_min', opts.snr_min, ...
        'ca_field', A.meta.ca_field, 'built', string(datetime('now')));

    % ---- falsifiable check that the mapping is the right one ---------------
    % An ROI-match file carries no session provenance, so applying one at the
    % wrong plane silently pairs unrelated cells and still returns tiny p-values.
    % Same check compare_group runs: real pairs sit on top of each other, shuffled
    % pairs are a FOV apart.
    c1 = A.roi_centroid(idx1, 1:2); c2 = Bs.roi_centroid(idx2, 1:2);
    d_real = sqrt(sum((c1 - c2).^2, 2));
    sh = randperm(n_matched);
    d_shuf = sqrt(sum((c1 - c2(sh, :)).^2, 2));
    R.match_check = struct('median_px', median(d_real), 'p90_px', prctile(d_real, 90), ...
        'median_shuffled_px', median(d_shuf), ...
        'corr_x', corr(c1(:,1), c2(:,1)), 'corr_y', corr(c1(:,2), c2(:,2)));
    R.match_check.pass = median(d_real) < 0.25 * median(d_shuf);
    if R.match_check.pass; mverdict = 'PASS'; else; mverdict = 'FAIL'; end
    fprintf(['ROI match: %d pairs, centroid distance %.1f px (90th %.1f) vs %.1f ' ...
             'shuffled, corr x %.4f y %.4f -> %s\n'], n_matched, ...
        R.match_check.median_px, R.match_check.p90_px, R.match_check.median_shuffled_px, ...
        R.match_check.corr_x, R.match_check.corr_y, mverdict);

    % ---- the null gate -----------------------------------------------------
    R.gate = struct('session1', A.gate, 'session2', Bs.gate, ...
        'both_have_rfs', A.gate.has_rfs && Bs.gate.has_rfs);
    fprintf('NULL GATE  %s: ratio %.3f p=%.2g | %s: ratio %.3f p=%.2g\n', ...
        A.meta.session_id, A.gate.median_ratio, A.gate.p_signrank, ...
        Bs.meta.session_id, Bs.gate.median_ratio, Bs.gate.p_signrank);
    if ~R.gate.both_have_rfs
        msg = sprintf(['checkers2 RFs do not clear the circular-shift null in at ' ...
            'least one session, so there is nothing to compare: any difference in ' ...
            'RF position, size or latency below would be a difference between two ' ...
            'noise fields.']);
        if opts.require_gate
            R.aborted = true; R.reason = msg;
            fprintf(2, '\nABORTED: %s\n', msg);
            return;
        else
            warning('compare_checkers_rf:noRFs', '%s', msg);
        end
    end
    R.aborted = false;

    % ---- per-cell properties, matched --------------------------------------
    f = intersect(fieldnames(A.props), fieldnames(Bs.props));
    for i = 1:numel(f)
        R.base.(f{i}) = A.props.(f{i})(idx1);
        R.drug.(f{i}) = Bs.props.(f{i})(idx2);
    end

    % Selection on the SUM of both sides' evidence, so it cannot favour either.
    R.good = (R.base.snr_db + R.drug.snr_db) / 2 >= opts.snr_min;
    fprintf('%d/%d matched cells clear %g dB (mean of the two conditions)\n', ...
        sum(R.good), n_matched, opts.snr_min);

    % ---- deltas + paired stats ---------------------------------------------
    R.delta = struct();
    R.delta.d_center_deg = sqrt((R.drug.az_deg - R.base.az_deg).^2 + ...
                                (R.drug.el_deg - R.base.el_deg).^2);
    R.delta.d_diameter_deg = R.drug.diameter_deg - R.base.diameter_deg;
    R.delta.d_peak_lag_s   = R.drug.peak_lag_s   - R.base.peak_lag_s;
    R.delta.d_ttp_s        = R.drug.ttp_s        - R.base.ttp_s;
    R.delta.gain_ratio     = abs(R.drug.peak_amp) ./ abs(R.base.peak_amp);
    R.delta.mod_index      = (abs(R.drug.peak_amp) - abs(R.base.peak_amp)) ./ ...
                             (abs(R.drug.peak_amp) + abs(R.base.peak_amp));

    keys = {'peak_amp_abs', 'diameter_deg', 'peak_lag_s', 'ttp_s', 'snr_db'};
    R.stats = struct();
    for i = 1:numel(keys)
        k = keys{i};
        if strcmp(k, 'peak_amp_abs')
            a = abs(R.base.peak_amp); b = abs(R.drug.peak_amp);
        else
            a = R.base.(k); b = R.drug.(k);
        end
        m = R.good & isfinite(a) & isfinite(b);
        st = struct('n', sum(m), 'median_base', median(a(m)), 'median_drug', median(b(m)), ...
                    'p', NaN);
        if sum(m) > 5 && any(a(m) ~= b(m)); st.p = signrank(a(m), b(m)); end
        R.stats.(k) = st;
        fprintf('  %-14s n=%4d  %8.4g -> %8.4g   p=%.3g\n', k, st.n, ...
            st.median_base, st.median_drug, st.p);
    end
end


% ======================================================================
function X = load_rf(x)
    if ischar(x) || isstring(x)
        assert(isfile(x), 'File not found: %s', x);
        L = load(x);
        fn = fieldnames(L);
        if isfield(L, 'R'); X = L.R; else; X = L.(fn{1}); end
    else
        X = x;
    end
    assert(isfield(X, 'props') && isfield(X, 'gate'), ...
        'not a checkers_rf_session output (missing .props/.gate)');
end
