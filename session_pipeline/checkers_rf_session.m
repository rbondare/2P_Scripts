function R = checkers_rf_session(preprocessed_file, cfg)
% CHECKERS_RF_SESSION  Receptive fields from the checkers2 white-noise block.
%
%   R = checkers_rf_session(preprocessed_file)
%   R = checkers_rf_session(preprocessed_file, cfg)
%
% checkers2 is the one stimulus organize_session cannot really "organize": it has
% trials==1 and no internal trial structure (one continuous 300 s block), so the
% block paradigm just stores the raw trace. The question it actually answers is a
% reverse-correlation one, which is what this does.
%
% ESTIMATOR. Calcium is far slower than the 5 Hz pattern rate, so this is not a
% spike-triggered average -- it is the stimulus/response covariance
%
%       RF(x, y, tau) = < [S(x,y,p) - <S(x,y)>] * [R(t_p + tau) - <R>] >
%
% which is the same quantity a spike-triggered average estimates
% (STA - <S> = Cov(S,R)/<R>) but written for a continuous response. Because S is
% +/-1 with unit variance, RF is directly a regression slope: dF/F per unit
% contrast. It is NOT normalised per cell, so amplitudes stay comparable between
% conditions -- which is the whole point when comparing baseline vs drug.
%
% <S(x,y)> is subtracted per element rather than assumed flat. Each pattern is a
% permutation of a fixed 306-bright/306-dark multiset, so the spatial mean of
% every frame is exactly 0, but the per-element temporal mean is only 0 to within
% 1/sqrt(n_patterns) (~0.03, up to ~0.09 for the most extreme of 612 elements).
%
% THE NULL IS COMPUTED, NOT ASSUMED. Every RF is calculated twice -- once with
% the real pattern order and once with the pattern sequence circularly shifted by
% cfg.null_shift. The shifted version preserves the stimulus statistics, the
% response, the lag axis and the sample count exactly and breaks only the pairing
% between stimulus and response, so it is this cell's true noise floor. snr_db is
% peak/peak against it, which means 0 dB is chance BY CONSTRUCTION -- already
% corrected for the peak being a maximum over ~10k (element x lag) samples, which
% a border or negative-lag noise estimate is not. R.gate carries the population
% verdict; if gate.has_rfs is false, nothing downstream of this is meaningful.
% Negative lags are kept as a second, independent null (stimulus that arrived
% after the response cannot have caused it), so cfg.lags must retain some.
%
% cfg fields (all optional):
%   .ca_type            1=Ca_dFF (default), 2=Ca_deconvolved, 3=Ca_F.
%                       dF/F is the pipeline default; deconvolved gives a much
%                       sharper temporal filter but its own kinetics assumptions.
%   .selected_plane_idx imaging plane, 1-based (default 1). Must be the plane the
%                       ROI-match file was built on -- see group_config.m.
%   .lags               lag axis in seconds (default -1:0.2:4)
%   .fs_range           plausible volume-rate bounds (default [2 30])
%
% Returns R with
%   .rf         [n_el_rows x n_el_cols x n_lags x n_rois] single, dF/F per contrast
%   .lags, .az_deg, .el_deg, .roi_index_global, .roi_centroid
%   .props      per-ROI table of spatial + temporal properties (see rf_properties)
%   .meta       provenance + everything needed to redo this
%
% Only the ON-SCREEN elements are used: the pattern carries a one-element border
% that is never displayed (it is the margin the sub-element jitter slides into).
%
% KNOWN LIMITATION -- the stimulus jitters its grid by up to one full element
% (8 deg) between patterns (subshift='pixel'; see reconstruct_checkers2). The
% jitter is independent of pattern content, so it does not bias RF position; it
% blurs the estimate by a uniform box one element wide. C.shift is carried
% through in .meta if you ever want to accumulate on a sub-element grid instead.

    if nargin < 2 || isempty(cfg); cfg = struct(); end
    if ~isfield(cfg, 'ca_type');            cfg.ca_type = 1;             end
    if ~isfield(cfg, 'selected_plane_idx'); cfg.selected_plane_idx = 1;  end
    if ~isfield(cfg, 'lags');               cfg.lags = -1:0.2:4;         end
    if ~isfield(cfg, 'fs_range');           cfg.fs_range = [2 30];       end
    if ~isfield(cfg, 'null_shift');         cfg.null_shift = 700;        end
    if ~isfield(cfg, 'keep_null');          cfg.keep_null = false;       end

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);                       % 2P_Scripts/
    addpath(root);                                % reconstruct_checkers2

    assert(isfile(preprocessed_file), 'File not found: %s', preprocessed_file);
    fprintf('\n=== checkers_rf_session: %s ===\n', preprocessed_file);
    B = load(preprocessed_file);

    % ---- calcium field + plane (same resolution rules as organize_session) --
    ca_fields = {'Ca_dFF', 'Ca_deconvolved', 'Ca_F'};
    ca_field  = ca_fields{cfg.ca_type};
    assert(isfield(B.CaData(1), ca_field), 'CaData(1).%s not present', ca_field);

    centroidZ = B.CaData(1).Ca_centroid_voxel(:, 3);
    planes    = unique(centroidZ);
    assert(cfg.selected_plane_idx <= numel(planes), ...
        'Requested plane %d but only %d present', cfg.selected_plane_idx, numel(planes));
    plane_z  = planes(cfg.selected_plane_idx);
    roi_idx  = find(centroidZ == plane_z);
    dff      = double(B.CaData(1).(ca_field)(roi_idx, :));

    if isvector(B.TimeCa)
        time_vec = B.TimeCa(:)'; time_row = 1;
    elseif size(B.TimeCa, 1) >= cfg.selected_plane_idx
        time_row = cfg.selected_plane_idx;            % this plane's own row
        time_vec = B.TimeCa(time_row, :);
    else
        time_row = 1; time_vec = B.TimeCa(1, :);
    end
    assert(numel(time_vec) == size(dff, 2), ...
        'TimeCa row %d has %d samples but calcium has %d frames', ...
        time_row, numel(time_vec), size(dff, 2));

    fs = 1 / median(diff(time_vec), 'omitnan');
    assert(fs >= cfg.fs_range(1) && fs <= cfg.fs_range(2), ...
        'volume rate %.2f Hz outside plausible range [%g %g]', fs, cfg.fs_range);

    % ---- the stimulus ------------------------------------------------------
    types = arrayfun(@(x) string(x.type), B.Stimuli);
    k = find(types == "checkers2", 1);
    assert(~isempty(k), 'no checkers2 block in %s', preprocessed_file);
    s = B.Stimuli(k);
    assert(s.RedFrameSynchronized, ...
        ['checkers2 RedFrameSynchronized is false: TimeStimulusFrame would be on ' ...
         'the un-synced stimulus-PC clock and every lag below would be wrong.']);

    C = reconstruct_checkers2(s);
    fprintf('  %d patterns at %.2f Hz, %d x %d visible elements (%.1f x %.1f deg)\n', ...
        C.meta.n_patterns, C.meta.pattern_change_rate, ...
        numel(C.vis_rows), numel(C.vis_cols), C.meta.deg_per_row, C.meta.deg_per_col);

    % ---- design matrix: visible elements only ------------------------------
    S = C.contrast(C.vis_rows, C.vis_cols, :);        % [nr x nc x n_pat]
    [nr, nc, n_pat] = size(S);
    S = reshape(S, nr*nc, n_pat);                     % [n_el x n_pat]
    S = double(S) - mean(double(S), 2);               % subtract <S> per element

    lags   = cfg.lags(:)';
    n_lag  = numel(lags);
    n_roi  = size(dff, 1);
    onset  = C.onset_t(:)';

    % Keep only patterns whose whole lag window is inside the recorded calcium.
    ok = (onset + min(lags)) >= time_vec(1) & (onset + max(lags)) <= time_vec(end);
    assert(any(ok), 'no pattern has its full lag window inside the calcium record');
    if ~all(ok)
        fprintf('  dropping %d/%d patterns whose lag window runs off the recording\n', ...
            sum(~ok), n_pat);
    end
    S = S(:, ok); onset_ok = onset(ok); n_use = sum(ok);

    fprintf('  %d ROIs in plane %d, %.2f Hz volume rate, %d lags %.1f..%.1f s\n', ...
        n_roi, cfg.selected_plane_idx, fs, n_lag, lags(1), lags(end));

    % ---- reverse correlation ----------------------------------------------
    % One interpolation + one GEMM per lag. dff' is [n_frames x n_roi] so interp1
    % resamples every ROI at once.
    %
    % Everything is computed twice: once with the real pattern order and once
    % with the pattern sequence circularly shifted. The shifted version keeps the
    % stimulus statistics, the response, the lag axis and the number of samples
    % byte for byte and breaks ONLY the pairing between them, so it is the exact
    % noise floor of this estimator for this cell. A per-cell SNR measured against
    % it cannot be inflated by slow drift, by bursty cells or by the fact that the
    % peak is a maximum over ~10k (element x lag) samples -- all of which a border-
    % or negative-lag-only noise estimate gets wrong.
    null_shift = round(cfg.null_shift);
    assert(mod(null_shift, n_use) > 20 || mod(n_use - null_shift, n_use) > 20, ...
        'cfg.null_shift = %d leaves the null nearly aligned with the real order', null_shift);
    Snull = circshift(S, [0 null_shift]);
    rf = zeros(nr*nc, n_lag, n_roi, 'single');
    rfn = zeros(nr*nc, n_lag, n_roi, 'single');
    dffT = dff';
    for i = 1:n_lag
        Rq = interp1(time_vec, dffT, onset_ok + lags(i), 'linear');   % [n_use x n_roi]
        Rq = Rq - mean(Rq, 1, 'omitnan');
        Rq(~isfinite(Rq)) = 0;
        rf(:, i, :)  = single((S     * Rq) / n_use);                  % [n_el x n_roi]
        rfn(:, i, :) = single((Snull * Rq) / n_use);
    end
    rf  = reshape(rf,  nr, nc, n_lag, n_roi);
    rfn = reshape(rfn, nr, nc, n_lag, n_roi);

    % ---- package -----------------------------------------------------------
    R = struct();
    R.rf     = rf;
    R.lags   = lags;
    R.az_deg = C.az_deg;
    R.el_deg = C.el_deg;
    R.roi_index_global = roi_idx(:)';
    R.roi_centroid     = B.CaData(1).Ca_centroid_voxel(roi_idx, :);
    R.props  = rf_properties(rf, rfn, lags, C.az_deg, C.el_deg);
    if cfg.keep_null; R.rf_null = rfn; end

    % Population gate. If the real peaks are not above the shifted-stimulus
    % peaks, this session has no recoverable receptive fields and nothing
    % downstream of here means anything.
    R.gate = struct();
    R.gate.null_shift    = null_shift;
    R.gate.peak_real     = R.props.peak_abs;
    R.gate.peak_null     = R.props.peak_abs_null;
    R.gate.median_ratio  = median(R.props.peak_abs ./ R.props.peak_abs_null, 'omitnan');
    R.gate.p_signrank    = signrank(R.props.peak_abs, R.props.peak_abs_null);
    R.gate.frac_above_3dB = mean(R.props.snr_db >= 3);
    R.gate.has_rfs = R.gate.median_ratio > 1.15 && R.gate.p_signrank < 0.05;
    if R.gate.has_rfs; verdict = 'RFs present'; else; verdict = 'NO RECOVERABLE RFs'; end
    fprintf(['  NULL GATE: real/null peak ratio %.3f (p=%.2g), %.1f%% of cells ' ...
             '>= 3 dB  ->  %s\n'], R.gate.median_ratio, R.gate.p_signrank, ...
             100*R.gate.frac_above_3dB, verdict);

    [~, sess, ~] = fileparts(preprocessed_file);
    R.meta = struct('session_id', erase(string(sess), "_preprocessed"), ...
        'source_file', string(preprocessed_file), 'ca_field', string(ca_field), ...
        'selected_plane_idx', cfg.selected_plane_idx, 'plane_z', plane_z, ...
        'time_row', time_row, 'fs_hz', fs, 'n_rois', n_roi, ...
        'n_patterns_used', n_use, 'n_patterns_total', n_pat, ...
        'stim_index', k, 'lags', lags, 'stim', C.meta, 'shift', C.shift, ...
        'built', string(datetime('now')));

    fprintf('  done: %d/%d ROIs above 3 dB vs the null\n', ...
        sum(R.props.snr_db >= 3), n_roi);
end


% ======================================================================
function P = rf_properties(rf, rfn, lags, az_deg, el_deg)
% Per-ROI spatial and temporal descriptors.
%
% SNR is measured against the circular-shift NULL (rfn), not against a border or
% the negative lags: the null is the identical statistic computed on the same
% cell with the stimulus/response pairing broken, so `peak_abs / peak_abs_null`
% is already corrected for the fact that the peak is a maximum over ~10k samples.
% 0 dB is chance by construction.

    [nr, nc, n_lag, n_roi] = size(rf);
    pre = lags < 0;
    assert(any(pre), ['cfg.lags has no negative entries: they are a second, ' ...
                      'independent null (stimulus after the response).']);
    post_idx = find(lags >= 0);

    z = @(n) nan(n, 1);
    P = struct('snr_db', z(n_roi), 'peak_abs', z(n_roi), 'peak_abs_null', z(n_roi), ...
               'peak_amp', z(n_roi), 'polarity', z(n_roi), ...
               'peak_lag_s', z(n_roi), 'az_deg', z(n_roi), 'el_deg', z(n_roi), ...
               'sigma_az_deg', z(n_roi), 'sigma_el_deg', z(n_roi), ...
               'diameter_deg', z(n_roi), 'noise', z(n_roi), ...
               'ttp_s', z(n_roi), 'fwhm_s', z(n_roi), 'biphasic_index', z(n_roi), ...
               'peak_row', z(n_roi), 'peak_col', z(n_roi));

    [AZ, EL] = meshgrid(az_deg, el_deg);
    for u = 1:n_roi
        A = double(rf(:, :, :, u));
        noise = std(reshape(double(rfn(:, :, post_idx, u)), [], 1), 0, 'omitnan');
        if ~isfinite(noise) || noise == 0; continue; end

        post = A(:, :, post_idx);
        [pk, li] = max(abs(post(:)));
        [ri, ci, ti] = ind2sub(size(post), li);
        val = post(ri, ci, ti);
        pkn = max(abs(reshape(double(rfn(:, :, post_idx, u)), [], 1)));

        P.noise(u)         = noise;
        P.peak_abs(u)      = pk;
        P.peak_abs_null(u) = pkn;
        P.snr_db(u)        = 20*log10(pk / pkn);
        P.peak_amp(u)   = val;
        P.polarity(u)   = sign(val);
        P.peak_lag_s(u) = lags(post_idx(ti));
        P.peak_row(u)   = ri;  P.peak_col(u) = ci;

        % --- spatial: second moments of the supra-threshold same-sign lobe
        M = post(:, :, ti) * sign(val);            % peak now positive
        M = M - median(M(:));
        m = M >= 0.5*max(M(:));                    % half-max footprint
        if nnz(m) >= 3
            w = M .* m; sw = sum(w(:));
            cx = sum(w(:).*AZ(:))/sw;  cy = sum(w(:).*EL(:))/sw;
            sx = sqrt(max(sum(w(:).*(AZ(:)-cx).^2)/sw, 0));
            sy = sqrt(max(sum(w(:).*(EL(:)-cy).^2)/sw, 0));
            P.az_deg(u) = cx;  P.el_deg(u) = cy;
            P.sigma_az_deg(u) = sx;  P.sigma_el_deg(u) = sy;
            P.diameter_deg(u) = 2*sqrt(sx*sy);
        else
            P.az_deg(u) = AZ(ri, ci);  P.el_deg(u) = EL(ri, ci);
        end

        % --- temporal: trace at the peak element (3x3 mean), full lag axis
        rr = max(ri-1,1):min(ri+1,nr);  cc = max(ci-1,1):min(ci+1,nc);
        tr = squeeze(mean(mean(A(rr, cc, :), 1), 2))' * sign(val);
        [tpk, tix] = max(tr);
        P.ttp_s(u) = lags(tix);
        half = tpk/2; j = tix;
        while j > 1 && tr(j) > half; j = j-1; end
        k2 = tix;
        while k2 < n_lag && tr(k2) > half; k2 = k2+1; end
        if k2 > j; P.fwhm_s(u) = lags(k2) - lags(j); end
        after = tr(tix:end);
        P.biphasic_index(u) = abs(min([min(after) 0])) / max(tpk, eps);
    end
end
