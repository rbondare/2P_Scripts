function [S, qc] = qc_session(S, B, cfg)
% QC_SESSION  Quality-control a canonical session struct (see PIPELINE.md).
%
% Two tiers:
%   Tier 1 (FATAL, data integrity)   -> assert(): halts this session.
%   Tier 2 (SCIENTIFIC VALIDITY)     -> recorded pass/fail flags; a failure
%                                       sets S.meta.qc_passed=false but does
%                                       NOT halt (so a batch can complete and
%                                       the QC figure can be inspected).
%
%   [S, qc] = qc_session(S, B, cfg)
%
% cfg fields used: .selected_plane_idx, .ca_type, .fs_range,
%   .nan_frac_tol (default 0.01). Returns S with S.meta.qc / qc_passed set,
% and the qc struct (qc.checks: struct array; qc.passed: overall logical).

    if nargin < 3 || isempty(cfg); cfg = struct(); end
    if ~isfield(cfg, 'selected_plane_idx'); cfg.selected_plane_idx = S.meta.selected_plane_idx; end
    if ~isfield(cfg, 'ca_type');            cfg.ca_type = S.meta.ca_type; end
    if ~isfield(cfg, 'fs_range');           cfg.fs_range = [2 30]; end
    if ~isfield(cfg, 'nan_frac_tol');       cfg.nan_frac_tol = 0.01; end

    fprintf('  QC: %s\n', S.meta.session_id);

    % Independently re-derive the raw plane data (QC must not trust the
    % organize outputs for integrity).
    ca_fields = {'Ca_dFF', 'Ca_deconvolved', 'Ca_F'};
    ca_field  = ca_fields{cfg.ca_type};
    dff_full  = B.CaData(1).(ca_field);
    centroidZ = B.CaData(1).Ca_centroid_voxel(:, 3);
    unique_planes = unique(centroidZ);
    plane_z = unique_planes(cfg.selected_plane_idx);
    dff = dff_full(centroidZ == plane_z, :);
    if isvector(B.TimeCa)
        time_vec = B.TimeCa(:)';
    elseif size(B.TimeCa,1) >= cfg.selected_plane_idx
        time_vec = B.TimeCa(cfg.selected_plane_idx, :);
    else
        time_vec = B.TimeCa(1, :);
    end
    t0 = time_vec(1); tN = time_vec(end);

    checks = struct('name', {}, 'tier', {}, 'pass', {}, 'detail', {});

    % ---------- TIER 1: FATAL integrity (assert) ----------
    add = @(n, t, p, d) add_check(n, t, p, d);

    ok = (numel(time_vec) == size(dff, 2));
    checks(end+1) = add('time_len_matches_frames', 1, ok, ...
        sprintf('time=%d frames=%d', numel(time_vec), size(dff,2)));
    assert(ok, 'QC FATAL [%s]: TimeCa length (%d) ~= dff frames (%d)', ...
        S.meta.session_id, numel(time_vec), size(dff,2));

    fs = S.meta.fs_hz;
    ok = fs >= cfg.fs_range(1) && fs <= cfg.fs_range(2);
    checks(end+1) = add('fs_plausible', 1, ok, sprintf('fs=%.3f Hz (allowed %g-%g)', fs, cfg.fs_range(1), cfg.fs_range(2)));
    assert(ok, 'QC FATAL [%s]: fs=%.3f Hz outside [%g %g]', S.meta.session_id, fs, cfg.fs_range(1), cfg.fs_range(2));

    nan_frac = mean(~isfinite(dff(:)));
    ok = nan_frac <= cfg.nan_frac_tol;
    checks(end+1) = add('raw_finite', 1, ok, sprintf('%.3f%% non-finite (tol %.1f%%)', 100*nan_frac, 100*cfg.nan_frac_tol));
    assert(ok, 'QC FATAL [%s]: %.2f%% non-finite dF/F exceeds tol %.2f%%', S.meta.session_id, 100*nan_frac, 100*cfg.nan_frac_tol);

    % TimeStimulusFrame monotonic + within recording span, per stimulus block
    fns = fieldnames(S.stimuli);
    for i = 1:numel(fns)
        E = S.stimuli.(fns{i});
        s = B.Stimuli(E.stim_index);
        tsf = s.TimeStimulusFrame(:);
        mono = all(diff(tsf) >= -1e-6);
        inspan = tsf(1) >= t0 - 1e-6 && max(tsf) <= tN + 1e-6;
        ok = mono && inspan;
        checks(end+1) = add(sprintf('tsf_valid_%s', fns{i}), 1, ok, ...
            sprintf('mono=%d inspan=%d [%.2f..%.2f] rec[%.2f..%.2f]', mono, inspan, tsf(1), max(tsf), t0, tN)); %#ok<AGROW>
        assert(ok, 'QC FATAL [%s/%s]: TimeStimulusFrame invalid (mono=%d inspan=%d)', ...
            S.meta.session_id, fns{i}, mono, inspan);
    end

    % ---------- TIER 2: scientific validity (flags) ----------
    for i = 1:numel(fns)
        E = S.stimuli.(fns{i});
        nm = fns{i};

        % red-frame sync
        rs = E.RedFrameSynchronized;
        ok = islogical(rs) && rs;
        checks(end+1) = add(sprintf('red_sync_%s', nm), 2, ok, ...
            sprintf('RedFrameSynchronized=%s', mat2str(rs))); %#ok<AGROW>

        switch E.paradigm
            case 'trial_onset'
                checks = onset_checks(checks, nm, E, t0, tN, B);
            case 'direction'
                checks = onset_checks(checks, nm, E, t0, tN, B);
                % expected subtrial count = n_dir * outer repeats
                s = B.Stimuli(E.stim_index);
                n_dir = numel(E.dir_values);
                exp_sub = n_dir * s.trials;
                found = size(E.subtrials, 3);
                ok = (found == exp_sub);
                checks(end+1) = add(sprintf('subtrial_count_%s', nm), 2, ok, ...
                    sprintf('found=%d expected=%d (n_dir=%d x trials=%d)', found, exp_sub, n_dir, s.trials)); %#ok<AGROW>
                % direction reconstruction verification (N/A when <2 repeat
                % cycles, e.g. grating trials=1 -- not a failure)
                ok = isfield(E,'verify') && E.verify.pass;
                if isfield(E,'verify') && ~E.verify.applicable
                    vdet = 'N/A (single repeat -- cannot cross-validate)';
                else
                    vdet = sprintf('r_same=%.3f r_shifted=%.3f', E.verify.r_same, E.verify.r_shifted);
                end
                checks(end+1) = add(sprintf('dir_reconstruction_%s', nm), 2, ok, vdet); %#ok<AGROW>
            case 'block'
                ok = size(E.block, 2) >= 2;
                checks(end+1) = add(sprintf('block_nonempty_%s', nm), 2, ok, ...
                    sprintf('%d frames, %.1fs', size(E.block,2), E.duration_sec)); %#ok<AGROW>
        end
    end

    % ---------- summarize ----------
    tier2 = checks(arrayfun(@(c) c.tier == 2, checks));
    passed = all(arrayfun(@(c) c.pass, tier2));
    qc = struct();
    qc.checks = checks;
    qc.passed = passed;
    qc.n_fail = sum(~arrayfun(@(c) c.pass, checks));

    S.meta.qc = qc;
    S.meta.qc_passed = passed;

    print_qc_table(S.meta.session_id, checks, passed);
end

% ======================================================================
function checks = onset_checks(checks, nm, E, t0, tN, B)
    onsets = E.onsets(:);
    pre = 0; post = 0;
    if isfield(E, 'peak_window'); end
    % coverage window is the analysis window used for this paradigm
    if isfield(E, 't_com') && ~isempty(E.t_com)
        pre = -E.t_com(1); post = E.t_com(end);
    end

    n_on = numel(onsets);
    in_range = all(onsets >= t0 - 1e-6 & onsets <= tN + 1e-6);
    covered  = all(onsets - pre >= t0 - 1e-6 & onsets + post <= tN + 1e-6);
    checks(end+1) = add_check(sprintf('onsets_in_range_%s', nm), 2, in_range, ...
        sprintf('%d onsets, all_in_rec=%d, window_covered=%d', n_on, in_range, covered)); %#ok<AGROW>

    % expected count for looming = Stimuli.trials
    if strcmp(E.type, 'looming')
        s = B.Stimuli(E.stim_index);
        ok = (n_on == s.trials);
        checks(end+1) = add_check(sprintf('onset_count_%s', nm), 2, ok, ...
            sprintf('found=%d expected(trials)=%d', n_on, s.trials)); %#ok<AGROW>
    end

    % flash timing structure: the flash count is RNG-determined (no fixed
    % metadata integer to assert), but the *structure* is falsifiable -- every
    % inter-flash gap must lie in [interval], and the last onset must sit
    % within the stimulus window. A wrong seed / wrong interval breaks this.
    if strcmp(E.type, 'sparse_local_global_flashes')
        s = B.Stimuli(E.stim_index);
        if isfield(s, 'specParams') && isfield(s.specParams, 'interval') ...
                && isfield(s, 'stimulus_trial_t') && isfield(s, 'ifi')
            rel  = onsets(:) - s.TimeStimulusFrame(1);      % onset times relative to stim start
            gaps = diff([0; rel]);                          % first gap measured from stim start
            iv   = s.specParams.interval;
            tol  = 2 * s.ifi(1);                            % 60 Hz display-frame quantization slack
            gaps_ok = all(gaps >= iv(1) - tol & gaps <= iv(2) + tol);
            last_ok = max(rel) <= s.stimulus_trial_t + tol;
            ok = gaps_ok && last_ok;
            checks(end+1) = add_check(sprintf('flash_timing_%s', nm), 2, ok, ...
                sprintf('gaps in [%g %g]=%d, last<=STT=%d (n=%d, gaps %.1f-%.1fs)', ...
                iv(1), iv(2), gaps_ok, last_ok, n_on, min(gaps), max(gaps))); %#ok<AGROW>
        end
    end

    % no all-NaN trial/subtrial slices
    A = [];
    if isfield(E, 'trials'); A = E.trials; elseif isfield(E, 'subtrials'); A = E.subtrials; end
    if ~isempty(A)
        empty_slices = squeeze(all(all(isnan(A), 1), 2));
        n_empty = sum(empty_slices);
        ok = (n_empty == 0);
        checks(end+1) = add_check(sprintf('no_empty_trials_%s', nm), 2, ok, ...
            sprintf('%d/%d all-NaN slices', n_empty, size(A,3))); %#ok<AGROW>
    end
end

% ======================================================================
function c = add_check(name, tier, pass, detail)
    c = struct('name', name, 'tier', tier, 'pass', logical(pass), 'detail', detail);
end

% ======================================================================
function print_qc_table(sid, checks, passed)
    fprintf('    +-- QC checks for %s --\n', sid);
    for i = 1:numel(checks)
        c = checks(i);
        tag = 'PASS'; if ~c.pass; tag = 'FAIL'; end
        fprintf('    | T%d %-4s %-28s %s\n', c.tier, tag, c.name, c.detail);
    end
    if passed
        fprintf('    +-- OVERALL: PASS\n');
    else
        fprintf('    +-- OVERALL: *** QC FAILED *** (inspect QC figure)\n');
    end
end
