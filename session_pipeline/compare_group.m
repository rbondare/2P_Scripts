function G = compare_group(pairs, opts)
% COMPARE_GROUP  Group layer over several matched baseline/drug recordings.
%
% One level above compare_sessions (see PIPELINE.md): it runs the validated
% pairwise comparison on every pair, verifies each ROI match, then aggregates
% within drug condition BOTH ways --
%
%   * POOLED     : every matched cell from every recording in one distribution
%                  (max power; cells are nested in sessions, so a single strong
%                  recording can carry it -- always read next to the per-pair
%                  numbers).
%   * PER-PAIR   : each recording reduced to one median, then summarised across
%                  recordings (pseudoreplication-safe, but n = #recordings).
%
% Re-derives NOTHING: metrics and traces come from compare_sessions, which
% itself only reads the cached *_organized.mat structs.
%
%   G = compare_group(pairs)
%   G = compare_group(pairs, opts)
%
% pairs : struct array, one element per baseline/drug recording pair:
%   .baseline_file    aggregated *_preprocessed.mat of the BASELINE session
%                     (= column 1 of the ROI-match mapping), or the organized
%                     *_organized.mat directly
%   .drug_file        same for the DRUG session (= column 2)
%   .roi_match_file   roiMatch_plane<P>_<date>_<t1>_<t2>_V1.mat
%   .condition        e.g. 'guanfacine' | 'yohimbine'
%   .plane            OPTIONAL 1-based plane; default inferred from the
%                     roiMatch filename (plane0 -> 1, plane1 -> 2, ...)
%   .plane_base       OPTIONAL 1-based plane for the baseline session, and
%   .plane_drug       for the drug session -- use these when the FOV shifted
%                     in z between sessions, so the match links DIFFERENT
%                     planes (0310: baseline plane2 <-> drug plane0). The
%                     roiMatch filename records only one number and cannot
%                     express that; the .match check below is what catches it.
%   .label            OPTIONAL short name; default '<date>_<t1>_<t2>'
%
% opts:
%   .organized_dir    default session_pipeline/organized_data
%   .verbose          default true
%   .match_check      default true -- matched-cell centroid distance vs a
%                     shuffled control. This is the falsifiable test that the
%                     ROI-match file really belongs to this session/plane pair
%                     (the roiMatch files store NO session provenance).
%
% Returns G with:
%   G.meta                     run info, conditions, pair count
%   G.pairs(k)                 label/condition/sessions/n_matched/match/R
%   G.cond.<condition>.<stim>  pooled + per_pair + across_pair + traces
%   G.contrast.<stim>          between-condition ranksum on pooled MI

    if nargin < 2 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'verbose');     opts.verbose = true;     end
    if ~isfield(opts, 'match_check'); opts.match_check = true; end
    here = fileparts(mfilename('fullpath'));
    if ~isfield(opts, 'organized_dir')
        opts.organized_dir = fullfile(here, 'organized_data');
    end
    addpath(here); addpath(fullfile(fileparts(here), 'analysis_functions'));

    n_pairs = numel(pairs);
    assert(n_pairs > 0, 'compare_group: no pairs given');

    G = struct();
    G.meta = struct('created', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
        'n_pairs', n_pairs, 'organized_dir', opts.organized_dir); %#ok<TNOW1,DATST>
    G.pairs = struct('label', {}, 'condition', {}, 'plane_base', {}, 'plane_drug', {}, ...
        'sess_base', {}, 'sess_drug', {}, 'org_base', {}, 'org_drug', {}, ...
        'roi_match_file', {}, 'n_matched', {}, 'match', {}, 'R', {});

    % ================= per-pair pairwise comparison =================
    for k = 1:n_pairs
        P = pairs(k);
        [pb_idx, pd_idx] = resolve_plane(P);
        [ob, sb] = resolve_org(P.baseline_file, pb_idx, opts.organized_dir);
        [od, sd] = resolve_org(P.drug_file,     pd_idx, opts.organized_dir);
        label = pair_label(P, sb, sd);

        if opts.verbose
            fprintf('\n[%d/%d] %-22s  %s (plane%d) -> %s (plane%d)\n', k, n_pairs, ...
                label, sb, pb_idx - 1, sd, pd_idx - 1);
        end

        S1 = load(ob); S2 = load(od);
        R  = compare_sessions(S1, S2, char(P.roi_match_file), struct('verbose', opts.verbose));

        m = struct('checked', false, 'pass', NaN, 'med_dist', NaN, ...
                   'med_dist_shuffled', NaN, 'corr_xy', NaN);
        if opts.match_check
            m = check_roi_match(S1, S2, char(P.roi_match_file));
            if opts.verbose
                fprintf('  ROI-match: median matched centroid dist %.1f px vs %.1f px shuffled, corr(x,y)=%.3f -> %s\n', ...
                    m.med_dist, m.med_dist_shuffled, m.corr_xy, tf(m.pass));
            end
            if ~isnan(m.pass) && ~m.pass
                warning(['compare_group: ROI match for %s looks wrong -- matched centroids are no ' ...
                    'closer than chance (%.0f px vs %.0f px shuffled). The mapping is not linking ' ...
                    'baseline plane%d to drug plane%d. Try other plane combinations ' ...
                    '(set pairs(k).plane_base / .plane_drug); the FOV can shift in z between sessions.'], ...
                    label, m.med_dist, m.med_dist_shuffled, pb_idx - 1, pd_idx - 1);
            end
        end

        G.pairs(end+1) = struct('label', label, 'condition', char(P.condition), ...
            'plane_base', pb_idx, 'plane_drug', pd_idx, 'sess_base', sb, 'sess_drug', sd, ...
            'org_base', ob, 'org_drug', od, 'roi_match_file', char(P.roi_match_file), ...
            'n_matched', R.meta.n_matched, 'match', m, 'R', R); %#ok<AGROW>
        clear S1 S2;
    end

    % ================= aggregate within condition =================
    conds = unique({G.pairs.condition}, 'stable');
    G.meta.conditions = conds;
    G.cond = struct();
    for c = 1:numel(conds)
        sel = find(strcmp({G.pairs.condition}, conds{c}));
        G.cond.(matlab.lang.makeValidName(conds{c})) = aggregate(G.pairs(sel), conds{c}, opts);
    end

    % ================= between-condition contrast =================
    G.contrast = struct();
    if numel(conds) == 2
        A = G.cond.(matlab.lang.makeValidName(conds{1}));
        Bc = G.cond.(matlab.lang.makeValidName(conds{2}));
        stims = intersect(fieldnames(A), fieldnames(Bc));
        stims = stims(~strcmp(stims, 'meta'));
        for i = 1:numel(stims)
            a = A.(stims{i}).pooled.mi_z; b = Bc.(stims{i}).pooled.mi_z;
            a = a(isfinite(a)); b = b(isfinite(b));
            p = NaN;
            if numel(a) >= 2 && numel(b) >= 2; p = ranksum(a, b); end
            G.contrast.(stims{i}) = struct( ...
                'stim_label', A.(stims{i}).stim_label, ...
                'cond1', conds{1}, 'cond2', conds{2}, ...
                'med_mi_1', median(a), 'med_mi_2', median(b), ...
                'n1', numel(a), 'n2', numel(b), 'p_ranksum', p);
        end
    end

    if opts.verbose; print_report(G); end
end

% ======================================================================
function A = aggregate(prs, cond_name, opts)
% Pool + per-pair summarise every stimulus present in >=1 pair of this
% condition.
    A = struct();
    A.meta = struct('condition', cond_name, 'n_pairs', numel(prs), ...
        'pair_labels', {{prs.label}}, 'n_matched_per_pair', [prs.n_matched]);

    stims = {};
    for k = 1:numel(prs)
        f = fieldnames(prs(k).R); f = f(~strcmp(f, 'meta'));
        stims = union(stims, f, 'stable');
    end

    for i = 1:numel(stims)
        nm = stims{i};
        have = find(arrayfun(@(p) isfield(p.R, nm), prs));
        if isempty(have); continue; end

        E = struct();
        E.stim_label   = prs(have(1)).R.(nm).stim_label;
        E.paradigm     = getf(prs(have(1)).R.(nm), 'paradigm', '');
        E.n_pairs      = numel(have);
        E.pair_labels  = {prs(have).label};

        % ---------- pooled over cells ----------
        pool = struct('z_base', [], 'z_drug', [], 'raw_base', [], 'raw_drug', [], ...
                      'mi_z', [], 'mi_raw', [], 'pair_id', []);
        pp = struct('label', {}, 'n', {}, 'med_z_base', {}, 'med_z_drug', {}, ...
                    'med_raw_base', {}, 'med_raw_drug', {}, 'med_mi_z', {}, ...
                    'med_mi_raw', {}, 'p_signrank', {});
        for j = 1:numel(have)
            r = prs(have(j)).R.(nm);
            pool.z_base   = [pool.z_base;   r.z_base(:)];
            pool.z_drug   = [pool.z_drug;   r.z_drug(:)];
            pool.raw_base = [pool.raw_base; r.raw_base(:)];
            pool.raw_drug = [pool.raw_drug; r.raw_drug(:)];
            pool.mi_z     = [pool.mi_z;     r.mi_z(:)];
            pool.mi_raw   = [pool.mi_raw;   r.mi_raw(:)];
            pool.pair_id  = [pool.pair_id;  repmat(j, numel(r.z_base), 1)];
            pp(end+1) = struct('label', prs(have(j)).label, 'n', numel(r.z_base), ...
                'med_z_base', median(r.z_base, 'omitnan'), 'med_z_drug', median(r.z_drug, 'omitnan'), ...
                'med_raw_base', median(r.raw_base, 'omitnan'), 'med_raw_drug', median(r.raw_drug, 'omitnan'), ...
                'med_mi_z', median(r.mi_z, 'omitnan'), 'med_mi_raw', median(r.mi_raw, 'omitnan'), ...
                'p_signrank', r.p_signrank); %#ok<AGROW>
        end
        pool.p_signrank_z   = paired_p(pool.z_base, pool.z_drug);
        pool.p_signrank_raw = paired_p(pool.raw_base, pool.raw_drug);
        pool.n_cells        = numel(pool.z_base);
        E.pooled   = pool;
        E.per_pair = pp;

        % ---------- across-pair (one value per recording) ----------
        mz = [pp.med_mi_z]; mr = [pp.med_mi_raw];
        E.across_pair = struct( ...
            'med_mi_z',        median(mz, 'omitnan'), ...
            'med_mi_raw',      median(mr, 'omitnan'), ...
            'mi_z_per_pair',   mz, ...
            'consistent_sign', all(sign(mz(~isnan(mz))) == sign(mz(find(~isnan(mz), 1)))) && any(~isnan(mz)), ...
            'p_signrank_pairs', pairs_p(mz), ...
            'n_pairs',         numel(mz), ...
            'underpowered',    numel(mz) < 6);   % signrank cannot reach p<0.05 below n=6

        % ---------- population traces averaged across recordings ----------
        E.traces = group_traces(prs(have), nm);

        A.(nm) = E;
    end
end

% ======================================================================
function T = group_traces(prs, nm)
% Mean-over-matched-cells trace per recording, then the grand mean across
% recordings. Only for paradigms compare_sessions attaches traces to
% (trial_onset / direction); block responses have session-specific lengths.
    T = struct('ok', false, 't_com', [], 'base', [], 'drug', [], ...
        'grand_base', [], 'grand_drug', [], 'sem_base', [], 'sem_drug', [], ...
        'pair_labels', {{}});
    r0 = prs(1).R.(nm);
    if ~isfield(r0, 't_com') || isempty(r0.t_com); return; end
    t = r0.t_com(:)';
    B = nan(numel(prs), numel(t)); D = B;
    for j = 1:numel(prs)
        r = prs(j).R.(nm);
        if ~isfield(r, 't_com') || isempty(r.t_com); continue; end
        b = mean(r.resp_base, 1, 'omitnan');
        d = mean(r.resp_drug, 1, 'omitnan');
        if numel(r.t_com) ~= numel(t) || any(abs(r.t_com(:)' - t) > 1e-9)
            b = interp1(r.t_com(:)', b, t, 'linear', NaN);
            d = interp1(r.t_com(:)', d, t, 'linear', NaN);
        end
        B(j, :) = b; D(j, :) = d;
    end
    n = size(B, 1);
    T.ok = true; T.t_com = t; T.base = B; T.drug = D;
    T.pair_labels = {prs.label};
    T.grand_base = mean(B, 1, 'omitnan');  T.grand_drug = mean(D, 1, 'omitnan');
    T.sem_base = std(B, 0, 1, 'omitnan') ./ sqrt(n);
    T.sem_drug = std(D, 0, 1, 'omitnan') ./ sqrt(n);
end

% ======================================================================
function m = check_roi_match(S1, S2, roi_match_file)
% Falsifiable test that this ROI-match file belongs to these two sessions at
% this plane: matched cells must sit at nearly the same FOV location (up to a
% rigid drift), i.e. far closer than randomly re-paired cells.
    m = struct('checked', false, 'pass', NaN, 'med_dist', NaN, ...
               'med_dist_shuffled', NaN, 'corr_xy', NaN, 'n', NaN);
    if ~isfield(S1.meta, 'roi_centroid') || ~isfield(S2.meta, 'roi_centroid')
        warning(['compare_group: organized struct has no meta.roi_centroid ' ...
                 '(cached before that field existed) -- ROI-match check skipped. ' ...
                 'Re-run run_organize to enable it.']);
        return;
    end
    M = load(roi_match_file);
    A = M.roiMatchData.allSessionMapping;
    i1 = A(:,1); i2 = A(:,2);
    ok = i1 >= 1 & i1 <= S1.meta.n_rois & i2 >= 1 & i2 <= S2.meta.n_rois;
    i1 = i1(ok); i2 = i2(ok);
    if numel(i1) < 10; return; end

    c1 = S1.meta.roi_centroid(i1, 1:2);
    c2 = S2.meta.roi_centroid(i2, 1:2);
    d  = sqrt(sum((c1 - c2).^2, 2));

    rng(0);                                   % reproducible shuffled control
    perm = randperm(numel(i2));
    ds = sqrt(sum((c1 - c2(perm, :)).^2, 2));

    cxy = mean([corr(c1(:,1), c2(:,1)), corr(c1(:,2), c2(:,2))]);
    m.checked = true; m.n = numel(i1);
    m.med_dist = median(d); m.med_dist_shuffled = median(ds); m.corr_xy = cxy;
    m.pass = (m.med_dist < 0.5 * m.med_dist_shuffled) && (cxy > 0.9);
end

% ======================================================================
function [pb, pd] = resolve_plane(P)
% Per-session plane. Explicit .plane_base/.plane_drug win (needed when z-drift
% made the best match a CROSS-plane one, e.g. 0310: baseline plane2 <-> drug
% plane0); then .plane for both; else the roiMatch filename.
    dflt = [];
    if isfield(P, 'plane') && ~isempty(P.plane)
        dflt = P.plane;
    else
        tok = regexp(char(P.roi_match_file), 'roiMatch_plane(\d+)_', 'tokens', 'once');
        if ~isempty(tok); dflt = str2double(tok{1}) + 1; end   % suite2p 0-based -> 1-based
    end
    pb = getf(P, 'plane_base', dflt);
    pd = getf(P, 'plane_drug', dflt);
    assert(~isempty(pb) && ~isempty(pd), ['compare_group: cannot determine planes for "%s" -- ' ...
        'set pairs(k).plane (or .plane_base/.plane_drug), or name the file roiMatch_plane<P>_...'], ...
        char(P.roi_match_file));
end

% ======================================================================
function [p, sid] = resolve_org(f, plane, org_dir)
% Accept either the aggregated *_preprocessed.mat (resolved to its cached
% organized struct) or an organized file directly.
    f = char(f);
    if endsWith(f, '_organized.mat')
        p = f; [~, sid] = fileparts(f);
        sid = regexprep(sid, '_plane\d+_organized$', '');
    else
        [p, sid] = organized_path(f, plane, org_dir);
    end
    assert(isfile(p), ['compare_group: no organized cache for %s at plane%d.\n' ...
        '  expected: %s\n  -> add it to run_organize.m''s session list and run it first.'], ...
        sid, plane - 1, p);
end

% ======================================================================
function label = pair_label(P, sb, sd)
    if isfield(P, 'label') && ~isempty(P.label); label = char(P.label); return; end
    tb = regexp(sb, '_(\d{6})_(\d{4})$', 'tokens', 'once');
    td = regexp(sd, '_(\d{6})_(\d{4})$', 'tokens', 'once');
    if ~isempty(tb) && ~isempty(td)
        label = sprintf('%s_%s_%s', tb{1}(3:end), tb{2}, td{2});
    else
        label = sprintf('%s_vs_%s', sb, sd);
    end
end

% ======================================================================
function p = paired_p(a, b)
    v = isfinite(a) & isfinite(b);
    if nnz(v) >= 2; p = signrank(a(v), b(v)); else; p = NaN; end
end

function p = pairs_p(m)
% Signrank of per-recording median MI against 0. With few recordings this
% cannot reach significance (n<6 -> min p = 0.0625..0.5); reported anyway,
% flagged by .underpowered.
    m = m(isfinite(m));
    if numel(m) >= 2; p = signrank(m); else; p = NaN; end
end

function v = getf(s, f, dflt)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = dflt; end
end

function s = tf(x)
    if isnan(x); s = 'n/a'; elseif x; s = 'PASS'; else; s = 'FAIL'; end
end

% ======================================================================
function print_report(G)
    fprintf('\n================ GROUP SUMMARY ================\n');
    for k = 1:numel(G.pairs)
        p = G.pairs(k);
        pl = sprintf('plane%d', p.plane_base - 1);
        if p.plane_drug ~= p.plane_base
            pl = sprintf('plane%d->%d', p.plane_base - 1, p.plane_drug - 1);
        end
        fprintf('  %-20s %-12s %-12s n=%4d matched  ROI-match %s (%.1f px vs %.0f shuffled)\n', ...
            p.label, p.condition, pl, p.n_matched, tf(p.match.pass), ...
            p.match.med_dist, p.match.med_dist_shuffled);
    end
    conds = fieldnames(G.cond);
    for c = 1:numel(conds)
        A = G.cond.(conds{c});
        fprintf('\n--- %s  (%d recordings, %s cells) ---\n', A.meta.condition, ...
            A.meta.n_pairs, mat2str(A.meta.n_matched_per_pair));
        stims = fieldnames(A); stims = stims(~strcmp(stims, 'meta'));
        fprintf('  %-28s %6s %8s %10s %10s %8s\n', 'stimulus', 'pairs', 'cells', ...
            'MI(pooled)', 'p(pooled)', 'consist');
        for i = 1:numel(stims)
            E = A.(stims{i});
            fprintf('  %-28s %6d %8d %+10.3f %10.3g %8s\n', E.stim_label, E.n_pairs, ...
                E.pooled.n_cells, median(E.pooled.mi_z, 'omitnan'), ...
                E.pooled.p_signrank_z, yn(E.across_pair.consistent_sign));
        end
    end
    if ~isempty(fieldnames(G.contrast))
        fprintf('\n--- between conditions (pooled MI, ranksum) ---\n');
        fs = fieldnames(G.contrast);
        for i = 1:numel(fs)
            C = G.contrast.(fs{i});
            fprintf('  %-28s %s %+0.3f  vs  %s %+0.3f   p=%.3g\n', C.stim_label, ...
                C.cond1, C.med_mi_1, C.cond2, C.med_mi_2, C.p_ranksum);
        end
    end
    fprintf('===============================================\n');
end

function s = yn(x)
    if x; s = 'yes'; else; s = 'no'; end
end
