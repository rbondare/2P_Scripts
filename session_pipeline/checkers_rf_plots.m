function fh = checkers_rf_plots(kind, varargin)
% CHECKERS_RF_PLOTS  Views of checkers2 receptive fields. Pure visualisation of
% what checkers_rf_session / compare_checkers_rf already computed, so a figure
% can never disagree with the stats (same rule as compare_plots.m).
%
%   checkers_rf_plots('diagnostic', A, Bs)   real vs circular-shift null. LOOK AT
%                                            THIS FIRST -- it is the only figure
%                                            that is meaningful when the null gate
%                                            has failed.
%   checkers_rf_plots('maps', A, Bs, R)      best matched cells, baseline vs drug
%   checkers_rf_plots('props', R)            paired property scatters
%
%   A, Bs  outputs of checkers_rf_session (baseline, drug)
%   R      output of compare_checkers_rf
%
% Optional name/value: 'save', <path>, 'n_show', <k>.

    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'save', '');
    addParameter(p, 'n_show', 8);
    args = varargin;
    isnv = find(cellfun(@(x) ischar(x)||isstring(x), args), 1);
    if isempty(isnv); pos = args; nv = {}; else; pos = args(1:isnv-1); nv = args(isnv:end); end
    parse(p, nv{:});
    opt = p.Results;

    switch lower(kind)
        case 'diagnostic'; fh = plot_diagnostic(pos{1}, pos{2});
        case 'maps';       fh = plot_maps(pos{1}, pos{2}, pos{3}, opt.n_show);
        case 'props';      fh = plot_props(pos{1});
        otherwise; error('unknown kind ''%s''', kind);
    end
    if ~isempty(opt.save)
        exportgraphics(fh, opt.save, 'Resolution', 200);
        fprintf('saved %s\n', opt.save);
    end
end


% ======================================================================
function fh = plot_diagnostic(A, Bs)
% Is there any stimulus-locked signal at all? Each panel compares the real
% reverse correlation against the identical computation with the pattern order
% circularly shifted. If the two distributions sit on top of each other there are
% no receptive fields, and everything else in this file is a picture of noise.
    fh = figure('Position', [80 80 1150 700], 'Color', 'w');
    S = {A, Bs}; lbl = {char(A.meta.session_id), char(Bs.meta.session_id)};

    for j = 1:2
        X = S{j};
        % --- peak real vs peak null, per cell
        subplot(2, 3, (j-1)*3 + 1);
        loglog(X.props.peak_abs_null, X.props.peak_abs, '.', 'MarkerSize', 3, ...
            'Color', [0.2 0.6 1]); hold on;
        lm = [min(X.props.peak_abs_null) max(X.props.peak_abs)];
        plot(lm, lm, 'k-', 'LineWidth', 0.8);
        axis square; box off;
        xlabel('peak |RF|, shifted null'); ylabel('peak |RF|, real');
        title(sprintf('%s\nratio %.3f, p=%.2g', lbl{j}, X.gate.median_ratio, ...
            X.gate.p_signrank), 'Interpreter', 'none', 'FontSize', 8);

        % --- SNR distribution against the 0 dB chance line
        subplot(2, 3, (j-1)*3 + 2);
        histogram(X.props.snr_db, 60, 'FaceColor', [0.2 0.6 1], 'EdgeColor', 'none');
        hold on; yl = ylim; plot([0 0], yl, 'k--', 'LineWidth', 1);
        plot([3 3], yl, 'r--', 'LineWidth', 1); box off;
        xlabel('SNR vs null (dB)'); ylabel('cells');
        title(sprintf('%.1f%% >= 3 dB  (0 dB = chance)', 100*X.gate.frac_above_3dB), ...
            'FontSize', 8);

        % --- the best cell's map, real and null, on ONE colour scale
        [~, u] = max(X.props.snr_db);
        li = find(X.lags >= 0, 1) + max(0, round(X.props.peak_lag_s(u)/median(diff(X.lags))));
        li = min(max(li, 1), numel(X.lags));
        M = double(X.rf(:, :, li, u));
        subplot(2, 3, (j-1)*3 + 3);
        v = max(abs(M(:)));
        imagesc(X.az_deg, X.el_deg, M, [-v v]); axis image; colormap(gca, redblue_local());
        set(gca, 'YDir', 'normal'); colorbar; box off;
        xlabel('azimuth (deg)'); ylabel('elevation (deg)');
        title(sprintf('best cell #%d, %.1f dB, lag %.1f s', u, X.props.snr_db(u), ...
            X.lags(li)), 'FontSize', 8);
    end
    sgtitle(['checkers2 RF diagnostic — real vs circular-shift null' newline ...
        'points on the unity line = no receptive fields'], 'FontSize', 10);
end


% ======================================================================
function fh = plot_maps(A, Bs, R, n_show)
% Matched cells, baseline and drug spatial maps at each cell's own peak lag,
% ranked by the SUM of both conditions' SNR (never by one side).
    score = (R.base.snr_db + R.drug.snr_db) / 2;
    [~, ord] = sort(score, 'descend', 'MissingPlacement', 'last');
    ord = ord(1:min(n_show, numel(ord)));

    fh = figure('Position', [60 60 160*numel(ord)+120 430], 'Color', 'w');
    for j = 1:numel(ord)
        u = ord(j);
        for row = 1:2
            if row == 1; X = A; k = R.meta.idx1(u); lag = R.base.peak_lag_s(u);
            else;        X = Bs; k = R.meta.idx2(u); lag = R.drug.peak_lag_s(u); end
            [~, li] = min(abs(X.lags - lag));
            M = double(X.rf(:, :, li, k));
            subplot(2, numel(ord), (row-1)*numel(ord) + j);
            v = max(abs(M(:)));
            imagesc(X.az_deg, X.el_deg, M, [-v v]); axis image;
            colormap(gca, redblue_local()); set(gca, 'YDir', 'normal');
            set(gca, 'XTick', [], 'YTick', []); box off;
            if row == 1
                title(sprintf('#%d  %.1f dB', u, score(u)), 'FontSize', 7);
            end
            if j == 1
                if row == 1; ylabel('baseline'); else; ylabel('drug'); end
            end
        end
    end
    sgtitle(sprintf('%s vs %s — matched-cell checkers2 RFs at each cell''s peak lag', ...
        char(R.meta.session1), char(R.meta.session2)), 'FontSize', 9, 'Interpreter', 'none');
end


% ======================================================================
function fh = plot_props(R)
% Paired baseline-vs-drug scatters for the gated cells.
    keys = {'peak_amp_abs', 'RF amplitude (dF/F per contrast)'; ...
            'diameter_deg', 'RF diameter (deg)'; ...
            'peak_lag_s',   'peak lag (s)'; ...
            'ttp_s',        'time to peak (s)'; ...
            'az_deg',       'RF azimuth (deg)'; ...
            'el_deg',       'RF elevation (deg)'};
    fh = figure('Position', [60 60 1000 620], 'Color', 'w');
    for i = 1:size(keys, 1)
        k = keys{i, 1};
        if strcmp(k, 'peak_amp_abs')
            a = abs(R.base.peak_amp); b = abs(R.drug.peak_amp);
        else
            a = R.base.(k); b = R.drug.(k);
        end
        m = R.good & isfinite(a) & isfinite(b);
        subplot(2, 3, i);
        plot(a(m), b(m), 'o', 'MarkerSize', 3, 'Color', [0.2 0.6 1], ...
            'MarkerFaceColor', [0.2 0.6 1]); hold on;
        if any(m)
            lm = [min([a(m); b(m)]) max([a(m); b(m)])];
            if diff(lm) > 0; plot(lm, lm, 'k-', 'LineWidth', 0.7); xlim(lm); ylim(lm); end
        end
        axis square; box off; xlabel('baseline'); ylabel('drug');
        pv = NaN;
        if isfield(R.stats, k); pv = R.stats.(k).p; end
        title(sprintf('%s\nn=%d, p=%.3g', keys{i,2}, sum(m), pv), 'FontSize', 8);
    end
    sgtitle(sprintf('%s vs %s — checkers2 RF properties, matched cells >= %g dB', ...
        char(R.meta.session1), char(R.meta.session2), R.meta.snr_min), ...
        'FontSize', 9, 'Interpreter', 'none');
end


% ======================================================================
function c = redblue_local()
    n = 256; h = round(n/2);
    r = [linspace(0.23, 1, h)'; linspace(1, 0.71, n-h)'];
    g = [linspace(0.30, 1, h)'; linspace(1, 0.09, n-h)'];
    b = [linspace(0.75, 1, h)'; linspace(1, 0.15, n-h)'];
    c = [r g b];
end
