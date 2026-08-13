function out = compare_plots(kind, R, varargin)
% COMPARE_PLOTS  Exploratory comparison figures from a compare_sessions R struct.
% Pure visualization: every value plotted is one compare_sessions already
% computed (z_base/z_drug/raw_base/raw_drug/mi_z/mi_raw). No re-analysis here,
% so a plot can never disagree with the stats -- see PIPELINE.md.
%
%   compare_plots('violins',    R)      baseline-vs-drug matched-cell distributions
%   compare_plots('modindex',   R)      modulation-index distribution per stimulus
%   compare_plots('histograms', R)      overlaid baseline/drug histograms
%   compare_plots('scatter',    R)      paired base-vs-drug per stimulus, unity line
%   T = compare_plots('summary', R)     table of medians + signrank p across stimuli
%
% Name-value options:
%   'metric'  'z' (default) | 'raw'     which response metric to show
%   'stimuli' cellstr subset (default: all common stimuli in R)
%   'save'    output PNG path (default '' -> show on screen)
%
% One view per call, one figure (or table) out -- so a notebook cell is one line.

    opt.metric = 'z'; opt.stimuli = {}; opt.save = '';
    for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end
    metric = lower(char(opt.metric));
    assert(any(strcmp(metric, {'z','raw'})), 'metric must be ''z'' or ''raw''');
    savep = char(opt.save);

    stims = fieldnames(R); stims = stims(~strcmp(stims, 'meta'));
    if ~isempty(opt.stimuli); stims = intersect(stims, opt.stimuli, 'stable'); end
    assert(~isempty(stims), 'No stimuli in R to plot');

    ttl = sprintf('%s  vs  %s   (%d matched cells)', ...
        R.meta.session1, R.meta.session2, R.meta.n_matched);

    out = [];
    switch lower(char(kind))
        case 'violins';    fig = plot_violins(R, stims, metric, ttl);
        case 'modindex';   fig = plot_modindex(R, stims, metric, ttl);
        case 'histograms'; fig = plot_histograms(R, stims, metric, ttl);
        case 'scatter';    fig = plot_scatter(R, stims, metric, ttl);
        case 'summary';    out = summary_table(R, stims, metric); return;
        otherwise; error('compare_plots: unknown kind ''%s''', char(kind));
    end

    if ~isempty(savep)
        try; exportgraphics(fig, savep, 'Resolution', 150); catch; print(fig, savep, '-dpng','-r150'); end
        fprintf('  -> %s\n', savep);
    end
end

% ======================================================================
function fig = plot_violins(R, stims, metric, ttl)
    n = numel(stims);
    fig = figure('Position', [80 80 max(460, 300*n) 470]);
    tl = tiledlayout(fig, 1, n, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, {'baseline vs drug -- matched-cell responses', ttl}, 'Interpreter', 'none', 'FontWeight', 'bold');
    cB = [0.30 0.45 0.75]; cD = [0.80 0.40 0.35];
    for i = 1:n
        res = R.(stims{i}); [b, d, ml] = getvals(res, metric);
        ax = nexttile(tl); hold(ax, 'on');
        violin_kde(ax, 1, b, 0.35, cB);
        violin_kde(ax, 2, d, 0.35, cD);
        plot(ax, [1 2], [median(b,'omitnan') median(d,'omitnan')], 'k-o', ...
            'LineWidth', 1.2, 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        yline(ax, 0, ':', 'Color', [.5 .5 .5]);
        set(ax, 'XTick', [1 2], 'XTickLabel', {'baseline','drug'}); xlim(ax, [0.4 2.6]); box(ax, 'off');
        ylim(ax, robust_lim([b; d], [1 99]));   % clip extreme z outliers for readability
        title(ax, sprintf('%s\nsignrank p=%.3g', res.stim_label, res.p_signrank), 'Interpreter', 'none', 'FontSize', 9);
        if i == 1; ylabel(ax, ml); end
    end
end

% ======================================================================
function fig = plot_modindex(R, stims, metric, ttl)
    n = numel(stims);
    fig = figure('Position', [80 80 max(460, 280*n) 470]);
    tl = tiledlayout(fig, 1, n, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, {'modulation index  (drug vs baseline, bounded [-1,1])', ttl}, 'Interpreter', 'none', 'FontWeight', 'bold');
    for i = 1:n
        res = R.(stims{i}); m = getmi(res, metric);
        ax = nexttile(tl); hold(ax, 'on');
        violin_kde(ax, 1, m, 0.35, [0.5 0.5 0.55]);
        if ~isempty(m)
            xj = 1 + (rand(numel(m),1)-0.5)*0.25;
            scatter(ax, xj, m, 6, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.25);
            plot(ax, 1, median(m,'omitnan'), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        end
        yline(ax, 0, 'k--');
        set(ax, 'XTick', 1, 'XTickLabel', {sprintf('n=%d', numel(m))}); xlim(ax, [0.4 1.6]); ylim(ax, [-1.05 1.05]); box(ax, 'off');
        title(ax, sprintf('%s\nmed MI=%+.3f, p=%.3g', res.stim_label, median(m,'omitnan'), res.p_signrank), 'Interpreter', 'none', 'FontSize', 9);
        if i == 1; ylabel(ax, sprintf('modulation index (%s)', metric)); end
    end
end

% ======================================================================
function fig = plot_histograms(R, stims, metric, ttl)
    n = numel(stims);
    fig = figure('Position', [80 80 max(460, 300*n) 470]);
    tl = tiledlayout(fig, 1, n, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, {'response distributions -- baseline vs drug', ttl}, 'Interpreter', 'none', 'FontWeight', 'bold');
    cB = [0.30 0.45 0.75]; cD = [0.80 0.40 0.35];
    for i = 1:n
        res = R.(stims{i}); [b, d, ml] = getvals(res, metric);
        ax = nexttile(tl); hold(ax, 'on');
        lim = robust_lim([b; d], [1 99]);        % clip extreme outliers for readability
        edges = linspace(lim(1), lim(2), 30);
        histogram(ax, b, edges, 'FaceColor', cB, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
        histogram(ax, d, edges, 'FaceColor', cD, 'FaceAlpha', 0.45, 'EdgeColor', 'none');
        xline(ax, median(b,'omitnan'), '-', 'Color', cB*0.7, 'LineWidth', 1.4);
        xline(ax, median(d,'omitnan'), '-', 'Color', cD*0.7, 'LineWidth', 1.4);
        xlim(ax, lim); box(ax, 'off'); xlabel(ax, ml);
        title(ax, sprintf('%s (p=%.3g)', res.stim_label, res.p_signrank), 'Interpreter', 'none', 'FontSize', 9);
        if i == 1; ylabel(ax, '# cells'); legend(ax, {'baseline','drug'}, 'Box', 'off', 'FontSize', 8); end
    end
end

% ======================================================================
function fig = plot_scatter(R, stims, metric, ttl)
    n = numel(stims);
    fig = figure('Position', [80 80 max(460, 300*n) 470]);
    tl = tiledlayout(fig, 1, n, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, {'paired matched-cell responses -- baseline vs drug', ttl}, 'Interpreter', 'none', 'FontWeight', 'bold');
    for i = 1:n
        res = R.(stims{i}); [b, d, ml] = getvals(res, metric);
        ax = nexttile(tl); hold(ax, 'on');
        scatter(ax, b, d, 10, [0.25 0.25 0.3], 'filled', 'MarkerFaceAlpha', 0.35);
        lim = robust_lim([b; d], [1 99]);        % clip extreme outliers for readability
        plot(ax, lim, lim, 'r--', 'LineWidth', 1); xlim(ax, lim); ylim(ax, lim);
        axis(ax, 'square'); box(ax, 'off');
        xlabel(ax, sprintf('baseline %s', ml));
        if i == 1; ylabel(ax, sprintf('drug %s', ml)); end
        title(ax, sprintf('%s (p=%.3g)', res.stim_label, res.p_signrank), 'Interpreter', 'none', 'FontSize', 9);
    end
end

% ======================================================================
function T = summary_table(R, stims, metric)
    n = numel(stims);
    Stimulus = strings(n,1); nMatched = zeros(n,1);
    MedBase = zeros(n,1); MedDrug = zeros(n,1); MedMI = zeros(n,1); p_signrank = zeros(n,1);
    for i = 1:n
        res = R.(stims{i}); [b, d] = getvals(res, metric); m = getmi(res, metric);
        Stimulus(i) = string(res.stim_label); nMatched(i) = res.n_matched;
        MedBase(i) = median(b,'omitnan'); MedDrug(i) = median(d,'omitnan');
        MedMI(i) = median(m,'omitnan'); p_signrank(i) = res.p_signrank;
    end
    T = table(Stimulus, nMatched, MedBase, MedDrug, MedMI, p_signrank);
    T.Properties.VariableNames = {'Stimulus','nMatched', ...
        ['MedBase_' metric], ['MedDrug_' metric], ['MedMI_' metric], 'p_signrank'};
    T.Properties.Description = sprintf('%s vs %s (metric=%s)', R.meta.session1, R.meta.session2, metric);
    disp(T);
end

% ======================================================================
function [b, d, ml] = getvals(res, metric)
    if strcmp(metric, 'z'); b = res.z_base; d = res.z_drug; ml = 'baseline-subtracted z';
    else;                   b = res.raw_base; d = res.raw_drug; ml = 'peak dF/F'; end
    b = b(:); d = d(:);
end

function m = getmi(res, metric)
    if strcmp(metric, 'z'); m = res.mi_z; else; m = res.mi_raw; end
    m = m(:); m = m(isfinite(m));
end

% robust_lim.m and violin_kde.m are shared files in this folder (used by
% group_plots too) -- deliberately not local copies.

