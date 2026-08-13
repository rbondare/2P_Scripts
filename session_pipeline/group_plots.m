function out = group_plots(kind, G, varargin)
% GROUP_PLOTS  Figures for a compare_group G struct (several matched
% baseline/drug recordings averaged within drug condition).
%
% Pure visualization -- every value drawn is one compare_group already
% computed. Rows of each figure = drug condition, columns = stimulus.
%
%   group_plots('traces',    G)   grand-mean population trace across recordings,
%                                 baseline vs drug (+-s.e.m. ACROSS recordings,
%                                 thin lines = individual recordings)
%   group_plots('superplot', G)   pooled matched-cell violins with each
%                                 recording's median overlaid (the honest view:
%                                 pooled distribution + per-recording spread)
%   group_plots('mi',        G)   pooled modulation-index violin + per-recording
%                                 medians
%   group_plots('forest',    G)   one row per recording: median MI with the
%                                 condition pooled median as a reference line
%   T = group_plots('summary', G) table: per condition x stimulus, n, pooled
%                                 medians, pooled p, per-recording consistency
%
% Name-value options:
%   'metric'    'z' (default) | 'raw'
%   'stimuli'   cellstr subset (default: all stimuli present)
%   'save'      output PNG path ('' -> screen)

    opt.metric = 'z'; opt.stimuli = {}; opt.save = '';
    for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end
    metric = lower(char(opt.metric));
    assert(any(strcmp(metric, {'z','raw'})), 'metric must be ''z'' or ''raw''');

    conds = fieldnames(G.cond);
    stims = {};
    for c = 1:numel(conds)
        f = fieldnames(G.cond.(conds{c})); f = f(~strcmp(f, 'meta'));
        stims = union(stims, f, 'stable');
    end
    if ~isempty(opt.stimuli); stims = intersect(stims, opt.stimuli, 'stable'); end
    assert(~isempty(stims), 'group_plots: no stimuli to plot');

    out = [];
    switch lower(char(kind))
        case 'traces';    fig = plot_traces(G, conds, stims);
        case 'superplot'; fig = plot_superplot(G, conds, stims, metric);
        case 'mi';        fig = plot_mi(G, conds, stims, metric);
        case 'forest';    fig = plot_forest(G, conds, stims, metric);
        case 'summary';   out = summary_table(G, conds, stims, metric); return;
        otherwise; error('group_plots: unknown kind ''%s''', char(kind));
    end

    if ~isempty(char(opt.save)); save_png(fig, char(opt.save)); end
end

% ======================================================================
function save_png(fig, savep)
% Hide the interactive axes toolbar first -- exportgraphics bakes it into the
% image (over a title, if you are unlucky) when it happens to be showing.
    for ax = findall(fig, 'Type', 'axes')'
        try
            ax.Toolbar.Visible = 'off';
        catch
        end
    end
    try; exportgraphics(fig, savep, 'Resolution', 150); catch; print(fig, savep, '-dpng','-r150'); end
    fprintf('  -> %s\n', savep);
end

% ======================================================================
function fig = plot_traces(G, conds, stims)
    % Only trial-aligned stimuli have a trace: block stimuli (spontaneous,
    % checkers2) have session-specific lengths and no onset, so they are
    % dropped here rather than left as empty columns.
    keep = false(1, numel(stims));
    for i = 1:numel(stims)
        for r = 1:numel(conds)
            A = G.cond.(conds{r});
            if isfield(A, stims{i}) && A.(stims{i}).traces.ok; keep(i) = true; end
        end
    end
    stims = stims(keep);
    assert(~isempty(stims), 'group_plots(traces): no trial-aligned stimuli (block only)');

    [fig, tl] = grid_fig(conds, stims, 330, 250);
    title(tl, {'population response of matched cells, averaged ACROSS recordings', ...
        'thick = grand mean over recordings (+-s.e.m. across recordings), thin = one recording'}, ...
        'FontWeight', 'bold');
    cB = [0.30 0.45 0.75]; cD = [0.80 0.40 0.35];
    for r = 1:numel(conds)
        A = G.cond.(conds{r});
        for i = 1:numel(stims)
            ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
            if ~isfield(A, stims{i}) || ~A.(stims{i}).traces.ok
                text(ax, 0.5, 0.5, 'not recorded', 'Units','normalized', ...
                    'HorizontalAlignment','center','FontSize',8,'Color',[.5 .5 .5]);
                set(ax, 'XTick', [], 'YTick', []);
            else
                E = A.(stims{i}); T = E.traces;
                for j = 1:size(T.base, 1)
                    plot(ax, T.t_com, T.base(j,:), '-', 'Color', [cB 0.35], 'LineWidth', 0.7);
                    plot(ax, T.t_com, T.drug(j,:), '-', 'Color', [cD 0.35], 'LineWidth', 0.7);
                end
                shade(ax, T.t_com, T.grand_base, T.sem_base, cB);
                shade(ax, T.t_com, T.grand_drug, T.sem_drug, cD);
                hb = plot(ax, T.t_com, T.grand_base, '-', 'Color', cB, 'LineWidth', 2);
                hd = plot(ax, T.t_com, T.grand_drug, '-', 'Color', cD, 'LineWidth', 2);
                xline(ax, 0, 'k:');
                xlim(ax, [T.t_com(1) T.t_com(end)]);
                title(ax, sprintf('%s  (%d rec)', E.stim_label, E.n_pairs), 'Interpreter','none','FontSize',9);
                if r == 1 && i == numel(stims)
                    legend(ax, [hb hd], {'baseline','drug'}, 'Box','off','FontSize',8, 'Location','best');
                end
            end
            if i == 1
                ylabel(ax, sprintf('%s\nmean dF/F', A.meta.condition), ...
                    'Interpreter','none', 'FontWeight','bold');
            end
            if r == numel(conds); xlabel(ax, 'time from onset (s)'); end
        end
    end
end

% ======================================================================
function fig = plot_superplot(G, conds, stims, metric)
    [fig, tl] = grid_fig(conds, stims, 300, 260);
    title(tl, {sprintf('matched-cell responses pooled across recordings (%s)', metric), ...
        'violin = all cells pooled; open circles = each recording''s median'}, 'FontWeight', 'bold');
    cB = [0.30 0.45 0.75]; cD = [0.80 0.40 0.35];
    for r = 1:numel(conds)
        A = G.cond.(conds{r});
        for i = 1:numel(stims)
            ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
            if ~isfield(A, stims{i}); axis(ax, 'off'); continue; end
            E = A.(stims{i}); [b, d, ml, p] = getvals(E, metric);
            violin_kde(ax, 1, b, 0.35, cB);
            violin_kde(ax, 2, d, 0.35, cD);
            for j = 1:numel(E.per_pair)
                pb = pairmed(E.per_pair(j), metric, 'base');
                pd = pairmed(E.per_pair(j), metric, 'drug');
                plot(ax, [1 2], [pb pd], '-o', 'Color', [0.25 0.25 0.25], ...
                    'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1);
            end
            yline(ax, 0, ':', 'Color', [.5 .5 .5]);
            set(ax, 'XTick', [1 2], 'XTickLabel', {'baseline','drug'});
            xlim(ax, [0.4 2.6]); ylim(ax, robust_lim([b; d], [1 99]));
            title(ax, sprintf('%s\n%d cells / %d rec, pooled p=%.3g', E.stim_label, ...
                E.pooled.n_cells, E.n_pairs, p), 'Interpreter','none','FontSize',8);
            if i == 1; ylabel(ax, sprintf('%s\n%s', A.meta.condition, ml), 'Interpreter','none'); end
        end
    end
end

% ======================================================================
function fig = plot_mi(G, conds, stims, metric)
    [fig, tl] = grid_fig(conds, stims, 260, 260);
    title(tl, {sprintf('modulation index (drug vs baseline, %s)', metric), ...
        'violin = pooled cells; red = pooled median; open circles = per-recording medians'}, ...
        'FontWeight', 'bold');
    for r = 1:numel(conds)
        A = G.cond.(conds{r});
        for i = 1:numel(stims)
            ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
            if ~isfield(A, stims{i}); axis(ax, 'off'); continue; end
            E = A.(stims{i}); m = getmi(E, metric);
            violin_kde(ax, 1, m, 0.35, [0.5 0.5 0.55]);
            if ~isempty(m)
                xj = 1 + (rand(numel(m),1)-0.5)*0.25;
                scatter(ax, xj, m, 5, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.18);
                plot(ax, 1, median(m,'omitnan'), 'ro', 'MarkerFaceColor','r','MarkerSize',7);
            end
            pm = permeds(E, metric);
            plot(ax, 1 + linspace(-0.18, 0.18, max(numel(pm),1)), pm, 'o', ...
                'Color', [0 0 0], 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1);
            yline(ax, 0, 'k--');
            set(ax, 'XTick', 1, 'XTickLabel', sprintf('n=%d cells', numel(m)));
            xlim(ax, [0.4 1.6]); ylim(ax, [-1.05 1.05]);
            title(ax, sprintf('%s\nMI=%+.3f (rec med %+.3f)', E.stim_label, ...
                median(m,'omitnan'), median(pm,'omitnan')), 'Interpreter','none','FontSize',8);
            if i == 1; ylabel(ax, sprintf('%s\nmodulation index', A.meta.condition), 'Interpreter','none'); end
        end
    end
end

% ======================================================================
function fig = plot_forest(G, conds, stims, metric)
    n = numel(stims);
    fig = figure('Position', [60 60 max(520, 300*n) 420]);
    tl = tiledlayout(fig, 1, n, 'Padding','compact','TileSpacing','compact');
    title(tl, {sprintf('per-recording modulation index (%s)', metric), ...
        'dot = that recording''s median MI; line = condition pooled median'}, 'FontWeight','bold');
    cc = lines(numel(conds));
    % Fixed row per recording, identical in every panel -- otherwise a panel a
    % recording is missing from (e.g. moving_bar) would silently shift rows
    % relative to the y-labels drawn only on the first panel.
    rowLab = {}; rowCond = [];
    for r = 1:numel(conds)
        A = G.cond.(conds{r});
        for j = 1:numel(A.meta.pair_labels)
            rowLab{end+1} = sprintf('%s (%s)', A.meta.pair_labels{j}, A.meta.condition); %#ok<AGROW>
            rowCond(end+1) = r; %#ok<AGROW>
        end
    end
    nrow = numel(rowLab);

    for i = 1:n
        ax = nexttile(tl); hold(ax, 'on'); box(ax, 'off');
        for r = 1:numel(conds)
            A = G.cond.(conds{r});
            if ~isfield(A, stims{i}); continue; end
            E = A.(stims{i}); pm = permeds(E, metric);
            rows = find(rowCond == r);
            plot(ax, repmat(median(getmi(E, metric), 'omitnan'), 1, 2), ...
                [min(rows)-0.4 max(rows)+0.4], '-', 'Color', cc(r,:), 'LineWidth', 2);
            for j = 1:numel(E.per_pair)
                y = find(strcmp(rowLab, sprintf('%s (%s)', E.per_pair(j).label, A.meta.condition)), 1);
                if isempty(y); continue; end
                plot(ax, pm(j), y, 'o', 'MarkerFaceColor', cc(r,:), 'MarkerEdgeColor','k','MarkerSize',7);
            end
        end
        xline(ax, 0, 'k--');
        set(ax, 'YTick', 1:nrow, 'YTickLabel', rowLab, 'TickLabelInterpreter','none','FontSize',8);
        ylim(ax, [0.3, nrow+0.7]); xlim(ax, [-1.05 1.05]);
        title(ax, firststim(G, conds, stims{i}), 'Interpreter','none','FontSize',9);
        xlabel(ax, 'modulation index');
        if i > 1; set(ax, 'YTickLabel', {}); end
    end
end

% ======================================================================
function T = summary_table(G, conds, stims, metric)
    Condition = strings(0,1); Stimulus = strings(0,1);
    nRec = []; nCells = []; MedBase = []; MedDrug = []; MedMI = [];
    p_pooled = []; MedMI_perRec = []; Consistent = strings(0,1); p_acrossRec = [];
    for r = 1:numel(conds)
        A = G.cond.(conds{r});
        for i = 1:numel(stims)
            if ~isfield(A, stims{i}); continue; end
            E = A.(stims{i}); [b, d, ~, p] = getvals(E, metric); m = getmi(E, metric);
            Condition(end+1,1) = string(A.meta.condition); %#ok<AGROW>
            Stimulus(end+1,1)  = string(E.stim_label); %#ok<AGROW>
            nRec(end+1,1)   = E.n_pairs; %#ok<AGROW>
            nCells(end+1,1) = E.pooled.n_cells; %#ok<AGROW>
            MedBase(end+1,1) = median(b,'omitnan'); %#ok<AGROW>
            MedDrug(end+1,1) = median(d,'omitnan'); %#ok<AGROW>
            MedMI(end+1,1)   = median(m,'omitnan'); %#ok<AGROW>
            p_pooled(end+1,1) = p; %#ok<AGROW>
            MedMI_perRec(end+1,1) = median(permeds(E, metric),'omitnan'); %#ok<AGROW>
            Consistent(end+1,1) = string(yn(E.across_pair.consistent_sign)); %#ok<AGROW>
            p_acrossRec(end+1,1) = E.across_pair.p_signrank_pairs; %#ok<AGROW>
        end
    end
    T = table(Condition, Stimulus, nRec, nCells, MedBase, MedDrug, MedMI, ...
        p_pooled, MedMI_perRec, Consistent, p_acrossRec);
    T.Properties.Description = sprintf('compare_group, metric=%s', metric);
    disp(T);
end

% ======================================================================
function [fig, tl] = grid_fig(conds, stims, w, h)
    nr = numel(conds); nc = numel(stims);
    fig = figure('Position', [40 40 max(560, w*nc) max(360, h*nr + 90)]);
    tl = tiledlayout(fig, nr, nc, 'Padding','compact','TileSpacing','compact');
end

function shade(ax, t, y, e, c)
    ok = isfinite(y) & isfinite(e);
    if nnz(ok) < 2; return; end
    t = t(ok); y = y(ok); e = e(ok);
    fill(ax, [t fliplr(t)], [y+e fliplr(y-e)], c, 'FaceAlpha', 0.18, 'EdgeColor','none');
end

function [b, d, ml, p] = getvals(E, metric)
    if strcmp(metric, 'z')
        b = E.pooled.z_base; d = E.pooled.z_drug; ml = 'baseline-subtracted z';
        p = E.pooled.p_signrank_z;
    else
        b = E.pooled.raw_base; d = E.pooled.raw_drug; ml = 'peak dF/F';
        p = E.pooled.p_signrank_raw;
    end
    b = b(:); d = d(:);
end

function m = getmi(E, metric)
    if strcmp(metric, 'z'); m = E.pooled.mi_z; else; m = E.pooled.mi_raw; end
    m = m(:); m = m(isfinite(m));
end

function pm = permeds(E, metric)
    if strcmp(metric, 'z'); pm = [E.per_pair.med_mi_z]; else; pm = [E.per_pair.med_mi_raw]; end
    pm = pm(:);
end

function v = pairmed(pp, metric, which)
    if strcmp(metric, 'z')
        if strcmp(which,'base'); v = pp.med_z_base; else; v = pp.med_z_drug; end
    else
        if strcmp(which,'base'); v = pp.med_raw_base; else; v = pp.med_raw_drug; end
    end
end

function lab = firststim(G, conds, stim)
    lab = stim;
    for r = 1:numel(conds)
        A = G.cond.(conds{r});
        if isfield(A, stim); lab = A.(stim).stim_label; return; end
    end
end

function s = yn(x)
    if x; s = 'yes'; else; s = 'no'; end
end
