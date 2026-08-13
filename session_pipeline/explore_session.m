function explore_session(src, out_dir)
% EXPLORE_SESSION  Single-session exploration from a cached *_organized.mat
% (no drug pair, no ROI match needed). One figure per stimulus showing the
% chunked/averaged response, tuning (direction stimuli), and a sorted
% heatmap. See PIPELINE.md.
%
%   explore_session('organized_data/AnimalRB19_260320_1444_organized.mat')
%   explore_session(S)                      % pass an in-memory struct
%   explore_session(src, out_dir)           % save PNGs instead of showing
%
% If out_dir is given, figures are saved there (invisible); otherwise shown.

    if ischar(src) || isstring(src)
        S = load(char(src));
    else
        S = src;
    end
    assert(isfield(S, 'stimuli'), 'Not an organized-session struct');

    save_mode = nargin >= 2 && ~isempty(out_dir);
    if save_mode && ~exist(out_dir, 'dir'); mkdir(out_dir); end
    vis = 'on'; if save_mode; vis = 'off'; end

    if ~S.meta.qc_passed
        fprintf('[explore_session] NOTE: %s did not pass QC -- interpret with care.\n', S.meta.session_id);
    end

    fns = fieldnames(S.stimuli);
    for i = 1:numel(fns)
        E = S.stimuli.(fns{i});
        fig = figure('Visible', vis, 'Position', [80 80 1100 420], ...
            'Name', sprintf('%s | %s', S.meta.session_id, E.type));
        switch E.paradigm
            case 'trial_onset'
                plot_trial_onset(fig, E, S);
            case 'direction'
                plot_direction(fig, E, S);
            case 'block'
                plot_block(fig, E, S);
        end
        if save_mode
            out_png = fullfile(out_dir, sprintf('%s_%s.png', S.meta.session_id, fns{i}));
            try; exportgraphics(fig, out_png, 'Resolution', 150); catch; print(fig, out_png, '-dpng', '-r150'); end
            close(fig);
            fprintf('  -> %s\n', out_png);
        end
    end
end

% ----------------------------------------------------------------------
function plot_trial_onset(fig, E, S)
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('%s  |  %s  |  %d trials', S.meta.session_id, E.type, E.n_trials), 'Interpreter', 'none');

    pop = mean(E.avg, 1, 'omitnan');
    sem = std(E.avg, 0, 1, 'omitnan') ./ sqrt(size(E.avg,1));
    ax = nexttile(tl); hold(ax, 'on');
    shade(ax, E.t_com, pop, sem, [0.2 0.2 0.2]);
    plot(ax, E.t_com, pop, 'k-', 'LineWidth', 2);
    xline(ax, 0, 'r--'); box(ax, 'off'); xlim(ax, [E.t_com(1) E.t_com(end)]);
    xlabel(ax, 'time from onset (s)'); ylabel(ax, 'pop. mean dF/F');
    title(ax, 'population mean \pm SEM (median-avg trials)');

    ax = nexttile(tl);
    heatmap_sorted(ax, E.t_com, E.avg, E.peak_window);
    xline(ax, 0, 'r--'); xlabel(ax, 'time from onset (s)'); ylabel(ax, 'ROI (sorted by peak)');
    title(ax, 'per-cell response');
end

% ----------------------------------------------------------------------
function plot_direction(fig, E, S)
    tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('%s  |  %s  |  %d directions', S.meta.session_id, E.type, numel(E.dir_values)), 'Interpreter', 'none');

    % tuning curve (population mean amplitude per direction)
    ax = nexttile(tl);
    tuning = mean(E.amp_by_dir, 1, 'omitnan');
    plot(ax, E.dir_values, tuning, 'ko-', 'LineWidth', 1.6, 'MarkerFaceColor', 'k');
    box(ax, 'off'); xlabel(ax, sprintf('%s (deg)', E.param_name)); ylabel(ax, 'peak dF/F');
    title(ax, 'tuning (pop. mean)');

    % per-direction population traces
    ax = nexttile(tl); hold(ax, 'on');
    cmap = parula(numel(E.dir_values));
    for d = 1:numel(E.dir_values)
        tr = mean(E.by_dir(:, :, d), 1, 'omitnan');
        plot(ax, E.t_com, tr, '-', 'Color', cmap(d,:), 'LineWidth', 1.2);
    end
    xline(ax, 0, 'r--'); box(ax, 'off'); xlim(ax, [E.t_com(1) E.t_com(end)]);
    xlabel(ax, 'time from onset (s)'); ylabel(ax, 'pop. mean dF/F');
    title(ax, 'per-direction traces');

    % preferred-direction heatmap
    ax = nexttile(tl);
    heatmap_sorted(ax, E.t_com, E.pref_resp, E.peak_window);
    xline(ax, 0, 'r--'); xlabel(ax, 'time from onset (s)'); ylabel(ax, 'ROI (sorted)');
    title(ax, 'preferred-direction response');
end

% ----------------------------------------------------------------------
function plot_block(fig, E, S)
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('%s  |  %s  |  %.0fs block', S.meta.session_id, E.type, E.duration_sec), 'Interpreter', 'none');

    ax = nexttile(tl);
    plot(ax, E.block_t, mean(E.block, 1, 'omitnan'), 'k-', 'LineWidth', 1.4);
    box(ax, 'off'); xlim(ax, [E.block_t(1) E.block_t(end)]);
    xlabel(ax, 'time (s)'); ylabel(ax, 'pop. mean dF/F'); title(ax, 'population mean');

    ax = nexttile(tl);
    heatmap_sorted(ax, E.block_t, E.block, [E.block_t(1) E.block_t(end)]);  % sort by whole-block mean
    xlabel(ax, 'time (s)'); ylabel(ax, 'ROI (sorted by mean)'); title(ax, 'per-cell activity');
end

% ----------------------------------------------------------------------
function shade(ax, x, y, e, c)
    xf = [x, fliplr(x)]; yf = [y + e, fliplr(y - e)];
    h = fill(ax, xf, yf, c, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    set(h, 'HandleVisibility', 'off');
end
