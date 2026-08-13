function fig = group_travelingwave(G, varargin)
% GROUP_TRAVELINGWAVE  Baseline vs drug traveling wave, per direction, with
% matched cells pooled across the recordings of one drug condition.
%
% The single-session version (walkthrough) showed one direction's repeats side
% by side to prove the latency sequence is reproducible. This is the paired
% version: for each direction, the SAME matched cells are drawn twice --
% baseline | drug -- in ONE cell order and ONE colour scale, so a change in
% the wave (amplitude, or the sequence breaking up) is read off directly.
%
%   group_travelingwave(G)
%   group_travelingwave(G, 'condition', 'yohimbine', 'save', 'wave.png')
%
% Options:
%   'condition'   which drug condition (default: first in G.cond)
%   'stimulus'    default 'moving_bar' (any 'direction'-paradigm stimulus)
%   'n_cells'     rows to show, most responsive first (default 300)
%   'sort_repeat' which baseline repeat defines the order (default 1)
%   'save'        PNG path ('' -> screen)
%
% NOT circular: cell selection AND latency sorting come from baseline repeat
% `sort_repeat` only; the baseline panel then shows the mean of the REMAINING
% baseline repeats. So a diagonal in the baseline panel is independent
% evidence of a real sequence, not an artefact of sorting the same data.
% Baseline and drug share one colour scale, so suppression shows as fading.

    opt.condition = ''; opt.stimulus = 'moving_bar'; opt.n_cells = 300;
    opt.sort_repeat = 1; opt.save = '';
    for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end

    conds = fieldnames(G.cond);
    if isempty(opt.condition); opt.condition = G.cond.(conds{1}).meta.condition; end
    cf = matlab.lang.makeValidName(opt.condition);
    assert(isfield(G.cond, cf), 'group_travelingwave: no condition ''%s'' in G', opt.condition);
    stim = matlab.lang.makeValidName(opt.stimulus);

    sel = find(strcmp({G.pairs.condition}, opt.condition));
    sel = sel(arrayfun(@(k) isfield(G.pairs(k).R, stim), sel));
    assert(~isempty(sel), 'group_travelingwave: no %s recordings have stimulus ''%s''', ...
        opt.condition, opt.stimulus);

    % ---- gather per direction, pooling cells over the condition's recordings ----
    t = []; dirs = []; Bsort = {}; Bshow = {}; Dshow = {}; labels = {};
    for k = sel
        P = G.pairs(k);
        S1 = load(P.org_base, 'stimuli'); S2 = load(P.org_drug, 'stimuli');
        E1 = S1.stimuli.(stim); E2 = S2.stimuli.(stim);
        assert(strcmp(E1.paradigm, 'direction'), ...
            'group_travelingwave: ''%s'' is paradigm ''%s'', needs ''direction''', opt.stimulus, E1.paradigm);
        i1 = P.R.meta.idx1; i2 = P.R.meta.idx2;

        if isempty(t)
            t = E1.t_com(:)';
            dirs = intersect(E1.dir_values, E2.dir_values);
            Bsort = cell(1, numel(dirs)); Bshow = Bsort; Dshow = Bsort;
        else
            assert(numel(E1.t_com) == numel(t), 'time grids differ across recordings');
            dirs = intersect(dirs, intersect(E1.dir_values, E2.dir_values));
        end

        for di = 1:numel(dirs)
            r1 = find(E1.direction_of == dirs(di));      % repeat subtrials, baseline
            r2 = find(E2.direction_of == dirs(di));      % ... drug
            if isempty(r1) || isempty(r2); continue; end
            sr = min(opt.sort_repeat, numel(r1));
            held = setdiff(1:numel(r1), sr);
            if isempty(held); held = sr; end             % single repeat -> no hold-out

            Bsort{di} = [Bsort{di}; E1.subtrials(i1, :, r1(sr))];
            Bshow{di} = [Bshow{di}; mean(E1.subtrials(i1, :, r1(held)), 3, 'omitnan')];
            Dshow{di} = [Dshow{di}; mean(E2.subtrials(i2, :, r2), 3, 'omitnan')];
        end
        labels{end+1} = sprintf('%s (n=%d)', P.label, numel(i1)); %#ok<AGROW>
        clear S1 S2;
    end

    % ---- one shared colour scale over everything drawn ----
    allv = cell2mat(cellfun(@(a) a(:), [Bshow(:); Dshow(:)], 'UniformOutput', false));
    lo = min(0, prctile(allv, 1)); hi = prctile(allv, 99);
    if ~(hi > lo); lo = min(allv); hi = max(allv); end

    pw = [0.5 4];                                        % moving_bar peak window
    nd = numel(dirs); nr = ceil(nd/4); nc = min(nd, 4);
    fig = figure('Position', [30 30 max(900, 330*nc) max(420, 330*nr)]);
    tl = tiledlayout(fig, nr, nc, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, {sprintf('%s  --  %s traveling wave, baseline | drug (same cells, same scale)', ...
        opt.condition, opt.stimulus), ...
        sprintf('%d recording(s): %s   |   top-%d cells, ordered by baseline repeat-%d latency', ...
        numel(sel), strjoin(labels, ', '), opt.n_cells, opt.sort_repeat)}, ...
        'Interpreter', 'none', 'FontWeight', 'bold');

    for di = 1:nd
        ax = nexttile(tl); hold(ax, 'on');
        Bs = Bsort{di}; Bh = Bshow{di}; Dh = Dshow{di};
        if isempty(Bs); axis(ax, 'off'); continue; end

        % selection + ordering from the sort repeat only
        pk = t >= pw(1) & t <= pw(2);
        amp = max(Bs(:, pk), [], 2, 'omitnan');
        [~, rank] = sort(amp, 'descend', 'MissingPlacement', 'last');
        keep = rank(1:min(opt.n_cells, numel(rank)));
        [~, pkidx] = max(Bs(keep, pk), [], 2, 'omitnan');
        [~, ord] = sort(pkidx, 'ascend');
        keep = keep(ord);

        % [baseline | drug] in one image, with a gap column between them
        gap = nan(numel(keep), 4);
        img = [Bh(keep, :), gap, Dh(keep, :)];
        imagesc(ax, 1:size(img, 2), 1:numel(keep), img, 'AlphaData', ~isnan(img));
        clim(ax, [lo hi]); colormap(ax, flipud(gray)); set(ax, 'YDir', 'normal');
        xlim(ax, [0.5 size(img,2)+0.5]); ylim(ax, [0.5 numel(keep)+0.5]);

        % onset markers + the divider
        n1 = numel(t);
        z0 = interp1(t, 1:n1, 0, 'linear', 'extrap');
        xline(ax, z0, 'b--', 'LineWidth', 1);
        xline(ax, n1 + 4 + z0, 'b--', 'LineWidth', 1);
        xline(ax, n1 + 2.5, 'r-', 'LineWidth', 2);
        tk = 0:2:floor(t(end));                       % seconds, both halves
        tx = interp1(t, 1:n1, tk, 'linear', 'extrap');
        set(ax, 'XTick', [tx, n1 + 4 + tx], 'XTickLabel', [ ...
            arrayfun(@(v) sprintf('%g', v), tk, 'UniformOutput', false), ...
            arrayfun(@(v) sprintf('%g', v), tk, 'UniformOutput', false)], 'FontSize', 7);
        text(ax, mean([1 n1]), -0.055, 'baseline', 'Units', 'data', 'FontSize', 8, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', [0.30 0.45 0.75]);
        text(ax, n1 + 4 + mean([1 n1]), -0.055, 'drug', 'Units', 'data', 'FontSize', 8, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'Color', [0.80 0.40 0.35]);
        title(ax, sprintf('%g\\circ  (%d cells)', dirs(di), numel(keep)), 'FontSize', 9);
        if mod(di-1, nc) == 0; ylabel(ax, 'cell (by baseline latency)'); end
        box(ax, 'off');
    end
    cb = colorbar(nexttile(tl, nd)); cb.Label.String = 'dF/F';

    if ~isempty(char(opt.save))
        savep = char(opt.save);
        for a = findall(fig, 'Type', 'axes')'      % keep the interactive toolbar
            try                                     % out of the exported image
                a.Toolbar.Visible = 'off';
            catch
            end
        end
        try
            exportgraphics(fig, savep, 'Resolution', 150);
        catch
            print(fig, savep, '-dpng', '-r150');
        end
        fprintf('  -> %s\n', savep);
    end
end
