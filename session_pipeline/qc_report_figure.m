function qc_report_figure(S, out_png)
% QC_REPORT_FIGURE  Render the per-session timing-alignment verification
% figure (see PIPELINE.md). One row per stimulus: population-mean response
% (left) + raw cells x time heatmap (right), onset marked. A CORRECT
% alignment shows the response band at a consistent latency; a timing bug
% smears it. Title is stamped PASS / *** QC FAILED ***.
%
%   qc_report_figure(S, out_png)

    fns = fieldnames(S.stimuli);
    n = numel(fns);
    if n == 0
        warning('qc_report_figure: no stimuli in %s', S.meta.session_id);
        return;
    end

    fig = figure('Visible', 'off', 'Position', [100 100 1000 max(260, 230*n)]);
    tl = tiledlayout(fig, n, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    passed = S.meta.qc_passed;
    if passed; status = 'PASS'; else; status = '*** QC FAILED ***'; end
    title(tl, sprintf('%s   |   plane %d   |   fs %.2f Hz   |   %d ROIs   |   QC: %s', ...
        S.meta.session_id, S.meta.selected_plane_idx, S.meta.fs_hz, S.meta.n_rois, status), ...
        'Interpreter', 'none', 'FontWeight', 'bold');

    for i = 1:n
        E = S.stimuli.(fns{i});
        switch E.paradigm
            case 'trial_onset'
                [x, mat, subtitle] = trial_onset_panel(E);
            case 'direction'
                [x, mat, subtitle] = direction_panel(E);
            case 'block'
                [x, mat, subtitle] = block_panel(E);
            otherwise
                continue;
        end
        pop = mean(mat, 1, 'omitnan');

        % ---- left: population-mean trace ----
        ax1 = nexttile(tl);
        plot(ax1, x, pop, 'k-', 'LineWidth', 1.6); hold(ax1, 'on');
        if ~strcmp(E.paradigm, 'block')
            xline(ax1, 0, 'r--', 'LineWidth', 1);
        end
        xlim(ax1, [x(1) x(end)]); box(ax1, 'off');
        ylabel(ax1, 'pop. mean dF/F');
        title(ax1, subtitle, 'Interpreter', 'none', 'FontSize', 9);
        if i == n; xlabel(ax1, 'time (s)'); end

        % ---- right: raw heatmap (unsorted -- shows the raw response band) ----
        ax2 = nexttile(tl);
        heatmap_sorted(ax2, x, mat);        % no sort_window => rows as-is
        hold(ax2, 'on');
        if ~strcmp(E.paradigm, 'block')
            xline(ax2, 0, 'r--', 'LineWidth', 1);
        end
        ylabel(ax2, 'ROI');
        title(ax2, sprintf('%s raw (cells x time)', E.type), 'Interpreter', 'none', 'FontSize', 9);
        if i == n; xlabel(ax2, 'time (s)'); end
    end

    try
        exportgraphics(fig, out_png, 'Resolution', 150);
    catch
        print(fig, out_png, '-dpng', '-r150');
    end
    close(fig);
    fprintf('    QC figure -> %s\n', out_png);
end

% ----------------------------------------------------------------------
function [x, mat, subtitle] = trial_onset_panel(E)
    x = E.t_com; mat = E.avg;
    rs = tern(islogical(E.RedFrameSynchronized) && E.RedFrameSynchronized, 'sync', 'NO-SYNC');
    subtitle = sprintf('%s | %d trials | %s | median-avg', E.type, E.n_trials, rs);
end

function [x, mat, subtitle] = direction_panel(E)
    x = E.t_com; mat = E.pref_resp;
    v = '';
    if isfield(E,'verify')
        if ~E.verify.applicable
            v = ' | recon N/A';
        else
            v = sprintf(' | recon %.2f/%.2f', E.verify.r_same, E.verify.r_shifted);
        end
    end
    subtitle = sprintf('%s | %d dirs | pref-dir collapse%s', E.type, numel(E.dir_values), v);
end

function [x, mat, subtitle] = block_panel(E)
    x = E.block_t; mat = E.block;
    subtitle = sprintf('%s | continuous block | %.0fs', E.type, E.duration_sec);
end

function s = tern(c, a, b)
    if c; s = a; else; s = b; end
end
