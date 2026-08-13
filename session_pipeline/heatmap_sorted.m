function heatmap_sorted(ax, x, mat, sort_window)
% HEATMAP_SORTED  Shared styled cells x time heatmap for the session pipeline.
% Centralizes the visual convention used by qc_report_figure / explore_session:
% flipud(gray) colormap, robust [min(0,p1) .. p99] color limits, colorbar,
% and YDir normal (ROI 1 at bottom).
%
%   heatmap_sorted(ax, x, mat)                 rows in given order (no sort)
%   heatmap_sorted(ax, x, mat, [lo hi])        rows sorted (desc) by mean
%                                              response in the x-window [lo hi]
%   heatmap_sorted(ax, x, mat, [])             same as no sort
%
% mat is [n_rois x numel(x)]. Overlays (onset xline, titles) are the caller's
% job -- this only draws + styles the image.

    if nargin < 4 || isempty(sort_window)
        ord = 1:size(mat, 1);
    else
        pk = x >= sort_window(1) & x <= sort_window(2);
        if ~any(pk); pk = true(1, numel(x)); end
        [~, ord] = sort(mean(mat(:, pk), 2, 'omitnan'), 'descend');
    end

    imagesc(ax, x, 1:size(mat, 1), mat(ord, :));
    hi = prctile(mat(:), 99);
    if ~isfinite(hi) || hi <= 0; hi = max(mat(:)); end
    lo = min(0, prctile(mat(:), 1));
    if isfinite(lo) && isfinite(hi) && hi > lo; clim(ax, [lo hi]); end
    colormap(ax, flipud(gray)); colorbar(ax); set(ax, 'YDir', 'normal');
end
