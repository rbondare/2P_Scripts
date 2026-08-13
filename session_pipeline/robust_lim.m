function lim = robust_lim(x, p)
% ROBUST_LIM  Axis limits from percentiles (default [1 99]) so a handful of
% extreme z-score outliers don't squash the informative part of a plot.
% Inference (median/IQR/signrank) is unaffected -- this is display only.
% Shared by compare_plots and group_plots.
    if nargin < 2 || isempty(p); p = [1 99]; end
    x = x(isfinite(x));
    if isempty(x); lim = [0 1]; return; end
    lo = prctile(x, p(1)); hi = prctile(x, p(2));
    if ~(hi > lo); lo = min(x); hi = max(x); end
    if ~(hi > lo); hi = lo + 1; end
    pad = 0.05 * (hi - lo);
    lim = [lo - pad, hi + pad];
end
