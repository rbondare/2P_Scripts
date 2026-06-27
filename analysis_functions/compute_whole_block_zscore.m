function [z_metric, Z_trace] = compute_whole_block_zscore(block_dff, match_idx)
    % For stimuli with no discrete onset/trial structure (spontaneous,
    % checkers2): there's no pre/post window to z-score a cell's own trace
    % against. Z-scoring a trace against its OWN mean/std and then
    % averaging is mathematically guaranteed to give ~0 for every cell (the
    % mean of a self-z-scored signal is its own, now-zero, mean) -- not a
    % usable metric, just floating-point noise around an exact zero.
    %
    % Instead: z-score each cell's mean dF/F over the block AGAINST THE
    % POPULATION of cells in that same session (same plane, same
    % condition): z = (cell_mean - population_mean) / population_std. This
    % reflects a cell's activity level relative to its own session's
    % population, in SD units -- comparable in spirit to the other
    % stimuli's z-score (a normalized, non-raw view of the same data) even
    % though it's population- rather than onset-relative. A uniform,
    % population-wide shift (every cell scaling/shifting together) will
    % still wash out here by construction -- that's exactly what the raw
    % dF/F metric (computed separately) is for; this one is sensitive to
    % CHANGES IN RELATIVE RANKING across cells instead.
    sel = block_dff(match_idx, :);
    cell_mean = mean(sel, 2, 'omitnan');             % [n_cells x 1]
    pop_mean = mean(cell_mean, 'omitnan');
    pop_std  = std(cell_mean, 'omitnan');
    z_metric = (cell_mean - pop_mean) / pop_std;     % [n_cells x 1]
    Z_trace  = (sel - pop_mean) / pop_std;           % time-resolved, same scaling, for the heatmap
end
