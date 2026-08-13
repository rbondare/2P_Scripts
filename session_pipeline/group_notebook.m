%% GROUP_NOTEBOOK  Several matched baseline/drug recordings, averaged per drug.
% Run cell-by-cell (Ctrl+Enter). Heavy lifting stays in compare_group (which
% only calls the validated compare_sessions on each pair) and group_plots.
% See PIPELINE.md.
%
% PREREQUISITE: every session below must be organized AT THE PLANE OF ITS
% ROI-MATCH FILE. run_organize.m's config already lists exactly these eight.

%% 1  CONFIG -- one entry per baseline/drug recording pair
clear; clc;
root = 'C:\Users\rbondarenko\projects\2P_Scripts';
addpath(fullfile(root, 'session_pipeline'));
addpath(fullfile(root, 'analysis_functions'));

% The cohort lives in group_config.m (single source of truth -- edit it to
% add/remove a recording pair). baseline_file = column 1 of the ROI-match
% mapping, drug_file = column 2; the plane comes from the roiMatch filename.
pairs = group_config();
disp(struct2table(pairs));

metric   = 'z';      % 'z' = baseline-subtracted z-score | 'raw' = peak dF/F
save_dir = '';       % '' = show on screen; or a folder to save PNGs

%% 2  COMPUTE  -- pairwise comparison of every pair, then pool per condition
% Prints, per pair: matched-cell count + the ROI-match sanity check
% (matched centroids must be far closer than shuffled -- this is what would
% catch a wrong plane or a mislabelled roiMatch file).
G = compare_group(pairs);

%% 3  TRACES  -- population response averaged ACROSS recordings
% Thick = grand mean over recordings +- s.e.m. across recordings;
% thin = each recording. This is the literal "average across recordings".
group_plots('traces', G, 'save', figpath(save_dir, 'traces', metric));

%% 4  SUPERPLOT  -- pooled cells + each recording's median
% Violin = every matched cell from every recording; open circles = one
% recording each. If the circles all move the same way, the pooled effect is
% real; if one recording carries it, you see that immediately.
group_plots('superplot', G, 'metric', metric, 'save', figpath(save_dir, 'superplot', metric));

%% 5  MODULATION INDEX  -- pooled MI + per-recording medians
group_plots('mi', G, 'metric', metric, 'save', figpath(save_dir, 'mi', metric));

%% 6  FOREST  -- per-recording median MI vs the condition pooled median
group_plots('forest', G, 'metric', metric, 'save', figpath(save_dir, 'forest', metric));

%% 7  SUMMARY TABLE  -- pooled medians/p AND per-recording consistency
T = group_plots('summary', G, 'metric', metric);

%% 7b TRAVELING WAVE  -- baseline | drug per direction, cells pooled per drug
% The paired version of the single-session traveling-wave figure: same cells,
% same colour scale, ordered by baseline latency, so both the amplitude change
% and whether the latency SEQUENCE survives are visible at once.
for c = fieldnames(G.cond)'
    cond = G.cond.(c{1}).meta.condition;
    group_travelingwave(G, 'condition', cond, ...
        'save', figpath(save_dir, ['travelingwave_' cond], metric));
end

%% 8  BETWEEN-DRUG CONTRAST  -- guanfacine vs yohimbine (pooled MI, ranksum)
disp(struct2table(structfun(@(c) c, G.contrast, 'UniformOutput', false)));

%% 9  DRILL DOWN  -- any single pair, with the existing pairwise figures
k = 1;                                   % index into G.pairs
fprintf('drill-down: %s (%s)\n', G.pairs(k).label, G.pairs(k).condition);
compare_plots('violins', G.pairs(k).R, 'metric', metric);

%% ---- local helper ----
function p = figpath(save_dir, name, metric)
    if isempty(save_dir); p = ''; return; end
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end
    p = fullfile(save_dir, sprintf('group_%s_%s.png', name, metric));
end
