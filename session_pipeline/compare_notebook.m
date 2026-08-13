%% COMPARE_NOTEBOOK  Exploratory comparison of two matched sessions.
% Run cell-by-cell (Ctrl+Enter). Each view is its own cell -- run only what
% you want. Heavy lifting stays in compare_sessions (metrics) and
% compare_plots (figures); this file is just orchestration. See PIPELINE.md.

%% 1  CONFIG -- the two sessions + ROI-match file + metric
clear; clc;
root = 'C:\Users\rbondarenko\projects\2P_Scripts';
addpath(fullfile(root, 'session_pipeline'));
addpath(fullfile(root, 'analysis_functions'));
od = fullfile(root, 'session_pipeline', 'organized_data');

% src1 = column 1 of the ROI-match mapping (baseline); src2 = column 2 (drug).
sess1     = fullfile(od, 'AnimalRB19_260320_1444_plane0_organized.mat');   % baseline
sess2     = fullfile(od, 'AnimalRB19_260320_1609_plane0_organized.mat');   % drug
roi_match = fullfile(root, 'roi_matching', 'roiMatch_plane0_0320_1444_1609_V1.mat');

metric = 'z';        % 'z' = baseline-subtracted z-score (default) | 'raw' = peak dF/F
save_dir = '';       % '' = show on screen; or a folder path to save PNGs

%% 2  LOAD + COMPUTE  -- matched-cell metrics for every common stimulus
R = compare_sessions(sess1, sess2, roi_match);
% R.<stim>.{raw_base,raw_drug,z_base,z_drug,mi_raw,mi_z,p_signrank,n_matched}

%% 3  PER-SESSION RECAP  -- the single-session figures you already have
explore_session(sess1);
explore_session(sess2);

%% 4  VIOLINS  -- baseline vs drug matched-cell response distributions, per stimulus
compare_plots('violins', R, 'metric', metric, 'save', figpath(save_dir, 'violins', metric));

%% 5  MODULATION INDEX  -- MI distribution per stimulus (median + signrank p)
compare_plots('modindex', R, 'metric', metric, 'save', figpath(save_dir, 'modindex', metric));

%% 6  HISTOGRAMS  -- overlaid baseline/drug response distributions
compare_plots('histograms', R, 'metric', metric, 'save', figpath(save_dir, 'histograms', metric));

%% 7  SCATTER  -- paired base-vs-drug per cell, unity line (per-cell shifts)
compare_plots('scatter', R, 'metric', metric, 'save', figpath(save_dir, 'scatter', metric));

%% 8  SUMMARY  -- one table: n_matched, medians, signrank p across stimuli
T = compare_plots('summary', R, 'metric', metric);

%% ---- local helper: build a save path only when save_dir is set ----
function p = figpath(save_dir, name, metric)
    if isempty(save_dir); p = ''; return; end
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end
    p = fullfile(save_dir, sprintf('compare_%s_%s.png', name, metric));
end
