%% RUN_ORGANIZE  Entry point: organize + QC + cache a list of sessions.
%
% Edit the CONFIG block, then run this script. For each aggregated file it:
%   1. builds the canonical struct   (organize_session)
%   2. runs QC                        (qc_session: Tier-1 asserts, Tier-2 flags)
%   3. renders the QC figure          (qc_report_figure)
%   4. saves <session>_organized.mat  (+ QC_<session>.png)
%
% Tier-1 (integrity) failures throw; each session is wrapped in try/catch so
% one bad file logs an ERROR and the batch continues. Tier-2 (scientific)
% failures don't halt but stamp qc_passed=false and mark the QC figure.
%
% See session_pipeline/PIPELINE.md for the full spec.

clear; clc;

%% ====================== CONFIG ======================
% One row per session:  { preprocessed_file , plane_idx }
% plane_idx is 1-BASED; it must match the plane of the ROI-match file the
% session will be compared with (roiMatch_plane0_* -> 1, plane1_* -> 2,
% plane2_* -> 3). Output is cached per session AND plane:
%   <session>_plane<K-1>_organized.mat
agg = "Z:\joeschgrp\Group Members\Rima\Aggregated\";
session_files = {
    agg + "AnimalRB19_260305_1327_preprocessed.mat", 1    % guanfacine pair, plane0
    agg + "AnimalRB19_260305_1409_preprocessed.mat", 1
    agg + "AnimalRB19_260320_1444_preprocessed.mat", 3    % guanfacine pair, plane2
    agg + "AnimalRB19_260320_1609_preprocessed.mat", 3
    agg + "AnimalRB19_260312_1115_preprocessed.mat", 2    % yohimbine pair, plane1
    agg + "AnimalRB19_260312_1243_preprocessed.mat", 2
    agg + "AnimalRB19_260310_1246_preprocessed.mat", 3    % yohimbine pair, matched
    agg + "AnimalRB19_260310_1341_preprocessed.mat", 1    % ACROSS planes (z-drift):
};                                                        % 1246 plane2 <-> 1341 plane0

cfg = struct();
cfg.ca_type            = 1;        % 1=Ca_dFF, 2=Ca_deconvolved, 3=Ca_F
cfg.fs_range           = [2 30];   % plausible sampling-rate bounds (Hz)
cfg.nan_frac_tol       = 0.01;     % max fraction of non-finite raw dF/F

skip_existing = true;              % don't re-organize a session already cached

out_dir = fullfile(fileparts(mfilename('fullpath')), 'organized_data');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% ====================== RUN ======================
here = fileparts(mfilename('fullpath'));
addpath(here);

summary = struct('session', {}, 'status', {}, 'qc_passed', {}, 'n_fail', {}, 'msg', {});
n_sess = size(session_files, 1);

for k = 1:n_sess
    f     = char(session_files{k, 1});
    plane = session_files{k, 2};
    fprintf('\n########## [%d/%d] %s  (plane idx %d) ##########\n', k, n_sess, f, plane);
    row = struct('session', sprintf('%s | plane%d', f, plane - 1), 'status', 'ERROR', ...
                 'qc_passed', false, 'n_fail', NaN, 'msg', '');

    out_mat = organized_path(f, plane, out_dir);
    if skip_existing && isfile(out_mat)
        C = load(out_mat, 'meta');
        row.status = 'CACHED'; row.qc_passed = C.meta.qc_passed;
        if isfield(C.meta, 'qc') && isfield(C.meta.qc, 'n_fail'); row.n_fail = C.meta.qc.n_fail; end
        fprintf('  cached -> %s (skipped)\n', out_mat);
        summary(end+1) = row; %#ok<SAGROW>
        continue;
    end

    try
        scfg = cfg; scfg.selected_plane_idx = plane;
        B = load(f);
        S = organize_session(f, scfg);
        [S, qc] = qc_session(S, B, scfg);

        [~, sid] = organized_path(f, plane, out_dir);
        out_png  = fullfile(out_dir, sprintf('QC_%s_plane%d.png', sid, plane - 1));
        save(out_mat, '-struct', 'S', '-v7.3');
        qc_report_figure(S, out_png);

        row.status    = 'OK';
        row.qc_passed = S.meta.qc_passed;
        row.n_fail    = qc.n_fail;
        fprintf('  saved -> %s\n', out_mat);
    catch ME
        row.msg = ME.message;
        fprintf(2, '  !! FAILED: %s\n', ME.message);
    end
    summary(end+1) = row; %#ok<SAGROW>
end

%% ====================== SUMMARY ======================
fprintf('\n================= RUN SUMMARY =================\n');
for k = 1:numel(summary)
    r = summary(k);
    parts = strsplit(r.session, ' | ');
    [~, sid] = fileparts(parts{1});
    if numel(parts) > 1; sid = sprintf('%s %s', sid, parts{2}); end
    if any(strcmp(r.status, {'OK', 'CACHED'}))
        qc_tag = 'QC-PASS'; if ~r.qc_passed; qc_tag = 'QC-FAIL'; end
        fprintf('  %-52s  %-6s  %s (%d checks failed)\n', sid, r.status, qc_tag, r.n_fail);
    else
        fprintf('  %-52s  %-6s  %s\n', sid, r.status, r.msg);
    end
end
fprintf('==============================================\n');
