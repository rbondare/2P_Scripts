% BUILD_RECDB  rebuild the RB22-24 recordings database (Z:\Group Members\Rima\Recordings_DB).
%   Stage A (MATLAB): export_session_meta  -> stimulus blocks + aggregated-file QC  (_build\m_sessions.csv, m_blocks.csv)
%   Stage B (Python): build_recdb.py       -> sessions.csv blocks.csv rois.csv axons.csv, FOV pictures, index.html
% Safe to re-run any time: the auto columns are regenerated, notes.csv / animals.csv (your hand-edited files) are never overwritten.
% Stage A re-reads a session only when its aggregated file or stim files changed (set force = true to redo everything).
force = false;          % true: re-read every session in stage A
snapshot = true;        % true: also copy the CSVs to axon analysis\recordings_db\snapshot (for git)

here = fileparts(mfilename('fullpath'));
addpath(here);
export_session_meta([], 'force', force);

py = 'C:\Users\rbondarenko\miniforge3\envs\suite2p_venv\python.exe';
args = '';  if snapshot; args = ' --snapshot'; end
cmd = sprintf('"%s" "%s"%s', py, fullfile(here, 'build_recdb.py'), args);
fprintf('running: %s\n', cmd);
status = system(cmd);
if status ~= 0; error('build_recdb:python', 'build_recdb.py failed (exit code %d)', status); end

% quick look
S = recdb_load('sessions');
fprintf('\n%d sessions: %s\n', height(S), strjoin(cellfun(@(s) sprintf('%s %d', s, sum(strcmp(S.status, s))), unique(S.status)', 'UniformOutput', false), ', '));
disp(S(:, {'session', 'status', 'task_type', 'stimuli', 'n_roi_accepted', 'pool'}));
