function results = run_aggregation_batch(recordings, varargin)
%RUN_AGGREGATION_BATCH Aggregate the named recordings, no prompts (for MATLAB -batch or a script)
%
%   run_aggregation_batch('AnimalRB24_261023_1538')
%   run_aggregation_batch({'AnimalRB24_261023_1538', 'AnimalRB24_261030_1744'})
%   run_aggregation_batch({'Z:\Group Members\Rima\DATA_2P\AnimalRB24_261030_1744'}, 'overwrite', true)
%   results = run_aggregation_batch(...)
%
%   Same steps as run_aggregation with config_source = 'manual', for exactly the recordings you list
%   (names, or folder paths of which only the last part is used):
%     1. every recording needs a suite2p folder in DATA_2P and stim*.mat in the stimulus folder
%     2. DLC import (import_video_behavior) for these recordings only; it archives their camera folder to
%        preCamera and needs <Camera>\<recording>\success.log, otherwise it is skipped. If a recording is
%        aggregated before its DLC results exist, run this again for it later (without 'overwrite'): the
%        video data is then added to the existing output
%     3. save_preprocessed for each recording; a recording that fails is reported and the others go on
%
%   Optional Name-Value Arguments:
%       'overwrite'  - false (default): recordings that already have an output are not aggregated again
%                      true: the existing output is deleted first (overwriting a v7.3 file in place makes it grow)
%       'import_dlc' - run the DLC import (default: true)
%       'ca_type', 'ca_format', 'cal_idx' - passed to save_preprocessed (defaults 1, 'suite2p', [])
%
%   Output (optional): struct array with name, status ('done', 'skipped' or 'failed'), message, output.
%   Without an output argument an error is raised at the end if any recording failed, so that
%   matlab -batch returns a non-zero exit code.
%
%   From a shell (the first run of a session reads all TIFF headers, ~30-40 min; later runs ~30 s):
%     matlab -batch "addpath('...\2P_Scripts\newAggregation'); run_aggregation_batch({'AnimalRB24_261023_1538'})"
%   The multi-acquisition question of save_preprocessed (non-monotonic frame timestamps) cannot be answered
%   in batch mode; such a recording fails and has to be aggregated from the MATLAB window.
%
%   See also: run_aggregation, save_preprocessed, import_video_behavior, get_experiment_list

%% Parse inputs
p = inputParser;
addRequired(p, 'recordings', @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(p, 'overwrite', false, @islogical);
addParameter(p, 'import_dlc', true, @islogical);
addParameter(p, 'ca_type', 1, @(x) isnumeric(x) && x >= 0 && x <= 3);
addParameter(p, 'ca_format', 'suite2p', @(x) ismember(x, {'auto', 'suite2p', 'caiman'}));
addParameter(p, 'cal_idx', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 1));
parse(p, recordings, varargin{:});
overwrite = p.Results.overwrite;

addpath(fileparts(mfilename('fullpath')));

names = cellfun(@recording_name, cellstr(recordings), 'UniformOutput', false);
names = unique(names(~cellfun(@isempty, names)), 'stable');
if isempty(names)
    error('run_aggregation_batch:NoRecordings', 'No recordings given');
end
n = numel(names);

results = struct('name', names(:), 'status', 'pending', 'message', '', 'output', '');
experiments = get_experiment_list('source', 'manual', 'recordings', names);

fprintf('=== Aggregating %d recording(s) ===\n', n);

%% Check what exists
todo = false(1, n);
for i = 1:n
    exp = experiments(i);
    out = fullfile(exp.config.OutputDir, [exp.name '_preprocessed.mat']);
    results(i).output = out;
    if ~exp.has_calcium
        results(i).status = 'failed';
        results(i).message = sprintf('no suite2p folder in %s', fullfile(exp.config.DATA_2P, exp.name));
    elseif ~exp.has_stimulus
        results(i).status = 'failed';
        results(i).message = sprintf('no stim*.mat in %s', fullfile(exp.config.StimBasePath, exp.animal_name, exp.name));
    elseif exp.is_processed && ~overwrite
        results(i).status = 'skipped';
        results(i).message = 'output exists (use ''overwrite'', true to aggregate again)';
    else
        todo(i) = true;
    end
    if exp.is_processed && overwrite && todo(i)
        delete(out);
        fprintf('Deleted old output: %s\n', out);
    end
end

%% Import DLC data (also for skipped recordings: adds the video data to an existing output that lacks it)
if p.Results.import_dlc
    importable = todo | strcmp({results.status}, 'skipped');
    animals = unique({experiments(importable).animal_name}, 'stable');
    for a = 1:numel(animals)
        in_animal = importable & strcmp({experiments.animal_name}, animals{a});
        fprintf('\n=== Importing DLC data: %s ===\n', strjoin({experiments(in_animal).name}, ', '));
        try
            import_video_behavior('config', experiments(find(in_animal, 1)).config, 'animals', animals(a), ...
                'sessions', {experiments(in_animal).name});
        catch ME
            fprintf('DLC import error: %s\n', ME.message);
        end
    end
end

%% Aggregate
for i = find(todo)
    exp = experiments(i);
    fprintf('\n[%d/%d] Aggregating: %s\n', i, n, exp.name);
    t0 = tic;
    try
        results(i).output = save_preprocessed(fullfile(exp.config.DATA_2P, exp.name), ...
            'ca_type', p.Results.ca_type, 'ca_format', p.Results.ca_format, 'config', exp.config, ...
            'cal_idx', p.Results.cal_idx, 'overwrite_intermediate', false);
        results(i).status = 'done';
        results(i).message = sprintf('%.1f min', toc(t0) / 60);
    catch ME
        results(i).status = 'failed';
        results(i).message = ME.message;
        if ~isempty(ME.stack)
            results(i).message = sprintf('%s (%s:%d)', ME.message, ME.stack(1).name, ME.stack(1).line);
        end
        fprintf('  ERROR: %s\n', results(i).message);
    end
end

%% Summary
fprintf('\n=== Summary ===\n');
for i = 1:n
    fprintf('  %-8s %-32s %s\n', results(i).status, results(i).name, results(i).message);
end
n_failed = nnz(strcmp({results.status}, 'failed'));
fprintf('Done: %d, skipped: %d, failed: %d\n', nnz(strcmp({results.status}, 'done')), ...
    nnz(strcmp({results.status}, 'skipped')), n_failed);

if nargout == 0
    if n_failed > 0
        error('run_aggregation_batch:Failed', '%d of %d recording(s) failed, see the summary above', n_failed, n);
    end
    clear results;
end

end


function name = recording_name(x)
%RECORDING_NAME Last part of a recording name or folder path
parts = regexp(x, '[\\/]', 'split');
parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    name = '';
else
    name = strtrim(parts{end});
end
end
