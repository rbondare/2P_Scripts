function T = recdb_load(name)
% RECDB_LOAD  read one table of the RB22-24 recordings database as a MATLAB table.
%   T = recdb_load('sessions')      % also 'blocks', 'rois', 'axons', 'animals', 'notes'
%   R = recdb_load('rois');  R = R(strcmp(R.session, 'AnimalRB24_261105_1423') & strcmp(R.accepted, 'yes'), :);
% Text columns come back as cellstr (blank = ''), number columns as double (blank = NaN); yes/no columns stay text ('yes' / 'no').
% Database folder: Z:\Group Members\Rima\Recordings_DB (environment variable RECDB_DIR overrides; falls back to the git snapshot).
here = fileparts(mfilename('fullpath'));
cands = {getenv('RECDB_DIR'), 'Z:\joeschgrp\Group Members\Rima\Recordings_DB', 'Z:\Group Members\Rima\Recordings_DB', fullfile(here, 'snapshot')};
d = '';
for i = 1:numel(cands)
    if ~isempty(cands{i}) && isfile(fullfile(cands{i}, 'sessions.csv')); d = cands{i}; break; end
end
if isempty(d); error('recdb_load:NoDatabase', 'no recordings database found (tried: %s)', strjoin(cands(~cellfun(@isempty, cands)), ', ')); end
opts = detectImportOptions(fullfile(d, [name '.csv']), 'TextType', 'char', 'VariableNamingRule', 'preserve');
opts = setvartype(opts, opts.VariableNames(strcmp(opts.VariableTypes, 'datetime') | strcmp(opts.VariableTypes, 'duration')), 'char');
for k = 1:numel(opts.VariableNames)
    if strcmp(opts.VariableTypes{k}, 'char') || strcmp(opts.VariableTypes{k}, 'string'); opts = setvaropts(opts, opts.VariableNames{k}, 'FillValue', ''); end
end
T = readtable(fullfile(d, [name '.csv']), opts);
end
