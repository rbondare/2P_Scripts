function [names, T] = recdb_select(varargin)
% RECDB_SELECT  pick recordings from the RB22-24 recordings database (sessions.csv).
%
%   names = recdb_select('animal', 'RB24', 'has_block', 'loom_expand', 'min_rois', 10, 'pool', true)
%   [names, T] = recdb_select('tags', {'task'}, 'exclude_tags', {'stale'}, 'date_from', '2026-10-27')
%
% Returns the session names (cellstr, sorted by animal / date) and the matching rows of sessions.csv (table).
% Filters, all optional and combined with AND (same names as recdb.py):
%   animal        'RB24' | {'RB23','RB24'}
%   status        'analysed' | 'aggregated' | 'suite2p done' | 'raw only' (cellstr allowed)
%   task_type     'open-loop visual' | 'lick task' | 'passive square' | 'spontaneous' | 'corridor (visuomotor)'
%   pool          true = pool=='yes' (usable in pooled statistics), false = pool=='no'
%   tags          all of these must be in all_tags         any_tags: at least one         exclude_tags: none of these
%   has_block     key in blocks.csv: flash fff bar loom_expand loom_recede square spont alleyCL;  min_trials: total completed trials >= this
%   min_rois      accepted ROIs >= this        min_minutes      recording length >= this
%   date_from / date_to   'yyyy-mm-dd'          aggregated  true | false
% Database folder: Z:\Group Members\Rima\Recordings_DB (environment variable RECDB_DIR overrides; falls back to the git snapshot).
% Use recdb_load('rois') / recdb_load('blocks') for the ROI and stimulus-block tables of the selected sessions.

p = inputParser;
for f = {'animal', 'status', 'task_type', 'tags', 'any_tags', 'exclude_tags', 'has_block', 'date_from', 'date_to'}; addParameter(p, f{1}, {}); end
addParameter(p, 'pool', []);  addParameter(p, 'min_trials', []);  addParameter(p, 'min_rois', []);  addParameter(p, 'min_minutes', []);  addParameter(p, 'aggregated', []);
parse(p, varargin{:});  q = p.Results;

T = recdb_load('sessions');
keep = true(height(T), 1);
tagsets = cellfun(@(s) split_(s), T.all_tags, 'UniformOutput', false);
if ~isempty(q.animal);       keep = keep & ismember(T.animal, cellstr(strrep(string(q.animal), 'Animal', ''))); end
if ~isempty(q.status);       keep = keep & ismember(T.status, cellstr(q.status)); end
if ~isempty(q.task_type);    keep = keep & strcmp(T.task_type, char(q.task_type)); end
if ~isempty(q.pool);         keep = keep & strcmp(T.pool, ternary_(q.pool, 'yes', 'no')); end
if ~isempty(q.tags);         need = cellstr(q.tags);        keep = keep & cellfun(@(t) all(ismember(need, t)), tagsets); end
if ~isempty(q.any_tags);     any_ = cellstr(q.any_tags);    keep = keep & cellfun(@(t) any(ismember(any_, t)), tagsets); end
if ~isempty(q.exclude_tags); ex = cellstr(q.exclude_tags);  keep = keep & ~cellfun(@(t) any(ismember(ex, t)), tagsets); end
if ~isempty(q.has_block)
    B = recdb_load('blocks');  B = B(strcmp(B.key, char(q.has_block)), :);
    tot = varfun(@(x) sum(x, 'omitnan'), B, 'InputVariables', 'n_trials', 'GroupingVariables', 'session');
    ok = ismember(T.session, tot.session);
    if ~isempty(q.min_trials); good = tot.session(tot.Fun_n_trials >= q.min_trials); ok = ismember(T.session, good); end
    keep = keep & ok;
end
if ~isempty(q.min_rois);     keep = keep & T.n_roi_accepted >= q.min_rois;  end
if ~isempty(q.min_minutes);  keep = keep & T.minutes >= q.min_minutes;  end
if ~isempty(q.date_from);    keep = keep & string(T.date) >= string(q.date_from);  end
if ~isempty(q.date_to);      keep = keep & string(T.date) <= string(q.date_to);  end
if ~isempty(q.aggregated);   keep = keep & strcmp(T.has_aggregated, ternary_(q.aggregated, 'yes', 'no')); end
T = sortrows(T(keep, :), {'animal', 'date', 'time'});
names = T.session(:)';
end

function c = split_(s)
if isempty(s) || (isnumeric(s) && isnan(s)); c = {}; else; c = cellstr(strsplit(char(s), ';')); c(cellfun(@isempty, c)) = []; end
end

function r = ternary_(c, a, b); if c; r = a; else; r = b; end; end
