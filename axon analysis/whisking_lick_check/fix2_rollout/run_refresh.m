function run_refresh(dry)
addpath('C:\Users\rbondarenko\projects\2P_Scripts\newAggregation');
f = dir('Z:\Group Members\Rima\Aggregated\AnimalRB2*_preprocessed.mat');
s = erase({f.name}, '_preprocessed.mat');
s = s(~cellfun(@isempty, regexp(s, '^AnimalRB2[234]_')));
fprintf('%d aggregated sessions\n', numel(s));
refresh_whisker_behavior(s, 'whisker_root', 'C:\Users\rbondarenko\DLC\work', 'dry_run', dry);
end
