function run_refresh_h264(dry)
%RUN_REFRESH_H264 Put the per-codec WhiskOn (H.264: on above 1.45 / off below 1.40, 2026-11-10) into the aggregated files and
%   DLC_filtered.mat of the H.264 sessions that are aggregated. Backups of the previous fields: C:\Users\rbondarenko\2P_cache\backup_whisk_v2\<s>_whisk_v2.mat
%   run_refresh_h264(true)  = dry run, nothing written
addpath('C:\Users\rbondarenko\projects\2P_Scripts\newAggregation');
T = readtable('C:\Users\rbondarenko\DLC\whisker_final_stage\stage_comparison.csv');
s = T.session(strcmp(T.codec, 'h264'));
have = cellfun(@(x) isfile(['Z:\Group Members\Rima\Aggregated\' x '_preprocessed.mat']) && isfile(['Z:\Group Members\Rima\DLCData\' x '_DLC_filtered.mat']), s);
fprintf('%d H.264 sessions, %d aggregated with a DLC_filtered file\n', numel(s), sum(have));
refresh_whisker_behavior(s(have), 'whisker_root', 'C:\Users\rbondarenko\DLC\work', 'backup_dir', 'C:\Users\rbondarenko\2P_cache\backup_whisk_v2', ...
    'backup_tag', 'v2', 'dry_run', dry);
end
