function verify_refresh_h264()
%VERIFY_REFRESH_H264 After run_refresh_h264(false): the aggregated files and DLC_filtered.mat of the H.264 sessions hold the new WhiskOn
%   and nothing else changed. Checks per session:
%     on==csv     Behav.WhiskOn equals whisk_on of the saved csv at the selected frames (NaN on artifact frames included)
%     others==csv WhiskIndex / WhiskAmp / WhiskArtifact / FacePadME / JawME equal the csv at the selected frames
%     same as old every whisking field except WhiskOn is identical to the backup of the previous version
%     new>=old    every frame that was on is still on (the lower levels can only add whisking)
%     nPlus       frames turned on by the new levels (aggregated rows)
%     levels      VideoParams.Whisker.whisk_off_below_index / whisk_on_above_index = 1.40 / 1.45
%     dlc         DLC_filtered Behav.WhiskOn == csv over all frames and its VideoParams.Whisker has the same levels
%     rest==dlc   every other Behav field equals DLC_filtered Behav(selected_frameIdx)
%     size        file size change vs before the refresh, MB (sizes recorded before: agg_sizes_before)
f = dir('Z:\Group Members\Rima\Aggregated\AnimalRB24_26*_preprocessed.mat');
s = erase({f.name}, '_preprocessed.mat');
T0 = readtable('C:\Users\rbondarenko\DLC\whisker_final_stage\stage_comparison.csv');
s = s(ismember(s, T0.session(strcmp(T0.codec, 'h264'))));
before = containers.Map();   % bytes before the refresh (listing saved by the rollout)
fid = fopen(fullfile(fileparts(mfilename('fullpath')), 'agg_sizes_before_h264_refresh.txt'));
while true
    ln = fgetl(fid); if ~ischar(ln); break; end
    p = strsplit(strtrim(ln)); before(erase(p{2}, '_preprocessed.mat')) = str2double(p{1});
end
fclose(fid);
whisk = {'WhiskIndex', 'WhiskAmp', 'WhiskArtifact', 'FacePadME', 'JawME'};
bad = 0;
fprintf('%-24s %6s %10s %9s %7s %8s %7s %5s %9s %8s\n', 'session', 'on==csv', 'others==csv', 'same as old', 'new>=old', 'nPlus', 'levels', 'dlc', 'rest==dlc', 'size MB');
for i = 1:numel(s)
    F = ['Z:\Group Members\Rima\Aggregated\' s{i} '_preprocessed.mat'];
    D = load(F, 'Behav', 'VideoParams');
    B = D.Behav; VP = D.VideoParams; sel = VP.selected_frameIdx;
    C = readtable(['C:\Users\rbondarenko\DLC\work\' s{i} '\whisker_results\' s{i} '_whisker.csv']);
    O = load(['C:\Users\rbondarenko\2P_cache\backup_whisk_v2\' s{i} '_whisk_v2.mat'], 'old');
    old = O.old.aggregated_Behav;
    okon = isequaln(B.WhiskOn, C.whisk_on(sel));
    csvmap = struct('WhiskIndex', 'whisk_index', 'WhiskAmp', 'whisk_amp', 'WhiskArtifact', 'artifact', 'FacePadME', 'me_pad', 'JawME', 'me_fan_down');
    okcsv = true; oksame = true;
    for q = 1:numel(whisk)
        okcsv = okcsv && isequaln(B.(whisk{q}), C.(csvmap.(whisk{q}))(sel));
        oksame = oksame && isequaln(B.(whisk{q}), old.(whisk{q}));
    end
    wasOn = old.WhiskOn > 0; isOn = B.WhiskOn > 0;
    oksup = ~any(wasOn & ~isOn);
    nplus = sum(~wasOn & isOn);
    oklev = isfield(VP.Whisker, 'whisk_on_above_index') && abs(VP.Whisker.whisk_on_above_index - 1.45) < 1e-9 && abs(VP.Whisker.whisk_off_below_index - 1.40) < 1e-9;
    L = load(['Z:\Group Members\Rima\DLCData\' s{i} '_DLC_filtered.mat'], 'Behav', 'VideoParams');
    okdlc = isequaln(L.Behav.WhiskOn, C.whisk_on) && isfield(L.VideoParams.Whisker, 'whisk_on_above_index') && abs(L.VideoParams.Whisker.whisk_on_above_index - 1.45) < 1e-9;
    % every other Behav field equals the DLC_filtered one at the selected frames
    fn = setdiff(fieldnames(B), [whisk, {'WhiskOn'}]);
    okrest = true; nchk = 0;
    for q = 1:numel(fn)
        if isfield(L.Behav, fn{q}) && size(L.Behav.(fn{q}), 1) == numel(L.Behav.FrameNum) && size(B.(fn{q}), 1) == numel(sel)
            nchk = nchk + 1;
            okrest = okrest && isequaln(B.(fn{q}), L.Behav.(fn{q})(sel, :));
        end
    end
    d = dir(F);
    dsz = (d.bytes - before(s{i})) / 1e6;
    ok = okon && okcsv && oksame && oksup && oklev && okdlc && okrest && abs(dsz) < 1;
    fprintf('%-24s %6d %10d %9d %7d %8d %7d %5d %9d (%d fields) %7.2f\n', s{i}, okon, okcsv, oksame, oksup, nplus, oklev, okdlc, okrest, nchk, dsz);
    bad = bad + ~ok;
end
fprintf('\nproblems: %d\n', bad);
end
