function verify_refresh()
f = dir('Z:\Group Members\Rima\Aggregated\AnimalRB2*_preprocessed.mat');
s = erase({f.name}, '_preprocessed.mat'); s = s(~cellfun(@isempty, regexp(s, '^AnimalRB2[234]_')));
bad = 0;
fprintf('%-24s %6s %8s %9s %9s %9s %9s %9s\n', 'session', 'jaw', 'box', 'idx==csv', 'jaw==csv', 'oldArt<=new', 'otherOK', 'size MB');
for i = 1:numel(s)
    F = ['Z:\Group Members\Rima\Aggregated\' s{i} '_preprocessed.mat'];
    D = load(F, 'Behav', 'VideoParams');
    B = D.Behav; VP = D.VideoParams;
    T = readtable(['C:\Users\rbondarenko\DLC\work\' s{i} '\whisker_results\' s{i} '_whisker.csv']);
    sel = VP.selected_frameIdx;
    okidx = isequaln(B.WhiskIndex, T.whisk_index(sel));
    okjaw = isequaln(B.JawME, T.me_fan_down(sel));
    okbox = isfield(VP.Whisker, 'whisker_box') && strcmp(VP.Whisker.whisker_box, 'front_top');
    O = load(['C:\Users\rbondarenko\2P_cache\backup_whisk_v1\' s{i} '_whisk_v1.mat'], 'old');
    oldart = O.old.aggregated_Behav.WhiskArtifact > 0;
    okart = all(B.WhiskArtifact(oldart) > 0);
    % an untouched field: pupil area equals the DLC_filtered one at the selected frames
    L = load(['Z:\Group Members\Rima\DLCData\' s{i} '_DLC_filtered.mat'], 'Behav');
    okother = isequaln(B.PupilArea, L.Behav.PupilArea(sel)) && isequaln(B.FacePadME, T.me_pad(sel));
    d = dir(F);
    fprintf('%-24s %6d %8d %9d %9d %9d %9d %9.1f\n', s{i}, okjaw, okbox, okidx, okjaw, okart, okother, d.bytes/1e6);
    bad = bad + ~(okidx && okjaw && okbox && okart && okother);
end
fprintf('\nproblems: %d\n', bad);
end
