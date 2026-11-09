function X = export_xs(session, out_dir)
% EXPORT_XS  standardised per-session export of an aggregated axon session for the cross-session analyses (state signal, stimulus responses)
%   export_xs('AnimalRB24_261106_1317')  ->  <out_dir>\<session>_xs.mat  (v7, readable from Python with scipy.io.loadmat)
% Content (everything on the 2P clock of TimeCa(1, :)):
%   tca, dFF_pipe / F / Fneu / centroid of the first populated CaData entry (accepted ROIs, rows = CaData rows)
%   camera behaviour on the camera clock: tcam, PupilArea, PupilBad, FacePadME, WhiskIndex, WhiskAmp, WhiskOn, EyeLidOpening, PupilAzElRaw (NaN vectors when missing)
%   ball readings: loco_t (s, 2P clock), loco_fwd, loco_side (mm per reading), loco_blk (LocomotionCal entry), loco_n_values / loco_n_trig / loco_sync_ok per entry
%   stimulus blocks: blk_* (type, key, span, interrupted) and trial onsets: on_t / on_key / on_blk (key = flash | fff | bar | loom_expand | loom_recede | square)
% No analysis here: it only restructures what the aggregation saved (same onset logic as axon_analysis_multi_V1).
if nargin < 2; out_dir = 'C:\Users\rbondarenko\2P_cache\xs'; end
if ~isfolder(out_dir); mkdir(out_dir); end
repo = 'C:\Users\rbondarenko\projects\2P_Scripts';  addpath(repo, fullfile(repo, 'analysis_functions'));
f = fullfile('Z:\Group Members\Rima\Aggregated', [session '_preprocessed.mat']);
m = whos('-file', f);  inf_ = {m.name};
want = {'CaData', 'TimeCa', 'Stimuli', 'Triggers', 'Behav', 'LocomotionCal'};
D = load(f, want{ismember(want, inf_)});
X = struct('session', session);
ca = D.CaData(find(~cellfun(@isempty, {D.CaData.Ca_dFF}), 1));
X.tca = double(D.TimeCa(1, :));
X.n_planes_timeca = size(D.TimeCa, 1);
X.dFF_pipe = single(ca.Ca_dFF);  X.F = single(ca.Ca_F);  X.Fneu = single(ca.Ca_F_neuropil);
X.centroid = double(ca.Ca_centroid_voxel);
if isfield(D, 'BestDFF'); X.has_bestdff = 1; else; X.has_bestdff = 0; end

%% camera behaviour
X.has_behav = 0;  X.tcam = [];
bf = {'PupilArea', 'PupilBad', 'FacePadME', 'WhiskIndex', 'WhiskAmp', 'WhiskOn', 'EyeLidOpening'};
if isfield(D, 'Behav') && isfield(D, 'Triggers') && isfield(D.Triggers, 'TimeCamera')
    tcam = double(D.Triggers.TimeCamera(:))';  X.tcam = tcam;  X.has_behav = 1;
    for q = 1:numel(bf)
        if isfield(D.Behav, bf{q}) && numel(D.Behav.(bf{q})) == numel(tcam); X.(bf{q}) = double(D.Behav.(bf{q})(:))'; else; X.(bf{q}) = nan(1, numel(tcam)); end
    end
    if isfield(D.Behav, 'PupilAzElRaw') && size(D.Behav.PupilAzElRaw, 1) == numel(tcam); X.PupilAzElRaw = double(D.Behav.PupilAzElRaw)'; else; X.PupilAzElRaw = nan(2, numel(tcam)); end
else
    for q = 1:numel(bf); X.(bf{q}) = []; end
    X.PupilAzElRaw = [];
end

%% ball readings (first readings pair with the first triggers: truncated open-loop logs keep only the first seconds of a block)
lt = [];  lf = [];  ls = [];  lb = [];  nv = [];  nt = [];  ok = [];
if isfield(D, 'LocomotionCal')
    for k = 1:numel(D.LocomotionCal)
        L = D.LocomotionCal(k);
        if ~isfield(L, 'Triggers_sec') || isempty(L.Triggers_sec) || ~isfield(L, 'Forward_mm') || isempty(L.Forward_mm); nv(end+1) = 0; nt(end+1) = 0; ok(end+1) = NaN; continue; end %#ok<AGROW>
        t = double(L.Triggers_sec(:));  fw = double(L.Forward_mm(:));  sd = double(L.SideRight_mm(:));  n = min([numel(t), numel(fw), numel(sd)]);
        lt = [lt; t(1:n)]; lf = [lf; fw(1:n)]; ls = [ls; sd(1:n)]; lb = [lb; k * ones(n, 1)]; %#ok<AGROW>
        nv(end+1) = numel(fw); nt(end+1) = numel(t); %#ok<AGROW>
        if isfield(L, 'SyncInfo') && isstruct(L.SyncInfo) && isfield(L.SyncInfo, 'ok'); ok(end+1) = double(L.SyncInfo.ok); else; ok(end+1) = NaN; end %#ok<AGROW>
    end
end
X.loco_t = lt(:)';  X.loco_fwd = lf(:)';  X.loco_side = ls(:)';  X.loco_blk = lb(:)';  X.loco_n_values = nv;  X.loco_n_trig = nt;  X.loco_sync_ok = ok;

%% stimulus blocks and trial onsets
blk_type = {};  blk_key = {};  blk_span = zeros(0, 2);  blk_int = [];  on_t = [];  on_key = {};  on_blk = [];  blk_loom_dir = [];
if isfield(D, 'Stimuli')
    for k = 1:numel(D.Stimuli)
        s = D.Stimuli(k);  key = 'other';  o = [];  ld = NaN;
        interrupted = isfield(s, 'interrupted') && ~isempty(s.interrupted) && any(s.interrupted(:));
        try
            switch s.type
                case 'sparse_local_global_flashes'; key = 'flash'; o = get_flash_onsets(s);
                case 'full_field_flash';            key = 'fff';   o = first_frames_(s);
                case 'moving_bar';                  key = 'bar';   o = first_frames_(s);
                case 'looming'
                    key = 'loom_expand';  o = first_frames_(s);  ld = 1;
                    if isfield(s, 'specParams') && isstruct(s.specParams) && isfield(s.specParams, 'loom_maxexpand') && ~isempty(s.specParams.loom_maxexpand)
                        ld = sign(s.specParams.loom_maxexpand(1));  if ld < 0; key = 'loom_recede'; end
                    end
                case 'training'
                    key = 'square';  Fi = s.Frameinfo;  n_w1 = size(Fi.frame_count_wait1, 1);  n_st = size(Fi.frame_count_stim, 1);
                    FT = s.TimeStimulusFrame;  FV = s.FrameValid & isfinite(FT);  o = nan(size(FT, 2), 1);
                    for t = 1:size(FT, 2); j = find(FV(n_w1 + (1:n_st), t), 1); if ~isempty(j); o(t) = FT(n_w1 + j, t); end; end
                case 'spontaneous'; key = 'spont';
                otherwise; key = ['other_' s.type];
            end
        catch ME
            warning('%s block %d (%s): onsets failed: %s', session, k, s.type, ME.message);  o = [];
        end
        o = double(o(isfinite(o)));  o = o(:);
        TF = s.TimeStimulusFrame;  span = [min(TF(isfinite(TF))) max(TF(isfinite(TF)))];  if isempty(span); span = [NaN NaN]; end
        blk_type{end+1} = s.type; blk_key{end+1} = key; blk_span(end+1, :) = span; blk_int(end+1) = interrupted; blk_loom_dir(end+1) = ld; %#ok<AGROW>
        on_t = [on_t; o]; on_key = [on_key; repmat({key}, numel(o), 1)]; on_blk = [on_blk; k * ones(numel(o), 1)]; %#ok<AGROW>
    end
end
X.blk_type = blk_type;  X.blk_key = blk_key;  X.blk_span = blk_span;  X.blk_interrupted = blk_int;  X.blk_loom_dir = blk_loom_dir;
X.on_t = on_t(:)';  X.on_key = on_key(:);  X.on_blk = on_blk(:)';
save(fullfile(out_dir, [session '_xs.mat']), '-struct', 'X', '-v7');
fprintf('%s: %d ROIs x %d frames (%.1f min, %.2f Hz) | behaviour %d | ball readings %d | blocks %d, onsets %s\n', session, size(X.dFF_pipe, 1), size(X.dFF_pipe, 2), ...
    (X.tca(end) - X.tca(1)) / 60, 1 / median(diff(X.tca)), X.has_behav, numel(X.loco_t), numel(blk_type), ...
    strjoin(cellfun(@(k) sprintf('%s x%d', k, sum(strcmp(on_key, k))), unique(on_key)', 'UniformOutput', false), ', '));
end

function o = first_frames_(s)
% trial onsets of full-field flash / looming / moving bar blocks: first frame of every trial (V3 section 8)
sfc = s.Frameinfo.subtrial_frame_count(~isnan(s.Frameinfo.frame_count));
o = s.TimeStimulusFrame(find(sfc == 1));
o = o(:);
end
