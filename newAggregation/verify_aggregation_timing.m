function R = verify_aggregation_timing(preprocessed_file, varargin)
%VERIFY_AGGREGATION_TIMING Independent timing / bookkeeping checks of an aggregated *_preprocessed.mat
%
%   R = verify_aggregation_timing(file)
%   R = verify_aggregation_timing(file, 'intermediate', '...\<rec>_S2Presult.mat', ...
%                                       'frames_csv', '...\*_frames.csv', ...
%                                       'shifts', [nVol x 1 lateral offset trace on the volume grid])
%
%   Checks (each prints PASS / WARN / FAIL; R.checks holds the same):
%     1  frame accounting   TimeCa columns == calcium trace columns == volumes in the tif headers
%     2  2P time base       frame spacing inside / between volumes; gaps between acquisitions against the tif header clocks
%     3  camera             trigger count vs camera frames.csv (hardware timestamps) vs DLC rows; the camera trigger times of the
%                           aggregated file against the camera's own host-clock times (constant offset = the two clocks agree,
%                           drift / steps = timing error); dropped camera frames
%     4  stimulus           sync status of every stimulus block / trial, stimulus times inside acquisitions
%     5  ball               locomotion time stamps (training sessions: LocomotionCal.SyncInfo)
%     6  end to end         optional: lag of the peak cross-correlation between a 2P motion trace ('shifts') and the video
%                           motion energy, in the aggregated time bases of both (a hardware-coupled pair must peak at ~0 lag)
%
%   Options: 'intermediate' (S2Presult.mat written by save_preprocessed; gives the untouched header clocks),
%            'frames_csv'   (new Basler *_frames.csv: frame, camera_timestamp_ticks, host_unix_seconds, queued),
%            'shifts'       (2P lateral offset per volume, same length as TimeCa columns), 'tz_hours' (default: autodetect).

p = inputParser;
addParameter(p, 'intermediate', '', @ischar);
addParameter(p, 'frames_csv', '', @ischar);
addParameter(p, 'shifts', [], @isnumeric);
addParameter(p, 'tz_hours', [], @isnumeric);
parse(p, varargin{:});
inter = p.Results.intermediate;
csvf = p.Results.frames_csv;
shifts = p.Results.shifts(:)';

R = struct('checks', {{}});
    function add(status, name, msg)
        R.checks(end + 1, :) = {status, name, msg};
        fprintf('  %-4s %-28s %s\n', status, name, msg);
    end

fprintf('=== Timing check of %s ===\n', preprocessed_file);
m = matfile(preprocessed_file);
vars = who('-file', preprocessed_file);
need = {'TimeCa', 'Triggers', 'CaData', 'SI_Info'};
for k = 1:numel(need)
    if ~ismember(need{k}, vars); error('verify_aggregation_timing:Missing', 'variable %s not in %s', need{k}, preprocessed_file); end
end
TimeCa = m.TimeCa;
Triggers = m.Triggers;
CaData = m.CaData;
SI = m.SI_Info;
ca = CaData(find(arrayfun(@(c) ~isempty(c.Ca_F), CaData), 1));
fs_frame = SI.hRoiManager.scanFrameRate;

%% 1 frame accounting
nrow = size(TimeCa, 1);
nvol = size(TimeCa, 2);
nF = size(ca.Ca_F, 2);
fprintf('1. frame accounting: TimeCa %d x %d, %d ROIs x %d samples\n', nrow, nvol, size(ca.Ca_F, 1), nF);
if nF == nvol
    add('PASS', 'trace length', sprintf('calcium traces have %d samples = TimeCa columns', nF));
else
    add('FAIL', 'trace length', sprintf('calcium traces %d samples vs TimeCa %d columns', nF, nvol));
end
H = [];
if ~isempty(inter) && exist(inter, 'file')
    S = load(inter, 'Headers');
    H = S.Headers;
    nfr = arrayfun(@(h) numel(h.frameTimestamps_sec), H);
    per = SI.hStackManager.numSlices + SI.hFastZ.numDiscardFlybackFrames;
    exp_vol = sum(floor(nfr / per));
    if exp_vol == nvol
        add('PASS', 'volumes vs headers', sprintf('%d tifs, %d frames, %d frames/volume -> %d volumes = TimeCa columns', numel(H), sum(nfr), per, exp_vol));
    else
        add('FAIL', 'volumes vs headers', sprintf('headers give %d volumes (%d frames / %d), TimeCa has %d', exp_vol, sum(nfr), per, nvol));
    end
    if any(mod(nfr, per) ~= 0)
        add('WARN', 'partial volumes', sprintf('%d tif(s) end inside a volume (frames %% %d != 0): timestamps after them may be shifted', nnz(mod(nfr, per)), per));
    end
end

%% 2 two-photon time base
fprintf('2. two-photon time base\n');
dwithin = diff(TimeCa, 1, 1);
tf = 1 / fs_frame;
add_status(max(abs(dwithin(:) - tf)) < 0.002, 'frame spacing in volume', ...
    sprintf('mean %.4f s (header %.4f s), max deviation %.2f ms', mean(dwithin(:)), tf, 1e3 * max(abs(dwithin(:) - tf))));
t1 = TimeCa(1, :);
dv = diff(t1);
vol_period = median(dv);
big = find(dv > 1.5 * vol_period);
add_status(all(diff(t1) > 0), 'monotonic', sprintf('TimeCa(1,:) strictly increasing; volume period %.5f s (header %.5f s)', vol_period, 1 / SI.hRoiManager.scanVolumeRate));
fprintf('     %d gaps > 1.5 volumes (acquisition boundaries); durations (s): %s\n', numel(big), mat2str(round(dv(big), 2)));
if ~isempty(H)
    k0 = find(arrayfun(@(h) ~isempty(h.epoch), H), 1);
    ep = arrayfun(@(h) datetime(h.epoch{1}), H(arrayfun(@(h) ~isempty(h.epoch), H)));
    first_ts = arrayfun(@(h) h.frameTimestamps_sec(1), H(arrayfun(@(h) ~isempty(h.epoch), H)));
    nfr = arrayfun(@(h) numel(h.frameTimestamps_sec), H(arrayfun(@(h) ~isempty(h.epoch), H)));
    vols = floor(nfr / (SI.hStackManager.numSlices + SI.hFastZ.numDiscardFlybackFrames));
    start_idx = [1, cumsum(vols(1:end-1)) + 1];
    exp_start = seconds(ep - ep(1)) + first_ts - first_ts(1);
    got_start = t1(start_idx) - t1(1);
    err = got_start - exp_start;
    add_status(max(abs(err)) < 0.01, 'acquisition starts', ...
        sprintf('first volume of each of %d tifs vs header clocks (epoch + frame timestamp): max |error| %.1f ms, median %.1f ms', numel(ep), 1e3 * max(abs(err)), 1e3 * median(abs(err))));
    overl = (t1(start_idx(2:end)) - (t1(start_idx(2:end) - 1) + vol_period));
    fprintf('     gap between consecutive acquisitions (first volume minus expected next): median %.2f s, range %.2f - %.2f s\n', median(overl), min(overl), max(overl));
end

%% 3 camera
fprintf('3. camera\n');
tc = Triggers.TimeCamera(:)';
nbeh = NaN;
if ismember('Behav', vars)
    B = m.Behav;
    if isfield(B, 'FrameNum'); nbeh = numel(B.FrameNum); end
end
fprintf('     camera triggers in the aggregated file: %d; Behav rows: %d\n', numel(tc), nbeh);
if ~isnan(nbeh)
    add_status(numel(tc) == nbeh, 'triggers == Behav rows', sprintf('%d vs %d', numel(tc), nbeh));
end
dtc = diff(tc);
fprintf('     trigger interval: median %.4f s (%.2f Hz), %d intervals > 25 ms (gaps), longest %.2f s\n', median(dtc), 1 / median(dtc), nnz(dtc > 0.025), max(dtc));
host = [];
if ~isempty(csvf) && exist(csvf, 'file')
    T = readtable(csvf);
    host = T.host_unix_seconds(:)';
    ticks = double(T.camera_timestamp_ticks(:))';
    fprintf('     frames.csv: %d frames, host interval median %.4f s, %d > 25 ms, longest %.2f s; queued==1: %d%%\n', ...
        numel(host), median(diff(host)), nnz(diff(host) > 0.025), max(diff(host)), round(100 * mean(T.queued == 1)));
    hw_dt = diff(ticks) * 1e-9;
    fprintf('     hardware ticks (ns): interval median %.5f s, max %.4f s\n', median(hw_dt), max(hw_dt));
    n_tail = 0;       % frames after the last trigger that match_camera_ticks_to_triggers left out
    if isfield(Triggers, 'CameraMatch'); n_tail = Triggers.CameraMatch.n_untimed_tail; end
    tail_note = '';
    if n_tail > 0; tail_note = sprintf(' (+ %d frames after the last trigger, left out)', n_tail); end
    add_status(abs(numel(host) - numel(tc) - n_tail) <= 1, 'frames.csv vs triggers', ...
        sprintf('%d camera frames vs %d triggers%s', numel(host), numel(tc), tail_note));
    if ismember('VideoParams', vars)
        VP = m.VideoParams;
        if isfield(VP, 'selected_frameIdx') && ~isempty(VP.selected_frameIdx)
            sel = VP.selected_frameIdx(:)';
            fprintf('     VideoParams.selected_frameIdx: %d frames, first %d, last %d, contiguous: %d\n', numel(sel), sel(1), sel(end), all(diff(sel) == 1));
        else
            sel = 1:numel(tc);
        end
    else
        sel = 1:numel(tc);
    end
    sel = sel(sel <= numel(host));
    n = min(numel(sel), numel(tc));
    tz = p.Results.tz_hours;
    clock0 = posixtime(datetime(SI.Clockstart, 'TimeZone', 'UTC'));    % header clock read as UTC
    off = host(sel(1:n)) - tc(1:n) - clock0;                          % = -tz*3600 + clock offset
    if isempty(tz); tz = -round(median(off) / 3600); end
    off = off + tz * 3600;
    fprintf('     host clock - (2P trigger clock): time zone assumed UTC%+d h, remaining offset median %.1f ms\n', tz, 1e3 * median(off));
    t_s = tc(1:n);
    pf = polyfit(t_s, off, 1);
    res = off - polyval(pf, t_s);
    host_ok = max(abs(off - median(off))) < 0.05;
    host_msg = sprintf('offset stays within %.1f ms of its median over %.1f min (drift %.2f ms/min, residual SD %.2f ms)', ...
        1e3 * max(abs(off - median(off))), (t_s(end) - t_s(1)) / 60, 1e3 * pf(1) * 60, 1e3 * std(res));
    if ~host_ok && isfield(Triggers, 'CameraMatch') && Triggers.CameraMatch.n_dropped > 0
        % the host time of a frame is taken when it leaves the camera queue and lags after dropped frames
        add('WARN', 'camera clock agreement', [host_msg '; host clock lags after dropped frames, see the hardware ticks check']);
    else
        add_status(host_ok, 'camera clock agreement', host_msg);
    end
    R.camera_offset_s = off;

    % The host clock above is read when a frame leaves the camera queue, so it lags when frames are dropped; the
    % hardware ticks are the camera's own clock: each video frame must sit at its own trigger (frame k <-> element k).
    n2 = min(numel(sel), numel(tc));
    tk2 = (ticks(sel(1:n2)) - ticks(sel(1))) * 1e-9;
    dd = tk2 - (tc(1:n2) - tc(1));
    pf2 = polyfit(tk2, dd, 1);
    r2 = dd - polyval(pf2, tk2);
    msg = sprintf('video frame vs its own trigger: residual SD %.3f ms, max %.2f ms, drift %.1f ppm', ...
        1e3 * std(r2), 1e3 * max(abs(r2)), 1e6 * pf2(1));
    if isfield(Triggers, 'CameraMatch')
        CMm = Triggers.CameraMatch;
        msg = sprintf('%s; %d dropped frames (%.2f %% of the triggers) matched with the camera time stamps', ...
            msg, CMm.n_dropped, 100 * CMm.n_dropped / CMm.n_triggers_raw);
    end
    add_status(std(r2) < 0.005 && max(abs(r2)) < 0.02, 'camera hardware ticks vs triggers', msg);
    R.camera_tick_residual_s = r2;
end

%% 4 stimulus
fprintf('4. stimulus\n');
if ismember('Stimuli', vars)
    St = m.Stimuli;
    for k = 1:numel(St)
        s = St(k);
        nm = '';
        if isfield(s, 'StimulusName'); nm = char(string(s.StimulusName)); end
        sy = NaN;
        if isfield(s, 'RedFrameSynchronized'); sy = s.RedFrameSynchronized; end
        line = sprintf('block %d %s: RedFrameSynchronized=%d', k, nm, sy);
        if isfield(s, 'RedFrameSynchronizedTrial')
            line = sprintf('%s, trials synced %d of %d', line, nnz(s.RedFrameSynchronizedTrial), numel(s.RedFrameSynchronizedTrial));
        end
        fprintf('     %s\n', line);
        if isfield(s, 'TimeStimulusFrame')
            tsf = s.TimeStimulusFrame(:);
            tsf = tsf(isfinite(tsf));
            if ~isempty(tsf)
                inside = tsf >= t1(1) - 1 & tsf <= TimeCa(end) + 1;
                add_status(all(inside), sprintf('block %d inside session', k), sprintf('stimulus frame times %.1f - %.1f s, 2P %.1f - %.1f s', min(tsf), max(tsf), t1(1), TimeCa(end)));
            end
        end
        if isfield(s, 'RedSyncInfo') && ~isempty(s.RedSyncInfo)
            ri = s.RedSyncInfo;
            if isfield(ri, 'n_pulses_lost')
                fprintf('       pulses lost per trial: median %g, max %g\n', median([ri.n_pulses_lost]), max([ri.n_pulses_lost]));
            end
        end
    end
else
    fprintf('     no Stimuli in the file\n');
end

%% 5 ball
fprintf('5. ball / locomotion\n');
if ismember('LocomotionCal', vars)
    LC = m.LocomotionCal;
    fn = fieldnames(LC);
    fprintf('     LocomotionCal fields: %s\n', strjoin(fn', ', '));
    if isfield(LC, 'SyncInfo')
        fprintf('     SyncInfo: %s\n', evalc('disp(LC.SyncInfo)'));
        if isfield(LC.SyncInfo, 'ok')
            if LC.SyncInfo.ok
                add('PASS', 'ball time stamps', sprintf('training_ball_data accepted the trigger/reading alignment (%s)', LC.SyncInfo.method));
            elseif isfield(LC.SyncInfo, 'n_readings') && abs(LC.SyncInfo.n_triggers - LC.SyncInfo.n_readings) <= 2
                add('WARN', 'ball time stamps', sprintf(['pairing not verified by the counter pattern, but %d readings for %d triggers: ' ...
                    'time error <= %d trigger interval(s) (~%.0f ms); check 6 (ball vs 2P) tests it end to end'], ...
                    LC.SyncInfo.n_readings, LC.SyncInfo.n_triggers, abs(LC.SyncInfo.n_triggers - LC.SyncInfo.n_readings), ...
                    1e3 * abs(LC.SyncInfo.n_triggers - LC.SyncInfo.n_readings) / 66));
            else
                add('FAIL', 'ball time stamps', sprintf('%d readings for %d triggers, pairing not verified: locomotion time stamps are unreliable', ...
                    LC.SyncInfo.n_readings, LC.SyncInfo.n_triggers));
            end
            if isfield(LC.SyncInfo, 'n_tail_readings') && LC.SyncInfo.n_tail_readings > 0
                fprintf('     %d readings came from the *_BallCamTail.mat file (aborted session)\n', LC.SyncInfo.n_tail_readings);
            end
        end
    end
else
    fprintf('     no LocomotionCal in the file\n');
end

%% 6 end to end: behaviours that are mechanically coupled must peak at ~0 lag in the aggregated time bases
% Signals are put on the volume grid (mean inside each volume window), smoothed over 1 s and cross-correlated over +-5 s;
% null = the same after circular shifts >= 60 s. A hardware-coupled pair (running -> brain displacement) is expected to
% peak within ~0.6 s; a larger lag is reported as WARN (it can be physiology, but it is where a timing error would show).
sigs = struct('name', {}, 't', {}, 'v', {});
if ~isempty(shifts)
    sh_t = t1(1:min(numel(t1), numel(shifts)));
    sigs(end + 1) = struct('name', '2P lateral offset', 't', sh_t, 'v', shifts(1:numel(sh_t)));
end
if ismember('LocomotionCal', vars)
    LC = m.LocomotionCal;
    if isfield(LC, 'Forward_mm') && isfield(LC, 'Triggers_sec')
        % one LocomotionCal per stimulus block in non-training sessions: put the blocks one after the other
        bt = []; bs = [];
        for kb = 1:numel(LC)
            if isempty(LC(kb).Triggers_sec) || isempty(LC(kb).Forward_mm); continue; end
            tb = LC(kb).Triggers_sec(:)';
            sb = hypot(double(LC(kb).Forward_mm(:)'), double(LC(kb).SideRight_mm(:)'));
            nb = min(numel(tb), numel(sb));
            bt = [bt, tb(1:nb)]; %#ok<AGROW>
            bs = [bs, sb(1:nb)]; %#ok<AGROW>
        end
        [bt, so] = sort(bt); bs = bs(so);
        n = min(numel(bt), numel(bs));
        sigs(end + 1) = struct('name', 'ball speed', 't', bt(1:n), 'v', bs(1:n));
        fprintf('     ball: %d readings with time stamps, covering %.0f - %.0f s (%.0f%% of the 2P session)\n', ...
            nnz(isfinite(bt(1:n))), min(bt(1:n)), max(bt(1:n)), 100 * (max(bt(1:n)) - min(bt(1:n))) / (TimeCa(end) - t1(1)));
    end
end
if ~isnan(nbeh) && isfield(B, 'FacePadME')
    nme = min(numel(tc), numel(B.FacePadME));
    sigs(end + 1) = struct('name', 'video face-pad motion', 't', tc(1:nme), 'v', double(B.FacePadME(1:nme))');
end
fprintf('6. end to end (%d signals)\n', numel(sigs));
edges = [t1, t1(end) + vol_period];
grid_sig = cell(1, numel(sigs));
for s = 1:numel(sigs)
    w = discretize(sigs(s).t, edges);
    ok = ~isnan(w) & isfinite(sigs(s).v);
    g = accumarray(w(ok)', sigs(s).v(ok)', [numel(t1) 1], @mean, NaN)';
    cov = isfinite(g);
    g = movmean(g, round(1 / vol_period), 'omitnan', 'Endpoints', 'shrink');
    grid_sig{s} = g;
    fprintf('     %-22s covers %.0f%% of the volumes\n', sigs(s).name, 100 * mean(cov));
end
maxlag = round(5 / vol_period);
for a = 1:numel(sigs)
    for b = a + 1:numel(sigs)
        x = grid_sig{a}; y = grid_sig{b};
        good = isfinite(x) & isfinite(y);
        if nnz(good) < 200; continue; end
        x = x - mean(x(good)); y = y - mean(y(good));
        x(~good) = 0; y(~good) = 0;
        [cc, lags] = xcorr(y, x, maxlag, 'coeff');          % positive lag: second signal later than the first
        [pk, ip] = max(cc);
        lag_s = lags(ip) * vol_period;
        rng(0);
        null = zeros(100, 1);
        for r = 1:100
            sh = round((60 + rand * (numel(x) * vol_period - 120)) / vol_period);
            null(r) = max(xcorr(y, circshift(x, sh), maxlag, 'coeff'));
        end
        nm = sprintf('%s vs %s', sigs(b).name, sigs(a).name);
        msg = sprintf('peak r %.2f at %+.2f s (%s later when positive), null 99th pct %.2f, %.0f%% of volumes', pk, lag_s, sigs(b).name, prctile(null, 99), 100 * mean(good));
        if pk <= prctile(null, 99)
            add('WARN', nm, ['no peak above the null: ' msg]);
        elseif abs(lag_s) <= 0.6
            add('PASS', nm, msg);
        else
            add('WARN', nm, msg);
        end
        R.xcorr.(matlab.lang.makeValidName(nm)) = struct('lags_s', lags * vol_period, 'cc', cc);
    end
end

fprintf('=== done: %d PASS, %d WARN, %d FAIL ===\n', nnz(strcmp(R.checks(:, 1), 'PASS')), nnz(strcmp(R.checks(:, 1), 'WARN')), nnz(strcmp(R.checks(:, 1), 'FAIL')));

    function add_status(ok, name, msg)
        if ok; add('PASS', name, msg); else; add('FAIL', name, msg); end
    end
end
