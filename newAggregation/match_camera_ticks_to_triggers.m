function [t, info] = match_camera_ticks_to_triggers(trig, ticks)
%MATCH_CAMERA_TICKS_TO_TRIGGERS Give every SAVED video frame the 2P time of the trigger it was taken at.
%
%   [t, info] = match_camera_ticks_to_triggers(trig, ticks)
%
%   trig   camera trigger times on the 2P clock (s), one per picture the 2P asked the camera for
%   ticks  camera hardware time stamp of every frame that reached the video (s, any origin; frames.csv column
%          camera_timestamp_ticks is in ns)
%   t      the trigger time of every video frame that has one (t = trig(info.idx), frames info.frame_idx,
%          normally all of them; frames after the last trigger are left out, info.n_untimed_tail)
%
%   Why: save_preprocessed used to pair video frame i with trigger i. That is only right if the camera saved a
%   frame for every trigger. When the camera (or the recording PC) drops frames the video is shorter than the
%   trigger log, every frame after a drop is paired with a trigger that is too early, and the trigger tail is
%   cut off (RB23_261112_1746: 1,341 dropped frames, video up to 26.8 s early; RB24_261112_1408: 11 % dropped).
%
%   How: the camera clock and the 2P clock both run continuously (also across the stops between acquisitions),
%   so a frame's 2P time is its tick time plus an offset that only drifts slowly (tens of ppm). Frames are matched
%   one after the other to the nearest LATER trigger (strictly increasing), the offset is tracked slowly from the
%   residuals, so a dropped frame simply skips a trigger. The first frame is assumed to belong to the first trigger.
%   If that is an integer number of periods wrong it shows at the first trigger after a gap of the trigger train
%   (the 2P was stopped, no trigger inside the gap): there the offset is re-anchored and the whole matching is
%   repeated. A session without any gap cannot reveal such an error (info.n_before_first / n_after_last tell how
%   many triggers are left over at either end). A second blind spot: if the first frame is misaligned AND the first
%   frame after that gap was also dropped, all frames are one trigger (20 ms) off without any warning.
%
%   info fields: ok (no re-anchoring needed in the final pass, <= 0.1 % of the frames further than a quarter
%   period from their trigger), idx, n_frames, n_triggers, n_dropped (triggers between the first and last matched
%   one without a frame), dropped_idx, n_before_first, n_after_last, residual_sd_ms, max_abs_residual_ms,
%   drift_ppm (clock drift tick -> 2P), n_bad, n_reanchor, reanchor (frame, correction_s, trigger), period_s.
%
%   See also: save_preprocessed, verify_aggregation_timing

trig = double(trig(:));
tk = double(ticks(:));
trig = trig(isfinite(trig));
nT = numel(trig);
nF = numel(tk);
if nT < 2 || nF < 2
    error('match_camera_ticks_to_triggers:TooShort', 'need at least 2 triggers and 2 frames');
end
if any(diff(trig) <= 0) || any(diff(tk) <= 0)
    error('match_camera_ticks_to_triggers:NotIncreasing', 'trigger times and tick times must be strictly increasing');
end

dT = diff(trig);
T = median(dT);
first_after_gap = false(nT, 1);
first_after_gap(find(dT > 1.5 * T) + 1) = true;   % trigger j+1 follows a stop of the 2P
tol = 0.25 * T;

off = trig(1) - tk(1);
for pass = 1:2
    [idx, e, reanchor] = match_pass(trig, tk, first_after_gap, T, tol, off);
    if pass == 1 && ~isempty(reanchor)
        off = off + reanchor(1, 2);      % the first frame was not taken at the first trigger: start again
    else
        break
    end
end

% Frames after the last trigger (the camera kept running after the 2P stopped) get no time: they are the tail of
% the video and are left out (info.frame_idx lists the frames that got a time; usually 1:numel(t)).
nv = nnz(idx > 0);
idx = idx(1:nv);
e = e(1:nv);
tkv = tk(1:nv);
t = trig(idx);

info = struct();
info.period_s = T;
info.n_frames = nF;
info.n_triggers = nT;
info.frame_idx = (1:nv)';
info.n_untimed_tail = nF - nv;
info.idx = idx;
info.reanchor = reanchor;
info.n_reanchor = size(reanchor, 1);
info.n_before_first = idx(1) - 1;
info.n_after_last = nT - idx(end);
info.dropped_idx = setdiff((idx(1):idx(end))', idx);
info.n_dropped = numel(info.dropped_idx);
d = t - tkv;                                   % clock offset tick -> 2P at every frame
pf = polyfit(tkv - tkv(1), d, 1);
r = d - polyval(pf, tkv - tkv(1));
info.drift_ppm = 1e6 * pf(1);
info.residual_sd_ms = 1e3 * std(r);
info.max_abs_residual_ms = 1e3 * max(abs(r));
info.n_bad = nnz(abs(e) > tol);
info.ok = info.n_reanchor == 0 && info.n_bad <= max(2, 0.001 * nF) && info.n_untimed_tail <= 0.01 * nF;
end


function [idx, e, reanchor] = match_pass(trig, tk, first_after_gap, T, tol, off)
% One sequential pass. idx(i) = trigger index of frame i (0 = no trigger left), e(i) = trigger minus predicted time.
nF = numel(tk);
nT = numel(trig);
idx = zeros(nF, 1);
e = zeros(nF, 1);
reanchor = zeros(0, 3);
j = 0;
for i = 1:nF
    jn = j + 1;
    if jn > nT
        break
    end
    p = tk(i) + off;
    while jn < nT && abs(trig(jn + 1) - p) < abs(trig(jn) - p)
        jn = jn + 1;
    end
    ei = trig(jn) - p;
    if first_after_gap(jn) && abs(ei) > tol && abs(ei) <= 20 * T
        % integer misalignment: only a trigger right after a gap can show it
        reanchor(end + 1, :) = [i, ei, jn]; %#ok<AGROW>
        off = off + ei;
        ei = 0;
    elseif abs(ei) < tol
        off = off + 0.02 * ei;                % follow the slow clock drift
    end
    idx(i) = jn;
    e(i) = ei;
    j = jn;
end
end
