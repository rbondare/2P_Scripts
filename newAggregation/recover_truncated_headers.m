function Headers = recover_truncated_headers(Headers, tif_files)
%RECOVER_TRUNCATED_HEADERS Rebuild the header of a one-frame TIFF whose ScanImage description was cut off
%
%   Headers = recover_truncated_headers(Headers, tif_files)
%
%   Headers    struct array as built by save_preprocessed (one element per TIFF; ImgInfo, SI, frameNumbers, ...)
%   tif_files  cell array of the full TIFF paths, same order ('' where there is no TIFF)
%
%   Problem: ScanImage stores the per-frame description (frame number, time stamp, aux trigger lists) in the TIFF
%   but cuts it at 2000 characters ('<output truncated>'). That happens to the first frame of a session when many
%   aux triggers were logged before the acquisition started (RB24_261104_1540: ~6 s of projector pulses before the
%   first frame). The ScanImage parser then throws and the header of that TIFF is empty, so
%     - its frame is missing from the frame time stamps while suite2p has a column for it, i.e. all calcium
%       traces would be shifted by one frame against the time stamps,
%     - its camera / ball triggers are missing, i.e. video and ball would be shifted against the 2P clock,
%     - save_preprocessed fails (frameTimestamps_sec([1 end]) of an empty vector).
%
%   What is done for each TIFF that has exactly one image, no parsed header, and is directly followed (consecutive
%   frame numbers, single plane) by an intact header:
%     1. frame number, time stamp, epoch etc. are read from the part of the description that survived
%        (parse_scanimage_frame_description),
%     2. the SI settings (Software tag, identical for every TIFF of a session) are taken from the next intact header,
%     3. camera (auxTrigger1) and ball (auxTrigger2) triggers that were cut off are reconstructed by extending the
%        regular trigger train of the following frames backwards to the start of the frame (window starts 1 ms
%        before the acquisition trigger or frame time stamp). Both trains run on the same clock as the frame, so they
%        are strictly periodic; a train is only extended if its first 25 intervals are regular (CV < 1%),
%     4. every other list that was cut off (projector auxTrigger0, auxTrigger3) is left EMPTY: a partial list would be
%        worse than none.
%   The element gets a field 'reconstructed' (struct: reason, description_chars, lists_cut, lists_empty, n_camera,
%   camera_s, n_ball, ball_s) so that the output file says what was rebuilt. No other element changes, and no field
%   is added at all when nothing needed recovering.
%
%   Check the reconstruction against an independent count: camera triggers against the rows of the camera's
%   frames.csv, ball triggers against the ball log (training_ball_data checks the camera counter pattern).
%
%   See also: parse_scanimage_frame_description, save_preprocessed

expected_scalars = {'frameNumbers', 'acquisitionNumbers', 'frameNumberAcquisition', 'frameTimestamps_sec', ...
    'acqTriggerTimestamps_sec', 'nextFileMarkerTimestamps_sec', 'endOfAcquisition', 'endOfAcquisitionMode', 'dcOverVoltage'};
needed = {'frameNumbers', 'frameTimestamps_sec', 'epoch'};

if ~isfield(Headers, 'SI'); return; end

for h = 1:numel(Headers)
    if ~isempty(Headers(h).SI); continue; end
    if h > numel(tif_files) || isempty(tif_files{h}); continue; end
    if ~isfield(Headers(h).ImgInfo, 'numImages') || Headers(h).ImgInfo.numImages ~= 1; continue; end
    if h == numel(Headers) || isempty(Headers(h+1).SI); continue; end   % needs the intact header that follows
    nxt = h + 1;

    try
        t = Tiff(tif_files{h}, 'r');
        desc = t.getTag('ImageDescription');
        close(t);
    catch
        continue;
    end
    [p, pinfo] = parse_scanimage_frame_description(desc);

    if ~all(isfield(p, needed))
        fprintf('  Header of %s is empty and could not be recovered (frame number / time stamp / epoch missing)\n', tif_files{h});
        continue;
    end
    try
        single_plane = Headers(nxt).SI.hStackManager.numSlices == 1;
    catch
        single_plane = false;
    end
    if ~single_plane || Headers(nxt).frameNumbers(1) ~= p.frameNumbers + 1 || Headers(nxt).frameTimestamps_sec(1) <= p.frameTimestamps_sec
        fprintf('  Header of %s is empty and the next TIFF does not follow it as a single-plane recording, not recovered\n', tif_files{h});
        continue;
    end

    % header in the same field order as a parsed one
    hdr = struct();
    hdr.SI = Headers(nxt).SI;
    for k = 1:numel(expected_scalars)
        if isfield(p, expected_scalars{k})
            hdr.(expected_scalars{k}) = p.(expected_scalars{k});
        else
            hdr.(expected_scalars{k}) = [];
        end
    end
    hdr.epoch = {p.epoch};

    lists_cut = {};
    lists_empty = {};
    for a = 0:3
        key = sprintf('auxTrigger%d', a);
        if isfield(p, key)
            hdr.(key) = {p.(key)};                 % survived: use as it is
        else
            hdr.(key) = {[]};
            lists_cut{end + 1} = key; %#ok<AGROW>
        end
    end
    hdr.I2CData = {[]};

    % camera / ball: extend the regular train of the following frames backwards
    t_start = p.frameTimestamps_sec;
    if isfield(p, 'acqTriggerTimestamps_sec'); t_start = min(t_start, p.acqTriggerTimestamps_sec); end
    t_start = t_start - 1e-3;
    n_rec = struct('auxTrigger1', 0, 'auxTrigger2', 0);
    t_rec = struct('auxTrigger1', [], 'auxTrigger2', []);
    for a = [1 2]
        key = sprintf('auxTrigger%d', a);
        if ~ismember(key, lists_cut); continue; end
        c = cellfun(@(x) reshape(x, 1, []), Headers(nxt).(key)(1:min(60, end)), 'UniformOutput', false);
        nxt_t = cat(2, c{:});
        nxt_t = nxt_t(1:min(25, numel(nxt_t)));
        if numel(nxt_t) < 6 || std(diff(nxt_t)) / mean(diff(nxt_t)) >= 0.01
            fprintf('  %s of %s: the trigger train that follows is missing or not regular, nothing reconstructed\n', key, tif_files{h});
            continue;
        end
        period = median(diff(nxt_t));
        tt = sort(nxt_t(1) - (1:ceil((nxt_t(1) - t_start) / period)) * period);
        tt = tt(tt >= t_start);
        hdr.(key) = {tt};
        n_rec.(key) = numel(tt);
        t_rec.(key) = tt;
    end
    for a = 0:3
        key = sprintf('auxTrigger%d', a);
        if ismember(key, lists_cut) && isempty(hdr.(key){1})
            lists_empty{end + 1} = key; %#ok<AGROW>
        end
    end

    fn = fieldnames(hdr);
    for k = 1:numel(fn)
        Headers(h).(fn{k}) = hdr.(fn{k});
    end
    Headers(h).reconstructed = struct( ...
        'reason', 'ScanImage cut the frame description at 2000 characters; scalars from the surviving text, SI from the next TIFF, cut camera/ball triggers extrapolated from the regular train that follows', ...
        'description_chars', numel(desc), 'truncated_marker', pinfo.truncated, ...
        'lists_cut', {lists_cut}, 'lists_empty', {lists_empty}, ...
        'n_camera', n_rec.auxTrigger1, 'camera_s', t_rec.auxTrigger1, ...
        'n_ball', n_rec.auxTrigger2, 'ball_s', t_rec.auxTrigger2);
    fprintf(['  Header of %s was cut off by ScanImage (%d characters): frame %d at %.4f s recovered, %d camera and %d ball ' ...
        'triggers reconstructed, lists left empty: %s\n'], tif_files{h}, numel(desc), p.frameNumbers, p.frameTimestamps_sec, ...
        n_rec.auxTrigger1, n_rec.auxTrigger2, strjoin(lists_empty, ', '));
end
end
