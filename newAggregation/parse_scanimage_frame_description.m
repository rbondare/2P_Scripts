function [h, info] = parse_scanimage_frame_description(desc)
%PARSE_SCANIMAGE_FRAME_DESCRIPTION Parse the per-frame ImageDescription of a ScanImage TIFF, tolerating truncation
%
%   [h, info] = parse_scanimage_frame_description(desc)
%
%   desc is the ImageDescription tag of one frame ("key = value" lines: frameNumbers, frameTimestamps_sec,
%   epoch, auxTrigger0..3, ...). ScanImage cuts that text off at 2000 characters and appends
%   '<output truncated>' when the aux trigger lists of the frame are long. This happens for the first frame of a
%   session when triggers (e.g. the projector) were logged for several seconds before the acquisition started. The
%   ScanImage library then throws while parsing and fast_get_scanimage_header_notMROI returns an empty header.
%
%   h      struct with every field of the description that is COMPLETE: numeric scalars as double, 'epoch' as a
%          1x6 vector, aux trigger lists as row vectors ([] for an empty list). A list or value that is cut off, and
%          everything after the cut, is left out.
%   info   struct: truncated (marker found), complete (names in h), incomplete (names that are cut off or follow the cut)
%
%   See also: recover_truncated_headers

h = struct();
info = struct('truncated', contains(desc, '<output truncated>'), 'complete', {{}}, 'incomplete', {{}});

lines = regexp(char(desc), '\r?\n', 'split');
cut = false;                                   % once a line is cut off, nothing after it can be trusted
for k = 1:numel(lines)
    tok = regexp(lines{k}, '^\s*(\w+)\s*=\s*(.*?)\s*$', 'tokens', 'once');
    if isempty(tok); continue; end
    key = tok{1};
    val = strrep(tok{2}, '<output truncated>', '');
    if cut || (info.truncated && k == numel(lines))
        info.incomplete{end + 1} = key; %#ok<AGROW>
        cut = true;
        continue;
    end
    if startsWith(val, '[')
        if ~endsWith(val, ']')                 % list cut off
            info.incomplete{end + 1} = key; %#ok<AGROW>
            cut = true;
            continue;
        end
        nums = sscanf(regexprep(val, '[\[\],]', ' '), '%f')';
        if isempty(nums); nums = []; end
        h.(key) = nums;
    elseif startsWith(val, '{')
        if ~endsWith(val, '}')
            info.incomplete{end + 1} = key; %#ok<AGROW>
            cut = true;
            continue;
        end
        h.(key) = {};                           % I2C data: not needed, only its presence
    else
        num = sscanf(val, '%f');
        if isscalar(num) && ~isempty(regexp(val, '^[-+0-9.eE]+$', 'once'))
            h.(key) = num;
        else
            info.incomplete{end + 1} = key; %#ok<AGROW>
            cut = true;
            continue;
        end
    end
    info.complete{end + 1} = key; %#ok<AGROW>
end
end
