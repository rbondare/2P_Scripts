function export_proxy_inputs(session, out_dir, ball_block)
% Small v7 .mat of everything locomotion_proxy_from_video.py needs from one aggregated session (the aggregated file is v7.3, which the
% python side does not read): camera and imaging time stamps, suite2p rigid shifts and, if ball_block is given, the ball speed of that block.
%   export_proxy_inputs('AnimalRB24_261023_1538', work_dir, 3)   reference session: block 3 (corridor) has the complete ball log
%   export_proxy_inputs('AnimalRB24_261023_1527', work_dir, 1)   target session: block 1 only holds the first readings (truncated log), used as a check
% Writes <out_dir>/<session>_proxyin.mat with tcam, timeca, xoff, yoff, and (ball_block) t_ball, fwd, side, n_read.
agg_dir = 'Z:\Group Members\Rima\Aggregated';
s2p_dir = fullfile('Z:\Group Members\Rima\DATA_2P', session, 'suite2p', 'plane0');
A = load(fullfile(agg_dir, [session '_preprocessed.mat']), 'Triggers', 'TimeCa', 'LocomotionCal');
o = load(fullfile(s2p_dir, 'Fall.mat'), 'ops');
out.tcam = double(A.Triggers.TimeCamera(:));
out.timeca = double(A.TimeCa(:));
out.xoff = double(o.ops.xoff(:));
out.yoff = double(o.ops.yoff(:));
if nargin >= 3 && ~isempty(ball_block)
    L = A.LocomotionCal(ball_block);
    n = min([numel(L.Forward_mm), numel(L.SideRight_mm), numel(L.Triggers_sec)]);
    out.n_read = size(L.Values, 1);
    out.t_ball = double(L.Triggers_sec(1:n));
    out.fwd = double(L.Forward_mm(1:n));
    out.side = double(L.SideRight_mm(1:n));
end
if ~isfolder(out_dir); mkdir(out_dir); end
save(fullfile(out_dir, [session '_proxyin.mat']), '-struct', 'out', '-v7');
fprintf('%s: %d camera triggers, %d imaging frames, %d motion samples%s\n', session, numel(out.tcam), numel(out.timeca), numel(out.xoff), ...
    ternary_txt(isfield(out, 'fwd'), sprintf(', ball block %d: %d readings of %d triggers', ball_block, out.n_read, numel(out.t_ball)), ''));
end

function t = ternary_txt(c, a, b)
if c; t = a; else; t = b; end
end
