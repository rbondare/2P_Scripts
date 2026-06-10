function [dFF, Ca_baseline] = compute_dff_from_fall(fall_path, save_path)
% Load fall.mat from suite2p plane0 and compute windowed percentile dFF.
% fall_path : path to fall.mat (e.g. '...suite2p/plane0/fall.mat')
% save_path : (optional) path to save output .mat (e.g. 'dff_plane0.mat')
%
% Returns:
%   dFF         : [nROI x nFrames] dF/F trace
%   Ca_baseline : [nROI x nFrames] estimated baseline

% --- parameters ---
neuropil_factor = 0.7;
window_base_s   = 60;   % seconds for baseline estimation window (prctile)
window_F_s      = 60;   % seconds for normalization window (movmedian)
prctile_val     = 8;    % percentile for baseline (low percentile ~= F0)
fs              = [];   % will be inferred from ops

% --- load ---
fall = load(fall_path, 'F', 'Fneu', 'ops', 'iscell');

F      = fall.F;      % [nROI x nFrames]
Fneu   = fall.Fneu;
ops    = fall.ops;
iscell = fall.iscell;

% keep only cells suite2p classified as ROIs
cell_idx = logical(iscell(:, 1));
F    = F(cell_idx, :);
Fneu = Fneu(cell_idx, :);

% build a TimeCa row vector [1 x nFrames] in seconds
nFrames = size(F, 2);
dt      = 1 / ops.fs;
TimeCa  = (0 : nFrames-1) * dt;
TimeCa  = repmat(TimeCa, 1, 1);   % keep as row; function expects [1 x T]

[dFF, Ca_baseline] = estimate_continuous_dff( ...
    TimeCa, F, Fneu, neuropil_factor, ...
    [window_base_s, window_F_s], ...
    'windowed', ...
    'prctile', prctile_val);

if nargin > 1 && ~isempty(save_path)
    save(save_path, 'dFF', 'Ca_baseline', 'TimeCa');
end

end
