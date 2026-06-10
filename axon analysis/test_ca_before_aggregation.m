compute_dff_from_fall( ...
    ('Z:\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\AnimalRB6_250302_1628\suite2p\plane0\Fall.mat')', ...
    'Z:\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\AnimalRB6_250302_1628\dff_plane0.mat');

%%

data1 = load('Z:\joeschgrp\Group Members\Rima\DATA_2P_OTHER\AnimalRB6\AnimalRB6_250302_1628\dff_plane0.mat');

%%
dFF = data1.dFF;
dFF_z = (dFF - mean(dFF, 2)) ./ std(dFF, 0, 2);
dFF_selected = dFF_z(:, :);
figure;
imagesc(dFF_selected)
title("Raw")
colormap('jet'); % Use a colormap like 'hot' or 'jet'
colorbar;  % Add a color scale


%% 

figure;
plot(data1.TimeCa, dFF(5,:), 'LineWidth', 1);
ylabel('dF/F');
ylim([-3 10]);
ax1 = gca; % Save axis handle for linking



%% ── SNR rescue figure ────────────────────────────────────────────────────
% Explores how much signal survives different denoising strategies.
% Uses dFF and TimeCa already in workspace.
% Adjust roi_idx and smooth_s to taste.

fs        = 10;          % Hz
roi_idx   = 1;           % which ROI to inspect in trace panels
smooth_s  = [1 3 5];     % smoothing windows to compare, in seconds

smooth_fr = round(smooth_s * fs);   % convert to frames
t         = data1.TimeCa;
raw       = dFF(roi_idx, :);
mean_all  = mean(dFF, 1);           % average across all ROIs

colors = [0.6 0.6 0.6;   % grey  – raw
          0.2 0.6 1.0;   % blue  – 1 s
          0.1 0.8 0.4;   % green – 3 s
          0.9 0.3 0.1];  % red   – 5 s

figure('Name', 'Axon SNR rescue', 'Position', [100 100 1400 700]);

% ── panel 1: single ROI raw vs smoothed ──────────────────────────────────
subplot(3,2,[1 2]);
hold on;
plot(t, raw, 'Color', colors(1,:), 'LineWidth', 0.8);
for k = 1:numel(smooth_fr)
    plot(t, movmean(raw, smooth_fr(k)), 'Color', colors(k+1,:), 'LineWidth', 1.5);
end
hold off;
ylabel('dF/F');
title(sprintf('ROI %d — raw vs movmean', roi_idx));
legend(['raw', arrayfun(@(s) sprintf('%d s avg', s), smooth_s, 'UniformOutput', false)], ...
       'Location', 'northeast');
xlim([t(1) t(end)]);

% ── panel 2: mean across all ROIs raw vs smoothed ────────────────────────
subplot(3,2,[3 4]);
hold on;
plot(t, mean_all, 'Color', colors(1,:), 'LineWidth', 0.8);
for k = 1:numel(smooth_fr)
    plot(t, movmean(mean_all, smooth_fr(k)), 'Color', colors(k+1,:), 'LineWidth', 1.5);
end
hold off;
ylabel('dF/F');
title(sprintf('Mean across all ROIs (n=%d)', size(dFF,1)));
legend(['raw', arrayfun(@(s) sprintf('%d s avg', s), smooth_s, 'UniformOutput', false)], ...
       'Location', 'northeast');
xlim([t(1) t(end)]);

% ── panel 3: heatmap raw ─────────────────────────────────────────────────
subplot(3,2,5);
imagesc(t, 1:size(dFF,1), dFF);
clim([-2 4]);
colormap('jet'); colorbar;
xlabel('Time (s)'); ylabel('ROI');
title('Raw dF/F — all ROIs');

% ── panel 4: heatmap smoothed with middle window ──────────────────────────
win_mid = smooth_fr(ceil(numel(smooth_fr)/2));
dFF_sm  = movmean(dFF, win_mid, 2);
subplot(3,2,6);
imagesc(t, 1:size(dFF_sm,1), dFF_sm);
clim([-2 4]);
colormap('jet'); colorbar;
xlabel('Time (s)'); ylabel('ROI');
title(sprintf('Smoothed %d s — all ROIs', smooth_s(ceil(numel(smooth_s)/2))));

sgtitle(sprintf('Axon calcium — SNR rescue (fs=%d Hz)', fs));

%%

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
fall = load(fall_path, 'F', 'Fneu', 'stat', 'iscell');

F      = fall.F;      % [nROI x nFrames]
Fneu   = fall.Fneu;
stat = fall.stat;
iscell = fall.iscell;

% keep only cells suite2p classified as ROIs
cell_idx = logical(iscell(:, 1));
F    = F(cell_idx, :);
Fneu = Fneu(cell_idx, :);

% build a TimeCa row vector [1 x nFrames] in seconds
nFrames = size(F, 2);
dt      = 1 / 10;
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