%% Diag_Flash_Reconstruction_Check.m
% Step back to basics: rebuild the sparse_local_global_flashes schedule
% directly from specParams + randomseed (same algorithm as
% reconstructFlashes.m, but kept inline here so every step is visible),
% then check -- on individual real flash events, not just the pooled
% average -- whether the calcium response actually starts at, before, or
% after the reconstructed onset time. Two onset-time conventions are
% compared:
%   "late"  = k*ifi      (current reconstructFlashes.m convention:
%                          Frametime(k) = k*ifi)
%   "early" = (k-1)*ifi  (candidate fix: frame k is the k-th rendered
%                          frame of the trial, occurring at elapsed
%                          (k-1)*ifi after the trial's first frame,
%                          since trialframecount resets to 0 and the
%                          first frame increments it to 1 -- see
%                          present_stimulus_addnew.m lines 1454/1489)
% Scratch/diagnostic script -- not part of the analysis pipeline.

baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat";
B = load(baseline_file);

flash_idx = find(strcmp({B.Stimuli.type}, 'sparse_local_global_flashes'), 1);
s = B.Stimuli(flash_idx);
block_start = s.TimeStimulusFrame(1);

ifi = s.ifi(1);
fprintf('ifi = %.6f s (%.2f Hz)\n', ifi, 1/ifi);
fprintf('specParams.interval = [%g %g], duration = %g\n', s.specParams.interval(1), s.specParams.interval(2), s.specParams.duration);

%% Rebuild schedule exactly as present_stimulus_addnew.m does
rng(s.randomseed);
N = 200 + ceil(s.stimulus_trial_t/ifi);
Object_in_frame = zeros(N, 1);
Frametime_sched = (1:N) * ifi;   % schedule's own internal label, k*ifi
t = 0;
while t <= Frametime_sched(end)
    t = t + rand(1)*diff(s.specParams.interval) + s.specParams.interval(1);
    Object_in_frame(Frametime_sched >= t & Frametime_sched < t + s.specParams.duration) = ...
        randi(size(s.specParams.types, 2), 1);
end
is_on = Object_in_frame > 0;
onset_mask  = diff([0; is_on]) == 1;
offset_mask = diff([is_on; 0]) == -1;
onset_k  = find(onset_mask);   % frame indices (k) where flash turns on
offset_k = find(offset_mask);

onset_late  = onset_k * ifi;          % current convention
onset_early = (onset_k - 1) * ifi;    % candidate fix
offset_late  = offset_k * ifi;
offset_early = (offset_k - 1) * ifi;

keep = onset_late <= s.stimulus_trial_t;
onset_late = onset_late(keep); onset_early = onset_early(keep);
offset_late = offset_late(keep); offset_early = offset_early(keep);
n_flash = numel(onset_late);
fprintf('n flashes reconstructed = %d\n', n_flash);

%% Anchor to real camera-clock time via nearest red pulse
PT = B.Triggers.TimeProjector;
[anchor_dt, anchor_loc] = min(abs(PT - block_start));
anchor = PT(anchor_loc);
fprintf('anchor = %.4f (block_start=%.4f, correction=%.4f)\n', anchor, block_start, anchor - block_start);

abs_onset_late  = anchor + onset_late;
abs_onset_early = anchor + onset_early;
abs_offset_late  = anchor + offset_late;
abs_offset_early = anchor + offset_early;

disp(table((1:n_flash)', onset_early, onset_late, offset_early-onset_early, ...
    'VariableNames', {'flash','onset_early_rel_s','onset_late_rel_s','duration_check_s'}));

%% Pull population mean dF/F (all ROIs, selected plane) for context
plane_idx = 1;
centroids = B.CaData(1).Ca_centroid_voxel;
roi_idx = find(round(centroids(:,3)) == plane_idx);
dff = B.CaData(1).Ca_dFF(roi_idx, :);
pop_mean = mean(dff, 1, 'omitnan');
time_vec = B.TimeCa(1, :);

%% Plot: individual flashes, raw population trace with both onset candidates marked
n_show = min(8, n_flash);
figure('Name', 'Individual flash events: raw population trace vs onset candidates', 'Position', [50 50 1500 900]);
for i = 1:n_show
    subplot(ceil(n_show/2), 2, i); hold on;
    t0 = abs_onset_late(i);
    win = time_vec > t0-2 & time_vec < t0+3;
    plot(time_vec(win) - t0, pop_mean(win), 'k-', 'LineWidth', 1);
    xline(0, 'r-', 'LineWidth', 1.5, 'DisplayName', 'onset (late = k*ifi)');
    xline(abs_onset_early(i)-t0, 'b--', 'LineWidth', 1.5, 'DisplayName', 'onset (early = (k-1)*ifi)');
    xline(abs_offset_late(i)-t0, 'r:', 'DisplayName', 'offset (late)');
    title(sprintf('Flash %d  (abs t=%.2fs)', i, t0), 'FontSize', 9);
    if i==1; legend('Location','best','FontSize',7); end
    xlabel('Time rel. to onset\_late (s)'); ylabel('dF/F (population mean)');
    set(gca,'Box','off');
end

%% Pooled/aligned average using each onset convention (median across all flashes)
pre = 1.5; post = 3; n_com = 100;
t_com = linspace(-pre, post, n_com);
Ca_late  = nan(n_com, n_flash);
Ca_early = nan(n_com, n_flash);
for i = 1:n_flash
    idxL = find(time_vec > abs_onset_late(i)-pre & time_vec < abs_onset_late(i)+post);
    if numel(idxL) >= 2
        Ca_late(:,i) = interp1(time_vec(idxL)-abs_onset_late(i), pop_mean(idxL), t_com, 'linear', 'extrap');
    end
    idxE = find(time_vec > abs_onset_early(i)-pre & time_vec < abs_onset_early(i)+post);
    if numel(idxE) >= 2
        Ca_early(:,i) = interp1(time_vec(idxE)-abs_onset_early(i), pop_mean(idxE), t_com, 'linear', 'extrap');
    end
end

figure('Name', 'Pooled flash-triggered average: onset_late vs onset_early', 'Position', [100 100 900 500]);
hold on;
plot(t_com, median(Ca_late, 2, 'omitnan'), 'r-', 'LineWidth', 2, 'DisplayName', 'aligned to onset\_late (k*ifi)');
plot(t_com, median(Ca_early, 2, 'omitnan'), 'b-', 'LineWidth', 2, 'DisplayName', 'aligned to onset\_early ((k-1)*ifi)');
xline(0, 'k--', 'HandleVisibility','off');
xlabel('Time from onset (s)'); ylabel('dF/F (population median across flashes)');
title('Flash-triggered average (all plane-1 ROIs, baseline) -- onset convention comparison', 'Interpreter','none');
legend('Location','best');
set(gca,'Box','off');
