%% Diag_RawRedFrames_Check.m
% Inspect RAW auxTrigger0 (projector/red-sync) pulses straight from the
% cached ScanImage *_header.mat files, BEFORE any of save_preprocessed.m's
% process_triggers() filtering (the unique-frame-number dedup and the
% 1/60s minimum-gap filter). Goal: see whether pulses look like they're
% being dropped, and understand the multi-file session structure, before
% trusting anything process_triggers() produces.
%
% Does not touch the multi-GB raw .tif files -- the *_header.mat caches
% already contain everything ScanImage extracted (auxTrigger0 included).

basedir = 'Z:\joeschgrp\Group Members\Rima\DATA_2P\AnimalRB19_260312_1115\';
files = dir(fullfile(basedir, '*_header.mat'));
[~, sidx] = sort({files.name});
files = files(sidx);
nF = numel(files);

survey = struct('name', {}, 'nFrames', {}, 'FT_start', {}, 'FT_end', {}, ...
    'epoch', {}, 'nPulses', {}, 'pulses', {});

fprintf('%-40s %10s %14s %14s %10s\n', 'File', 'nFrames', 'FT_start(s)', 'FT_end(s)', 'nPulses');
for f = 1:nF
    S = load(fullfile(basedir, files(f).name), 'header');
    h = S.header;

    pulses = cat(2, h.auxTrigger0{:});
    pulses = sort(pulses(:)');

    FT = h.frameTimestamps_sec;

    survey(f).name     = files(f).name;
    survey(f).nFrames  = numel(h.frameNumbers);
    survey(f).FT_start = FT(1);
    survey(f).FT_end   = FT(end);
    survey(f).epoch    = h.epoch{1};
    survey(f).nPulses  = numel(pulses);
    survey(f).pulses   = pulses;

    fprintf('%-40s %10d %14.3f %14.3f %10d\n', files(f).name, survey(f).nFrames, FT(1), FT(end), numel(pulses));
end

%% Check exactly the same monotonicity condition process_frame_timestamps() checks
FT_all = arrayfun(@(s) [s.FT_start, s.FT_end], survey, 'UniformOutput', false);
fprintf('\n--- Cross-file frameTimestamps_sec continuity (file N end vs file N+1 start) ---\n');
nonmonotonic_found = false;
for f = 1:nF-1
    gap = survey(f+1).FT_start - survey(f).FT_end;
    flag = '';
    if gap < 0
        flag = '  <-- NON-MONOTONIC (reset/new acquisition)';
        nonmonotonic_found = true;
    end
    fprintf('  %s -> %s : gap = %+.3f s%s\n', survey(f).name, survey(f+1).name, gap, flag);
end
fprintf('multi-acquisition stitching required (per process_frame_timestamps logic): %d\n', nonmonotonic_found);

%% Compute epoch-based dt exactly as process_frame_timestamps() does, for stitching the raw pulses
epoch_dt = arrayfun(@(s) datetime(s.epoch), survey);
dt = seconds(epoch_dt - epoch_dt(1));
fprintf('\n--- epoch-based dt per file (s, relative to file 1 epoch) ---\n');
for f = 1:nF
    fprintf('  %s : epoch=%s  dt=%.3f s\n', survey(f).name, datestr(epoch_dt(f)), dt(f));
end

%% Stitch raw pulses across files using epoch-based dt, keep track of source file
PT_raw = [];
file_id = [];
for f = 1:nF
    PT_raw = [PT_raw, survey(f).pulses + dt(f)]; %#ok<AGROW>
    file_id = [file_id, repmat(f, 1, numel(survey(f).pulses))]; %#ok<AGROW>
end
[PT_raw, ord] = sort(PT_raw);
file_id = file_id(ord);

% Apply the SAME 1/60s dedup filter process_triggers() applies, for comparison
dPT = [inf, diff(PT_raw)];
keep = dPT >= 1/60;
PT_filtered = PT_raw(keep);

fprintf('\nTotal raw pulses (all files, stitched): %d\n', numel(PT_raw));
fprintf('After 1/60s dedup filter: %d (removed %d, %.1f%%)\n', ...
    numel(PT_filtered), numel(PT_raw) - numel(PT_filtered), ...
    100 * (numel(PT_raw) - numel(PT_filtered)) / numel(PT_raw));

save(fullfile(basedir, 'redframe_raw_survey.mat'), 'survey', 'dt', 'PT_raw', 'file_id', 'PT_filtered');

%% Figure 1: per-file structural overview
short_labels = arrayfun(@(f) sprintf('f%02d', f), 1:nF, 'UniformOutput', false);
figure('Name', 'Per-file structure', 'Position', [100 100 1000 450]);
subplot(2,1,1);
bar([survey.nFrames]);
set(gca, 'XTick', 1:nF, 'XTickLabel', short_labels);
ylabel('# imaging frames'); title('Frames per file', 'Interpreter', 'none');
subplot(2,1,2);
bar([survey.nPulses]);
set(gca, 'XTick', 1:nF, 'XTickLabel', short_labels);
ylabel('# raw auxTrigger0 pulses'); title('Raw red-sync pulses per file (before any filtering)', 'Interpreter', 'none');
saveas(gcf, fullfile(basedir, 'diag_fig1_perfile_overview.png'));

%% Figure 2: raster of stitched raw pulses across the whole session, colored by source file
figure('Name', 'Raw pulse raster', 'Position', [100 100 1200 400]);
gscatter(PT_raw, file_id, file_id, [], '.', 6);
xlabel('Stitched session time (s)'); ylabel('Source file index');
title('Raw auxTrigger0 pulses across session (pre-filter)', 'Interpreter', 'none');
legend('off');
saveas(gcf, fullfile(basedir, 'diag_fig2_raw_raster.png'));

%% Figure 3: inter-pulse-interval histograms, raw vs after 1/60s filter
% Self-calibrate bin range from the actual data (log-spaced) instead of guessing.
d_raw = diff(PT_raw);
d_raw = d_raw(d_raw > 0);
d_filt = diff(PT_filtered);
d_filt = d_filt(d_filt > 0);

fprintf('\n--- Raw inter-pulse interval stats (s) ---\n');
fprintf('  min=%.5f  1st pct=%.5f  median=%.5f  99th pct=%.5f  max=%.5f\n', ...
    min(d_raw), prctile(d_raw,1), median(d_raw), prctile(d_raw,99), max(d_raw));
fprintf('--- Filtered (post 1/60s) inter-pulse interval stats (s) ---\n');
fprintf('  min=%.5f  1st pct=%.5f  median=%.5f  99th pct=%.5f  max=%.5f\n', ...
    min(d_filt), prctile(d_filt,1), median(d_filt), prctile(d_filt,99), max(d_filt));

edges = logspace(log10(min(d_raw)), log10(max(d_raw)), 150);
figure('Name', 'Inter-pulse intervals', 'Position', [100 100 1100 420]);
subplot(1,2,1);
histogram(d_raw, edges);
set(gca, 'YScale', 'log', 'XScale', 'log');
xline(1/60, 'r--', '1/60s filter cutoff', 'LineWidth', 1.2);
xlabel('Inter-pulse interval (s, log scale)'); ylabel('Count (log scale)');
title('RAW (before 1/60s filter)');
subplot(1,2,2);
histogram(d_filt, edges);
set(gca, 'YScale', 'log', 'XScale', 'log');
xline(1/60, 'r--', '1/60s filter cutoff', 'LineWidth', 1.2);
xlabel('Inter-pulse interval (s, log scale)'); ylabel('Count (log scale)');
title('AFTER 1/60s dedup filter');
saveas(gcf, fullfile(basedir, 'diag_fig3_intervals.png'));

%% Figure 4: zoomed single-block raster -- actual pulse-by-pulse cadence, raw vs filtered
% Pick the first dense block (lowest file index with > 100 raw pulses) so we can
% see real tick-by-tick spacing instead of a solid-looking line.
dense_files = find([survey.nPulses] > 100);
fb = dense_files(1);
block_mask_raw  = file_id == fb;
block_t_raw     = PT_raw(block_mask_raw);

% recompute the filtered subset for this block directly from block_t_raw
block_dPT = [inf, diff(block_t_raw)];
block_t_filt = block_t_raw(block_dPT >= 1/60);

win = block_t_raw(1) + [0 5]; % first 5 seconds of this block
in_win_raw  = block_t_raw  >= win(1) & block_t_raw  <= win(2);
in_win_filt = block_t_filt >= win(1) & block_t_filt <= win(2);

figure('Name', 'Zoomed single-block cadence', 'Position', [100 100 1100 420]);
subplot(2,1,1);
stem(block_t_raw(in_win_raw), ones(sum(in_win_raw),1), 'Marker', '.', 'BaseValue', 0);
xlim(win); ylim([0 1.2]);
xlabel('Time (s)'); title(sprintf('%s -- RAW pulses, first 5s of block (file %s)', short_labels{fb}, survey(fb).name), 'Interpreter', 'none');
subplot(2,1,2);
stem(block_t_filt(in_win_filt), ones(sum(in_win_filt),1), 'Marker', '.', 'BaseValue', 0, 'Color', 'r');
xlim(win); ylim([0 1.2]);
xlabel('Time (s)'); title('Same window, AFTER 1/60s dedup filter');
saveas(gcf, fullfile(basedir, 'diag_fig4_zoomed_block.png'));

fprintf('\nFigures saved to %s\n', basedir);
fprintf('  diag_fig1_perfile_overview.png\n  diag_fig2_raw_raster.png\n  diag_fig3_intervals.png\n');
