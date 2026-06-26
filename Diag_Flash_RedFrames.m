%% Diag_Flash_RedFrames.m
% Diagnostic: inspect what raw trigger data is actually available, then
% visualize red-sync pulses (kept vs filtered) around the
% sparse_local_global_flashes block, for both baseline and drug
% recordings. Scratch/diagnostic script -- not part of the analysis
% pipeline.

baseline_file = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat";
drug_file     = "Z:\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat";

B = load(baseline_file);
D = load(drug_file);

fprintf('=== Top-level fields in B ===\n');
disp(fieldnames(B));
fprintf('=== Fields in B.Triggers ===\n');
disp(fieldnames(B.Triggers));

%% Block bounds for all stimuli, both recordings
for lbl = {'B','D'}
    S = eval(lbl{1});
    fprintf('\n=== All stimulus blocks (%s) ===\n', lbl{1});
    for i = 1:numel(S.Stimuli)
        s = S.Stimuli(i);
        fprintf('%2d  type=%-30s  red_sync_nth_frame=%-4g  trials=%-4g  start=%.3f  end=%.3f  RedFrameSync=%d\n', ...
            i, s.type, s.red_sync_nth_frame, s.trials, s.TimeStimulusFrame(1), s.TimeStimulusFrame(end), s.RedFrameSynchronized);
    end
end

%% Visualize red pulses (final, filtered TimeProjector) around the flash block
for lbl = {'B','D'}
    S = eval(lbl{1});
    flash_idx = find(strcmp({S.Stimuli.type}, 'sparse_local_global_flashes'), 1);
    s = S.Stimuli(flash_idx);
    block_start = s.TimeStimulusFrame(1);
    block_end   = s.TimeStimulusFrame(end);

    PT = S.Triggers.TimeProjector;
    pad = 15;
    win_mask = PT > block_start-pad & PT < block_end+pad;
    PT_win = PT(win_mask);
    dPT_win = diff(PT_win);

    fprintf('\n=== %s: flash block [%.3f, %.3f], red pulses within +/-%ds ===\n', lbl{1}, block_start, block_end, pad);
    fprintf('n pulses in window = %d\n', numel(PT_win));
    [anchor_dt, anchor_loc] = min(abs(PT - block_start));
    anchor = PT(anchor_loc);
    fprintf('anchor (nearest to block_start) = %.4f  (dt=%.4f, %s)\n', anchor, anchor_dt, ...
        merge(anchor>block_start,'AFTER block_start','BEFORE block_start'));

    figure('Name', sprintf('%s: red pulses near flash block', lbl{1}), 'Position', [100 100 1100 500]);
    subplot(2,1,1); hold on;
    stem(PT_win, ones(size(PT_win)), 'k.', 'LineWidth', 1);
    xline(block_start, 'b-', 'LineWidth', 1.5, 'DisplayName', 'block\_start (TimeStimulusFrame(1))');
    xline(block_end, 'b--', 'LineWidth', 1.5, 'DisplayName', 'block\_end');
    xline(anchor, 'r-', 'LineWidth', 1.5, 'DisplayName', 'chosen anchor pulse');
    legend('Location', 'best');
    ylim([0 1.5]);
    title(sprintf('%s: final (filtered) red pulses, +/-%ds around flash block', lbl{1}, pad), 'Interpreter', 'none');
    xlabel('Time (s, camera clock)');
    set(gca,'Box','off');

    subplot(2,1,2);
    stem(PT_win(2:end), dPT_win, 'k.');
    hold on; yline(1/60, 'r--', '1/60s filter threshold');
    xline(block_start, 'b-'); xline(block_end, 'b--');
    xlabel('Time (s)'); ylabel('\Delta t to previous pulse (s)');
    title('Inter-pulse interval (post-filter) -- anything still under 1/60s would mean filter logic itself is suspect', 'Interpreter','none');
    set(gca,'Box','off');
end

function out = merge(cond, a, b)
    if cond; out = a; else; out = b; end
end
