function [Ca_trials, t_com] = extract_looming_trials(Stimuli, dff, TimeCa, Triggers)
    % Slices and aligns each looming trial in one recording to a common
    % time grid, using the per-trial red sync frame as the true onset.

    Ca_trials = [];
    t_com     = [];

    loom_idx = find(strcmp({Stimuli.type}, 'looming'), 1);
    if isempty(loom_idx)
        return;
    end

    trial_t    = Stimuli(loom_idx).stimulus_trial_t;
    if isfield(Stimuli(loom_idx), 'RedFrameSynchronized') && ~Stimuli(loom_idx).RedFrameSynchronized
        warning('Looming: RedFrameSynchronized=false -- timing fell back to PC clock, expect an offset');
    end
    time_start = Stimuli(loom_idx).TimeStimulusFrame(1);
    time_end   = Stimuli(loom_idx).TimeStimulusFrame(end);

    % Inclusive bounds: TimeStimulusFrame(1) lands exactly on the trial-1
    % pulse, so a strict > would silently drop it.
    PT          = Triggers.TimeProjector;
    loom_onsets = PT(PT >= time_start & PT <= time_end)';
    n_trials    = numel(loom_onsets);

    pre_window_sec  = 2;                        % adjustable
    post_window_sec = trial_t - pre_window_sec; % rest of the trial cycle
    n_com           = 100;
    t_com           = linspace(-pre_window_sec, post_window_sec, n_com);

    n_rois    = size(dff, 1);
    Ca_trials = nan(n_rois, n_com, n_trials);
    time_vec  = TimeCa(1, :);

    for k = 1:n_trials
        idx = find(time_vec > loom_onsets(k) - pre_window_sec & ...
                   time_vec < loom_onsets(k) + post_window_sec);
        if numel(idx) < 2; continue; end
        t_rel = time_vec(idx) - loom_onsets(k);
        Ca_trials(:, :, k) = interp1(t_rel(:), dff(:, idx)', t_com(:), 'linear', 'extrap')';
    end
end
