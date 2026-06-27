function [Ca_subtrials, t_com, direction_of] = extract_direction_trials(...
        Stimuli, dff, TimeCa, stype, orientation_field, pre_window_sec, post_window_sec, n_com)
    % Slices every individual subtrial (one per orientation presentation)
    % to a common time grid, anchored on that subtrial's own onset
    % (Frameinfo.subtrial_frame_count==1, same mechanism already
    % validated for looming's per-trial onsets), and labels each subtrial
    % with its reconstructed realized direction.
    Ca_subtrials = []; t_com = []; direction_of = [];

    idx = find(strcmp({Stimuli.type}, stype), 1);
    if isempty(idx); return; end
    s = Stimuli(idx);
    if isfield(s, 'RedFrameSynchronized') && ~s.RedFrameSynchronized
        warning('%s: RedFrameSynchronized=false -- timing fell back to PC clock, expect an offset', stype);
    end

    order = reconstruct_direction_order(s, orientation_field);
    n_dir = numel(order);

    sfc = s.Frameinfo.subtrial_frame_count(~isnan(s.Frameinfo.frame_count));
    subtrial_starts = find(sfc == 1);
    n_subtrials = numel(subtrial_starts);
    if n_subtrials == 0; return; end

    pos_of       = mod((0:n_subtrials-1), n_dir) + 1;  % 1..n_dir, repeating
    direction_of = order(pos_of);
    onsets       = s.TimeStimulusFrame(subtrial_starts);

    t_com    = linspace(-pre_window_sec, post_window_sec, n_com);
    n_rois   = size(dff, 1);
    Ca_subtrials = nan(n_rois, n_com, n_subtrials);
    time_vec = TimeCa(1, :);
    for k = 1:n_subtrials
        fidx = find(time_vec >= onsets(k)-pre_window_sec & time_vec <= onsets(k)+post_window_sec);
        if numel(fidx) < 2; continue; end
        t_rel = time_vec(fidx) - onsets(k);
        Ca_subtrials(:, :, k) = interp1(t_rel(:), dff(:, fidx)', t_com(:), 'linear', 'extrap')';
    end
end
