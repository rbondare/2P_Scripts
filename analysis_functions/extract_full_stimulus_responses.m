function [responses, frame_ranges, time_ranges, properties] = extract_full_stimulus_responses(Stimuli, dff, stim_type, TimeCa, roi_indices)
    % Extract full dF/F responses for all presentations of a stimulus type
    %
    % Inputs:
    %   Stimuli: stimulus array from preprocessed file
    %   dff: calcium data (n_neurons x n_frames, full dataset)
    %   stim_type: stimulus type to extract (string)
    %   TimeCa: time matrix from preprocessed file (2xn_frames, row 1 = time)
    %   roi_indices: (optional) row indices to subset dff to specific plane (default: use all ROIs)
    %
    % Returns:
    %   responses: cell array, each element is (n_neurons x stim_duration_frames)
    %   frame_ranges: (n_presentations x 2) [frame_start, frame_end]
    %   time_ranges: (n_presentations x 2) [time_start, time_end]
    %   properties: struct array with stimulus metadata

    % If roi_indices provided, subset to those ROIs only
    if nargin >= 5 && ~isempty(roi_indices)
        dff = dff(roi_indices, :);
    end

    n_frames_max = size(dff, 2);

    % Extract time vector from TimeCa (row 1 contains time values)
    if size(TimeCa, 1) < 1 || ~ismatrix(TimeCa)
        error('TimeCa must be at least 1xn_frames (preferably 2xn_frames)');
    end
    if isvector(TimeCa)
        time_vector = TimeCa(:)';
    else
        time_vector = TimeCa(1, :);  % First row contains time values
    end

    if length(time_vector) ~= n_frames_max
        error('TimeCa row 1 length (%d) must match dff columns (%d)', length(time_vector), n_frames_max);
    end

    responses = {};
    frame_ranges = [];
    time_ranges = [];
    properties = {};

    pres_count = 0;

    for i = 1:numel(Stimuli)
        if ~isfield(Stimuli(i), 'type') || ~strcmp(Stimuli(i).type, stim_type)
            continue;
        end

        if ~isfield(Stimuli(i), 'TimeStimulusFrame') || isempty(Stimuli(i).TimeStimulusFrame)
            continue;
        end

        % Extract start time from FIRST element of TimeStimulusFrame
        stim_time_values = Stimuli(i).TimeStimulusFrame;
        time_start = stim_time_values(1);  % Use FIRST value as start

        % Calculate end time: start + (duration x num_trials), but never
        % less than the raw TimeStimulusFrame span. For grating/moving_bar,
        % "trials" means outer REPEAT count and stimulus_trial_t is ONE
        % SUBTRIAL's duration -- stimulus_trial_t*trials silently omits the
        % per-direction multiplier and truncates the block to ~1/n_dir of
        % its real length (confirmed against this dataset: moving_bar's
        % stimulus_trial_t*trials=21s vs a true ~168s span; grating's
        % 10s vs a true ~160s span). TimeStimulusFrame's own span is always
        % the real, complete block regardless of what "trials" means for a
        % given stimulus type, so take whichever estimate is larger.
        if isfield(Stimuli(i), 'stimulus_trial_t') && isfield(Stimuli(i), 'trials')
            stim_total_time = Stimuli(i).stimulus_trial_t * Stimuli(i).trials;
            time_end = max(time_start + stim_total_time, max(stim_time_values(:)));
        else
            % Fallback if metadata missing
            time_end = max(stim_time_values(:));
        end

        % Inclusive bounds: TimeStimulusFrame(1) can land exactly on a real
        % red-sync pulse (RedFrameSynchronized=true), so a strict > would
        % occasionally drop the boundary frame -- see
        % project_aggregation_pipeline_fixes memory.
        frame_idx = find(time_vector >= time_start & time_vector <= time_end);

        % Validate
        if isempty(frame_idx) || length(frame_idx) < 2
            continue;
        end

        pres_count = pres_count + 1;

        % Extract response for all neurons using found frame indices
        st_frame = frame_idx(1);
        en_frame = frame_idx(end);
        response_matrix = dff(:, st_frame:en_frame);
        responses{pres_count, 1} = response_matrix;

        % Store frame indices (for reference)
        frame_ranges = [frame_ranges; st_frame, en_frame];

        % Store actual time ranges
        time_ranges = [time_ranges; time_start, time_end];

        % Store properties
        props = struct();
        props.stimulus_index = i;
        props.frame_start = st_frame;
        props.frame_end = en_frame;
        props.time_start_sec = time_start;
        props.time_end_sec = time_end;
        props.duration_frames = en_frame - st_frame + 1;
        props.duration_sec = time_end - time_start;
        props.n_frames_matched = length(frame_idx);

        % Extract metadata
        if isfield(Stimuli(i), 'stimulus_trial_t')
            props.stimulus_trial_t = Stimuli(i).stimulus_trial_t;
        end
        if isfield(Stimuli(i), 'trials')
            props.trials = Stimuli(i).trials;
        end
        if isfield(Stimuli(i), 'specParams') && isstruct(Stimuli(i).specParams)
            props.spec_params = Stimuli(i).specParams;
        end

        properties{pres_count, 1} = props;
    end
end
