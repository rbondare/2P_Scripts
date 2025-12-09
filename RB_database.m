%% --- Configuration ---
rootDir  = 'Z:\Group Members\Rima\Stimulus';   % main Stimulus directory
saveFile = fullfile(rootDir, 'RB_Opto_Database.mat');  % database save path

% --- Load existing database if it exists ---
if isfile(saveFile)
    load(saveFile, 'DB');
    fprintf('Loaded existing database (%d entries)\n', height(DB));
else
    % Initialize new table with correct types
    DB = table('Size', [0 11], ...
        'VariableTypes', {'string','string','string','string','string', ...
                          'string','double','double','double','double','double'}, ...
        'VariableNames', {'MouseID','Date','SessionName','FilePath','condition', ...
                          'LocationMode','Level','waitperiod','trials','StimSize','StimulusContrast'});
end

%% --- Loop through each animal folder ---
mice = dir(fullfile(rootDir, 'AnimalRB*')); % all folders starting with RB

for m = 1:numel(mice)
    mouseID = mice(m).name;
    trainPath = fullfile(rootDir, mouseID, 'Training');

    if ~isfolder(trainPath)
        warning('No Training folder found for %s', mouseID);
        continue;
    end

    % Find all training session files
    files = dir(fullfile(trainPath, 'stim_001__training_*.mat'));
    for f = 1:numel(files)
        fname = files(f).name;
        fpath = fullfile(trainPath, fname);

        % --- Extract date and time from file name ---
        tokens = regexp(fname, 'training_(\d{8})_(\d{6})', 'tokens');
        if isempty(tokens)
            warning('Could not parse date/time from %s', fname);
            continue;
        end
        dateStr = tokens{1}{1};
        timeStr = tokens{1}{2};
        sessionName = sprintf('stim_001__training_%s_%s', dateStr, timeStr);

        % --- Skip if this session already exists in the database ---
        if any(strcmp(DB.SessionName, sessionName))
            fprintf('Skipping existing session: %s\n', sessionName);
            continue;
        end

        fprintf('Processing: %s ...\n', fname);

        % --- Load the stimulus file (read-only) ---
        S = load(fpath, '-mat');  

        % --- Determine condition ---
        if isfield(S, 'options') && isfield(S.options, 'opto_trial')
            if any(S.options.opto_trial == 1)
                condition = "opto";
            else
                condition = "control";
            end
        else
            condition = "other";
        end

        % --- Extract stimulus mode ---
        if isfield(S, 'locationMode')
            LocationMode = string(S.locationMode);
        elseif isfield(S, 'options') && isfield(S.options, 'specParams') && ...
               isfield(S.options.specParams, 'locationMode')
            LocationMode = string(S.options.specParams.locationMode);
        else
            LocationMode = "unknown";
        end
        % --- Extract Level ---
        if isfield(S, 'level')
            level = double(S.level);
        else
            level = NaN;
        end

        % --- Extract waitperiod ---
        if isfield(S,'subsessionwait')
            waitperiod = double(S.subsessionwait);
        else
            waitperiod = NaN;
        end

        % --- Extract number of trials ---
        if isfield(S, 'i')
            trials = double(S.i);
        else
            trials = NaN;
        end

        % --- Extract Stimulus size and contrast ---
        if isfield(S, 'options') && isfield(S.options, 'specParams') && isfield(S.options.specParams, 'StimSize')
            stimSize = double(S.options.specParams.StimSize);
        else
            stimSize = NaN;
        end 
       
        if isfield(S, 'options') && isfield(S.options, 'specParams') && isfield(S.options.specParams, 'StimulusContrast')
            stimContrast = double(S.options.specParams.StimulusContrast);
        else
            stimContrast = NaN;
        end

        % --- Add row to table (matching column order exactly) ---
        newRow = {string(mouseID), ...      % MouseID
                  string(dateStr), ...       % Date
                  string(sessionName), ...   % SessionName
                  string(fpath), ...         % FilePath
                  condition, ...             % condition (string)
                  LocationMode, ...          % LocationMode (string)
                  level, ...                 % Level
                  waitperiod, ...            % waitperiod
                  trials, ...                % trials
                  stimSize, ...              % StimSize
                  stimContrast};             % StimulusContrast

        DB = [DB; newRow];
    end
end

%% --- Save database ---
save(saveFile, 'DB');
fprintf('Database saved successfully! (%d total sessions)\n', height(DB));
