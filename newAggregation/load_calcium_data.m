function [Fall, S2Pdirs, multicolor] = load_calcium_data(basedir, varargin)
%LOAD_CALCIUM_DATA Load calcium imaging data from various formats
%
%   [Fall, S2Pdirs, multicolor] = load_calcium_data(basedir)
%   [Fall, S2Pdirs, multicolor] = load_calcium_data(basedir, 'format', 'auto')
%
%   Abstracts calcium data loading to support multiple formats:
%   - suite2p: Traditional Fall.mat and .npy files
%   - caiman: HDF5 format (.h5/.hdf5)
%   - none: Returns empty structures for behavior-only processing
%   - auto: Auto-detect format from files present
%
%   Input:
%       basedir - Base directory of the recording (contains suite2p/ folder)
%
%   Optional Name-Value Arguments:
%       'format' - Calcium data format:
%           'auto'    - Auto-detect from files present (default)
%           'suite2p' - Load from suite2p Fall.mat and .npy files
%           'caiman'  - Load from Caiman HDF5 files
%           'none'    - Return empty structures (for behavior-only)
%
%       'SI_Info' - ScanImage info structure (needed for volume calculations)
%           If not provided, will try to extract from headers
%
%   Output:
%       Fall       - Suite2p-style data structure(s), one per plane
%                    For 'none' format, returns empty []
%       S2Pdirs    - Directory listing of plane folders
%                    For 'none' format, returns empty struct
%       multicolor - Boolean indicating if multicolor imaging
%                    For 'none' format, returns false
%
%   Example:
%       % Auto-detect and load calcium data
%       [Fall, S2Pdirs, multicolor] = load_calcium_data('D:\DATA_2P\recording\');
%
%       % Force behavior-only (no calcium)
%       [Fall, S2Pdirs, multicolor] = load_calcium_data('D:\DATA_2P\recording\', 'format', 'none');
%
%       % Load from Caiman HDF5
%       [Fall, S2Pdirs, multicolor] = load_calcium_data('D:\DATA_2P\recording\', 'format', 'caiman');
%
%   See also: save_preprocessed, get_user_config

% Parse inputs
p = inputParser;
addRequired(p, 'basedir', @ischar);
addParameter(p, 'format', 'auto', @(x) ismember(x, {'auto', 'suite2p', 'caiman', 'none'}));
addParameter(p, 'SI_Info', [], @(x) isempty(x) || isstruct(x));
parse(p, basedir, varargin{:});

format = p.Results.format;
SI_Info = p.Results.SI_Info;

% Ensure basedir ends with filesep
if ~endsWith(basedir, filesep)
    basedir = [basedir filesep];
end

% Auto-detect format if needed
if strcmp(format, 'auto')
    format = detect_calcium_format(basedir);
    fprintf('Auto-detected calcium format: %s\n', format);
end

% Load based on format
switch format
    case 'none'
        [Fall, S2Pdirs, multicolor] = load_empty_calcium();

    case 'suite2p'
        [Fall, S2Pdirs, multicolor] = load_suite2p_calcium(basedir, SI_Info);

    case 'caiman'
        [Fall, S2Pdirs, multicolor] = load_caiman_calcium(basedir);

    otherwise
        error('load_calcium_data:UnknownFormat', 'Unknown format: %s', format);
end

end


function format = detect_calcium_format(basedir)
%DETECT_CALCIUM_FORMAT Auto-detect calcium data format from files present

suite2p_dir = fullfile(basedir, 'suite2p');
caiman_files = [dir(fullfile(basedir, '*.h5')); dir(fullfile(basedir, '*.hdf5'))];

if exist(suite2p_dir, 'dir') && ~isempty(dir(fullfile(suite2p_dir, 'plane*')))
    format = 'suite2p';
elseif ~isempty(caiman_files)
    format = 'caiman';
else
    format = 'none';
    fprintf('No calcium data detected in: %s\n', basedir);
end

end


function [Fall, S2Pdirs, multicolor] = load_empty_calcium()
%LOAD_EMPTY_CALCIUM Return empty structures for behavior-only processing

Fall = [];
S2Pdirs = struct('name', {}, 'folder', {}, 'date', {}, 'bytes', {}, 'isdir', {}, 'datenum', {});
multicolor = false;

fprintf('Using empty calcium data (behavior-only mode)\n');

end


function [Fall, S2Pdirs, multicolor] = load_suite2p_calcium(basedir, SI_Info)
%LOAD_SUITE2P_CALCIUM Load calcium data from suite2p output

suite2p_dir = [basedir 'suite2p' filesep];

% Get plane directories
S2Pdirs = dir([suite2p_dir 'plane*']);

if isempty(S2Pdirs)
    Fall = [];
    multicolor = false;
    return;
end

% Sort plane directories
temp = struct2cell(S2Pdirs);
temp = temp(1,:)';
temp = cellfun(@(x) str2double(x(6:end)), temp);
[~, sort_idx] = sort(temp);
S2Pdirs = S2Pdirs(sort_idx);

% Check for combined and projection data
combined = exist([suite2p_dir 'combined'], 'dir') == 7;
proj = exist([suite2p_dir 'proj'], 'dir') == 7;

% Determine if multicolor from SI_Info
if ~isempty(SI_Info) && isfield(SI_Info, 'hChannels') && isfield(SI_Info.hChannels, 'channelsActive')
    multicolor = numel(SI_Info.hChannels.channelsActive) > 1;
else
    multicolor = false;
end

fprintf('Loading suite2p data from %s\n', suite2p_dir);
fprintf('  Planes: %d, Combined: %d, Projection: %d, Multicolor: %d\n', ...
    numel(S2Pdirs), combined, proj, multicolor);

% Load combined iscell if available
if combined
    iscellComb = readNPY([S2Pdirs(1).folder filesep 'combined' filesep 'iscell.npy']);
    if multicolor && exist([S2Pdirs(1).folder filesep 'combined' filesep 'redcell.npy'], 'file')
        redcellComb = readNPY([S2Pdirs(1).folder filesep 'combined' filesep 'redcell.npy']);
    else
        redcellComb = [];
    end
else
    iscellComb = [];
    redcellComb = [];
end

% Load Fall data
if proj
    % Load from projection folder
    Fall = load([S2Pdirs(1).folder filesep 'proj' filesep 'suite2p' filesep 'plane0' filesep 'Fall.mat']);
else
    % Load from first plane
    Fall = load([S2Pdirs(1).folder filesep S2Pdirs(1).name filesep 'Fall.mat']);

    % Apply combined iscell if available
    if combined
        Fall.iscell = iscellComb(1:size(Fall.iscell, 1), :);
        iscellComb(1:size(Fall.iscell, 1), :) = [];
    end

    % Handle redcell field
    if isfield(Fall, 'redcell') && combined
        Fall = rmfield(Fall, 'redcell');
    end
    if multicolor && combined && ~isempty(redcellComb)
        Fall.redcell = redcellComb(1:size(Fall.iscell, 1), :);
        redcellComb(1:size(Fall.iscell, 1), :) = [];
    end

    % Load additional planes if needed
    if numel(S2Pdirs) > 1 && ~isempty(SI_Info)
        if ~isfield(SI_Info.hFastZ, 'numFramesPerVolume')
            if isfield(SI_Info.hStackManager, 'numFramesPerVolumeWithFlyback')
                SI_Info.hFastZ.numFramesPerVolume = SI_Info.hStackManager.numFramesPerVolumeWithFlyback;
            else
                warning('Plane number not found in SI_Info, loading only first plane');
                return;
            end
        end

        num_planes = SI_Info.hFastZ.numFramesPerVolume - SI_Info.hFastZ.numDiscardFlybackFrames;
        for f = 2:min(numel(S2Pdirs), num_planes)
            temp = load([S2Pdirs(f).folder filesep S2Pdirs(f).name filesep 'Fall.mat']);

            % Handle multicolor field consistency
            if multicolor && ~isfield(temp, 'redcell')
                temp.redcell = [];
            elseif ~multicolor && isfield(temp, 'redcell')
                temp = rmfield(temp, 'redcell');
            end

            Fall(f) = temp;

            % Apply combined iscell
            if combined
                Fall(f).iscell = iscellComb(1:size(Fall(f).iscell, 1), :);
                iscellComb(1:size(Fall(f).iscell, 1), :) = [];
                if multicolor && ~isempty(redcellComb)
                    Fall(f).redcell = redcellComb(1:size(Fall(f).iscell, 1), :);
                    redcellComb(1:size(Fall(f).iscell, 1), :) = [];
                end
            end
        end
    end
end

end


function [Fall, S2Pdirs, multicolor] = load_caiman_calcium(basedir)
%LOAD_CAIMAN_CALCIUM Load calcium data from Caiman HDF5 output
%
%   Converts Caiman HDF5 format to suite2p-compatible structure

% Find HDF5 files
h5_files = [dir(fullfile(basedir, '*.h5')); dir(fullfile(basedir, '*.hdf5'))];

if isempty(h5_files)
    error('load_calcium_data:NoCaimanFiles', 'No Caiman HDF5 files found in: %s', basedir);
end

% Use first HDF5 file found
h5_file = fullfile(h5_files(1).folder, h5_files(1).name);
fprintf('Loading Caiman data from: %s\n', h5_file);

% Create placeholder S2Pdirs
S2Pdirs = struct('name', 'caiman', 'folder', basedir, 'date', '', 'bytes', 0, 'isdir', true, 'datenum', 0);
multicolor = false;

% Read Caiman data and convert to suite2p format
try
    % Standard Caiman HDF5 structure
    % Paths may vary based on Caiman version

    % Try common Caiman output paths
    if h5_path_exists(h5_file, '/estimates/C')
        C = h5read(h5_file, '/estimates/C');           % Denoised temporal components
        S = h5read(h5_file, '/estimates/S');           % Deconvolved spikes
        A = h5read(h5_file, '/estimates/A');           % Spatial components
    elseif h5_path_exists(h5_file, '/C')
        C = h5read(h5_file, '/C');
        S = h5read(h5_file, '/S');
        A = h5read(h5_file, '/A');
    else
        error('load_calcium_data:CaimanStructure', ...
            'Unrecognized Caiman HDF5 structure. Expected /estimates/C or /C');
    end

    % Convert to suite2p-like structure
    num_cells = size(C, 1);
    num_frames = size(C, 2);

    Fall = struct();
    Fall.F = C;                              % Fluorescence traces
    Fall.Fneu = zeros(size(C));              % Neuropil (Caiman handles differently)
    Fall.spks = S;                           % Deconvolved activity
    Fall.iscell = ones(num_cells, 2);        % Mark all as cells [iscell, probability]
    Fall.iscell(:, 2) = 1;                   % Probability = 1

    % Create cell stats (spatial footprints)
    Fall.stat = cell(num_cells, 1);
    for n = 1:num_cells
        Fall.stat{n} = struct();
        % Extract spatial footprint
        spatial = full(A(:, n));
        if ~isempty(spatial)
            [row, col] = ind2sub(sqrt(numel(spatial)) * [1 1], find(spatial > 0));
            Fall.stat{n}.ypix = row;
            Fall.stat{n}.xpix = col;
            Fall.stat{n}.med = [median(col), median(row)];
        else
            Fall.stat{n}.ypix = [];
            Fall.stat{n}.xpix = [];
            Fall.stat{n}.med = [NaN, NaN];
        end
    end

    % Try to read ops-like information
    Fall.ops = struct();
    Fall.ops.nframes = num_frames;
    Fall.ops.nchannels = 1;

    fprintf('Loaded Caiman data: %d cells, %d frames\n', num_cells, num_frames);

catch ME
    error('load_calcium_data:CaimanReadError', ...
        'Error reading Caiman file: %s\n%s', h5_file, ME.message);
end

end


function exists = h5_path_exists(filename, datapath)
%H5_PATH_EXISTS Check if a path exists in HDF5 file
try
    h5info(filename, datapath);
    exists = true;
catch
    exists = false;
end
end
