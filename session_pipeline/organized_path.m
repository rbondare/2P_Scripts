function [p, sid] = organized_path(preprocessed_file, plane_idx, out_dir)
% ORGANIZED_PATH  Canonical cache path for one session x plane.
%
%   p = organized_path(preprocessed_file, plane_idx)
%   p = organized_path(preprocessed_file, plane_idx, out_dir)
%
% plane_idx is 1-BASED (as cfg.selected_plane_idx), but the filename carries
% the 0-BASED suite2p plane number so it reads the same way the ROI-match
% files do:  cfg.selected_plane_idx = 3  ->  ..._plane2_organized.mat
% matching  roiMatch_plane2_0320_1444_1609_V1.mat.
%
% Single source of truth for the naming convention -- run_organize writes it,
% compare_group resolves it.

    if nargin < 3 || isempty(out_dir)
        out_dir = fullfile(fileparts(mfilename('fullpath')), 'organized_data');
    end
    [~, sid] = fileparts(char(preprocessed_file));
    sid = regexprep(sid, '_preprocessed$', '');
    p = fullfile(out_dir, sprintf('%s_plane%d_organized.mat', sid, plane_idx - 1));
end
