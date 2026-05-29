%% DEBUG: Inspect data structure from preprocessed files

rec1_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat';
rec2_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat';

fprintf('\n========== LOADING REC1 ==========\n');
rec1 = load(rec1_file);

fprintf('\nRec1 fields:\n');
disp(fieldnames(rec1));

fprintf('\nRec1.CaData(1) fields:\n');
disp(fieldnames(rec1.CaData(1)));

fprintf('\nRec1.CaData(1).Ca_centroid_voxel:\n');
fprintf('  Class: %s\n', class(rec1.CaData(1).Ca_centroid_voxel));
fprintf('  Size: %s\n', mat2str(size(rec1.CaData(1).Ca_centroid_voxel)));
fprintf('  First 5 rows:\n');
disp(rec1.CaData(1).Ca_centroid_voxel(1:min(5, size(rec1.CaData(1).Ca_centroid_voxel,1)), :));

fprintf('\nRec1.CaData(1).Ca_centroid_voxel(:, 3) (Z coordinates):\n');
centroidZ = rec1.CaData(1).Ca_centroid_voxel(:, 3);
fprintf('  Class: %s\n', class(centroidZ));
fprintf('  Unique values: %s\n', mat2str(unique(centroidZ)));

fprintf('\nRec1.CaData(1).Ca_dFF:\n');
fprintf('  Class: %s\n', class(rec1.CaData(1).Ca_dFF));
fprintf('  Size: %s\n', mat2str(size(rec1.CaData(1).Ca_dFF)));

fprintf('\nRec1.Stimuli:\n');
fprintf('  Number of stimuli: %d\n', length(rec1.Stimuli));
fprintf('  Stimulus types: %s\n', strjoin({rec1.Stimuli(:).type}, ', '));

fprintf('\nRec1.Stimuli(1) fields:\n');
stim_fields = fieldnames(rec1.Stimuli(1));
for i = 1:length(stim_fields)
    field = stim_fields{i};
    val = rec1.Stimuli(1).(field);
    if isstruct(val)
        fprintf('  %s: [struct]\n', field);
    elseif ismatrix(val) && (numel(val) > 10 || size(val,1) > 3)
        fprintf('  %s: [%s] (size: %s)\n', field, class(val), mat2str(size(val)));
    else
        fprintf('  %s: %s\n', field, mat2str(val));
    end
end

fprintf('\nRec1.TimeCa:\n');
fprintf('  Class: %s\n', class(rec1.TimeCa));
fprintf('  Size: %s\n', mat2str(size(rec1.TimeCa)));

fprintf('\n========== DONE ==========\n');
