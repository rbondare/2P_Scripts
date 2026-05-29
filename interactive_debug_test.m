%% INTERACTIVE TEST: Debug organize_stimulus_by_parameters
% Run this directly in MATLAB editor for real-time feedback

clear all; close all;

rec1_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1115_preprocessed.mat';
rec2_file = 'Z:\group\joeschgrp\Group Members\Rima\Aggregated\AnimalRB19_260312_1243_preprocessed.mat';

fprintf('\n========== STEP 1: Load first file only ==========\n');
fprintf('Loading: %s\n', rec1_file);
rec1 = load(rec1_file);
fprintf('  ✓ Loaded successfully\n\n');

fprintf('========== STEP 2: Inspect centroid data ==========\n');
centroid = rec1.CaData(1).Ca_centroid_voxel;
fprintf('Centroid shape: %s\n', mat2str(size(centroid)));
fprintf('Centroid data type: %s\n\n', class(centroid));

fprintf('Centroid first 5 rows:\n');
disp(centroid(1:5, :));

fprintf('========== STEP 3: Extract Z coordinates ==========\n');
centroidZ = centroid(:, 3);
fprintf('Column 3 data type: %s\n', class(centroidZ));
fprintf('Column 3 first 10 values: %s\n\n', mat2str(centroidZ(1:10)'));

fprintf('========== STEP 4: Get unique planes ==========\n');
unique_planes = unique(centroidZ);
unique_planes = sort(unique_planes);
fprintf('Unique planes found: %s\n', mat2str(unique_planes));
fprintf('Number of unique planes: %d\n\n', length(unique_planes));

fprintf('========== STEP 5: Try indexing ==========\n');
try
    plane_1 = unique_planes(1);
    fprintf('Plane at index 1: %s (class: %s)\n', mat2str(plane_1), class(plane_1));
    
    plane_2 = unique_planes(2);
    fprintf('Plane at index 2: %s (class: %s)\n', mat2str(plane_2), class(plane_2));
    fprintf('  ✓ Indexing works\n\n');
catch ME
    fprintf('  ✗ ERROR: %s\n\n', ME.message);
end

fprintf('========== STEP 6: Count neurons per plane ==========\n');
for i = 1:length(unique_planes)
    plane_z = unique_planes(i);
    count = sum(centroidZ == plane_z);
    fprintf('Plane %d (Z=%.1f): %d neurons\n', i, plane_z, count);
end

fprintf('\n========== STEP 7: Extract dFF for plane 1 ==========\n');
dff_full = rec1.CaData(1).Ca_dFF;
centroidZ_num = double(centroidZ);
plane_1_z = unique_planes(1);
roi_indices = find(centroidZ_num == plane_1_z);
dff_plane1 = dff_full(roi_indices, :);
fprintf('dFF for plane 1: %s (expected: %d neurons × time points)\n', ...
    mat2str(size(dff_plane1)), length(roi_indices));

fprintf('\n========== SUCCESS ========== \n');
fprintf('Data structure analysis complete. Check the output above for issues.\n');
