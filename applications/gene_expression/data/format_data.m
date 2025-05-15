% put genetics data in mat file

addpath('./raw');

opts = detectImportOptions('GSE10072_series_matrix.csv');
vars = opts.SelectedVariableNames;
opts.SelectedVariableNames = vars(2:end);
G = readmatrix('GSE10072_series_matrix.csv',opts);
save('gene.mat', 'G');

opts.SelectedVariableNames = vars(1);
probes = readtable('GSE10072_series_matrix.csv',opts);
save('probes.mat', 'probes');

opts = detectImportOptions('tumor_idx.csv');
tumor_idx = readmatrix('tumor_idx.csv',opts);
save('tumor_idx.mat', 'tumor_idx'); 