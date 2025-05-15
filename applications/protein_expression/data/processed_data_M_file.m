% Create MAT file for processed dataset

%% read in processed data 

opts = detectImportOptions('Preprocessed_Data.csv');
disp([opts.VariableNames' opts.VariableTypes'])
opts = setvartype(opts,{'Mouse_ID'},'char');

pro = readtable('Preprocessed_Data.csv',opts);

M = pro{:,3:79};  % 1080x77
proteins = pro.Properties.VariableNames(3:79);
classes = pro{:,83};
mouse_ids = pro{:,1};
trial_num = pro{:,2};

save('proteins_processed.mat', 'M', 'proteins', 'classes', 'mouse_ids', 'trial_num');