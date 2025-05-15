% ------------------------------------
% Control Experiment Driver
% SOM Control Experiment - extension of that in PLOS One paper by Higuera, 
% Gardiner, and Cios.
% ----------------------------------

%% add paths

addpath('../../ls_cur', '../../sf_cur', '../../deim_cur', '../../qr_cur', ...
    './data', './helper_fcns', './figures');

%% Set-up Variables ---------------------

% set SOM size and seed
SOM_size_x = 7;
SOM_size_y = 7;
num_nodes = SOM_size_x*SOM_size_y;

% filter data for specific experiment: 'c' - control group
subgroup = 'c';

% Control experiment comparisons 
class_A = ["c-CS-s", "c-CS-m", "c-CS-m", "c-CS-s", "c-SC-m", "c-CS-m"];
class_B = ["c-SC-s", "c-SC-m", "c-SC-s", "c-SC-m", "c-SC-s", "c-CS-s"];

% Read in processed data + filter for specific experiment -----------

load proteins_processed.mat

temp = cat(1,classes{1:end});  % 1080 x 6 character array
subgroup_idx = (temp(:,1) == subgroup);

P = M(subgroup_idx, :);
classes_filt = classes(subgroup_idx);
mouse_ids_filt = mouse_ids(subgroup_idx);

%% -------------- All 77 proteins experiment -----------------

% Select an SOM seed
[~,seed] = SOM_select(SOM_size_x, SOM_size_y, P, classes_filt);

% Using defaults for SOM and training - no info given in paper other than hextop
rng(seed); % SOM is randomized - use seed found in SOM_select()
net = selforgmap([SOM_size_x SOM_size_y], 100, 3, 'hextop', 'linkdist');  % other than size - defaults.  
net = train(net,P');  % observations need to be columns, not rows, P': 77x570
y = net(P'); % get SOM node for each observation
 
% results
R = org_SOM_results(classes_filt, y, num_nodes);
[R_new] = cluster_SOM(hextop([SOM_size_x SOM_size_y]), mouse_ids_filt, y, R, num_nodes);
figure(1);
visualize_SOM(R_new, SOM_size_x, SOM_size_y)
print('./figures/all_proteins.eps', '-depsc')
[num_mixed_CS_nodes, num_mixed_CS_meas] = calc_measures(R,"c-CS-m","c-CS-s")
[num_mixed_SC_nodes, num_mixed_SC_meas] = calc_measures(R,"c-SC-m","c-SC-s")

% ------------ Feature Selection Experiments -----------------------

%% Wilcoxon rank sum test
disp('----- Wilcoxon rank sum test -----')
control_exp(true, '', "./figures/wilcoxon/", SOM_size_x, SOM_size_y, ...
    num_nodes, classes_filt, mouse_ids_filt, P, net, R_new, class_A, class_B);

%% SF CUR 
disp('----- SF CUR -----')
control_exp(false, 'sf', "./figures/sf_cur/", SOM_size_x, SOM_size_y, ...
    num_nodes, classes_filt, mouse_ids_filt, P, net, R_new, class_A, class_B);

%% LS-D CUR 
disp('----- LS-D CUR -----')
control_exp(false, 'ls-d', "./figures/ls-d_cur/", SOM_size_x, SOM_size_y, ...
    num_nodes, classes_filt, mouse_ids_filt, P, net, R_new, class_A, class_B);

%% DEIM CUR 
disp('----- DEIM CUR -----')
control_exp(false, 'deim', "./figures/deim_cur/", SOM_size_x, SOM_size_y, ...
    num_nodes, classes_filt, mouse_ids_filt, P, net, R_new, class_A, class_B);

%% QR CUR 
disp('----- QR CUR -----')
control_exp(false, 'qr', "./figures/qr_cur/", SOM_size_x, SOM_size_y, ...
    num_nodes, classes_filt, mouse_ids_filt, P, net, R_new, class_A, class_B);

