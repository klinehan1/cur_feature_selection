function [num_mixed_nodes, num_measurements] = calc_measures(R, class_A, class_B)

% --------------------------------------------------
% Calculate the number of mixed nodes and measurements in those nodes for
% the given classes.
%
%  Input:
%   R: table of SOM results created in org_SOM_results()
%   class_A, class_B: classes
%
%  Output:
%   num_mixed_nodes: number of mixed (CS or SC - depends on class_A and class_B) nodes
%   num_measurements: number of measurements in mixed nodes
% -----------------------------------------------

% filter R for nodes with class A or B
class_idx = (R.(class_A) > 0) | (R.(class_B) > 0);      
idx = (R.num_classes > 1);

num_mixed_nodes = sum(class_idx & idx);
num_measurements = sum(R.total_observations(class_idx & idx));