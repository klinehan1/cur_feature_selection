function [err, num_mixed, num_mixed_meas] = calc_SOM_error(N, net, P, y, R)

% ---------------------------------------------------
% Calculate average quantization error for SOM, the number of mixed class
% neurons, and the number of measurements in mixed class neurons
%
% Input:
%   N: number of neurons in the SOM
%   net: SOM 
%   P: input data matrix to SOM
%   y: SOM results - num_nodes x num_observations, binary data, 
%       y(i,j) = 1 if observation j is in node i. 
%   R: table of SOM results created in org_SOM_results() & cluster_SOM()
%
% Output:
%   err: average quantization error of SOM
%   num_mixed: number of mixed class neurons
%   num_mixed_meas: number of measurements in mixed class neurons
% --------------------------------------------------

    % weight vector of node i is in net.IW{1}(i,:), 49x77
    % 570 observations, node for observation j in y(:,j), data in P(j,:)
    
    node_errs = zeros(N,1);

    % calculate error in each neuron
    for i=1:N
        weight_vec = net.IW{1}(i,:);
        data_obs_idx = find(y(i,:) == 1);
        node_data = P(data_obs_idx,:);
    
        for j=1:length(data_obs_idx)
            node_errs(i) = node_errs(i) + norm(node_data(j,:) - weight_vec);
        end
    end
    
    [num_obs,~] = size(P);
    err = sum(node_errs)/num_obs;
    
    % number of mixed class neurons
    
    num_mixed_idx = (R.("num_classes") > 1);
    num_mixed = sum(num_mixed_idx);
    
    % number of measurements in mixed class neurons
    
    num_mixed_meas = sum(R.("total_observations")(num_mixed_idx));

end