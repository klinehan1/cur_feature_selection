function [T, seed] = SOM_select(SOM_size_x, SOM_size_y, P, classes_filt)

%--------------------------------------------------------------------------
%  Function to run ten SOMs and choose the "best" 
%  - Results include average quantization error, the number of mixed class
%   neurons, and the total number of measurements in mixed class neurons 
%
%  Input:
%    SOM_size_x, SOM_size_y: SOM size
%    P: data for SOM to cluster
%    classes_filt: vector of mice classes 
%
%  Output: 
%    T: table of results for 10 SOMs
%    seed: seed for "best" SOM
%--------------------------------------------------------------------------

%% Set-up

err = zeros(10,1);
num_mixed = zeros(10,1);
num_mixed_meas = zeros(10,1);
num_nodes = SOM_size_x*SOM_size_y;

%% Run 10 SOMs, organize results + find clusters + calc error for each

for iter=1:10

    % Run SOM  -----------------------
    
    % Using defaults for SOM and training - no info given in paper other than
    % hextop

    rng(iter); % SOM is randomized
    net = selforgmap([SOM_size_x SOM_size_y], 100, 3, 'hextop', 'linkdist');  % other than size - defaults.  
    net = train(net,P');  % observations need to be columns, not rows, P': 77x570
    y = net(P');

    % Organize SOM results ---------------------------

    R = org_SOM_results(classes_filt, y, num_nodes);

    % Find Error for SOM -----------------------------
    
    [err(iter), num_mixed(iter), num_mixed_meas(iter)] = calc_SOM_error(num_nodes, net, P, y, R);

    % find seed with minimum # of mixed nodes.  If a tie use the min of #
    % mixed measurements, then err if another tie.
    
    seed = find(num_mixed == min(num_mixed));
    if length(seed) > 1
        seed = find(num_mixed_meas == min(num_mixed_meas(seed)));

        if length(seed) > 1
            seed = find(err == min(err(seed)));
        end
    end

end

T =  array2table([(1:10)', err, num_mixed, num_mixed_meas], 'VariableNames', ...
    ['Seed', 'Error', "Mixed Class Neurons", "Measures in Mixed Class Neurons"]);
