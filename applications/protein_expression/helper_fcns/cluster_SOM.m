function [R_new] = cluster_SOM(pos, mouse_ids, y, R, num_nodes)

%--------------------------------------------------------------------------
%  Find "class-specific" clusters in SOM 
%
%  Input:
%    pos: SOM topology 
%    mouse_ids: vector of mouse ids
%    y: SOM results - num_nodes x num_observations, binary data, 
%       y(i,j) = 1 if observation j is in node i. 
%    R: table of SOM results created in org_SOM_results()  
%    num_nodes = number of neurons in SOM
%
%  Output: R_new, table of results - R appended with columns for clustering 
%  results and visualization, i.e., cluster, edge color, edge width   
%--------------------------------------------------------------------------

    R.("cluster") = strings(num_nodes,1);
    R.("edge_color") = strings(num_nodes,1);
    R.("edge_color")(:) = "black";
    R.("edge_width") = 0.5*ones(num_nodes,1);
    
    d = linkdist(pos);
    
    % find clusters in SOM
    
    for i = 1:num_nodes
    
        % only one class in node
        if R.("num_classes")(i) == 1  
    
            c = R.("majority_class")(i);
    
            % look at neighbors
            neighbors = find(d(:,i) == 1);
            
            adj_nbhd = (R.("num_classes")(neighbors) == 1 & ...
                R.("majority_class")(neighbors) == c);
    
            if sum(adj_nbhd) > 0
                R.("cluster")(i) = c;
            else
                % singleton nodes
                mouse_idx = find(y(i,:));
                mice = mouse_ids(mouse_idx);
                T = tabulate(mice);
                
                if sum([T{:,2}] >= 12) > 0   % 12 is the threshold for clustering in the paper
                    R.("cluster")(i) = c;
                end
            end
        end 
    end
    
    R.("edge_color")(R.cluster == "c-CS-s") = "#009e73";
    R.("edge_color")(R.cluster == "c-CS-m") = "#C70039";
    R.("edge_color")(R.cluster == "c-SC-s") = "black";
    R.("edge_color")(R.cluster == "c-SC-m") = "#65350f";
    R.("edge_width")(R.cluster ~= "") = 3;
    
    R_new = R;

end


