function [R] = org_SOM_results(classes, y, num_nodes)

%--------------------------------------------------------------------------
%  Organize SOM results for analysis
%
%  Input:
%    classes: vector of mice classes 
%    y: SOM results - num_nodes x num_observations, binary data, 
%       y(i,j) = 1 if observation j is in node i. 
%    num_nodes = number of neurons in SOM
%
%  Output: R, table of results - each row corresponds to a node, columns
%  include counts for each class, number of observations, number of 
%  classes, majority and minority class, color of node for visualization   
%--------------------------------------------------------------------------

    class_names = unique(classes);
    R = array2table(zeros(num_nodes,length(class_names)), 'VariableNames', class_names);
    
    % find frequency of each class in each node
    for i=1:num_nodes
        idx = find(y(i,:));  % gives indices of observations in node i
        t = tabulate(classes(idx));  % cell array of value_counts for classes in node i
        tbl = cell2table(t,'VariableNames', {'Value','Count','Percent'});
    
        for j=1:height(tbl) % nrow(tbl)
            c = tbl.Value{j};
            R.(c)(i) = tbl.Count(j);
        end
    end
    
    % calculate total number of observations in each node
    R.("total_observations") = zeros(num_nodes,1);
    
    for i=1:num_nodes
        R.total_observations(i) = sum(R{i,1:length(class_names)});
    end
    
    R.("num_classes") = zeros(num_nodes,1);
    R.("majority_class") = strings(num_nodes,1);
    R.("minority_class") = strings(num_nodes,1);
    R.("face_color") = strings(num_nodes,1);
    R.("face_color")(:) = "none";
    
    % find number of classes, majority and minority class for each node
    for i=1:num_nodes
        %disp(i)
    
        [res, ix] = sort(R{i,1:length(class_names)}, 'descend');
        gt_zero_idx = (res > 0);
        col_num = ix(gt_zero_idx);
    
        R.num_classes(i) = length(col_num);

        %if length(col_num) > 2
        %    disp('MORE THAN 2 CLASSES IN A NODE')
        %end
    
        if length(col_num) >= 1
    
            R.majority_class(i) = R.Properties.VariableNames{col_num(1)};
        
            if length(col_num) > 1  % >1 class present in node i 

                s = "";
                for j=2:length(col_num)
                    if j > 2
                        s = s + "|";
                    end
                    s = s + R.Properties.VariableNames{col_num(j)};
                end
                R.minority_class(i) = s;
            end
        end
    
    end
    
    R.("face_color")(R.majority_class == "c-CS-s") = "#b8e186";
    R.("face_color")(R.majority_class == "c-CS-m") = "#f0e442";
    R.("face_color")(R.majority_class == "c-SC-s") = "#d55e00";
    R.("face_color")(R.majority_class == "c-SC-m") = "#e69f00";

end