function visualize_SOM(R, x_size, y_size)

%--------------------------------------------------------------------------
%  Visualize SOM with clusters using a hex grid with boxes 
%
%  Input:
%    R: table of SOM results created in org_SOM_results() and cluster_SOM()
%    x_size: number of boxes in grid in x dimension
%    y_size: number of boxes in grid in y dimension   
%
%  Output: R_new, table of results - R appended with columns for clustering 
%  results and visualization, i.e., cluster, edge color, edge width   
%--------------------------------------------------------------------------
    
    x_pts=[0.1 0.9 0.9 0.1]; % x-coordinates of the vertices
    y_pts=[0.1 0.1 0.9 0.9]; % y-coordinates of the vertices
    
    for j=0:(y_size-1)
        for i=0:(x_size-1)
            k = (i+1) + (x_size)*j; % map i,j index to 1 to num_nodes

            if strlength(R.minority_class(k)) <= 6
                s = {R.majority_class(k), R.minority_class(k), R.total_observations(k)};
                font_s = 7;
            else % more than 1 minority class
                min_class_s = {};
                temp = strsplit(R.minority_class(k), "|");
                for a=1:length(temp)
                    min_class_s{a} = temp(a);
                end
                s = {R.majority_class(k), min_class_s{1:end}, R.total_observations(k)};
                if length(temp) <= 2
                    font_s = 5.5;
                else
                    font_s = 4.5;
                end
            end

            if mod(j,2) == 0
                patch('XData', x_pts+i, 'YData', y_pts+j, 'FaceColor', R.face_color(k), ...
                      'EdgeColor', R.edge_color(k), 'LineWidth', R.edge_width(k)) % make a square with lower right corner at [i,j]
                text(0.5+i, 0.5+j, s, "HorizontalAlignment", "center", "Color", "black", "FontSize", font_s)
            else
                patch('XData', x_pts+0.5+i, 'YData', y_pts+j, 'FaceColor', R.face_color(k), ... 
                      'EdgeColor', R.edge_color(k), 'LineWidth', R.edge_width(k)) % make a square with lower right corner at [i,j]
                text(1+i, 0.5+j, s, "HorizontalAlignment", "center", "Color", "black", "FontSize", font_s)
            end
    
            hold on
        end
    end
    
    axis off

end


