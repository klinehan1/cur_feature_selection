function [protein_idx] = wilcoxon_feature_sel(class_A_dat, class_B_dat)

%--------------------------------------------------------------------------
%  Feature selection using the Wilcoxon Rank Sum Test.  This is a test
%  between two classes, checks for equal medians. Will find discriminant 
%  proteins between two classes.  
%
%  Input:
%   class_A_dat: input data for class A (dim 77)
%   class_B_dat: input data for class B (dim 77)
%
%  Output: 
%   protein_idx: columns (proteins) selected from the data  
%--------------------------------------------------------------------------

    % Wilcoxon test for each protein
    
    protein_idx = [];
    j = 0;
    
    for i=1:77   % 77 can be hardcoded - the number of proteins is a constant
    
        % the more different the two vectors, the lower p is. default - two-sided test
        p = ranksum(class_A_dat(:,i), class_B_dat(:,i));  
    
        if p < 0.05
            j = j + 1;
            protein_idx(j) = i;
        end
    end

end