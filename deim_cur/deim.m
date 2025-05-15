function [p] = deim(V)

% -----------------------------------------------
% DEIM technique as in Sorenson and Embree
%
% Input:
% - V: data matrix
%
% Output: 
% - p: indices of chosen rows.
%------------------------------------------------

    [~,k] = size(V);

    v = V(:,1);
    [~,p_idx] = max(abs(v));
    p = [p_idx];

    for j = 2:k
        v = V(:,j);
        c = V(p, 1:j-1)\v(p);
        r = v - V(:, 1:j-1)*c;
        [~,p_idx] = max(abs(r));
        p = [p; p_idx];
    end

end