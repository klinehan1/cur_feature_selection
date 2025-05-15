function [prox, tau] = prox_op(x, alpha)

% --------------------------------------------------------------------
% Algorithm to find the proximal operator of the L-infinity norm from 
% our arXiv paper: https://arxiv.org/abs/2408.11211
%   prox(x) = argmin_y [ (1/2)||y-x||_2^2 + alpha*||y||_inf ] 
% 
% Input:
%   x: data (vector)
%   alpha: constant on penalty term
%    
% Output: 
%   prox: prox(x)
%   tau: threshold to compute prox(x)
% --------------------------------------------------------------------

    if alpha >= norm(x,1) 
        prox = zeros(size(x));
        tau = 0;
    else     
        % permute x to be in decreasing order in abs value
        m = length(x);
        s = sort(abs(x),'descend');
        s(end+1) = 0;

        % find value for minimizer    
        s_sum = 0;
        i = 1;
        while i <= m  % len(x) = m
            s_i = s(i);
            s_sum = s_sum + s_i;
            
            % check for repeated elements
            j = 1;
            while (i+j <= m) && (s(i+j) == s_i)        
                s_sum = s_sum + s_i;
                j = j+1;       
            end
            i = i + (j-1); 

            t0 = (s_sum - alpha)/i;  % minimizer

            if (t0 <= s(i)) && (t0 > s(i+1)) 
                tau = t0;
                break;
            end

            i = i+1;
        end

        % compute proximal operator
        prox = x;
        idx = (abs(x) > tau);
        prox(idx) = sign(x(idx))*tau;
        
    end
end