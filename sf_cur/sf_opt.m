function [W, exitflag] = sf_opt(matrix, m, n, XXtX, XXt, XtX, lamda, mu, max_iter)

%--------------------------------------------------------------
% Surrogate Functional Optimization 
%
%  Input:
%    matrix: matrix to solve for ('C' or 'R')
%    m: number of rows of X, the original data matrix
%    n: number of columns of X (matrix == 'C'), or number of columns of C (matrix == 'R')
%    XXtX: X*X'*X when solving for C, X*X'*C when solving for R
%    XXt: X*X' 
%    XtX: X'*X when solving for C, C'*C when solving for R
%    lamda: scalar
%    mu: scalar
%    max_iter: maximum number of iterations for surrogate functional
%      optimizations (good default = 20)
%
%  Output: 
%    W 
%    exitflag: flag for number of iterations completed
%--------------------------------------------------------------

% Minimize J_hat using surrogate functions

    W = zeros(n,m);

    k = 0; d = 100; 
    while (d > 0.00001) && (k < max_iter)  % d: diff in W and W_old
    
        k = k+1;
        W_old = W;  % W_old is Z
       
        L = XXtX + mu*W_old' - XXt*(W_old'*XtX); % (mxn)+(mxn)-(mxm)(mxn)(nxn)

        if matrix == 'C'

            % form rows of W using iteration - vectors for prox op are
            % in columns of L  
            for i = 1:n  
                [W(i,:), ~] = prox_op((1/mu)*L(:,i)', lamda/(2*mu));
            end
        
        else % matrix == 'R'
        
            % form columns of W using iteration - vectors for prox op
            % are in rows of L
            for j = 1:m                
                [W(:,j), ~] = prox_op((1/mu)*L(j,:)', lamda/(2*mu));
            end

        end

        d = norm(W-W_old, 'fro')^2;
         
    end

    if k == max_iter
        exitflag = 'reached maximum number of iterations';
    else
        exitflag = sprintf('convergence criterion achieved in %d iterations', k);
    end 

end

