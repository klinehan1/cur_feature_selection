function [C,U,R,cols,rows, failed] = sf_cur(X, nr, nc, max_iter)

%--------------------------------------------------------------
% CUR using convex optimization and surrogate functions
% For this implementation, we have prioritized speed over storage.
%
%  Input:
%    X: data matrix to be approximated 
%    nr: number of rows for R
%    nc: number of cols for C
%    tol: used in determining which cols/rows are all zeros
%    max_iter: maximum number of iterations for surrogate functional
%       optimizations (good default = 20)
%
%  Output: 
%    C, U, R
%    cols, rows: indices of chosen columns and rows
%    failed: true if bisection cannot reach nc columns within 200 loops (i.e., if X has any identical columns) 
%--------------------------------------------------------------

failed = false;
[m,n] = size(X);

% Find matrix C ----------------------------------------------
%   solving: ||X - XWX||_F^2 + lambda_C*||W||_1,inf   
%       use surrogate functions to decouple the problem   
%   loop using bisection on lambda to gaurantee number of columns

% critical lambda - bisection on lambda from 0 to critical lambda
if m >= n
    XtX = X'*X;
    M = full(XtX*X');  % reshape(X(:),m,n) = X    (nxm)(mxn)(nxm) = mn^2
else  % m < n
    XXt = full(X*X');
    M = X'*(XXt);  % reshape(X(:),m,n) = X    (nxm)(mxn)(nxm) = m^2n
end

lamda_crit = 2*norm(M, Inf);  % matrix infinity norm: max abs row sum

% loop through until desired number of columns chosen using bisection.

lamda_max = lamda_crit;  % choose 0 columns
lamda_min = 0;   % choose all columns
c = 0;

sigmaX = svds(X,1);
mu = (sigmaX^4)*1.001;  

% premultiply matrices needed in optimization 
if m >= n
    XXt = full(X*X');
else
    XtX = X'*X;
end
XXtX = M';
clear M;

stop_count = 0;
while c ~= nc    
    
    stop_count = stop_count + 1;
    %lamda_min
    %lamda_max
    
    lamda = (lamda_max + lamda_min)/2;
    
    % optimization 
    [W, exitflag] = sf_opt('C', m, n, XXtX, XXt, XtX, lamda, mu, max_iter);
    %exitflag

    % find nonzero rows in W
    s = max(abs(W),[],2);
    cols = (s ~= 0);  % can be exactly zero rather than a tolerance, thm 3.5 in prox op paper
    c = sum(cols);

    % update lamda for next iteration
    if c > nc
       lamda_min = lamda; 
    elseif c < nc
       lamda_max = lamda; 
    end

    if (stop_count == 200)
        C = []; U = []; R = [];
        cols = []; rows = [];
        failed = true;
        disp('FAILED')
        return;
    end
    
end

C = X(:,cols);
cols = find(cols);


% Find matrix R ----------------------------------------------
%   solving: ||X - CWX||_F^2 + lambda_R*||W||_inf,1   
%       use surrogate functions to decouple the problem   
%   loop using bisection on lambda to gaurantee number of rows

% critical lambda - bisection on lambda from 0 to critical lambda

if m >= c 
   N = (C'*X)*X';  % reshape(X(:),m,n) = X   
else % m < c 
   N = C'*(XXt);  % reshape(X(:),m,n) = X
end

lamda_crit = 2*norm(N, 1);  % matrix one norm: max abs col sum

% loop through until desired number of rows chosen using bisection.

lamda_max = lamda_crit;  % choose 0 rows
lamda_min = 0;   % choose all rows
r = 0;

sigmaC = svds(C,1);
mu = (sigmaX^2*sigmaC^2)*1.001;

% premultiply matrices needed in optimization
% XXt already computed
CtC = C'*C;
XXtC = full(N');
clear N;

while r ~= nr    
    
    %lamda_min
    %lamda_max
    
    lamda = (lamda_max + lamda_min)/2;
    
    % optimization     
    [W, exitflag] = sf_opt('R', m, c, XXtC, XXt, CtC, lamda, mu, max_iter);
    %exitflag
    
    % find nonzero cols in W
    s = max(abs(W),[],1);
    rows = (s ~= 0);
    r = sum(rows);

    % update lamda for next iteration
    if r > nr
       lamda_min = lamda; 
    elseif r < nr
       lamda_max = lamda; 
    end

end

R = X(rows,:);   
rows = find(rows);

% Compute U -------------------------------------

% fine for the dense case, but not sparse case
% U = pinv(C)*X*pinv(R);

U = pinv(full(C))*X*pinv(full(R));  % (cxm)(mxn)(nxr)

end








