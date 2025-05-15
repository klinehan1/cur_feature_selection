function [C,U,R, col_idx, row_idx] = ls_cur(A, k, nr, nc, sample)

% -----------------------------------------------
% CUR from Mahoney and Drineas, 2009
%
% Input:
% - A: data matrix
% - k: rank to use in leverage score calculation
% - nr, nc: number of rows and columns to select
% - sample: true - use sampling, false - return cols/rows with highest
% leverage scores
%
% Output: 
% - C,U,R, 
% - col_idx and row_idx: indices of chosen columns and rows.
%------------------------------------------------

[m,n] = size(A);

% calculate column and row sampling probabilities (q, p) ---------

[Ua,~,Va] = svds(A,k);

q = zeros(n,1);

for j = 1:n
    q(j) = (1/k)*(norm(Va(j,:))^2); 
end


p = zeros(m,1);

for i = 1:m
    p(i) = (1/k)*(norm(Ua(i,:))^2); 
end

if sample == true
    % form C and R using Bernoulli trials --------------------
    
    rc = rand(n,1);
    adj_q = nc*q;
    col_idx = (rc <= adj_q);
    C = A(:,col_idx);
    col_idx = find(col_idx);
    
    rr = rand(m,1);
    adj_p = nr*p;
    row_idx = (rr <= adj_p);
    R = A(row_idx,:);
    row_idx = find(row_idx);

    % sample with replacement ----------------------
    %C = datasample(full(A),nc,2,'Replace',false,'Weights',q);
    %col_idx = sort(datasample(1:n,nc,'Replace', false, 'Weights',q));
    %C = A(:,col_idx);
    
    %R = datasample(full(A),nr,1,'Replace',false,'Weights',p);
    %row_idx = sort(datasample(1:m,nr,'Replace', false, 'Weights',p));
    %R = A(row_idx,:);

else % sample == false
    % choose cols/rows with largest leverage scores ----------
    [~,col_idx] = sort(q, "descend");
    col_idx = col_idx(1:nc);
    C = A(:, col_idx);

    [~,row_idx] = sort(p, "descend");
    row_idx = row_idx(1:nr);
    R = A(row_idx, :);

end

% Compute U ---------------------------

% fine for dense case
%U = pinv(C)*A*pinv(R);

U = pinv(full(C))*A*pinv(full(R));

end