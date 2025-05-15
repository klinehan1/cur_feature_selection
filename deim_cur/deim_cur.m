function [C,U,R,cols,rows] = deim_cur(A,k)

% ----------------------------------------
% CUR using DEIM - Sorenson and Embree
%
% Input:
% - A: data matrix
% - k = # of cols to choose = # of rows to choose
%
% Output:
% - C,U,R
% - cols, rows: indices of chosen rows/columns
% ---------------------------------------

[V,~,W] = svds(A,k);   % A ~ VSW'

rows = deim(V);
cols = deim(W);

C = A(:,cols);
R = A(rows,:);

U = pinv(full(C))*A*pinv(full(R));

end