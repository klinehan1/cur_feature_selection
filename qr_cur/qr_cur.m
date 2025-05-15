function [C,U,R, q, p] = qr_cur(A, nr, nc)

% -----------------------------------------------
% CUR from Pivoted QR
% modification from DEIM CUR, Sorenson and Embree (original: Stewart)
%
% Input:
% - A: data matrix
% - nr, nc: number of rows and columns to select
%
% Output: 
% - C,U,R, 
% - q, p: indices of chosen columns and rows.
%------------------------------------------------

    [~,~,q] = qr(full(A), "vector");
    q = q(1:nc);
    C = A(:,q);

    [~,~,p] = qr(full(C'), "vector");  
    p = p(1:nr);    
    R = A(p,:);
    
    U = pinv(full(C))*A*pinv(full(R));

end