function [protein_idx_aic, protein_idx_bic] = cur_feature_sel(class_A_dat, ...
    class_B_dat, cur_alg)

%--------------------------------------------------------------------------
%  Feature selection using CUR algorithms.  Will find discriminant proteins between 
%  two classes, and the number of discriminant proteins using AIC and BIC.  
%
% CUR algorithms: 
% - SF (ours)  
% - LS-D (Mahoney and Drineas)
% - DEIM (Sorenson and Embree)  
% - QR (Stewart)
%
%  Input:
%   class_A_dat: input data for class A (dim 77)
%   class_B_dat: input data for class B (dim 77)
%   cur_alg: CUR algorithm to use ('sf', 'ls-d', 'deim', or 'qr')
%
%  Output: 
%   protein_idx_aic: columns (proteins) selected from the data using AIC
%   protein_idx_bic: columns (proteins) selected from the data using BIC
%--------------------------------------------------------------------------

    % Let the matrix D contain all pairwise differences of the data matrices
    % for class A and B.  

    sz_A = size(class_A_dat,1);
    sz_B = size(class_B_dat,1);
    
    D = zeros(sz_A*sz_B ,77);
    
    for i = 1:sz_A
    
        A_node_wgt = class_A_dat(i,:); 
    
        for j = 1:sz_B    
            D(j+(i-1)*sz_B,:) = (A_node_wgt - class_B_dat(j,:));
        end
    end
    
    % run CUR on D
    [m,n] = size(D);
    
    if strcmp(cur_alg,'ls-d') | strcmp(cur_alg,'deim') 
        % look at spectrum for LS-D + precompute for DEIM
        [~,S,V] = svd(D);
        s = diag(S);
    end

    rss = zeros(77,1);
    aic = zeros(77,1);
    bic = zeros(77,1);

    % precompute for QR
    if strcmp(cur_alg,'qr') == 1
        [~,~,q] = qr(D, "vector");
    end

    for i=1:77

        %disp(i)
        failed = false;
        if strcmp(cur_alg,'sf') == 1
            [C,U,R,~,~,failed] = sf_cur(D, m, i, 20);
    
        elseif strcmp(cur_alg,'ls-d') == 1
            [C,U,R,~,~] = ls_cur(D, 2, m, i, false); 
    
        elseif strcmp(cur_alg,'deim') == 1
            %[C,U,R, ~,~] = deim_cur(D, i);
            cols = deim(V(:,1:i));
            C = D(:,cols);
            R = D;
            U = pinv(full(C))*D*pinv(full(R));

        else % 'qr'
            %[C,U,R, ~,~] = qr_cur(X, nc_nr(i), nc_nr(i));
            C = D(:, q(1:i) );     
            R = D;
            U = pinv(full(C))*D*pinv(full(R));
        end 

        % compute RSS, AIC, BIC 
        % Let k:# entries in U + noise variance + column/row encoding, n:# entries in D
        % AIC = 2k + nlog(RSS/n), BIC = klog(n) + nlog(RSS/n)
        
        if ~failed
            rss(i) = norm(D-C*U*R, 'fro')^2; 
            aic(i) = 2*(m*i+3) + (m*n)*log(rss(i)/(m*n));
            bic(i) = (m*i+3)*log(m*n) + (m*n)*log(rss(i)/(m*n)); 
        end

    end

    %figure(5);
    %subplot(3,1,1);
    %semilogy(1:77, rss, '-o', 'LineWidth', 2);
    %xlabel('Number of Columns')
    %ylabel('RSS')
    
    %subplot(3,1,2);
    %plot(1:77, aic, '-o', 'LineWidth', 2);
    %xlabel('Number of Columns')
    %ylabel('AIC')
    
    %subplot(3,1,3);
    %plot(1:77, bic, '-o', 'LineWidth', 2);
    %xlabel('Number of Columns')
    %ylabel('BIC')

    % select number of columns based on AIC and BIC
    [~,num_cols_aic] = min(aic);
    [~,num_cols_bic] = min(bic);

    % Run CUR for given number of columns - only need C
    if strcmp(cur_alg,'sf') == 1
        [~,~,~,protein_idx_aic,~,~] = sf_cur(D, m, num_cols_aic, 20);
        [~,~,~,protein_idx_bic,~,~] = sf_cur(D, m, num_cols_bic, 20);

    elseif strcmp(cur_alg,'ls-d') == 1
        [~,~,~,protein_idx_aic,~] = ls_cur(D, 2, m, num_cols_aic, false);
        [~,~,~,protein_idx_bic,~] = ls_cur(D, 2, m, num_cols_bic, false);

    elseif strcmp(cur_alg,'deim') == 1
        %[~,~,~,protein_idx,~] = deim_cur(D, i);
        protein_idx_aic = deim(V(:,1:num_cols_aic));
        protein_idx_bic = deim(V(:,1:num_cols_bic));

    else % 'qr'
        %[~,~,~,protein_idx,~] = qr_cur(X, nc_nr(i), nc_nr(i));
        protein_idx_aic = q(1:num_cols_aic);
        protein_idx_bic = q(1:num_cols_bic);

    end 

end