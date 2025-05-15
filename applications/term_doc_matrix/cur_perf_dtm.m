% Compare performance of CUR algorithms and SVD on term-document matrix. 
%   1. Relative Error: ||A-CUR||_F/||A||_F
%   2. Time for matrix approximation
%
% CUR algorithms: 
% - SF (ours)  
% - LS-R, LS-D (Mahoney and Drineas)
% - DEIM (Sorenson and Embree)  
% - QR (Sorenson and Embree modification, original: Stewart)
% - SVD
%--------------------------------------------------------

%% add paths

addpath('../../ls_cur', '../../sf_cur', '../../deim_cur', '../../qr_cur');
disp("paths added")

%% load data

load newsgroups.mat
[m,n] = size(X);

rng(1);  % setting for LS randomized CUR
disp("data loaded")

%% 1. Find Relative Error of each CUR approximation for ncol = nrow 

nc_nr = [200:200:2200 2295];
a = length(nc_nr);

sf_err = zeros(a,1);
ls_r_err = zeros(a,1);
ls_r_sd = zeros(a,1);
ls_d_err = zeros(a,1);
deim_err = zeros(a,1);
qr_err = zeros(a,1);
svd_err = zeros(a,1);
X_fro = norm(X,'fro');

sf_t = zeros(a,1);
ls_r_t = zeros(a,1);
ls_r_t_sd = zeros(a,1);
ls_d_t = zeros(a,1);
deim_t = zeros(a,1);
qr_t = zeros(a,1);
svd_t = zeros(a,1);

%[~,~,q] = qr(X, "vector");
%[W,S,V] = svds(X, nc_nr(end));

for i=1:a

    disp(i)

    % SF
    tic;
    [C,U,R, ~,~] = cur_sf(X, nc_nr(i), nc_nr(i), 20);
    sf_t(i) = toc;
    sf_err(i) = norm(X-C*U*R,'fro')/X_fro;
    
    % LS-R - sampling, avg over 5 runs
    temp = zeros(5,1);
    temp_t = zeros(5,1);
    for j = 1:5
        tic;
        [C,U,R, ~,~] = ls_cur(X, 10, nc_nr(i), nc_nr(i), true);
        temp_t(j) = toc;
        temp(j) = norm(X-C*U*R,'fro')/X_fro;
    end
    ls_r_err(i) = mean(temp);
    ls_r_sd(i) = std(temp);
    ls_r_t(i) = mean(temp_t);
    ls_r_t_sd(i) = std(temp_t);
    disp('LS-R done')
    
    % LS-D - highest leverage scores
    tic;
    [C,U,R, ~,~] = ls_cur(X, 10, nc_nr(i), nc_nr(i), false);
    ls_d_t(i) = toc;
    ls_d_err(i) = norm(X-C*U*R,'fro')/X_fro;
    disp('LS-D done')

    % DEIM
    tic;
    [C,U,R, ~,~] = deim_cur(X, nc_nr(i));
    deim_t(i) = toc;
    deim_err(i) = norm(X-C*U*R,'fro')/X_fro;
    disp('DEIM done')

    % QR 
    tic;
    [C,U,R, ~,~] = qr_cur(X, nc_nr(i), nc_nr(i));
    qr_t(i) = toc;
    % C = X(:, q(1:nc_nr(i)) );
    % [~,~,p] = qr(C', "vector");      
    % R = X( p(1:nc_nr(i)), :);
    % U = pinv(full(C))*X*pinv(full(R));   
    qr_err(i) = norm(X-C*U*R,'fro')/X_fro;
    disp('QR done')

    % SVD
    tic;
    [U,S,V] = svds(X, nc_nr(i));
    svd_t(i) = toc;
    svd_err(i) = norm(X-U*S*V','fro')/X_fro;
    %svd_err(i) = norm(X-W(:,1:nc_nr(i))*S(1:nc_nr(i),1:nc_nr(i))*V(:,1:nc_nr(i))','fro')/X_fro;
    disp('SVD done')
end

%% save errors and times

save("accuracy_results.mat", "sf_err", "ls_r_err", "ls_r_sd", ...
    "ls_d_err", "deim_err", "qr_err", "svd_err", "nc_nr"); 

save("time_results.mat", "sf_t", "ls_r_t", "ls_r_t_sd", ...
    "ls_d_t", "deim_t", "qr_t", "svd_t", "nc_nr");

