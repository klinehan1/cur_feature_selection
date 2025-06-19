% Compare performance of CUR algorithms on the Genetics dataset. 
% This experiment is an extension of that in Section 6.3 of 
% Sorenson and Embree.
%   1. Relative Error: ||A-CUR||_F/||A||_F
%   2. Time for matrix approximation
%   3. Classification task performance
%
% CUR algorithms: 
% - SF (ours)  
% - LS-R, LS-D (Mahoney and Drineas)
% - DEIM (Sorenson and Embree)  
% - QR (Sorenson and Embree modification, original: Stewart)
% - SVD
%--------------------------------------------------------

%% add paths

addpath('../../ls_cur', '../../sf_cur', '../../deim_cur', '../../qr_cur', ...
    './data');

%% load data

load gene.mat

% center rows of matrix 
X = G - mean(G,2);
clear G;
X = X';  % patients x probes
[m,n] = size(X);

rng(1);  % setting for LS randomized CUR

%% 1. Find Relative Error of each CUR approximation for ncol = nrow 

nc_nr = [5:5:105 107];
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
    [C,U,R, ~,~] = sf_cur(X, nc_nr(i), nc_nr(i), 20);
    sf_t(i) = toc;
    sf_err(i) = norm(X-C*U*R,'fro')/X_fro;
    
    % LS-R - sampling, avg over 5 runs
    temp = zeros(5,1);
    temp_t = zeros(5,1);
    for j = 1:5
        tic;
        [C,U,R, ~,~] = ls_cur(X, 2, nc_nr(i), nc_nr(i), true);
        temp_t(j) = toc;
        temp(j) = norm(X-C*U*R,'fro')/X_fro;
    end
    ls_r_err(i) = mean(temp);
    ls_r_sd(i) = std(temp);
    ls_r_t(i) = mean(temp_t);
    ls_r_t_sd(i) = std(temp_t);

    % LS-D - highest leverage scores
    tic;
    [C,U,R, ~,~] = ls_cur(X, 2, nc_nr(i), nc_nr(i), false);
    ls_d_t(i) = toc;
    ls_d_err(i) = norm(X-C*U*R,'fro')/X_fro;

    % DEIM
    tic
    [C,U,R, ~,~] = deim_cur(X, nc_nr(i));
    deim_t(i) = toc;
    deim_err(i) = norm(X-C*U*R,'fro')/X_fro;

    % QR 
    tic;
    [C,U,R, ~,~] = qr_cur(X, nc_nr(i), nc_nr(i));
    qr_t(i) = toc;
    % C = X(:, q(1:nc_nr(i)) );
    % [~,~,p] = qr(C', "vector");      
    % R = X( p(1:nc_nr(i)), :);
    % U = pinv(full(C))*X*pinv(full(R));   
    qr_err(i) = norm(X-C*U*R,'fro')/X_fro;

    % SVD
    tic;
    [U,S,V] = svds(X, nc_nr(i));
    svd_t(i) = toc;
    svd_err(i) = norm(X-U*S*V','fro')/X_fro;
    %svd_err(i) = norm(X-W(:,1:nc_nr(i))*S(1:nc_nr(i),1:nc_nr(i))*V(:,1:nc_nr(i))','fro')/X_fro;

end

%% 

plot(nc_nr,sf_err, '-', 'Marker', ".", 'MarkerSize', 17, 'LineWidth',2);
hold on; grid on;
errorbar(nc_nr, ls_r_err, ls_r_sd, '-', 'LineWidth',2);
plot(nc_nr,ls_d_err, '-d', 'LineWidth',2);
plot(nc_nr,deim_err, '-o', 'LineWidth',2);
plot(nc_nr,qr_err, '-^', 'LineWidth',2);
plot(nc_nr,svd_err, '-s', 'LineWidth',2);
legend('SF','LS-R','LS-D','DEIM','QR', 'SVD');
title('CUR Accuracy on Gene Expression Data');
xlabel("Number of Selected Columns/Rows");
ylabel("Relative Error");
ax = gca;
ax.Children = ax.Children([6 1 2 3 4 5]);

print('gene_exp_acc.eps', '-depsc')

%%

semilogy(nc_nr,sf_t, '-', 'Marker', ".", 'MarkerSize', 17, 'LineWidth',2);
hold on; grid on;
errorbar(nc_nr, ls_r_t, ls_r_t_sd, '-', 'LineWidth',2);
semilogy(nc_nr,ls_d_t, '-d', 'LineWidth',2);
semilogy(nc_nr,deim_t, '-o', 'LineWidth',2);
semilogy(nc_nr,qr_t, '-^', 'LineWidth',2);
semilogy(nc_nr,svd_t, '-s', 'LineWidth',2);
legend('SF','LS-R','LS-D','DEIM','QR', 'SVD', 'Location', 'best','NumColumns',2);
title('CUR Time on Gene Expression Data');
xlabel("Number of Selected Columns/Rows");
ylabel("Time (seconds)");
ax = gca;
ax.Children = ax.Children([6 1 2 3 4 5]);

print('gene_exp_time.eps', '-depsc')

%% 2. Probe selection performance - top 15 probes

load tumor_idx.mat 
load probes.mat

well_idx = setdiff(1:107,tumor_idx);

for i=1:5

    if i==1
        % SF 
        [~,~,~, cidx,~] = sf_cur(X, 15, 15, 20);
        disp('SF')
        f = "sf";
    elseif i==2
        % LS-D - highest leverage scores
        [~,~,~, cidx,~] = ls_cur(X, 2, 15, 15, false);
        disp('LS')
        f = "ls";
    elseif i==3
        % DEIM
        [~,~,~, cidx,~] = deim_cur(X, 15);
        disp('DEIM')
        f = "deim";
    elseif i==4
        % QR 
        [~,~,~, cidx,~] = qr_cur(X, 15, 15);
        disp('QR')
        f = "qr";
    else % i = 5
        % SVD - find proteins that correlate to first 2 principal comp's
        [U,S,V] = svds(X,2);
        P=U*S;
        C = corr(X,P);

        [VAL, IDX] = maxk(C,15,'ComparisonMethod','abs');
        [~,idx] = maxk(VAL(:),15,'ComparisonMethod','abs');
        cidx = IDX(idx);
        disp('SVD')
        f = "svd";
    end

    % probe names
    probe_names = convertCharsToStrings(probes(cidx,:).ID_REF);
    
    % count number of well patients with values > 1
    well_large_val = sum(X(well_idx,cidx)>1, 1)';
    
    % count number of sick patients with values > 1
    sick_large_val = sum(X(tumor_idx,cidx)>1, 1)';
    
    T =  array2table([probe_names, well_large_val, sick_large_val], 'VariableNames', ...
            ["Probe", "Well", "Sick"]);
    disp(T)
    writetable(T, "classify_"+f+".txt", 'Delimiter', '|');

end


%% 3. Probe selection performance - top 5, 10, ..., 100 probes

load tumor_idx.mat 

well_idx = setdiff(1:107,tumor_idx);

num_probes = [5:5:100];
a = length(num_probes);

sf_diffmed = zeros(a,1);
ls_d_diffmed = zeros(a,1);
deim_diffmed = zeros(a,1);
qr_diffmed = zeros(a,1);
svd_diffmed = zeros(a,1);

sf_diffmean = zeros(a,1);
ls_d_diffmean = zeros(a,1);
deim_diffmean = zeros(a,1);
qr_diffmean = zeros(a,1);
svd_diffmean = zeros(a,1);

sf_diffstd = zeros(a,1);
ls_d_diffstd = zeros(a,1);
deim_diffstd = zeros(a,1);
qr_diffstd = zeros(a,1);
svd_diffstd = zeros(a,1);

sf_ls_overlap = zeros(a,1);

for i=1:a

    disp(i)

    % SF 
    [~,~,~, cidx,~] = sf_cur(X, 15, num_probes(i), 20);
    well_large_val = sum(X(well_idx,cidx)>1, 1)';
    sick_large_val = sum(X(tumor_idx,cidx)>1, 1)';
    sf_diffmed(i) = median(abs(well_large_val - sick_large_val));
    sf_diffmean(i) = mean(abs(well_large_val - sick_large_val));
    sf_diffstd(i) = std(abs(well_large_val - sick_large_val));
    sf_cidx = cidx;

    % LS-D - highest leverage scores
    [~,~,~, cidx,~] = ls_cur(X, 2, 15, num_probes(i), false);
    well_large_val = sum(X(well_idx,cidx)>1, 1)';
    sick_large_val = sum(X(tumor_idx,cidx)>1, 1)';
    ls_d_diffmed(i) = median(abs(well_large_val - sick_large_val));
    ls_d_diffmean(i) = mean(abs(well_large_val - sick_large_val));
    ls_d_diffstd(i) = std(abs(well_large_val - sick_large_val));
    ls_cidx = cidx;

    % overlap between SF and LS selected probes
    sf_ls_overlap(i) = length(intersect(sf_cidx, ls_cidx));

    % DEIM
    [~,~,~, cidx,~] = deim_cur(X, num_probes(i));
    well_large_val = sum(X(well_idx,cidx)>1, 1)';
    sick_large_val = sum(X(tumor_idx,cidx)>1, 1)';
    deim_diffmed(i) = median(abs(well_large_val - sick_large_val));
    deim_diffmean(i) = mean(abs(well_large_val - sick_large_val));
    deim_diffstd(i) = std(abs(well_large_val - sick_large_val));

    % QR 
    [~,~,~, cidx,~] = qr_cur(X, 15, num_probes(i));
    well_large_val = sum(X(well_idx,cidx)>1, 1)';
    sick_large_val = sum(X(tumor_idx,cidx)>1, 1)';
    qr_diffmed(i) = median(abs(well_large_val - sick_large_val));
    qr_diffmean(i) = mean(abs(well_large_val - sick_large_val));
    qr_diffstd(i) = std(abs(well_large_val - sick_large_val));

    % SVD - find proteins that correlate to first 2 principal comp's
    [U,S,V] = svds(X,2);
    P=U*S;
    C = corr(X,P);

    [VAL, IDX] = maxk(C, num_probes(i),'ComparisonMethod','abs');
    [~,idx] = maxk(VAL(:), num_probes(i),'ComparisonMethod','abs');
    cidx = IDX(idx);
    well_large_val = sum(X(well_idx,cidx)>1, 1)';
    sick_large_val = sum(X(tumor_idx,cidx)>1, 1)';
    svd_diffmed(i) = median(abs(well_large_val - sick_large_val));
    svd_diffmean(i) = mean(abs(well_large_val - sick_large_val));
    svd_diffstd(i) = std(abs(well_large_val - sick_large_val));

end

%%

T =  array2table([num_probes', sf_diffmed, round(sf_diffmean,2), ...
    ls_d_diffmed, round(ls_d_diffmean,2), deim_diffmed, round(deim_diffmean,2), ...
    qr_diffmed, round(qr_diffmean,2), svd_diffmed, round(svd_diffmean,2)], 'VariableNames', ...
            ["Number of Probes", "SF-mdn", "SF-mean", ...
            "LS-D-mdn", "LSD-D-mean", "DEIM-mdn", "DEIM-mean", ...
            "QR-med", "QR-mean", "SVD-med", "SVD-mean"]);
disp(T)
writetable(T, "probe_sel_results.txt", 'Delimiter', '|');

T2 =  array2table([num_probes', round(sf_diffstd,2), ...
    round(ls_d_diffstd,2), round(deim_diffstd,2), ...
    round(qr_diffstd,2), round(svd_diffstd,2)], 'VariableNames', ...
            ["Number of Probes", "SF-std", ...
            "LS-D-std", "DEIM-std", ...
            "QR-std", "SVD-std"]);
disp(T2)
writetable(T2, "probe_sel_std.txt", 'Delimiter', '|');

num_probes' - sf_ls_overlap
