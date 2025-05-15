function control_exp(wrs, cur_alg, fig_path, SOM_size_x, ...
    SOM_size_y, num_nodes, classes_filt, mouse_ids_filt, P, net, R_new, ...
    class_A, class_B)

% --------------------------------------------
% SOM Control Experiment - extension of that in PLOS One paper by Higuera, 
% Gardiner, and Cios.
%  
%   1. SOM on all proteins + clustering 
%   2. Feature selection using either Wilcoxon rank sum test or CUR
%   3. SOM using selected proteins
%   4. Compare results of 1 and 3
%
% Each run will provide results for a CUR algorithm using AIC and BIC model 
% criterion selection or the Wilcoxon rank sum test.  
% 
% Input:
%  fig_path: e.g., './figures/wilcoxon/', './figures/sf_cur/'
%  wrs: true (use Wilcoxon rank sum test) or false (use CUR) 
%  cur_alg: 'sf', 'ls-d', 'deim', or 'qr' ('' if wrs=True)
%  SOM_size_x, SOM_size_y, num_nodes, classes_filt, mouse_ids_filt, P, 
%    net, R_new, class_A, class_B: all set in control_exp_driver.m
%-----------------------------------------

% -------- Feature Selection between pairwise class comparisons ------------

% get weights - weight vector of node i is in net.IW{1}(i,:)
W = net.IW{1};  % num_nodes x 77

for i=1:6

    % filter data
    class_A_nodes = find(R_new.cluster == class_A(i));
    class_B_nodes = find(R_new.cluster == class_B(i));
    class_A_weights = W(class_A_nodes, :);
    class_B_weights = W(class_B_nodes, :);
    
    if wrs
        % Wilcoxon rank sum test feature selection
        [sel_protein_idx] = wilcoxon_feature_sel(class_A_weights, class_B_weights);
    
        if i==1
            CSs_SCs_w = sel_protein_idx;
        elseif i==2
            CSm_SCm_w = sel_protein_idx;
        elseif i==3
            CSm_SCs_w = sel_protein_idx;
        elseif i==4
            CSs_SCm_w = sel_protein_idx;
        elseif i==5
            SCm_SCs_w = sel_protein_idx;
        else % i==6
            CSm_CSs_w = sel_protein_idx;
        end
    
    else 
        % CUR feature selection
        [sel_protein_idx_aic, sel_protein_idx_bic] = cur_feature_sel(class_A_weights, class_B_weights, cur_alg);
    
        if i==1
            CSs_SCs_aic = sel_protein_idx_aic;
            CSs_SCs_bic = sel_protein_idx_bic;
        elseif i==2
            CSm_SCm_aic = sel_protein_idx_aic;
            CSm_SCm_bic = sel_protein_idx_bic;
        elseif i==3
            CSm_SCs_aic = sel_protein_idx_aic;
            CSm_SCs_bic = sel_protein_idx_bic;
        elseif i==4
            CSs_SCm_aic = sel_protein_idx_aic;
            CSs_SCm_bic = sel_protein_idx_bic;
        elseif i==5
            SCm_SCs_aic = sel_protein_idx_aic;
            SCm_SCs_bic = sel_protein_idx_bic;
        else % i==6
            CSm_CSs_aic = sel_protein_idx_aic;
            CSm_CSs_bic = sel_protein_idx_bic;
        end
       
    end

    %proteins(sel_protein_idx)

end

for i=1:2

    if wrs
        CSs_SCs = CSs_SCs_w; CSm_SCm = CSm_SCm_w; CSm_SCs = CSm_SCs_w; 
        CSs_SCm = CSs_SCm_w; SCm_SCs = SCm_SCs_w; CSm_CSs = CSm_CSs_w;
        fp_cont = '';
    else % cur
        if i==1 % aic
            CSs_SCs = CSs_SCs_aic; CSm_SCm = CSm_SCm_aic; CSm_SCs = CSm_SCs_aic; 
            CSs_SCm = CSs_SCm_aic; SCm_SCs = SCm_SCs_aic; CSm_CSs = CSm_CSs_aic;
            fp_cont = '_aic';
            disp('---- AIC ----')
        else % bic
            CSs_SCs = CSs_SCs_bic; CSm_SCm = CSm_SCm_bic; CSm_SCs = CSm_SCs_bic; 
            CSs_SCm = CSs_SCm_bic; SCm_SCs = SCm_SCs_bic; CSm_CSs = CSm_CSs_bic;
            fp_cont = '_bic';
            disp('---- BIC ----')
        end
    end

    % find intersection of four control experiments 
    common_prot_idx = intersect( intersect(CSs_SCm,CSs_SCs), ...
        intersect(CSm_SCm,CSm_SCs));
    
    %proteins(common_prot_idx)
    
    % ------------- CS experiment --------------
    disp('--- CS ---')

    % Select SOM using selected proteins 
    combined_prot_idx2 = union(common_prot_idx, CSm_CSs);
    num_disc_prots = length(combined_prot_idx2)
    P_fs = P(:,combined_prot_idx2);
    [~, seed] = SOM_select(SOM_size_x, SOM_size_y, P_fs, classes_filt);
    
    % SOM using selected proteins
    rng(seed);
    net_fs = selforgmap([SOM_size_x SOM_size_y], 100, 3, 'hextop', 'linkdist');  % other than size - defaults.  
    net_fs = train(net_fs,P_fs');  % observations need to be columns, not rows, P': ?x570
    y_fs = net_fs(P_fs'); % get SOM node for each observation
    
    % SOM results 
    R_fs = org_SOM_results(classes_filt, y_fs, num_nodes);
    [R_fs_new] = cluster_SOM(hextop([SOM_size_x SOM_size_y]), mouse_ids_filt, y_fs, R_fs, num_nodes);
    figure(2*(i-1)+1);
    visualize_SOM(R_fs_new, SOM_size_x, SOM_size_y)
    print(fig_path + 'CS' + fp_cont + '.eps', '-depsc')
    [num_mixed_CS_nodes, num_mixed_CS_meas] = calc_measures(R_fs,"c-CS-m","c-CS-s")

 % -------------- SC experiment -----------------
     disp('--- SC ---')

    % Select SOM using selected proteins 
    combined_prot_idx1 = union(common_prot_idx, SCm_SCs);
    num_disc_prots = length(combined_prot_idx1)
    P_fs = P(:,combined_prot_idx1);
    [~, seed] = SOM_select(SOM_size_x, SOM_size_y, P_fs, classes_filt);
    
    % SOM using selected proteins  
    rng(seed);
    net_fs = selforgmap([SOM_size_x SOM_size_y], 100, 3, 'hextop', 'linkdist');  % other than size - defaults.  
    net_fs = train(net_fs,P_fs');  % observations need to be columns, not rows, P': ?x570
    y_fs = net_fs(P_fs'); % get SOM node for each observation
    
    % SOM results 
    R_fs = org_SOM_results(classes_filt, y_fs, num_nodes);
    [R_fs_new] = cluster_SOM(hextop([SOM_size_x SOM_size_y]), mouse_ids_filt, y_fs, R_fs, num_nodes);
    figure(2*(i-1)+2);
    visualize_SOM(R_fs_new, SOM_size_x, SOM_size_y)
    print(fig_path + 'SC' + fp_cont + '.eps', '-depsc')
    [num_mixed_SC_nodes, num_mixed_SC_meas] = calc_measures(R_fs,"c-SC-m","c-SC-s")    

    if wrs
        break
    end
end

