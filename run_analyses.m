
%%% Seguin et al. 2025, medRxiv (https://www.medrxiv.org/content/10.1101/2025.02.10.25322034v2)

% This script runs the main analyses of the manuscript. 
% It computes analyses and generates the figures of the manuscript.
% The script uses publicly available data from Weigand et al 2018 (Biol Psychiatry), named Cohort I in our study. 
% Data from Cohort II is not publicly available and therefore analyses using this patient group are not reproduced here. 
% We perform analyses on two normative connectomes constructed from publicly available datasets from 
% Griffa et al 2019 (Lausanne connectome) and Mansour et al 2021 (Schaefer connectome). 
% We use functions from the BCT (Rubinov & Sporns 2010) to compute network communication measures.

%%% Caio Seguin, 22 Oct 2025

%% Load data

clc; clear; close all;

fprintf('Load data...\n');

addpath data/
addpath func/
addpath func_plot/
addpath external/
addpath pre_computed_results/

%%% Load TMS coordinates and clinical improvement for cohort I (Weigand et al 2018)
load tms_data;

%%% Load connectome datasets
load parc_data;

fprintf('Done.\n');

%% Compute network communication models

fprintf('Compute network communication models...\n');

parc_vec = {'lau_l5', 'sch_800'};
n_parc_vec = length(parc_vec);

for k_parc = 1:n_parc_vec

    parc_lbl = parc_vec{k_parc};

    %%% Input matrices
    ED = parc_data.(parc_lbl).ED; % Euclidean distance
    W = parc_data.(parc_lbl).SC; % weighted SC
    A = double(~~W); % binary SC
    L = -log10(W./(max(W(:)) + min(W(W>0)))); L(L == inf) = 0; % connection lengths
    S = diag(sum(W,1)); % node strength

    %%% Communication matrices
    % Weighted communicability (normalised following Crofts et al)
    CMY_wei = expm(S^(-1/2)*W*S^(-1/2));
    % Binary communicability
    CMY_bin = expm(A);
    % Weighted shortest path length (SPL_wei_hop is used in the main text)
    [SPL_wei, SPL_wei_hop, SPL_wei_Pmat] = distance_wei_floyd(L);
    % Binary shortest path length
    SPL_bin = distance_wei_floyd(A);
    % Search information weighted
    SI_wei = search_information(W,L);
    SI_bin = search_information(A,A);
    % Navigation path length
    [~, NPL_bin, NPL_wei] = navigate_bct(L, ED);

    % Clean diagonals
    n = length(W);
    CMY_wei = CMY_wei.*(~eye(n));
    CMY_bin = CMY_bin.*(~eye(n));
    SPL_wei = SPL_wei.*(~eye(n));
    SPL_wei_hop = SPL_wei_hop.*(~eye(n));
    SPL_bin = SPL_bin.*(~eye(n));
    SI_wei = SI_wei.*(~eye(n));
    SI_bin = SI_bin.*(~eye(n));
    NPL_wei = NPL_wei.*(~eye(n));
    NPL_bin = NPL_bin.*(~eye(n));

    parc_data.(parc_lbl).CMY_wei = CMY_wei;
    parc_data.(parc_lbl).CMY_bin = CMY_bin;
    parc_data.(parc_lbl).SPL_wei = SPL_wei;
    parc_data.(parc_lbl).SPL_wei_hop = SPL_wei_hop;
    parc_data.(parc_lbl).SPL_wei_Pmat = SPL_wei_Pmat;
    parc_data.(parc_lbl).SPL_bin = SPL_bin;
    parc_data.(parc_lbl).SI_wei = SI_wei;
    parc_data.(parc_lbl).SI_bin = SI_bin;
    parc_data.(parc_lbl).NPL_wei = NPL_wei;
    parc_data.(parc_lbl).NPL_bin = NPL_bin;

end

fprintf('Done.\n');

%% Model communication from DLPFC TMS sites to the SGC and correlate it to treatment response 

fprintf('Model communication from DLPFC TMS sites to the SGC and correlate it to treatment response...\n');

sgc_mni = [6, 16, -10]; % from https://doi.org/10.1016/j.biopsych.2018.12.002 
r_sphere_mm = 10;

corr_type = 'Spearman';

cohort_vec = {'cohort_1'};
n_cohort_vec = length(cohort_vec);

parc_vec = fields(parc_data);
n_parc_vec = length(parc_vec);

stim_sphere_wei_vec = {'uni', 'em1'}; % kappa = 0, -1
n_stim_sphere_wei_vec = length(stim_sphere_wei_vec);

X_vec = {'SC', 'FC', 'ED', 'CMY_wei', 'CMY_bin', 'SPL_wei', 'SPL_wei_hop', 'SPL_bin', 'SI_wei', 'SI_bin', 'NPL_wei', 'NPL_bin'};
n_X_vec = length(X_vec);

for k_cohort = 1:n_cohort_vec

    for k_parc = 1:n_parc_vec

        cohort_lbl = cohort_vec{k_cohort};
        parc_lbl = parc_vec{k_parc};       

        T_vxl2mni = parc_data.(parc_lbl).T_vxl2mni;
        vxl_size = parc_data.(parc_lbl).vxl_size;
        V = parc_data.(parc_lbl).V;
        stim_xyz = mni2cor(tms_data.(cohort_lbl).tms_coords_mni_adjusted, T_vxl2mni);
        sgc_xyz = mni2cor(sgc_mni, T_vxl2mni);
        imp = tms_data.(cohort_lbl).improvement;
        m = size(stim_xyz, 1);

        % Compute stimulation blocks
        r_stim = r_sphere_mm./vxl_size;
        n_sphere_vxl = length(sphere_vxl(stim_xyz(1,:), r_stim, V));
        stim_block = zeros(m, n_sphere_vxl);
        stim_vxl_idx = zeros(m, n_sphere_vxl);
        stim_sphere_dis = zeros(m, n_sphere_vxl);
        for i = 1:m
            [stim_block(i,:), stim_vxl_idx(i,:), aux] = sphere_vxl(stim_xyz(i,:), r_stim, V);
            aux = aux.*vxl_size;
            stim_sphere_dis(i,:) = aux;
        end

        % Compute target (SGC) block
        r_target = r_sphere_mm./vxl_size;
        target_block = sphere_vxl(sgc_xyz, r_target, V);

        % Exclude subcortical regions from target block
        switch parc_lbl(1:3)
            case 'lau'
                target_block((target_block >= 502 & target_block <= 508) | (target_block >= 1008 & target_block <= 1014)) = 0;
            case 'sch'
                target_block(target_block >= 801) = 0;
            otherwise
                error('Unexpected parcellation');
        end

        tms_model.(parc_lbl).(cohort_lbl).r_sphere_mm = r_sphere_mm;
        tms_model.(parc_lbl).(cohort_lbl).stim_xyz = stim_xyz;
        tms_model.(parc_lbl).(cohort_lbl).sgc_xyz = sgc_xyz;
        tms_model.(parc_lbl).(cohort_lbl).stim_block = stim_block;
        tms_model.(parc_lbl).(cohort_lbl).target_block = target_block;
        tms_model.(parc_lbl).(cohort_lbl).stim_sphere_dis = stim_sphere_dis;


        for k_stim_sphere_wei = 1:n_stim_sphere_wei_vec

            stim_sphere_wei_lbl = stim_sphere_wei_vec{k_stim_sphere_wei};
            stim_sphere_wei = exp(stim_sphere_wei_lbl_to_kappa(stim_sphere_wei_lbl).* stim_sphere_dis);

            % Compute communication from stimulation to targets
            for k_X = 1:n_X_vec
    
                X_lbl = X_vec{k_X};
                X = parc_data.(parc_lbl).(X_lbl);
                x_uni = zeros(m,1);
                x_d2 = zeros(m,1);
    
                for i = 1:m
                    x_uni(i) = block_wei_avg(X, stim_block(i,:), target_block, stim_sphere_wei(i,:));
                end
    
                tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).x = x_uni;
    
                [tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).corr_imp.rho, ...
                 tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).corr_imp.p] = ...
                    corr(tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).x, imp, 'Type', corr_type);
    
                % Leave-one-out correlations
                aux_rho_loo_uni = zeros(m,1);
                aux_p_loo_uni = zeros(m,1);       
                for i = 1:m
                    aux_x_uni = tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).x;
                    aux_imp = imp;
                    aux_x_uni(i) = [];
                    aux_imp(i) = [];
                    [aux_rho_loo_uni(i), aux_p_loo_uni(i)] = corr(aux_x_uni, aux_imp, 'Type', corr_type);
                end
    
                tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).corr_imp_loo.rho = aux_rho_loo_uni;
                tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).corr_imp_loo.p = aux_p_loo_uni;
    
            end

        end

    end

end

fprintf('Done.\n');

%% Model communication from DLPFC TMS sites to the SGC + accumbens + caudate parallel circuit, and correlate it to treatment response 

fprintf('Model communication from DLPFC TMS sites to the parallel circuit and correlate it to treatment response...\n');

sgc_mni = [6, 16, -10]; % from https://doi.org/10.1016/j.biopsych.2018.12.002 
r_sphere_mm = 10;

corr_type = 'Spearman';

cohort_vec = {'cohort_1'};
n_cohort_vec = length(cohort_vec);

parc_vec = fields(parc_data);
n_parc_vec = length(parc_vec);

stim_sphere_wei_vec = {'em1'};
n_stim_sphere_wei_vec = length(stim_sphere_wei_vec);

X_vec = {'SPL_wei_hop'};
n_X_vec = length(X_vec);

for k_cohort = 1:n_cohort_vec

    for k_parc = 1:n_parc_vec

        cohort_lbl = cohort_vec{k_cohort};
        parc_lbl = parc_vec{k_parc};  

        parc_labels = parc_data.(parc_lbl).labels;
        T_vxl2mni = parc_data.(parc_lbl).T_vxl2mni;
        vxl_size = parc_data.(parc_lbl).vxl_size;
        V = parc_data.(parc_lbl).V;
        stim_xyz = mni2cor(tms_data.(cohort_lbl).tms_coords_mni_adjusted, T_vxl2mni);
        sgc_xyz = mni2cor(sgc_mni, T_vxl2mni);
        imp = tms_data.(cohort_lbl).improvement;
        m = size(stim_xyz, 1);

        % Compute stimulation blocks
        r_stim = r_sphere_mm./vxl_size;
        n_sphere_vxl = length(sphere_vxl(stim_xyz(1,:), r_stim, V));
        stim_block = zeros(m, n_sphere_vxl);
        stim_vxl_idx = zeros(m, n_sphere_vxl);
        stim_sphere_dis = zeros(m, n_sphere_vxl);
        for i = 1:m
            [stim_block(i,:), stim_vxl_idx(i,:), aux] = sphere_vxl(stim_xyz(i,:), r_stim, V);
            aux = aux.*vxl_size;
            stim_sphere_dis(i,:) = aux;
        end

        % Compute target (SGC) block
        r_target = r_sphere_mm./vxl_size;
        target_block = sphere_vxl(sgc_xyz, r_target, V);

        switch parc_lbl(1:3)
            
            case 'lau'
                
                % Exclude subcortical regions from target block
                target_block((target_block >= 502 & target_block <= 508) | (target_block >= 1008 & target_block <= 1014)) = 0;

                % Add circuit regions
                target_block = [target_block; sphere_vxl(mni2cor([-6, 16, -10], T_vxl2mni), r_target, V)]; % add left SGC
                aux = nnz(target_block);
                target_block = [target_block; strsea(parc_labels, 'lh_caudate')*ones(aux,1)];  
                target_block = [target_block; strsea(parc_labels, 'lh_accumbensarea')*ones(aux,1)];
                target_block = [target_block; strsea(parc_labels, 'rh_caudate')*ones(aux,1)];  
                target_block = [target_block; strsea(parc_labels, 'rh_accumbensarea')*ones(aux,1)];                            

            case 'sch'

                % Exclude subcortical regions from target block
                target_block(target_block >= 801) = 0;

                % Add circuit regions
                target_block = [target_block; sphere_vxl(mni2cor([-6, 16, -10], T_vxl2mni), r_target, V)]; % add left SGC
                aux = nnz(target_block);
                target_block = [target_block; strsea(parc_labels, 'CAUDATE_LEFT')*ones(aux,1)];  
                target_block = [target_block; strsea(parc_labels, 'ACCUMBENS_LEFT')*ones(aux,1)];
                target_block = [target_block; strsea(parc_labels, 'CAUDATE_RIGHT')*ones(aux,1)];  
                target_block = [target_block; strsea(parc_labels, 'ACCUMBENS_RIGHT')*ones(aux,1)];

            otherwise
                error('Unexpected parcellation');
        end

        tms_model_circuit.(parc_lbl).(cohort_lbl).r_sphere_mm = r_sphere_mm;
        tms_model_circuit.(parc_lbl).(cohort_lbl).stim_xyz = stim_xyz;
        tms_model_circuit.(parc_lbl).(cohort_lbl).sgc_xyz = sgc_xyz;
        tms_model_circuit.(parc_lbl).(cohort_lbl).stim_block = stim_block;
        tms_model_circuit.(parc_lbl).(cohort_lbl).target_block = target_block;
        tms_model_circuit.(parc_lbl).(cohort_lbl).stim_sphere_dis = stim_sphere_dis;


        for k_stim_sphere_wei = 1:n_stim_sphere_wei_vec

            stim_sphere_wei_lbl = stim_sphere_wei_vec{k_stim_sphere_wei};
            stim_sphere_wei = exp(stim_sphere_wei_lbl_to_kappa(stim_sphere_wei_lbl).* stim_sphere_dis);

            % Compute communication from stimulation to targets
            for k_X = 1:n_X_vec
    
                X_lbl = X_vec{k_X};
                X = parc_data.(parc_lbl).(X_lbl);
                x_uni = zeros(m,1);
                x_d2 = zeros(m,1);
    
                for i = 1:m
                    x_uni(i) = block_wei_avg(X, stim_block(i,:), target_block, stim_sphere_wei(i,:));
                end
    
                tms_model_circuit.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).x = x_uni;
    
                [tms_model_circuit.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).corr_imp.rho, ...
                 tms_model_circuit.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).corr_imp.p] = ...
                    corr(tms_model_circuit.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).x, imp, 'Type', corr_type);
    
                % Leave-one-out correlations
                aux_rho_loo_uni = zeros(m,1);
                aux_p_loo_uni = zeros(m,1);       
                for i = 1:m
                    aux_x_uni = tms_model_circuit.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).x;
                    aux_imp = imp;
                    aux_x_uni(i) = [];
                    aux_imp(i) = [];
                    [aux_rho_loo_uni(i), aux_p_loo_uni(i)] = corr(aux_x_uni, aux_imp, 'Type', corr_type);
                end
    
                tms_model_circuit.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).corr_imp_loo.rho = aux_rho_loo_uni;
                tms_model_circuit.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).corr_imp_loo.p = aux_p_loo_uni;

            end

        end

    end

end

fprintf('Done.\n');

%% Compute R-maps between WM hops from the DLPFC and treatment response 

%%% Load pre-computed results
fprintf('Load pre-computed R-maps between WM hops from the DLPFC and treatment response...\n');
load('rmap_SPL_imp.mat');

%%% Uncomment below to run analyses. Adjust number of parfor workers as
%%% desired. Run time is approx 20 minutes using 10 parallel workers.

% fprintf('Compute R-maps between WM hops from the DLPFC and treatment response...\n');
% 
% cohort_vec = {'cohort_1'};
% n_cohort_vec = length(cohort_vec);
% 
% parc_vec = fields(parc_data);
% n_parc_vec = length(parc_vec);
% 
% stim_sphere_wei_vec = {'em1'};
% n_stim_sphere_wei_vec = length(stim_sphere_wei_vec);
% 
% X_vec = {'SPL_wei_hop'};
% n_X_vec = length(X_vec);
% 
% vxl_skip = 1;
% 
% for k_cohort = 1:n_cohort_vec
% 
%     for k_parc = 1:n_parc_vec
% 
%         cohort_lbl = cohort_vec{k_cohort};
%         parc_lbl = parc_vec{k_parc};       
% 
%         stim_block = tms_model.(parc_lbl).(cohort_lbl).stim_block;
%         stim_sphere_dis = tms_model.(parc_lbl).(cohort_lbl).stim_sphere_dis;
%         n_sphere_vxl = length(stim_sphere_dis);
%         V = parc_data.(parc_lbl).V;
%         imp = tms_data.(cohort_lbl).improvement;
%         m = size(imp, 1);
% 
%         V_idx = find(V);
%         [vxl_x, vxl_y, vxl_z] = ind2sub(size(V), V_idx);
%         n_vxl = size(vxl_x,1);
%         vxl_rho_x_imp = nan(n_vxl, 1);
% 
%         for k_stim_sphere_wei = 1:n_stim_sphere_wei_vec
% 
%             stim_sphere_wei_lbl = stim_sphere_wei_vec{k_stim_sphere_wei};
%             stim_sphere_wei = exp(stim_sphere_wei_lbl_to_kappa(stim_sphere_wei_lbl).* stim_sphere_dis);
% 
%             for k_X = 1:n_X_vec
%     
%                 X_lbl = X_vec{k_X};
%                 X = parc_data.(parc_lbl).(X_lbl);
%     
%                 % for k_vxl = 1:vxl_skip:n_vxl
%                 parfor (k_vxl = 1:n_vxl, 10)
%         
%                     if mod(k_vxl,100) == 0
%                         fprintf('%s (%d/%d), %s (%d/%d), %s (%d/%d), (%d/%d) %.2f%% complete (skip = %d)\n', ...
%                             cohort_lbl, k_cohort, n_cohort_vec, ...
%                             parc_lbl, k_parc, n_parc_vec, ...
%                             X_lbl, k_X, n_X_vec, ...
%                             k_vxl, n_vxl, 100*k_vxl/n_vxl, vxl_skip);
%                     end
%                 
%                     vxl_target_block = sphere_vxl([vxl_x(k_vxl), vxl_y(k_vxl), vxl_z(k_vxl)], r_target, V);
% 
%                     % Exclude voxels from target block
%                     region_idx = V(vxl_x(k_vxl), vxl_y(k_vxl), vxl_z(k_vxl));
%                     switch parc_lbl(1:3)
%                         case 'lau'
%                             ctx_idx = [1:501, 509:1007];
%                             sub_ctx_idx = [502:508, 1008:1014];
%                         case 'sch'
%                             ctx_idx = 1:800;
%                             sub_ctx_idx = 801:814;
%                     end
%                     if isempty(find(ctx_idx == region_idx, 1)) % sphere centred at subcortical voxel
%                         vxl_target_block(ismember(vxl_target_block, ctx_idx)) = 0; % exclude cortex from target
%                     else % sphere centred at cortical voxel
%                         vxl_target_block(ismember(vxl_target_block, sub_ctx_idx)) = 0; % exclude subcortex from target
%                     end
% 
%                     x_stim_target = zeros(m,1);
%                     
%                     for i = 1:m
%                         x_stim_target(i) = block_wei_avg(X, stim_block(i,:), vxl_target_block, stim_sphere_wei(i,:));
%                     end
%                     
%                     vxl_rho_x_imp(k_vxl) = corr(x_stim_target, imp, 'Type', 'S');
%                 
%                 end
% 
%                 V_rho_x_imp = zeros(size(V));
%                 V_rho_x_imp(V ~= 0) = vxl_rho_x_imp;
% 
%                 rmap.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).vxl = vxl_rho_x_imp;
%                 rmap.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).V = V_rho_x_imp;
%                 rmap.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).vxl_skip = vxl_skip;
%           
%             end
% 
%         end
% 
%     end
% 
% end
% 
% %%% Compute exclusion struct
% for k_parc = 1:n_parc_vec
% 
%     parc_lbl = parc_vec{k_parc};
% 
%     V = parc_data.(parc_lbl).V;
%     V_idx = find(V);
%     [vxl_x, vxl_y, vxl_z] = ind2sub(size(V), V_idx);
%     n_vxl = size(vxl_x,1);
% 
%     parc_labels = parc_data.(parc_lbl).labels;
% 
%     % combine exc regions across both cohorts
%     exc_region_idx = [];
%     for k_cohort = 1:n_cohort_vec
% 
%         cohort_lbl = cohort_vec{k_cohort};
% 
%         [~, ~, ~, sgc_mask] = sphere_vxl(tms_model.(parc_lbl).(cohort_lbl).sgc_xyz, r_target, parc_data.(parc_lbl).V);
%         sgc_mask_vxl = find(sgc_mask);  
% 
%         stim_block = tms_model.(parc_lbl).(cohort_lbl).stim_block;
%         stim_block_labels = parc_labels(unique(stim_block(stim_block>0)));
%         exc_labels = cell(length(stim_block_labels),1);
%         for i = 1:length(stim_block_labels)
%             switch parc_lbl
%                 case 'lau_l5'
%                     aux = strsplit(stim_block_labels{i}, '_');
%                     exc_labels{i} = aux{2};
%                 case 'sch_800'
%                     aux = strsplit(stim_block_labels{i}, '_');
%                     exc_labels{i} = strjoin(aux(2:end-1), '_');
%                 otherwise
%                     error('Invalid parc_lbl');
%             end
%         end
%         exc_labels = unique(exc_labels);
%         exc_region_idx = [exc_region_idx; strsea(parc_labels, exc_labels)];
% 
%     end
%     exc_region_idx = unique(exc_region_idx); % remove regions repated across cohorts
% 
%     dlpfc_mask = zeros(size(V));
%     for i = 1:length(exc_region_idx)
%         dlpfc_mask(V == exc_region_idx(i)) = 1;
%     end
%     dlpfc_mask_vxl = find(dlpfc_mask);
%     exc_mask_vxl = union(dlpfc_mask_vxl, sgc_mask_vxl);
%     
%     vxl_exc_vec = [];
%     parfor (k_vxl = 1:n_vxl, 10)
% 
%         if mod(k_vxl,100) == 0
%             fprintf('%s (%d/%d), (%d/%d) %.2f%% complete\n', ...
%                 parc_lbl, k_parc, n_parc_vec, ...
%                 k_vxl, n_vxl, 100*k_vxl/n_vxl);
%         end
%     
%         [~, ~, ~, V_target_block] = sphere_vxl([vxl_x(k_vxl), vxl_y(k_vxl), vxl_z(k_vxl)], r_target, V);
% 
%         vxl_V_target_block = find(V_target_block);
% 
%         if any(ismember(vxl_V_target_block, exc_mask_vxl))
%             vxl_exc_vec = [vxl_exc_vec; k_vxl];
%         end
%     
%     end
%     
%     exc_vxl_struct.(parc_lbl).vxl_exc_vec = vxl_exc_vec;
% 
% end    
% 
% %%% Create new rmap struct without the excluded voxels
% rmap_exc_vxl = rmap;
% 
% for k_cohort = 1:n_cohort_vec
%     for k_parc = 1:n_parc_vec
%         for k_stim_sphere_wei = 1:n_stim_sphere_wei_vec
%             for k_X = 1:n_X_vec
% 
%                 cohort_lbl = cohort_vec{k_cohort};
%                 parc_lbl = parc_vec{k_parc};
%                 stim_sphere_wei_lbl = stim_sphere_wei_vec{k_stim_sphere_wei};
%                 X_lbl = X_vec{k_X};
% 
%                 switch parc_lbl(1:3)
%                     case 'lau'
%                         vxl_exc_vec = exc_vxl_struct.lau_l5.vxl_exc_vec;
%                         V = parc_data.lau_l5.V;
%                     case 'sch'
%                         vxl_exc_vec = exc_vxl_struct.sch_800.vxl_exc_vec;
%                         V = parc_data.sch_800.V;
%                     otherwise
%                         error('Unexpected parcellation');
%                 end
% 
%                 vxl_rho_x_imp = rmap_exc_vxl.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).vxl;
%                 % vxl_rho_x_imp(vxl_exc_vec) = [];
%                 vxl_rho_x_imp(vxl_exc_vec) = nan;
%                 
%                 rmap_exc_vxl.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).vxl = vxl_rho_x_imp;
% 
%                 V_rho_x_imp_exc_vxl = zeros(size(V));
%                 V_rho_x_imp_exc_vxl(V ~= 0) = vxl_rho_x_imp;
%                 rmap_exc_vxl.(parc_lbl).(cohort_lbl).(X_lbl).(stim_sphere_wei_lbl).V = V_rho_x_imp_exc_vxl;
% 
%             end
%         end
%     end
% end

fprintf('Done.\n');

%% Compute DLPFC maps

% Load MNI brain template edges
brain_edges = niftiread('MNI152_T1_2mm_edges.nii.gz');
% Load DLPFC mask
dlpfc_mask = niftiread('DLPFCcombo20mm.nii');

parc_vec = fields(parc_data);
n_parc_vec = length(parc_vec);

X_vec = {'FC', 'SPL_wei_hop'};
n_X_vec = length(X_vec);

for k_parc = 1:n_parc_vec

    parc_lbl = parc_vec{k_parc};

    target_block = tms_model.(parc_lbl).cohort_1.target_block;
    V = parc_data.(parc_lbl).V;
    se = strel('sphere', 1);
    brain_edges_dil = imdilate(brain_edges, se.Neighborhood);
    tms_mask = dlpfc_mask & brain_edges_dil & ~~V;
    tms_mask_idx = find(tms_mask);
    n_tms_mask_idx = length(tms_mask_idx);

    for k_X = 1:n_X_vec          

        X_lbl = X_vec{k_X};
        X = parc_data.(parc_lbl).(X_lbl);

        V_x_sgc = nan(size(tms_mask));
        
        for i = 1:n_tms_mask_idx

            fprintf('%s (%d/%d), %s (%d/%d), (%d/%d) %.2f%% complete\n', ...
            parc_lbl, k_parc, n_parc_vec, ...
            X_lbl, k_X, n_X_vec, ...
            i, n_tms_mask_idx, 100*i/n_tms_mask_idx);  

            [x,y,z] = ind2sub(size(V_x_sgc), tms_mask_idx(i));

            [stim_sphere, ~, aux] = sphere_vxl([x,y,z], r_stim, V);
            stim_sphere_dis = aux.*vxl_size;
            n_sphere_vxl = length(stim_sphere_dis);

            stim_sphere_wei = exp(0.* stim_sphere_dis);

            V_x_sgc(x,y,z) = block_wei_avg(X, stim_sphere, target_block, stim_sphere_wei);
        
        end

        dlpfc_maps.(parc_lbl).(X_lbl).V = V_x_sgc;
        dlpfc_maps.(parc_lbl).(X_lbl).tms_mask = tms_mask;

    end

end

%% Plot Figure 1B,C,D,E

% --- Result parameters --- %
    
parc_lbl = 'lau_l5';
stim_sphere_wei_lbl = 'em1';
com_lbl = 'SPL_wei_hop';

parc_xyz = parc_data.(parc_lbl).coords_xyz;

% --- Figure parameters --- %

font_size = 18;

pc = 1;

clear fig panel;

fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
fig.Position = [1 1 1349 947];
set(fig, 'AutoResizeChildren', 'off');

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 1: Cohort 1
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

cohort_lbl = 'cohort_1';

stim_xyz = tms_model.(parc_lbl).(cohort_lbl).stim_xyz; % tms sites
sgc_xyz = tms_model.(parc_lbl).(cohort_lbl).sgc_xyz; % sgc
imp = tms_data.(cohort_lbl).improvement; % clinical outcome
spl_stim_target = tms_model.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).x; % communication
vxl_rho = rmap_exc_vxl.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).vxl;  % rmap vector
vxl_rho(isnan(vxl_rho)) = [];
rho_rh_sgc = tms_model.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).corr_imp.rho; % RH SGC corr
rho_lh_sgc = rmap.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).V(48,71,31); % LH SGC corr
rmap_rh_sgc_perc = 100*nnz(vxl_rho <= rho_rh_sgc)./numel(vxl_rho); % RH SGC Rmap position
rmap_lh_sgc_perc = 100*nnz(vxl_rho <= rho_lh_sgc)./numel(vxl_rho); % LH SGC Rmap position

offset = 0;

% ------------------------------------------------------------------------------ %
% Cohort 1, brain plot + TMS coordiantes + SGC
% ------------------------------------------------------------------------------ %

% Brain coordinates 1
panel{pc} = subplot(2,5,1+offset);
pc = pc + 1;
scatter(parc_xyz(:,1), parc_xyz(:,2), 10, 'black', 'MarkerEdgeAlpha', 0.5); hold on;
scatter(stim_xyz(:,1), stim_xyz(:,2), 15, 'filled', 'blue'); hold on;
scatter(sgc_xyz(:,1), sgc_xyz(:,2), 30, 'filled', 'red'); hold on;
% scatter(sgc_xyz(:,1), sgc_xyz(:,2), 250, 'red', 'LineWidth', 1); hold on;
set(gca, 'XDir', 'reverse');
axis equal;
axis off;

% Brain coordinates 2
panel{pc} = subplot(2,5,2+offset);
pc = pc + 1;
scatter(parc_xyz(:,1), parc_xyz(:,3), 10, 'black', 'MarkerEdgeAlpha', 0.5); hold on;
scatter(stim_xyz(:,1), stim_xyz(:,3), 15, 'filled', 'blue'); hold on;
scatter(sgc_xyz(:,1), sgc_xyz(:,3), 30, 'filled', 'red'); hold on;
% scatter(sgc_xyz(:,1), sgc_xyz(:,3), 250, 'red', 'LineWidth', 1); hold on;
set(gca, 'XDir', 'reverse');
axis equal;
axis off;

% ------------------------------------------------------------------------------ %
% Cohort 1, scatter com vs improvement
% ------------------------------------------------------------------------------ %

panel{pc} = subplot(2,5,3+offset);
pc = pc + 1;
sct = scatter(spl_stim_target, imp, 'o', 'filled');
[rho_rh_sgc, p_rh_sgc] = corr(spl_stim_target, imp, 'Type', 'S');
sct.SizeData = 75;
sct.MarkerEdgeColor = 'k';
sct.LineWidth = 0.75;
xlabel({'WM hops from DLPFC to SGC'});
ylabel('Treatment response (%)');
set(gca, 'FontSize', font_size);
title(sprintf('r = %.2f, p = %.4f', rho_rh_sgc, p_rh_sgc/2), 'FontWeight', 'normal', 'FontSize', font_size+1);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
grid minor;
adjust_scatter_limits(panel{pc-1}, 0.1);

% ------------------------------------------------------------------------------ %
% Cohort 1, Rmap cortical projection
% ------------------------------------------------------------------------------ %

panel{pc} = subplot(2,5,4+offset);
pc = pc + 1;
img_path = 'pre_computed_results/img/rmap_fig_1d.png';
I = imread(img_path);
imshow(I(1:2400,:,:));
title({'WM hops from DLPFC vs. treatment response'}, 'FontWeight', 'normal', 'FontSize', font_size+1);

% ------------------------------------------------------------------------------ %
% Cohort 1, Rmap histogram + RH SGC position
% ------------------------------------------------------------------------------ %

% Histogram
panel{pc} = subplot(2,5,5+offset);
pc = pc + 1;
n_bins = 18;
histogram(vxl_rho, n_bins, 'EdgeColor', 'none'); hold on;
histogram(vxl_rho, n_bins, 'DisplayStyle', 'stairs', 'LineWidth', 1, 'EdgeColor', 'k'); hold on;
histogram(vxl_rho, n_bins, 'DisplayStyle', 'bar', 'LineWidth', 0.25, 'EdgeColor', 'k', 'FaceColor', 'none'); hold on;
xlabel('Correlation');
ylabel('Number of voxels');
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
set(gca, 'YTick', [0 4000 8000 12000]);
% adjust_histogram_limits(panel{pc-1}, 0.1);
xlim([-1, 1]);
ylim([0,13000]);
set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1]);
set(gca, 'XTickLabelRotation', 0);
panel{pc-1}.YAxis.Exponent = 3;

% Add LH SGC arrow
text(panel{pc-1},  0.35, 11000, ...
    {'LH SGC', sprintf('r = %.2f', rho_lh_sgc), sprintf('top %.0f%%', rmap_lh_sgc_perc)}, ...
    'Color', 'k', 'FontSize', font_size-4, 'HorizontalAlignment','left');
line(panel{pc-1}, [rho_lh_sgc, rho_lh_sgc], [0 7000], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 0.75);
line(panel{pc-1}, [rho_lh_sgc, 0.6], [7000 7000], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 0.75);
line(panel{pc-1}, [0.6, 0.6], [7000 9000], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 0.75);
scatter(panel{pc-1}, 0.6, 9000, '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(panel{pc-1}, rho_lh_sgc, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

% Add RH SGC arrow
text(panel{pc-1}, -0.9, 11000, ...
    {'RH SGC', sprintf('r = %.2f', rho_rh_sgc), sprintf('top %.0f%%', rmap_rh_sgc_perc)}, ...
    'Color', 'k', 'FontSize', font_size-4, 'HorizontalAlignment','left');
line(panel{pc-1}, [rho_rh_sgc, rho_rh_sgc], [0 7000], 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.75);
line(panel{pc-1}, [rho_rh_sgc, -0.65], [7000 7000], 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.75);
line(panel{pc-1}, [-0.65, -0.65], [7000 9000], 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.75);
scatter(panel{pc-1}, -0.65, 9000, '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
scatter(panel{pc-1}, rho_rh_sgc, 0, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 3: Adjustments
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

% [left bottom width height]

% for i = 1:10    
%     panel{i}.PlotBoxAspectRatio = [1 1 1];    
% end

% scatters
panel{3}.PlotBoxAspectRatio = [1 1 1];     

% hisotgrams
panel{5}.PlotBoxAspectRatio = [1 1 1];    

% left
% panel{1}.Position(1) = 0;
% panel{2}.Position(1) = 0.125;

position_offset_1 = [0.025 0.05 0 0];
panel{1}.Position = [0.0250    0.6400    0.1000    0.3412] + position_offset_1;
panel{2}.Position = [0.0250    0.4450    0.1000    0.3412] + position_offset_1;
panel{3}.Position = [0.2000    0.5838    0.1600    0.3412] + position_offset_1;
panel{4}.Position = [0.2200    0.5100    0.6000    0.3750] + position_offset_1;
panel{5}.Position = [0.6800    0.5838    0.1600    0.3412] + position_offset_1; % panel{5}.Position = [0.735    0.5838    0.1600    0.3412] + position_offset_1;

axes(panel{1});
text(80,107, 'Cohort I (N = 25)', 'FontWeight', 'normal', 'FontSize', font_size+3);

axes(panel{2});
text(80,0, {'\color{blue} \bullet \color{black}  DLPFC TMS sites','\color{red} \bullet \color{black} SGC seed'}, ...
    'Interpreter', 'tex', 'FontWeight', 'normal', 'FontSize', font_size);

%% Plot Figure 2A,B,C

% --- Result parameters --- %   

parc_lbl = 'lau_l5';
stim_sphere_wei_lbl = 'em1';

% --- Figure parameters --- %

clear path_table;

font_size = 19;
font_size_big = 22;

top_edge_method = 'cumshare'; top_edge_val = 1.1;

max_label_merger = 10;
merge_label_contribution_thr = 0.05;

% mediator_size_multiplier = 150;
mediator_size_multiplier = 100;

fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
% fig.Position = [1 1 1385 803];
fig.Position = [1 33 1473 803];
set(fig, 'AutoResizeChildren', 'off');

% --- Regular expressions to fix labels --- %

switch parc_lbl(1:3)
    case 'lau'
        process_parc_labels = @(s) regexprep(...
            regexprep(strtrim(strrep(strrep(strrep(...
            regexprep(regexprep(s, '^_|_$', ''), '__', ', '), ... % Replace "__" with Unicode right arrow
            '_', ' '), ...
            'lh', 'LH'), ...
            'rh', 'RH')), ...
            '([a-z])([A-Z])', '$1 $2'), ... % Add spaces between lowercase and uppercase transitions
            ' (\d+)', ' ($1)');         
    case 'sch'
        process_parc_labels = @(s) regexprep(regexprep(regexprep(regexprep(regexprep(regexprep(regexprep(regexprep(regexprep( ...
            regexprep(regexprep(regexprep(regexprep(regexprep(regexprep(regexprep( ...
            regexprep(regexprep(s,'7Networks_LH_','L '),'7Networks_RH_','R '),'Default',''), ...
            'Cont',''),'SomMot',''),'DorsAttn',''),'SalVentAttn',''),'Limbic',''), ...
            'Vis',''),'__',', '),'_',' '),'CIFTI_STRUCTURE_',''),'LEFT','L'),'RIGHT','R'), ...
            'PFCdPFCm','PFCdm'), 'CIFTI STRUCTURE',''), '  ',' '), ...
            '\<(\d+)\>','($1)');

end    

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 1: Cohort 1
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

cohort_lbl = 'cohort_1';

% Get info from tms model and parc
stim_block = tms_model.(parc_lbl).(cohort_lbl).stim_block;
target_block = tms_model.(parc_lbl).(cohort_lbl).target_block;
parc_labels = parc_data.(parc_lbl).labels;
parc_xyz = parc_data.(parc_lbl).coords_xyz;
m = size(stim_block,1);

SPL_h = parc_data.(parc_lbl).SPL_wei_hop;
SPL_Pmat = parc_data.(parc_lbl).SPL_wei_Pmat;
n = size(SPL_h,1);

% Stim sphere wei info
stim_sphere_dis = tms_model.(parc_lbl).(cohort_lbl).stim_sphere_dis;
stim_sphere_wei = exp(stim_sphere_wei_lbl_to_kappa(stim_sphere_wei_lbl).*stim_sphere_dis);

% Initialise tables and variables
stim_target_SP = table('Size', [0, 2], 'VariableTypes', {'double', 'double'}, 'VariableNames', {'path', 'count'});
stim_target_SP_str_labels = table('Size', [0, 2], 'VariableTypes', {'string', 'double'}, 'VariableNames', {'path', 'count'});
stim_target_SP_str_labels_inter = table('Size', [0, 2], 'VariableTypes', {'string', 'double'}, 'VariableNames', {'path', 'count'});
stim_target_SP_hops = [];

tab_stim = tabulate(flat(stim_block));
stim_block_unique = tab_stim(2:end, 1);
stim_block_count = tab_stim(2:end, 2);
tab_target = tabulate(flat(target_block));
target_block_unique = tab_target(2:end, 1);
target_block_count  = tab_target(2:end, 2);

% ==== weighted source multiplicities (replaces the tabulate on stim) ====
% stim_weights must be m x 485, aligned with stim_block (phi(d) per voxel).
flatStim = flat(stim_block);
flatW    = flat(stim_sphere_wei); 
mask     = flatStim > 0;
flatStim = flatStim(mask);
flatW    = flatW(mask);

% Sum weights by parcel label (1..n). 'n' is #parcels (size(SPL_h,1)).
stim_weight_sum = accumarray(flatStim, flatW, [n 1], @sum, 0);

% Keep only parcels with nonzero total weight
stim_block_unique_wei = find(stim_weight_sum > 0);
stim_block_count_wei  = stim_weight_sum(stim_block_unique_wei);  % NOTE: now "counts" == weighted sums
% ============================================================================

% Count paths
for i = 1:length(stim_block_unique)
    for j = 1:length(target_block_unique)

        stim_idx = stim_block_unique(i);
        target_idx = target_block_unique(j);

        n_stim_taget = stim_block_count_wei(i)*target_block_count(j);

        path = retrieve_shortest_path(stim_idx, target_idx, SPL_h, SPL_Pmat);
        n_hops = length(path);

        path_string = [];
        path_string_labels = [];
        path_string_inter = [];
        path_string_labels_inter = [];

        for k = 1:n_hops
            
            path_string = [path_string, '_', num2str(path(k)), '_'];
            path_string_labels = [path_string_labels, '_', parc_labels{path(k)}, '_'];

            if k ~= 1 && k ~= n_hops % intermdiate regions
                path_string_inter = [path_string_inter, '_', num2str(path(k)), '_'];
                path_string_labels_inter = [path_string_labels_inter, '_', parc_labels{path(k)}, '_'];
            end

        end

        stim_target_SP = add_to_path_table(stim_target_SP, path, n_stim_taget);
        stim_target_SP_str_labels = add_to_path_table(stim_target_SP_str_labels, path_string_labels, n_stim_taget);
        if n_hops > 2
            stim_target_SP_str_labels_inter = add_to_path_table(stim_target_SP_str_labels_inter, path_string_labels_inter, n_stim_taget);
        end

        stim_target_SP_hops = [stim_target_SP_hops; repmat(n_hops, stim_block_count(i)*target_block_count(j), 1)];

    end
end

% Normalise and sort counts
stim_target_SP_str_labels_inter.percentage = stim_target_SP_str_labels_inter.count./sum(stim_target_SP_str_labels_inter.count);
stim_target_SP_str_labels_inter = sortrows(stim_target_SP_str_labels_inter, 'count', 'descend');

stim_target_SP_str_labels_inter_top = merge_paths_by_region_names_thr_adj(...
    stim_target_SP_str_labels_inter, merge_label_contribution_thr, true, max_label_merger, ...
    full(build_label_adjacency_from_volume(parc_data.(parc_lbl).V)), parc_labels, true);

% Store path table
path_table.(parc_lbl).(cohort_lbl).stim_target_SP_str_labels_inter = stim_target_SP_str_labels_inter;
path_table.(parc_lbl).(cohort_lbl).stim_target_SP_str_labels_inter_top = stim_target_SP_str_labels_inter_top;

% Edge usage (top edges) -- used for plotting
stim_target_SP_top = filter_paths_with_coverage_ids(stim_target_SP, top_edge_method, top_edge_val);
E_top = zeros(n);
for i = 1:size(stim_target_SP_top,1)
    path = stim_target_SP_top.path{i};
    n_hops = length(path);
    for k = 1:(n_hops-1)
        E_top(path(k), path(k+1)) = E_top(path(k), path(k+1)) + stim_target_SP_top.count(i);
    end        
end

% Edge usage (all edges) -- used to compute mediators
E_all = zeros(n);
for i = 1:size(stim_target_SP,1)
    path = stim_target_SP.path{i};
    n_hops = length(path);
    for k = 1:(n_hops-1)
        E_all(path(k), path(k+1)) = E_all(path(k), path(k+1)) + stim_target_SP.count(i);
    end        
end

e_all = sum(E_all,2);
e_all(stim_block_unique) = 0;
e_all(target_block_unique) = 0;
warning('Excluding stim and target regions from mediators list');
mediator_regions = table(parc_labels, e_all, e_all/sum(e_all), 'VariableNames', {'region', 'count', 'percentage'});

% Construct graph obj for plotting
G_com = digraph(E_top);
stim_node_size = 15.*ones(n,1);
target_node_size = 15.*ones(n,1);
mediator_node_size = mediator_size_multiplier.*(mediator_regions.count./max(mediator_regions.count));
mediator_node_size(sum(E_top,2) == 0) = 0; % remove from graph plot but keep in bar plot
idx_stim_nodes = unique(stim_block(stim_block > 0)); % find stim nodes
idx_target_nodes = unique(target_block(target_block > 0)); % find target nodes
idx_mediator_nodes = mediator_node_size > 0;
stim_node_size(idx_stim_nodes == 0) = nan;
target_node_size(idx_target_nodes == 0) = nan;
mediator_node_size(idx_mediator_nodes == 0) = nan;
regular_node_size = 10*ones(n,1);
regular_node_size((stim_node_size + target_node_size + mediator_node_size) > 0) = nan;

%%% Plots

pc = 1;
offset = 0;
clear panel;

% Brain 1
panel{pc} = subplot(2,4,1+offset); pc = pc + 1;
g_com = plot(G_com, 'XData', parc_xyz(:,1), 'YData', parc_xyz(:,2), 'LineWidth', 2); hold on;
g_com.Marker = 'none';
g_com.EdgeAlpha = 0.75;
g_com.EdgeCData = log10(G_com.Edges.Weight);   
g_com.LineWidth = 1.5;
colormap(othercolor('BuPu3'));
cb_1 = colorbar('Location', 'south');
cb_1.Label.String = 'Edge contribution';
cb_1.Ticks = [-2, 3];
cb_1.TickLabels = {'10^{-2}', '10^{3}'};
cb_1.AxisLocation = 'out';
scatter(parc_xyz(:,1), parc_xyz(:,2), regular_node_size, 'black', 'MarkerEdgeAlpha', 0.6); hold on;
scatter(parc_xyz(:,1), parc_xyz(:,2), mediator_node_size, ...
    'MarkerFaceColor', [0.4940 0.1840 0.5560], ...
    'MarkerEdgeColor', 'k', ...
    'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5); hold on;
scatter(parc_xyz(idx_stim_nodes,1), parc_xyz(idx_stim_nodes,2), 15, ...
    'MarkerFaceColor', 'blue', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1); hold on;
scatter(parc_xyz(idx_target_nodes,1), parc_xyz(idx_target_nodes,2), 15, ...
    'MarkerFaceColor', 'red', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1); hold on;
axis equal;
axis off;

% Brain 2
panel{pc} = subplot(2,4,2+offset); pc = pc + 1;
g_com = plot(G_com, 'XData', parc_xyz(:,2), 'YData', parc_xyz(:,3), 'LineWidth', 2); hold on;
g_com.Marker = 'none';
g_com.EdgeAlpha = 0.75;
g_com.EdgeCData = log10(G_com.Edges.Weight);  
g_com.LineWidth = 1.5;
colormap(othercolor('BuPu3'));
scatter(parc_xyz(:,2), parc_xyz(:,3), regular_node_size, 'black', 'MarkerEdgeAlpha', 0.6); hold on;
scatter(parc_xyz(:,2), parc_xyz(:,3), mediator_node_size, ...
    'MarkerFaceColor', [0.4940 0.1840 0.5560], ...
    'MarkerEdgeColor', 'k', ...
    'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1.5); hold on;
scatter(parc_xyz(idx_stim_nodes,2), parc_xyz(idx_stim_nodes,3), 15, ...
    'MarkerFaceColor', 'blue', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1); hold on;
scatter(parc_xyz(idx_target_nodes,2), parc_xyz(idx_target_nodes,3), 15, ...
    'MarkerFaceColor', 'red', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerEdgeAlpha', 1, ...
    'LineWidth', 1); hold on;
axis equal;
axis off;

% Histogram
panel{pc} = subplot(2,4,3+offset); pc = pc + 1;
histogram(stim_target_SP_hops-1, 'LineWidth', 1);
xlabel('Hops');
ylabel('Number of paths');
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
adjust_histogram_limits(panel{pc-1}, 0.1);

% % Bar plot
panel{pc} = subplot(2,4,4+offset); pc = pc + 1;
n_top = min(5, size(stim_target_SP_str_labels_inter_top,1));
% Top mediators
[~, sidx] = sort(mediator_regions.count, 'descend');
top_mediators_val = flipud(100*mediator_regions.percentage(sidx(1:n_top)));
top_mediators_lbl = flipud(parc_labels(sidx(1:n_top)));
% Top paths
top_paths_val = flipud(100*stim_target_SP_str_labels_inter_top.percentage(1:n_top));
top_paths_lbl = flipud(stim_target_SP_str_labels_inter_top.path(1:n_top));
% Bar plot
createBarFigure_tex(top_mediators_val, top_paths_val, process_parc_labels(top_mediators_lbl), ...
    process_parc_labels(top_paths_lbl), [10, 20, 30], font_size-1, font_size-6);

% Mediator and path contributions
top_mediator_contribution_1 = round(sum(top_mediators_val));
top_path_contribution_1 = round(sum(top_paths_val));

%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% %%% ADJUSTMENTS %%% %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%

% [left bottom width height]

panel{1}.XDir = 'reverse';
panel{2}.XDir = 'reverse';
panel{5}.XDir = 'reverse';
panel{6}.XDir = 'reverse';

for i = 1:4
    panel{i}.PlotBoxAspectRatio = [1 1 1];    
end

% Bar plots
panel{4}.PlotBoxAspectRatio = [2.25 1 1];

% Histograms
panel{3}.PlotBoxAspectRatio = [1 2.5 1];

panel{1}.Position = [-0.175   0.6500    0.5000    0.32];
panel{2}.Position = [0.025  0.6500    0.4000    0.32];
panel{3}.Position = [0.3200    0.7162    0.1023    0.2000];
panel{4}.Position = [0.485   0.6325    0.200    0.3412];

axes(panel{1});
text(0, panel{1}.YLim(2)*0.95, 'Cohort I', ...
    'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'FontSize', font_size_big-1)

axes(panel{4}); text(-10, 7.3, {'Top 5 regions', sprintf('(%d%% contribution)', top_mediator_contribution_1)}, ...
    'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'FontSize', font_size);
axes(panel{4}); text(10, 7.3, {'Top 5 paths', sprintf('(%d%% contribution)', top_path_contribution_1)}, ...
    'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'FontSize', font_size);

for i = 1:4
    panel{i}.FontSize = font_size;
end

cb_1.Position = [0.1280    0.6734    0.0600    0.0075];
cb_1.Label.Position = [0.8396    4.4077         0];

cb_1.FontSize = font_size-3;



%% Plot Figure 3

% --- Result parameters --- %

com_lbl = 'SPL_wei_hop';

img_x_idx = 600:2330;
img_y_idx = 550:2675;

% --- Figure parameters --- %

font_size = 16;

clear fig panel;

fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
fig.Position = [1 1 1349 947];
set(fig, 'AutoResizeChildren', 'off');

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 1: LAU
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

parc_lbl = 'lau_l5';

V_spl_sgc = dlpfc_maps.(parc_lbl).(com_lbl).V;
V_fc_sgc = dlpfc_maps.(parc_lbl).FC.V;
tms_mask = dlpfc_maps.(parc_lbl).FC.tms_mask;
cb_limit_1_min = min(V_spl_sgc(:));
cb_limit_1_max = max(V_spl_sgc(:));
cb_limit_2_min = -max(abs(V_fc_sgc(:))); % should be centred around 0
cb_limit_2_max = max(abs(V_fc_sgc(:)));

panel{1} = subplot(2,3,1);
I = imread('pre_computed_results/img/dlpfc_map_lau_l5_spl_fig_3a.png');
imshow(I(img_x_idx, img_y_idx,:));
title({'Lausanne normative connectome', 'DLPFC-SGC WM hops'}, ...
    'FontWeight', 'normal', 'FontSize', font_size+1);

panel{2} = subplot(2,3,2);
I = imread('pre_computed_results/img/dlpfc_map_lau_l5_fc_fig_3b.png');
imshow(I(img_x_idx, img_y_idx,:));
title('DLPFC-SGC FC', 'FontWeight', 'normal', 'FontSize', font_size+1);

panel{3} = subplot(2,3,3);
sct = scatter(V_spl_sgc(tms_mask==1),V_fc_sgc(tms_mask==1));
[rho,p] = corr(V_spl_sgc(tms_mask==1),V_fc_sgc(tms_mask==1), 'Type', 'S');
% sct.SizeData = 75;
% sct.MarkerEdgeColor = 'k';
% sct.LineWidth = 0.75;
xlabel('DLPFC-SGC WM hops');
ylabel('DLPFC-SGC FC');
set(gca, 'FontSize', font_size);
title(sprintf('r = %.2f, p < 10^{-12}', rho), 'FontWeight', 'normal', 'FontSize', font_size+1);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
grid minor;
axis square;
adjust_scatter_limits(panel{3}, 0.1);

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 2: SCH
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

parc_lbl = 'sch_800';

V_spl_sgc = dlpfc_maps.(parc_lbl).(com_lbl).V;
V_fc_sgc = dlpfc_maps.(parc_lbl).FC.V;
tms_mask = dlpfc_maps.(parc_lbl).FC.tms_mask;
cb_limit_3_min = min(V_spl_sgc(:));
cb_limit_3_max = max(V_spl_sgc(:));
cb_limit_4_min = -max(abs(V_fc_sgc(:))); % should be centred around 0
cb_limit_4_max = max(abs(V_fc_sgc(:)));

panel{4} = subplot(2,3,4);
I = imread('pre_computed_results/img/dlpfc_map_sch_800_spl_fig_3d.png');    
imshow(I(img_x_idx, img_y_idx,:));
title({'Schaefer normative connectome', 'DLPFC-SGC WM hops'}, ...
    'FontWeight', 'normal', 'FontSize', font_size+1);

panel{5} = subplot(2,3,5);
I = imread('pre_computed_results/img/dlpfc_map_sch_800_fc_fig_3e.png');
imshow(I(img_x_idx, img_y_idx,:));
title('DLPFC-SGC FC', 'FontWeight', 'normal', 'FontSize', font_size+1);

panel{6} = subplot(2,3,6);
sct = scatter(V_spl_sgc(tms_mask==1),V_fc_sgc(tms_mask==1));
[rho,p] = corr(V_spl_sgc(tms_mask==1),V_fc_sgc(tms_mask==1), 'Type', 'S');
% sct.SizeData = 75;
% sct.MarkerEdgeColor = 'k';
% sct.LineWidth = 0.75;
xlabel('DLPFC-SGC WM hops');
ylabel('DLPFC-SGC FC');
set(gca, 'FontSize', font_size);
title(sprintf('r = %.2f, p < 10^{-12}', rho), 'FontWeight', 'normal', 'FontSize', font_size+1);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
grid minor;
axis square
adjust_scatter_limits(panel{6}, 0.1);

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 4: adjustments
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

axis(panel{1}, 'tight');
axis(panel{2}, 'tight');
axis(panel{4}, 'tight');
axis(panel{5}, 'tight');

panel{3}.Position = [0.6916    0.5838    0.2200    0.3412];
panel{6}.Position = [0.6916    0.1100    0.2200    0.3412];
%% Plot Fig S 1A,B,C

% --- Result parameters --- %
    
parc_lbl = 'sch_800';
stim_sphere_wei_lbl = 'em1';
com_lbl = 'SPL_wei_hop';

parc_xyz = parc_data.(parc_lbl).coords_xyz;

% --- Figure parameters --- %

font_size = 18;
font_size_title = 22;

pc = 1;

clear fig panel;

fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
fig.Position = [1 1 1349 947];
set(fig, 'AutoResizeChildren', 'off');

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 1: Cohort 1
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

cohort_lbl = 'cohort_1';

stim_xyz = tms_model.(parc_lbl).(cohort_lbl).stim_xyz; % tms sites
sgc_xyz = tms_model.(parc_lbl).(cohort_lbl).sgc_xyz; % sgc
imp = tms_data.(cohort_lbl).improvement; % clinical outcome
spl_stim_target = tms_model.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).x; % communication
vxl_rho = rmap_exc_vxl.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).vxl;  % rmap vector
vxl_rho(isnan(vxl_rho)) = [];
rho_rh_sgc = tms_model.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).corr_imp.rho; % RH SGC corr
rho_lh_sgc = rmap.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).V(48,71,31); % LH SGC corr
rmap_rh_sgc_perc = 100*nnz(vxl_rho <= rho_rh_sgc)./numel(vxl_rho); % RH SGC Rmap position
rmap_lh_sgc_perc = 100*nnz(vxl_rho <= rho_lh_sgc)./numel(vxl_rho); % LH SGC Rmap position

offset = 0;

% ------------------------------------------------------------------------------ %
% Cohort 1, brain plot + TMS coordiantes + SGC
% ------------------------------------------------------------------------------ %

% Brain coordinates 1
panel{pc} = subplot(2,5,1+offset);
pc = pc + 1;
scatter(parc_xyz(:,1), parc_xyz(:,2), 10, 'black', 'MarkerEdgeAlpha', 0.5); hold on;
scatter(stim_xyz(:,1), stim_xyz(:,2), 15, 'filled', 'blue'); hold on;
scatter(sgc_xyz(:,1), sgc_xyz(:,2), 30, 'filled', 'red'); hold on;
% scatter(sgc_xyz(:,1), sgc_xyz(:,2), 250, 'red', 'LineWidth', 1); hold on;
set(gca, 'XDir', 'reverse');
axis equal;
axis off;

% Brain coordinates 2
panel{pc} = subplot(2,5,2+offset);
pc = pc + 1;
scatter(parc_xyz(:,1), parc_xyz(:,3), 10, 'black', 'MarkerEdgeAlpha', 0.5); hold on;
scatter(stim_xyz(:,1), stim_xyz(:,3), 15, 'filled', 'blue'); hold on;
scatter(sgc_xyz(:,1), sgc_xyz(:,3), 30, 'filled', 'red'); hold on;
% scatter(sgc_xyz(:,1), sgc_xyz(:,3), 250, 'red', 'LineWidth', 1); hold on;
set(gca, 'XDir', 'reverse');
axis equal;
axis off;

% ------------------------------------------------------------------------------ %
% Cohort 1, scatter com vs improvement
% ------------------------------------------------------------------------------ %

panel{pc} = subplot(2,5,3+offset);
pc = pc + 1;
sct = scatter(spl_stim_target, imp, 'o', 'filled');
[rho_rh_sgc, p_rh_sgc] = corr(spl_stim_target, imp, 'Type', 'S');
sct.SizeData = 75;
sct.MarkerEdgeColor = 'k';
sct.LineWidth = 0.75;
xlabel({'WM hops from DLPFC to SGC'});
ylabel('Treatment response (%)');
set(gca, 'FontSize', font_size);
title(sprintf('r = %.2f, p = %.4f', rho_rh_sgc, p_rh_sgc/2), 'FontWeight', 'normal', 'FontSize', font_size+1);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
grid minor;
adjust_scatter_limits(panel{pc-1}, 0.1);

% ------------------------------------------------------------------------------ %
% Cohort 1, Rmap cortical projection
% ------------------------------------------------------------------------------ %

panel{pc} = subplot(2,5,4+offset);
pc = pc + 1;
img_path = 'pre_computed_results/img/rmap_fig_1d.png';
I = imread(img_path);
imshow(I(1:2400,:,:));
title({'WM hops from DLPFC vs. treatment response'}, 'FontWeight', 'normal', 'FontSize', font_size+1);

% ------------------------------------------------------------------------------ %
% Cohort 1, Rmap histogram + RH SGC position
% ------------------------------------------------------------------------------ %

% Histogram
panel{pc} = subplot(2,5,5+offset);
pc = pc + 1;
n_bins = 18;
histogram(vxl_rho, n_bins, 'EdgeColor', 'none'); hold on;
histogram(vxl_rho, n_bins, 'DisplayStyle', 'stairs', 'LineWidth', 1, 'EdgeColor', 'k'); hold on;
histogram(vxl_rho, n_bins, 'DisplayStyle', 'bar', 'LineWidth', 0.25, 'EdgeColor', 'k', 'FaceColor', 'none'); hold on;
xlabel('Correlation');
ylabel('Number of voxels');
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
set(gca, 'YTick', [0 4000 8000 12000]);
% adjust_histogram_limits(panel{pc-1}, 0.1);
xlim([-1, 1]);
ylim([0,13000]);
set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1]);
set(gca, 'XTickLabelRotation', 0);
panel{pc-1}.YAxis.Exponent = 3;

% Add LH SGC arrow
text(panel{pc-1},  0.35, 11000, ...
    {'LH SGC', sprintf('r = %.2f', rho_lh_sgc), sprintf('top %.0f%%', rmap_lh_sgc_perc)}, ...
    'Color', 'k', 'FontSize', font_size-4, 'HorizontalAlignment','left');
line(panel{pc-1}, [rho_lh_sgc, rho_lh_sgc], [0 7000], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 0.75);
line(panel{pc-1}, [rho_lh_sgc, 0.6], [7000 7000], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 0.75);
line(panel{pc-1}, [0.6, 0.6], [7000 9000], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 0.75);
scatter(panel{pc-1}, 0.6, 9000, '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(panel{pc-1}, rho_lh_sgc, 0, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

% Add RH SGC arrow
text(panel{pc-1}, -0.9, 11000, ...
    {'RH SGC', sprintf('r = %.2f', rho_rh_sgc), sprintf('top %.0f%%', rmap_rh_sgc_perc)}, ...
    'Color', 'k', 'FontSize', font_size-4, 'HorizontalAlignment','left');
line(panel{pc-1}, [rho_rh_sgc, rho_rh_sgc], [0 7000], 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.75);
line(panel{pc-1}, [rho_rh_sgc, -0.65], [7000 7000], 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.75);
line(panel{pc-1}, [-0.65, -0.65], [7000 9000], 'LineStyle', '-', 'Color', 'r', 'LineWidth', 0.75);
scatter(panel{pc-1}, -0.65, 9000, '^', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
scatter(panel{pc-1}, rho_rh_sgc, 0, 'filled', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 3: Adjustments
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

% [left bottom width height]

% for i = 1:10    
%     panel{i}.PlotBoxAspectRatio = [1 1 1];    
% end

% scatters
panel{3}.PlotBoxAspectRatio = [1 1 1];     

% hisotgrams
panel{5}.PlotBoxAspectRatio = [1 1 1];    

% left
% panel{1}.Position(1) = 0;
% panel{2}.Position(1) = 0.125;

position_offset_1 = [0.025 0.05 0 0];
panel{1}.Position = [0.0250    0.6400    0.1000    0.3412] + position_offset_1;
panel{2}.Position = [0.0250    0.4450    0.1000    0.3412] + position_offset_1;
panel{3}.Position = [0.2000    0.5838    0.1600    0.3412] + position_offset_1;
panel{4}.Position = [0.2200    0.5100    0.6000    0.3750] + position_offset_1;
panel{5}.Position = [0.6800    0.5838    0.1600    0.3412] + position_offset_1; % panel{5}.Position = [0.735    0.5838    0.1600    0.3412] + position_offset_1;

axes(panel{1});
text(80,107, 'Cohort I (N = 25)', 'FontWeight', 'normal', 'FontSize', font_size+3);

axes(panel{2});
text(80,0, {'\color{blue} \bullet \color{black}  DLPFC TMS sites','\color{red} \bullet \color{black} SGC seed'}, ...
    'Interpreter', 'tex', 'FontWeight', 'normal', 'FontSize', font_size);

%% Plot Fig S 2A,C

stim_sphere_wei_lbl = 'em1';

font_size = 14;

clear fig panel;
fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
fig.Position = [1 33 1091 833]; % [1 1 1349 947];
set(fig, 'AutoResizeChildren', 'off');

parc_lbl = 'lau_l5';

panel{1} = subplot(2,2,1);
loo_p = tms_model.(parc_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp_loo.p/2;
rho_rh_sgc = tms_model.(parc_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp.rho;
histogram(tms_model.(parc_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp_loo.rho, 12, 'LineWidth', 1);
adjust_histogram_limits;
title(sprintf('Lausanne connectome, Cohort I, max p = %.4f', max(loo_p)), ...
    'FontWeight', 'normal', 'FontSize', 13, 'Interpreter', 'none');
ylabel('Number of correlations');
xlabel({'Leave-one-out correlations between', 'treatment response and DLPFC-SGC WM hops'});
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
line_aux = line([rho_rh_sgc, rho_rh_sgc], ylim, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);
legend(line_aux, {'All patients'}, 'Box', 'off', 'Location', 'northwest', 'FontSize', font_size-2);
axis square;

parc_lbl = 'sch_800';

panel{3} = subplot(2,2,3);
loo_p = tms_model.(parc_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp_loo.p/2;
rho_rh_sgc = tms_model.(parc_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp.rho;
histogram(tms_model.(parc_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp_loo.rho, 12, 'LineWidth', 1);
adjust_histogram_limits;
title(sprintf('Schaefer connectome, Cohort I, max p = %.4f', max(loo_p)), ...
    'FontWeight', 'normal', 'FontSize', 13, 'Interpreter', 'none');
ylabel('Number of correlations');
xlabel({'Leave-one-out correlations between', 'treatment response and DLPFC-SGC WM hops'});
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
line([rho_rh_sgc, rho_rh_sgc], ylim, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);
axis square;

%% Plot Fig S 3A,C

% --- Result parameters --- %

com_stim_sphere_wei_com_lbl = 'em1';
com_stim_sphere_wei_fc_lbl = 'uni';

X_vec = {'SC', 'FC', 'ED', 'CMY_wei', 'CMY_bin', 'SPL_wei', 'SPL_wei_hop', 'SPL_bin', 'SI_wei', 'SI_bin', 'NPL_wei', 'NPL_bin'};
n_X_vec = length(X_vec);

% --- Figure parameters --- %

font_size = 14;

clear fig panel;

fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
fig.Position = [1 183 1048 683];
set(fig, 'AutoResizeChildren', 'off');

X_colors = (othercolor('Accent4', n_X_vec));
non_sig_color = [0.75 0.75 0.75];
non_sig_p_thr = 0.05;

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 1: LAU
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

parc_lbl = 'lau_l5';

% ------------------------------------------------------------------------------ %
% Cohort 1
% ------------------------------------------------------------------------------ %

cohort_lbl = 'cohort_1';

panel{1} = subplot(2,2,1);
rho_vec = zeros(n_X_vec,1);
p_vec = zeros(n_X_vec,1);
for k_X = 1:n_X_vec
    X_lbl = X_vec{k_X};
    if strcmp(X_lbl, 'FC')
        rho_vec(k_X) = ...
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(com_stim_sphere_wei_fc_lbl).corr_imp.rho;
        p_vec(k_X) = ...
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(com_stim_sphere_wei_fc_lbl).corr_imp.p./2; % one-tailed test
    else
        rho_vec(k_X) = ...
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(com_stim_sphere_wei_com_lbl).corr_imp.rho;
        p_vec(k_X) = ...
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(com_stim_sphere_wei_com_lbl).corr_imp.p./2; % one-tailed test        
    end
end
% Sort corr vector
tmp_X_colors = X_colors;
tmp_X_colors(p_vec > non_sig_p_thr,:) = repmat(non_sig_color, nnz(p_vec > non_sig_p_thr), 1);
abs_rho_vec = abs(rho_vec);
[abs_rho_vec_sidx, sidx] = sort(abs_rho_vec, 'descend', 'MissingPlacement', 'last');
rho_vec_sidx = rho_vec(sidx);

% % Correction for multiple comparisons
% % Uncomment below (requires the Bioinformatics Toolbox)
% fprintf('*** %s %s %s ***\n', cohort_lbl, parc_lbl, stim_sphere_wei_lbl);
% % FDR
% p_vec_corrected = mafdr(p_vec(2:end),'BHFDR', true); % correct for comparisons across network communication measures
% p_table = table(X_vec(2:end)', p_vec(2:end), p_vec_corrected, 'VariableNames', {'Measure', 'p-value', 'FDR BH p-value'});        
% p_table = sortrows(p_table, 'p-value', 'ascend');
% display(p_table)

% Bar plot
X_vec_sidx = X_vec(sidx);
bar(abs_rho_vec_sidx, 'FaceColor', 'flat', 'CData', tmp_X_colors(sidx,:)); hold on;
set(gca, 'YLim', [0, 0.7]);
set(gca, 'XTickLabel', strrep(X_vec_sidx, '_', ' '));
set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'XTickLabelRotation', 90);
ylabel({'Correlation to', 'treatment response (|r|)'});
set(gca, 'FontSize', font_size);
for k_X = 1:n_X_vec
    if rho_vec_sidx(k_X) > 0
        scatter(k_X, abs_rho_vec_sidx(k_X) + 0.025, 'Marker', '+', 'MarkerEdgeColor', 'r', 'LineWidth', 1); hold on;
    else
        scatter(k_X, abs_rho_vec_sidx(k_X) + 0.025, 'Marker', 'o',  'MarkerEdgeColor', 'b', 'LineWidth', 1); hold on;
    end
    if strcmp(X_vec_sidx(k_X), 'SPL_wei_hop')
        scatter(k_X, abs_rho_vec_sidx(k_X) + 0.07, 'Marker', 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1); hold on;
    end
end
title('Lausanne connectome, Cohort I', 'FontWeight', 'normal', 'Interpreter', 'none');
set(gca,'PlotBoxAspectRatio', [2 1 1]);

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 2: SCH
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

parc_lbl = 'sch_800';

% ------------------------------------------------------------------------------ %
% Cohort 1
% ------------------------------------------------------------------------------ %

cohort_lbl = 'cohort_1';

panel{3} = subplot(2,2,3);
rho_vec = zeros(n_X_vec,1);
p_vec = zeros(n_X_vec,1);
for k_X = 1:n_X_vec
    X_lbl = X_vec{k_X};
    if strcmp(X_lbl, 'FC')
        rho_vec(k_X) = ...
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(com_stim_sphere_wei_fc_lbl).corr_imp.rho;
        p_vec(k_X) = ...
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(com_stim_sphere_wei_fc_lbl).corr_imp.p./2; % one-tailed test
    else
        rho_vec(k_X) = ...
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(com_stim_sphere_wei_com_lbl).corr_imp.rho;
        p_vec(k_X) = ...
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).(com_stim_sphere_wei_com_lbl).corr_imp.p./2; % one-tailed test        
    end
end
% Sort corr vector
tmp_X_colors = X_colors;
tmp_X_colors(p_vec > non_sig_p_thr,:) = repmat(non_sig_color, nnz(p_vec > non_sig_p_thr), 1);
abs_rho_vec = abs(rho_vec);
[abs_rho_vec_sidx, sidx] = sort(abs_rho_vec, 'descend', 'MissingPlacement', 'last');
rho_vec_sidx = rho_vec(sidx);

% % Correction for multiple comparisons
% % Uncomment below (requires the Bioinformatics Toolbox)
% fprintf('*** %s %s %s ***\n', cohort_lbl, parc_lbl, stim_sphere_wei_lbl);
% % FDR
% p_vec_corrected = mafdr(p_vec(2:end),'BHFDR', true); % correct for comparisons across network communication measures
% p_table = table(X_vec(2:end)', p_vec(2:end), p_vec_corrected, 'VariableNames', {'Measure', 'p-value', 'FDR BH p-value'});        
% p_table = sortrows(p_table, 'p-value', 'ascend');
% display(p_table)

% Bar plot
X_vec_sidx = X_vec(sidx);
bar(abs_rho_vec_sidx, 'FaceColor', 'flat', 'CData', tmp_X_colors(sidx,:)); hold on;
set(gca, 'YLim', [0, 0.7]);
set(gca, 'XTickLabel', strrep(X_vec_sidx, '_', ' '));
set(gca, 'TickLabelInterpreter', 'none');
set(gca, 'XTickLabelRotation', 90);
ylabel({'Correlation to', 'treatment response (|r|)'});
set(gca, 'FontSize', font_size);
for k_X = 1:n_X_vec
    if rho_vec_sidx(k_X) > 0
        scatter(k_X, abs_rho_vec_sidx(k_X) + 0.025, 'Marker', '+', 'MarkerEdgeColor', 'r', 'LineWidth', 1); hold on;
    else
        scatter(k_X, abs_rho_vec_sidx(k_X) + 0.025, 'Marker', 'o',  'MarkerEdgeColor', 'b', 'LineWidth', 1); hold on;
    end
    if strcmp(X_vec_sidx(k_X), 'SPL_wei_hop')
        scatter(k_X, abs_rho_vec_sidx(k_X) + 0.07, 'Marker', 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1); hold on;
    end
end
title('Schaefer connectome, Cohort I', 'FontWeight', 'normal', 'Interpreter', 'none');
set(gca,'PlotBoxAspectRatio', [2 1 1]);

%% Plot Fig S 4A,C

%%% Pre-computed null analysis results
% Computation of surrogate connectomes is time consuming.
% Code to compute surrogate connectomes is from Betzel and Bassett 2018 and is available at https://www.brainnetworkslab.com/s/fcn_match_length_degree_distribution.m
load('pre_computed_results/null_model_SPL_imp.mat');

%%% Plots

stim_sphere_wei_lbl = 'em1';
null_type_lbl = 'fcn_match_length_degree_distribution_dir';

font_size = 14;

clear fig panel;
fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
fig.Position = [1 33 1091 833];
set(fig, 'AutoResizeChildren', 'off');

parc_lbl = 'lau_l5';

panel{1} = subplot(2,2,1);
rho_rh_sgc = tms_model.(parc_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp.rho;
rho_rh_sgc_null = tms_model_null.(parc_lbl).(null_type_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp_madrs.rho;
histogram(rho_rh_sgc_null, 12, 'LineWidth', 1);
adjust_histogram_limits;
title(sprintf('Lausanne connectome, Cohort I, p = %.4f', ...
    nnz(rho_rh_sgc_null  < rho_rh_sgc)/length(rho_rh_sgc_null)), ...
    'FontWeight', 'normal', 'FontSize', font_size, 'Interpreter', 'none');
ylabel('Number of correlations');
xlabel({'Correlations between treatment response and', 'DLPFC-SGC WM hops in surrogate connectomes'});
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
line_aux = line([rho_rh_sgc, rho_rh_sgc], ylim, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);
legend(line_aux, {'Empirical connectome'}, 'Box', 'off', 'Location', 'northeast', 'FontSize', font_size-2);
axis square;

parc_lbl = 'sch_800';

panel{3} = subplot(2,2,3);
rho_rh_sgc = tms_model.(parc_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp.rho;
rho_rh_sgc_null = tms_model_null.(parc_lbl).(null_type_lbl).cohort_1.SPL_wei_hop.(stim_sphere_wei_lbl).corr_imp_madrs.rho;
histogram(rho_rh_sgc_null, 12, 'LineWidth', 1);
adjust_histogram_limits;
title(sprintf('Schaefer connectome, Cohort I, p = %.4f', ...
    nnz(rho_rh_sgc_null  < rho_rh_sgc)/length(rho_rh_sgc_null)), ...
    'FontWeight', 'normal', 'FontSize', font_size, 'Interpreter', 'none');
ylabel('Number of correlations');
xlabel({'Correlations between treatment response and', 'DLPFC-SGC WM hops in surrogate connectomes'});
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
line([rho_rh_sgc, rho_rh_sgc], ylim, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);
axis square;

%% Plot Fig S 5A,C

% --- Result parameters --- %
    
stim_sphere_wei_lbl = 'em1';
com_lbl = 'SPL_wei_hop';

% --- Figure parameters --- %

font_size = 14;

marker_size = 25;

pc = 1;

clear fig panel;

fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
fig.Position = [1 33 1384 833]; % [1 33 1091 833];
set(fig, 'AutoResizeChildren', 'off');

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 1: LAU, cohort I
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

parc_lbl = 'lau_l5';
cohort_lbl = 'cohort_1';

panel{1} = subplot(2,5,1);
x = tms_model.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).x;
y = tms_data.(cohort_lbl).improvement;
[rho, p] = corr(x, y, 'Type', 'S');
sct = scatter(x, y, 'o', 'filled');
adjust_scatter_limits;
sct.SizeData = marker_size;
sct.MarkerEdgeColor = 'k';
sct.LineWidth = 0.75;
xlabel({'RH SGC'});
ylabel('Treatment response (%)');
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
axis square;
grid minor;
title(sprintf('Lausanne connectome, Cohort I\nr = %.2f, p = %.4f', rho, p/2), ...
    'FontWeight', 'normal', 'FontSize', font_size+1);

panel{2} = subplot(2,5,2);
x = tms_model_circuit.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).x;
y = tms_data.(cohort_lbl).improvement;
[rho, p] = corr(x, y, 'Type', 'S');
sct = scatter(x, y, 'o', 'filled');
adjust_scatter_limits;
sct.SizeData = marker_size;
sct.MarkerEdgeColor = 'k';
sct.LineWidth = 0.75;
xlabel({'SGC-caudate-accumbens', 'bilateral circuit'});
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
axis square;
grid minor;
title(sprintf('r = %.2f, p = %.4f', rho, p/2), ...
    'FontWeight', 'normal', 'FontSize', font_size+1);

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 3: SCH, cohort I
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

parc_lbl = 'sch_800';
cohort_lbl = 'cohort_1';

panel{5} = subplot(2,5,6);
x = tms_model.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).x;
y = tms_data.(cohort_lbl).improvement;
[rho, p] = corr(x, y, 'Type', 'S');
sct = scatter(x, y, 'o', 'filled');
adjust_scatter_limits;
sct.SizeData = marker_size;
sct.MarkerEdgeColor = 'k';
sct.LineWidth = 0.75;
xlabel({'RH SGC'});
ylabel('Treatment response (%)');
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
axis square;
grid minor;
title(sprintf('Schaefer connectome, Cohort I\nr = %.2f, p = %.4f', rho, p/2), ...
    'FontWeight', 'normal', 'FontSize', font_size+1);

panel{6} = subplot(2,5,7);
x = tms_model_circuit.(parc_lbl).(cohort_lbl).(com_lbl).(stim_sphere_wei_lbl).x;
y = tms_data.(cohort_lbl).improvement;
[rho, p] = corr(x, y, 'Type', 'S');
sct = scatter(x, y, 'o', 'filled');
adjust_scatter_limits;
sct.SizeData = marker_size;
sct.MarkerEdgeColor = 'k';
sct.LineWidth = 0.75;
xlabel({'SGC-caudate-accumbens', 'bilateral circuit'});
set(gca, 'FontSize', font_size);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
axis square;
grid minor;
title(sprintf('r = %.2f, p = %.4f', rho, p/2), ...
    'FontWeight', 'normal', 'FontSize', font_size+1);
%% Plot Fig S 16A,C

% --- Result parameters --- %

com_lbl = 'SPL_wei_hop';

% --- Figure parameters --- %

font_size = 14;

clear fig panel;

fig = figure;
fig.Renderer = 'painters';
fig.WindowStyle = 'normal';
fig.PaperPositionMode = 'auto';
fig.PaperOrientation = 'landscape';
fig.Position = [1 33 1091 833];
set(fig, 'AutoResizeChildren', 'off');

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 1: LAU
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

parc_lbl = 'lau_l5';

V_spl_sgc = dlpfc_maps.(parc_lbl).(com_lbl).V;
[~, aux] = min(V_spl_sgc(:));
[spl_opt_x, spl_opt_y, spl_opt_z] = ind2sub(size(V_spl_sgc), aux);
spl_opt_mni = cor2mni([spl_opt_x, spl_opt_y, spl_opt_z], parc_data.(parc_lbl).T_vxl2mni);

% ------------------------------------------------------------------------------ %
% Cohort 1
% ------------------------------------------------------------------------------ %

cohort_lbl = 'cohort_1';
panel{1} = subplot(2,2,1);
stim_xyz = tms_model.(parc_lbl).(cohort_lbl).stim_xyz;
imp = tms_data.(cohort_lbl).improvement;
D_spl_opt = squareform(pdist([[spl_opt_x, spl_opt_y, spl_opt_z]; stim_xyz]));
sct = scatter(D_spl_opt(2:end,1), imp, 'o', 'filled');
adjust_scatter_limits;
[rho,p] = corr(D_spl_opt(2:end,1), imp, 'Type', 'S');
sct.SizeData = 75;
sct.MarkerEdgeColor = 'k';
sct.LineWidth = 0.75;
xlabel({'Distance to model-based DLPFC target', sprintf('(MNI [%d %d %d])', ...
    spl_opt_mni(1), spl_opt_mni(2), spl_opt_mni(3))});
ylabel('Improvement');
set(gca, 'FontSize', font_size);
title(sprintf('Lausanne connectome, Cohort 1, r = %.2f, p = %.4f', rho, p/2), ...
    'FontWeight', 'normal', 'FontSize', font_size+1);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
grid minor;
axis square;

% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %
% Part 2: SCH
% ------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------ %

parc_lbl = 'sch_800';

V_spl_sgc = dlpfc_maps.(parc_lbl).(com_lbl).V;
[~, aux] = min(V_spl_sgc(:));
[spl_opt_x, spl_opt_y, spl_opt_z] = ind2sub(size(V_spl_sgc), aux);
spl_opt_mni = cor2mni([spl_opt_x, spl_opt_y, spl_opt_z], parc_data.(parc_lbl).T_vxl2mni);

% ------------------------------------------------------------------------------ %
% Cohort 1
% ------------------------------------------------------------------------------ %

cohort_lbl = 'cohort_1';
panel{3} = subplot(2,2,3);
stim_xyz = tms_model.(parc_lbl).(cohort_lbl).stim_xyz;
imp = tms_data.(cohort_lbl).improvement;
D_spl_opt = squareform(pdist([[spl_opt_x, spl_opt_y, spl_opt_z]; stim_xyz]));
sct = scatter(D_spl_opt(2:end,1), imp, 'o', 'filled');
adjust_scatter_limits;
[rho,p] = corr(D_spl_opt(2:end,1), imp, 'Type', 'S');
sct.SizeData = 75;
sct.MarkerEdgeColor = 'k';
sct.LineWidth = 0.75;
xlabel({'Distance to model-based DLPFC target', sprintf('(MNI [%d %d %d])', ...
    spl_opt_mni(1), spl_opt_mni(2), spl_opt_mni(3))});
ylabel('Improvement');
set(gca, 'FontSize', font_size);
title(sprintf('Schaefer connectome, Cohort 1, r = %.2f, p = %.4f', rho, p/2), ...
    'FontWeight', 'normal', 'FontSize', font_size+1);
set(gca, 'LineWidth', 1);
set(gca, 'Box', 'on');
set(gca, 'Color', [1 1 1]*0.975);
grid minor;
axis square;
