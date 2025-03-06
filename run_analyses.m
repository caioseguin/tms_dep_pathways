
%%% Neuroanatomical pathways of TMS therapy for depression
%%% Seguin et al., https://doi.org/10.1101/2025.02.10.25322034

%%% This script runs the main analyses of the manuscript using publicly
%%% available data from Weigand et al 2018, named Cohort I in our study.
%%% Data from Cohort II is not publicly available and therefore analyses
%%% using this patient group are not reproduced here.
%%% We perform analyses for two publicly available normative connectome
%%% datasets from Mansour et al 2021 and Griffa et al 2019.
%%% We use functions from the BCT (Rubinov & Sporns 2010) to compute
%%% network communication measures.

%%% Caio Seguin, 6 March 2025

clc; clear; close all;

fprintf('Load data...\n');

addpath data/
addpath func/
addpath pre_computed_results/

%%% Load TMS coordinates and clinical improvement for cohort I (Weigand et al 2018)
load tms_data;

%%% Load normative connectome datasets
load parc_data;

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
    L = log_transform(W); % connection lengths
    S = diag(sum(W,1));

    %%% Communication matrices
    % Weighted CMY
    CMY_wei = expm(S^(-1/2)*W*S^(-1/2));
    % Binary CMY
    CMY_bin = expm(A);
    % Weighted SPL
    [SPL_wei, SPL_wei_hop] = distance_wei_floyd(L);
    % Binary SPL
    SPL_bin = distance_wei_floyd(A);
    % Search information weighted
    SI_wei = search_information(W,L);
    SI_bin = search_information(A,A);
    % Navigation
    [~, NPL_bin, NPL_wei] = navigate_bct(L, ED);

    parc_data.(parc_lbl).CMY_wei = CMY_wei;
    parc_data.(parc_lbl).CMY_bin = CMY_bin;
    parc_data.(parc_lbl).SPL_wei = SPL_wei;
    parc_data.(parc_lbl).SPL_wei_hop = SPL_wei_hop;
    parc_data.(parc_lbl).SPL_bin = SPL_bin;
    parc_data.(parc_lbl).SI_wei = SI_wei;
    parc_data.(parc_lbl).SI_bin = SI_bin;
    parc_data.(parc_lbl).NPL_wei = NPL_wei;
    parc_data.(parc_lbl).NPL_bin = NPL_bin;

end

%% Model communication from DLPFC TMS sites to the SGC and correlate it to treatment response 

fprintf('Model communication from DLPFC TMS sites to the SGC and correlate it to treatment response...\n');

sgc_mni = [6, 16, -10]; % from https://doi.org/10.1016/j.biopsych.2018.12.002 
r_sphere_mm = 10;

corr_type = 'Spearman';

cohort_vec = {'cohort_1'};
n_cohort_vec = length(cohort_vec);

parc_vec = {'lau_l5', 'sch_800'};
n_parc_vec = length(parc_vec);

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
            aux(aux == 0) = 1;
            stim_sphere_dis(i,:) = aux;
        end

        % Compute target (SGC) block
        r_target = r_sphere_mm./vxl_size;
        target_block = sphere_vxl(sgc_xyz, r_target, V);

        % Exclude subcortical regions from target block
        switch parc_lbl
            case 'lau_l5'
                target_block((target_block >= 502 & target_block <= 508) | (target_block >= 1008 & target_block <= 1014)) = 0;
            case 'sch_800'
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

        uni_stim_sphere_wei = ones(m, n_sphere_vxl);

        % Compute communication from stimulation to targets
        for k_X = 1:n_X_vec

            X_lbl = X_vec{k_X};
            X = parc_data.(parc_lbl).(X_lbl);
            x_uni = zeros(m,1);
            x_d2 = zeros(m,1);

            for i = 1:m
                x_uni(i) = block_wei_avg(X, stim_block(i,:), target_block, uni_stim_sphere_wei(i,:));
            end

            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.x = x_uni;

            [tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.corr_imp_madrs.rho, ...
             tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.corr_imp_madrs.p] = ...
                corr(tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.x, imp, 'Type', corr_type);

            % Leave-one-out correlations
            aux_rho_loo_uni = zeros(m,1);
            aux_p_loo_uni = zeros(m,1);       
            for i = 1:m
                aux_x_uni = tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.x;
                aux_imp = imp;
                aux_x_uni(i) = [];
                aux_imp(i) = [];
                [aux_rho_loo_uni(i), aux_p_loo_uni(i)] = corr(aux_x_uni, aux_imp, 'Type', corr_type);
            end

            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.corr_imp_madrs_loo.rho = aux_rho_loo_uni;
            tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.corr_imp_madrs_loo.p = aux_p_loo_uni;

        end

    end

end


%% Correlation between DLPFC-SGC WM hops and treatment response (Cohort I, Fig 1C)

fprintf('Correlation between DLPFC-SGC WM hops and treatment response (Cohort I, Fig 1C)...\n');

figure;

my_scatter(tms_model.lau_l5.cohort_1.SPL_wei_hop.stim_sphere_uni.x, tms_data.cohort_1.improvement, ...
    'Cohort I, Lausanne connectome (Fig 1C)', 'WM hops from DLPFC to SGC', 'Treatment reponse (%)');

%% Histogram of the correlation between WM hops from the DLPFC and treatment response (Cohort I, Fig 1E) 

fprintf('Histogram of the correlation between WM hops from the DLPFC and treatment response (Cohort I, Fig 1E)...\n');

%%% Run time is approximately 45 mins. Uncomment section below to run or load pre-computed results

% V = parc_data.lau_l5.V;
% imp = tms_data.cohort_1.improvement;
% SPL = parc_data.lau_l5.SPL_wei_hop;
% stim_block = tms_model.lau_l5.cohort_1.stim_block;
% 
% V_idx = find(V);
% [vxl_x, vxl_y, vxl_z] = ind2sub(size(V), V_idx);
% 
% n_vxl = size(vxl_x,1);
% 
% vxl_rho_spl_imp = zeros(n_vxl, 1);
% 
% for k_vxl = 1:n_vxl
% 
%     fprintf('(%d/%d) %.2f%% complete\n', k_vxl, n_vxl, 100*k_vxl/n_vxl);
% 
%     vxl_target_block = sphere_vxl([vxl_x(k_vxl), vxl_y(k_vxl), vxl_z(k_vxl)], r_target, V);
% 
%     spl_stim_target = zeros(m,1);
%     
%     for i = 1:m
%         spl_stim_target(i) = block_wei_avg(SPL, stim_block(i,:), vxl_target_block);
%     end
%     
%     vxl_rho_spl_imp(k_vxl) = corr(spl_stim_target, imp, 'Type', 'S');
% 
% end


%%% Load pre-computed results
load('voxel_rho_spl_imp_cohort_1_laussane')

V = parc_data.lau_l5.V;
V_rho_spl_imp = zeros(size(V));
V_rho_spl_imp(V ~= 0) = vxl_rho_spl_imp;
rho_lh_sgc = V_rho_spl_imp(48,71,31);
rho_rh_sgc = -0.55908;

% Histogram
figure;
histogram(vxl_rho_spl_imp, 15);
xlabel('Correlation');
ylabel('Number of voxels');
set(gca, 'FontSize', 14);
xlim([-1, 1]);
ylim([0,14000]);
line([rho_lh_sgc, rho_lh_sgc], [0,14000], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5);
line([rho_rh_sgc, rho_rh_sgc], [0,14000], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);
title(sprintf('Fig 1E (Cohort I). Red: RH SGC r = %.2f (top %.2f%%). Black: LH SGC r = %.2f (top %.2f%%)', ...
    rho_rh_sgc, 100*nnz(vxl_rho_spl_imp <= rho_rh_sgc)./numel(vxl_rho_spl_imp), ...
    rho_lh_sgc, 100*nnz(vxl_rho_spl_imp <= rho_lh_sgc)./numel(vxl_rho_spl_imp)), ...
    'FontWeight', 'normal');
set(gca, 'XTick', [-1, -0.5, 0, 0.5, 1]);
set(gca, 'XTickLabelRotation', 0);

%% Identify neuroanatomical pathways from the DLPFC to the SGC (Cohort I, Figs 2B,C)

fprintf('Identify neuroanatomical pathways from the DLPFC to the SGC (Cohort I, Figs 2B,C)...\n');

L = log_transform(parc_data.lau_l5.SC);
[SPL, SPL_h, SPL_Pmat] = distance_wei_floyd(L);
stim_block = tms_model.lau_l5.cohort_1.stim_block;
target_block = tms_model.lau_l5.cohort_1.target_block;
parc_labels = parc_data.lau_l5.labels;
n = length(SPL);

stim_target_SP = table('Size', [0, 2], 'VariableTypes', {'double', 'double'}, 'VariableNames', {'path', 'count'});
stim_target_SP_str_labels = table('Size', [0, 2], 'VariableTypes', {'string', 'double'}, 'VariableNames', {'path', 'count'});
stim_target_SP_str_labels_inter = table('Size', [0, 2], 'VariableTypes', {'string', 'double'}, 'VariableNames', {'path', 'count'});

stim_target_SP_hops = [];

mediator_matrix = zeros(n,n-1);

E_sp = zeros(n);

tab_stim = tabulate(flat(stim_block));
stim_block_unique = tab_stim(2:end, 1);
stim_block_count = tab_stim(2:end, 2);

tab_target = tabulate(flat(target_block));
target_block_unique = tab_target(2:end, 1);
target_block_count = tab_target(2:end, 2);

for i = 1:length(stim_block_unique)
    for j = 1:length(target_block_unique)

        stim_idx = stim_block_unique(i);
        target_idx = target_block_unique(j);

        n_stim_taget = stim_block_count(i)*target_block_count(j);

        path = retrieve_shortest_path(stim_idx, target_idx, SPL_h, SPL_Pmat);
        n_hops = length(path);

        path_string = [];
        path_string_labels = [];
        path_string_inter = [];
        path_string_labels_inter = [];

        for k = 1:(n_hops-1)
            E_sp(path(k), path(k+1)) = E_sp(path(k), path(k+1)) + n_stim_taget;
            E_sp(path(k+1), path(k)) = E_sp(path(k+1), path(k)) + n_stim_taget;
        end

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

        stim_target_SP_hops = [stim_target_SP_hops; repmat(n_hops, n_stim_taget, 1)];

    end
end

stim_target_SP.percentage = stim_target_SP.count./sum(stim_target_SP.count);
stim_target_SP_str_labels.percentage = stim_target_SP_str_labels.count./sum(stim_target_SP_str_labels.count);
stim_target_SP_str_labels_inter.percentage = stim_target_SP_str_labels_inter.count./sum(stim_target_SP_str_labels_inter.count);

stim_target_SP = sortrows(stim_target_SP, 'count', 'descend');
stim_target_SP_str_labels = sortrows(stim_target_SP_str_labels, 'count', 'descend');
stim_target_SP_str_labels_inter = sortrows(stim_target_SP_str_labels_inter, 'count', 'descend');

mediator_regions = table(parc_labels, sum(mediator_matrix,2), sum(mediator_matrix,2)./sum(mediator_matrix(:)), 'VariableNames', {'region', 'count', 'percentage'});
mediator_regions = sortrows(mediator_regions, 'count', 'descend');

figure;
histogram(stim_target_SP_hops-1);
title('Fig 2B (Cohort I)', 'FontWeight', 'normal');
xlabel('Hops');
ylabel('Number of paths');
set(gca, 'FontSize', 14);

fprintf('Top regions (Fig 2C, purple bars):\n');
display(mediator_regions);

fprintf('*** Top paths (Fig 2C, blue bars) ***\n');
display(stim_target_SP_str_labels_inter);


%% Correlation between DLPFC maps: WM hops vs FC (Cohort I, Fig 3C,F)

fprintf('Correlation between DLPFC maps: WM hops vs FC (Cohort I, Fig 3C,F)...\n');

% Load MNI brain template edges
brain_edges = niftiread('MNI152_T1_2mm_edges.nii.gz');
% Load DLPFC mask
dlpfc_mask = niftiread('DLPFCcombo20mm.nii');

%%% Normative connectome 1 (Lausanne parcellation)

V = parc_data.lau_l5.V;
SPL = parc_data.lau_l5.SPL_wei_hop;
FC = parc_data.lau_l5.FC;
target_block = tms_model.lau_l5.cohort_1.target_block;

se = strel('sphere', 1);
brain_edges_dil = imdilate(brain_edges, se.Neighborhood);
tms_mask = dlpfc_mask & brain_edges_dil & ~~V;
tms_mask_idx = find(tms_mask);

V_spl_sgc = inf(size(tms_mask));
V_cmy_sgc = zeros(size(tms_mask));
V_fc_sgc = zeros(size(tms_mask));

for i = 1:length(tms_mask_idx)
    [x,y,z] = ind2sub(size(V_spl_sgc), tms_mask_idx(i));
    V_spl_sgc(x,y,z) = block_wei_avg(SPL, sphere_vxl([x,y,z], r_stim, V), target_block);
    V_fc_sgc(x,y,z) = block_wei_avg(FC, sphere_vxl([x,y,z], r_stim, V), target_block);
end

figure;
my_scatter(V_spl_sgc(tms_mask==1), V_fc_sgc(tms_mask==1), ...
    'Fig 3C. Normative connectome 1 (Lausanne parcellation)', 'DLPFC-SGC WM hops', 'DLPFC-SGC FC');

%%% Normative connectome 2 (Schaefer parcellation)

V = parc_data.sch_800.V;
SPL = parc_data.sch_800.SPL_wei_hop;
FC = parc_data.sch_800.FC;
target_block = tms_model.sch_800.cohort_1.target_block;

se = strel('sphere', 1);
brain_edges_dil = imdilate(brain_edges, se.Neighborhood);
tms_mask = dlpfc_mask & brain_edges_dil & ~~V;
tms_mask_idx = find(tms_mask);

V_spl_sgc = inf(size(tms_mask));
V_cmy_sgc = zeros(size(tms_mask));
V_fc_sgc = zeros(size(tms_mask));

for i = 1:length(tms_mask_idx)
    [x,y,z] = ind2sub(size(V_spl_sgc), tms_mask_idx(i));
    V_spl_sgc(x,y,z) = block_wei_avg(SPL, sphere_vxl([x,y,z], r_stim, V), target_block);
    V_fc_sgc(x,y,z) = block_wei_avg(FC, sphere_vxl([x,y,z], r_stim, V), target_block);
end

figure;
my_scatter(V_spl_sgc(tms_mask==1), V_fc_sgc(tms_mask==1), ...
    'Fig 3F. Normative connectome 2 (Schaefer parcellation)', 'DLPFC-SGC WM hops', 'DLPFC-SGC FC');

%% Leave-one-out outlier sensitivity tests (Cohort I, Fig S1A)

fprintf('Leave-one-out outlier sensitivity tests (Cohort I, Fig S1A)...\n');

rho_rh_sgc = -0.55908;

% Histogram
figure;
histogram(tms_model.lau_l5.cohort_1.SPL_wei_hop.stim_sphere_uni.corr_imp_madrs_loo.rho, 12, 'LineWidth', 1);
title(sprintf('Fig S1A. Cohort I, max p = %.4f', max(tms_model.lau_l5.cohort_1.SPL_wei_hop.stim_sphere_uni.corr_imp_madrs_loo.p)), ...
    'FontWeight', 'normal', 'FontSize', 14);
ylabel('Number of correlations');
xlabel({'Leave-one-out correlations between', 'treatment response and DLPFC-SGC WM hops'});
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
line([rho_rh_sgc, rho_rh_sgc], ylim, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);
axis square;

%% Surrogate connectome analyses (Cohort I, Fig S2 A)

fprintf('Surrogate connectome analyses (Cohort I, Fig S2 A)...\n');

%%% Pre-computed null analysis results
% Computation of surrogate connectomes is time consuming.
% Code to compute surrogate connectomes is from Betzel and Bassett 2018 and is available at https://www.brainnetworkslab.com/s/fcn_match_length_degree_distribution.m
load('SC_null_rho_imp_cohort_1_lausanne.mat');

rho_rh_sgc = -0.55908;

% Histogram
figure;
histogram(null_struct.rho_imp.spl, 12, 'LineWidth', 1);
xlabel({'Correlations between treatment response and', 'DLPFC-SGC WM hops in surrogate connectomes'});
ylabel('Number of surrogates');
set(gca, 'FontSize', 14);
title(sprintf('Fig S2A: Cohort I, p = %.4f', nnz((null_struct.rho_imp.spl)  <= (rho_rh_sgc))/length(null_struct.rho_imp.spl)), ...
    'FontWeight', 'normal', 'FontSize', 14);
set(gca, 'LineWidth', 1);
set(gca, 'Color', [1 1 1]*0.975);
line([rho_rh_sgc, rho_rh_sgc], ylim, 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);
axis square;

%% Correlation between treatment response and DLPFC-SGC communication using different network measures (Cohort I, Fig S3A,C)

fprintf('Correlation between treatment response and DLPFC-SGC communication using different network measures (Cohort I, Fig S3A,C)...\n');

font_size = 14;

n_X_vec = length(X_vec);

X_colors = (othercolor('Accent4', n_X_vec));
non_sig_color = [0.75 0.75 0.75];
non_sig_p_thr = 0.05;

title_vec = {'Fig S3A: Cohort I (USA), Normative connectome 1 (Lausanne parcellation)'; ''; ' Fig S3C: Cohort I (USA), Normative connectome 2 (Schaefer parcellation)'};

fig = figure;

subplot_counter = 1;
panel = {};

for k_parc = 1:n_parc_vec

    for k_cohort = 1:n_cohort_vec

        cohort_lbl = cohort_vec{k_cohort};
        parc_lbl = parc_vec{k_parc}; 
        rho_vec = zeros(n_X_vec,1);
        p_vec = zeros(n_X_vec,1);

        for k_X = 1:n_X_vec
            
            X_lbl = X_vec{k_X};
            rho_vec(k_X) = tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.corr_imp_madrs.rho;
            p_vec(k_X) = tms_model.(parc_lbl).(cohort_lbl).(X_lbl).stim_sphere_uni.corr_imp_madrs.p./2; % one-tailed test

        end

        panel{subplot_counter} = subplot(2,2,subplot_counter);
        tmp_X_colors = X_colors;
        tmp_X_colors(p_vec > non_sig_p_thr,:) = repmat(non_sig_color, nnz(p_vec > non_sig_p_thr), 1);
        
        abs_rho_vec = abs(rho_vec);
        [abs_rho_vec_sidx, sidx] = sort(abs_rho_vec, 'descend', 'MissingPlacement', 'last');
        rho_vec_sidx = rho_vec(sidx);

        % Correction for multiple comparisons
        fprintf('*** %s %s ***\n', cohort_lbl, parc_lbl);

        % FDR
        p_vec_corrected = mafdr(p_vec(4:end),'BHFDR', true); % correct for comparisons across network communication measures
        p_table = table(X_vec(4:end)', p_vec(4:end), p_vec_corrected, 'VariableNames', {'Measure', 'p-value', 'FDR BH p-value'});
        p_table = sortrows(p_table, 'p-value', 'ascend');
        display(p_table)
        

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
                scatter(k_X, abs_rho_vec_sidx(k_X) + 0.025, 'Marker', '_',  'MarkerEdgeColor', 'b', 'LineWidth', 1); hold on;
            end
            if strcmp(X_vec_sidx(k_X), 'SPL_wei_hop')
                scatter(k_X, abs_rho_vec_sidx(k_X) + 0.07, 'Marker', 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1); hold on;
            end
        end
        title(title_vec{subplot_counter}, 'FontWeight', 'normal', 'Interpreter', 'none');
        set(gca,'PlotBoxAspectRatio', [2 1 1]);
        
        subplot_counter = subplot_counter + 2;

    end
end