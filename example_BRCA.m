%% Example using BRCA

addpath(genpath(pwd))
load('all_data')   
% CCC_set is the target cluster number for 22 cancer types

load('clsnumber_set.mat')
% all_clin includes clinical information for the patients from 22 cancer types
% all_exp is the rna data for all the patinets
% all_mirna is the mirna data for all the patinets
% all_cna is the cna data for all the patinets
% all_pat is the index vector of patients indicating corresponding cancer types (See the first column of 
% all_clin for the original cancer name)

%%% Perform clustering analysis

%% Consider a breast cancer (BRCA)
iii=1;       
% iii is the cancer index from the 22 cancer types of the index and iii=1 represents BRCA

%% choose patients belonging to the group of cancer with index iii and do data processing
ind_set=find(all_pat==iii);  ttt1=ind_set; n=length(ttt1);
sel_exp=(all_exp(:,ttt1));sel_cna=(all_cna(:,ttt1));sel_mirna=(all_mirna(:,ttt1));sel_clin=all_clin(ttt1,:); 

%% survival time and status
surv_stat=sel_clin(:,end-4); 
surv_stat=strrep(table2array(surv_stat),'Alive','1'); surv_stat=strrep(surv_stat,'Dead','0'); 
surv_stat=str2double(surv_stat); 
surv_time=sel_clin(:,end-3);  surv_time=table2array(surv_time);  surv_time=str2double(surv_time);

% choose patients having information about survival time and status
ind_set2=setdiff(1:n,union(find(isnan(surv_time)), find(isnan(surv_stat))));  n=length(ind_set2);
surv_time=surv_time(ind_set2);  surv_stat=surv_stat(ind_set2);

%% write omics data sets and clinical information
sel_exp=(sel_exp(:,ind_set2)); sel_cna=(sel_cna(:,ind_set2)); sel_mirna=(sel_mirna(:,ind_set2));   
sel_clin=sel_clin(ind_set2,:); 

% this step is optional
sel_exp2 = reg_pca(sel_exp',min(n,100));
sel_mirna2 = reg_pca(sel_mirna',min(n,100));
sel_cna2 = reg_pca(sel_cna',min(n,100));

%% Generate similarity matrices for clustering analysis
K=3; gg=ones(1,K); sigma_set=1:0.25:2;  k_set=10:2:30;
data_set3={sel_exp2,sel_mirna2,sel_cna2};  
[Wfc0s_euc_near_n]=generate_sim_matrices2(K,data_set3, gg, 0,sigma_set,k_set);

%% Target clustering number 
CCC=4;

%% this is our clustering algorithm 
% regularization parameters are specified as suggested in the manuscript:
K=length(data_set3); gg=ones(1,K); c=0.1; rho=2; lam=0.001; mu=1; eta=1;   

% Solve the optimization problem
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update2_2HW(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, ...
Wfc0s_euc_near_n);  

% learned weight for each omic data  % ck123 records the weight of each data
ck123=ck;
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;

% Incorportating learned weights for clustering analysis
V_tot=[]; 
for dd=1:K;  
    V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; 
end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));

% run K-means on the final output
% Clus_ind_wd123 is the inferred cluster label
Clus_ind_wd123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations