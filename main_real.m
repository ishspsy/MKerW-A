%% This file runs diverse clustering methods for 22 cancer types and perform survival analysis.

addpath(genpath(pwd))
% The following file can be downloaded from 
% https://www.dropbox.com/s/v22fx0j2gnpeta6/all_data.mat?dl=0
load('all_data')   
% all_clin includes clinical information for the patients from 22 cancer types
% all_exp is the rna data for all the patinets
% all_mirna is the mirna data for all the patinets
% all_cna is the cna data for all the patinets
% all_pat is the index vector of patients indicating corresponding cancer types (See the first column of all_clin for the original cancer name)
load('clsnumber_set.mat')
% CCC_set is the target cluster number for 22 cancer types

% The followings are indices of each of the 22 cancer types:

% iii = 1 BRCA
% iii = 2 STAD
% iii = 3 LUAD
% iii = 4 LUSC
% iii = 5 COAD
% iii = 6 HNSC
% iii = 7 KIRC
% iii = 8 BLCA
% iii = 9 UVM
% iii = 10 PRAD
% iii = 11 SARC
% iii = 12 KIRP
% iii = 13 LIHC
% iii = 14 PAAD
% iii = 15 ESCA
% iii = 16 LGG
% iii = 17 MESO
% iii = 18 UCEC
% iii = 19 THCA
% iii = 20 READ
% iii = 21 OV
% iii = 22 CESC


%%% Part 1: Perform clustering analysis
iii_ind = 0;

for iii= setdiff(1:33, [10 11 12 16 20 24 25 27 29 31 33])   %% consider 22 cancer types
% iii_ind is the cancer index
iii_ind = iii_ind + 1

clearvars -except all_exp all_cna all_mirna all_clin all_pat CCC_set iii iii_ind

%% choose patients belonging to the group of cancer with index iii and do data processing
ind_set=find(all_pat==iii);  ttt1=ind_set; n=length(ttt1);

sel_exp=(all_exp(:,ttt1));   
sel_cna=(all_cna(:,ttt1));   
sel_mirna=(all_mirna(:,ttt1));   
sel_clin=all_clin(ttt1,:); 

surv_stat=sel_clin(:,end-4);  
surv_stat=strrep(table2array(surv_stat),'Alive','1');
surv_stat=strrep(surv_stat,'Dead','0');   
surv_stat=str2double(surv_stat); 
surv_time=sel_clin(:,end-3);  surv_time=table2array(surv_time);  surv_time=str2double(surv_time);

ind_set2=setdiff(1:n,union(find(isnan(surv_time)), find(isnan(surv_stat))));  n=length(ind_set2);
surv_time=surv_time(ind_set2);  surv_stat=surv_stat(ind_set2);

sel_exp=(sel_exp(:,ind_set2));   
sel_cna=(sel_cna(:,ind_set2));   
sel_mirna=(sel_mirna(:,ind_set2));   
sel_clin=sel_clin(ind_set2,:); 

% this step is optional
sel_exp2 = reg_pca(sel_exp',min(n,100));
sel_mirna2 = reg_pca(sel_mirna',min(n,100));
sel_cna2 = reg_pca(sel_cna',min(n,100));

% use the one of the following data sets for clustering analysis depending on the methods.
data_set3={sel_exp2,sel_mirna2,sel_cna2};  
data_set33={double(sel_exp2),double(sel_mirna2),double(sel_cna2)}; 
in_Xi= [sel_exp2,sel_mirna2,sel_cna2]; 



%% Generate similarity matrices for clustering analysis
K=3; gg=ones(1,K); sigma_set=1:0.25:2;  k_set=10:2:30;
[Wfc0s_euc_near_n]=generate_sim_matrices(K,data_set3, gg, 0,sigma_set,k_set);

%% Target clustering number for cancer iii
CCC=CCC_set(iii);

%% Arbitrary clustering labels
ini_labs=1:CCC;  ini_labs=repmat(ini_labs,[1, round(n/CCC)]); 
if length(ini_labs)<n;  ini_labs=[ini_labs, CCC*ones(1,n-length(ini_labs))];  
elseif length(ini_labs)>n; ini_labs=ini_labs(1:n);
end


%% Consensus clustering
U = {'U_H','std',[]};   
r = 100;
w = ones(r,1); % the weight of each partitioning
rep = 10;      % the number of ECC runs
maxIter = 40;
minThres = 1e-5;
utilFlag = 0;

% using rna for clustering
IDX1 = BasicCluster_RPS(data_set33{1},r,CCC,'sqEuclidean',1); %Generate basic partitions
[pi_sumbest,pi_index1,pi_converge,pi_utility,cons1] = RunECC(IDX1,CCC,U,w,rep,maxIter,minThres,utilFlag); % run KCC for consensus clustering

% using mirna for clustering
IDX2 = BasicCluster_RPS(data_set33{2},r,CCC,'sqEuclidean',1); %Generate basic partitions
[pi_sumbest,pi_index2,pi_converge,pi_utility,cons2] = RunECC(IDX2,CCC,U,w,rep,maxIter,minThres,utilFlag); % run KCC for consensus clustering

% using cna for clustering
IDX3 = BasicCluster_RPS(data_set33{3},r,CCC,'sqEuclidean',1); %Generate basic partitions
[pi_sumbest,pi_index3,pi_converge,pi_utility,cons2] = RunECC(IDX3,CCC,U,w,rep,maxIter,minThres,utilFlag); % run KCC for consensus clustering

IDXX123=[IDX1,IDX2,IDX3];
[pi_sumbest,pi_index123,pi_converge,pi_utility,cons3] = RunECC(IDXX123,CCC,U,ones(3*r,1),rep,maxIter,minThres,utilFlag); 


%% spectral_centroid_multiview (multi-view method)
num_views=length(data_set33); sigma_c=10*ones(1,num_views); 
[Fc Pc Rc nmi_cm2 avgent_cm AR_cm  idx_cm123] = spectral_centroid_multiview(data_set33,num_views,CCC,sigma_c, 0.01*ones(1,num_views),ini_labs',30);


%% spectral_pairwise_multview  (multi-view method)
num_views=length(data_set33); sigma_c=10*ones(1,num_views); 
[Fc Pc Rc nmi_pm2 avgent_cm AR_cm  idx_pm123] = spectral_pairwise_multview(data_set33,num_views,CCC,sigma_c, 0.01,ini_labs',30);



%% this is the k-menas 

% using rna for clustering
idxxk1=litekmeans(double(data_set3{1}),CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

% using mirna for clustering
idxxk2=litekmeans(double(data_set3{2}),CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

% using cna for clustering
idxxk3=litekmeans(double(data_set3{3}),CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

% using all omics data sets
idxxk123=litekmeans(double(in_Xi),CCC,'Replicates',50); %% run k-means on embeddings to get cell populations


%% this is the Spectral clustering

K=length(data_set3); gg=ones(1,K);

% Average of 55 similarity matrices
Wfc0s_euc_near_n_average=cell(1,K);
for i=1:K
    ini_average=0;
    for iiaverage=1:55
    Wfc0s_euc_near_n_average{i}{1}{1}=ini_average+Wfc0s_euc_near_n{i}{1}{iiaverage}/55;
    end
end

% using rna for clustering
[V_tot1, temp, evs]=eig1(double(Wfc0s_euc_near_n_average{1}{1}{1}), CCC); 
Clus_ind_spec1=litekmeans(V_tot1,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

% using mirna for clustering
[V_tot2, temp, evs]=eig1(double(Wfc0s_euc_near_n_average{2}{1}{1}), CCC); 
Clus_ind_spec2=litekmeans(V_tot2,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

% using cna for clustering
[V_tot3, temp, evs]=eig1(double(Wfc0s_euc_near_n_average{3}{1}{1}), CCC); 
Clus_ind_spec3=litekmeans(V_tot3,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

% using all omics data sets
Clus_ind_spec123=litekmeans([V_tot1,V_tot2,V_tot3],CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations


%% this is our method without similarity learning
K=length(data_set3); gg=ones(1,K);  c=0.1; rho=2;  lam=0.001; mu=1; eta=1; 

% Solve the optimization problem
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update(CCC, c,rho, n, K, 1,1, gg, lam, mu, eta, Wfc0s_euc_near_n_average); 
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set]; 

% Incorporate obtained target matrices for clustering analysis
V_tot=[]; 
for dd=1:K;  
    V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; 
end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wok123=litekmeans(V_tot,CCC,'Replicates',50);  


%% this is our algorithm with learned weight
K=length(data_set3); gg=ones(1,K); c=0.1; rho=2;  lam=0.001; mu=1; eta=1;   

% Solve the optimization problem
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n);  

% learned weight for each omic data
ck123=ck;
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;

% Incorportating learned weights for clustering analysis
V_tot=[]; 
for dd=1:K;  
    V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; 
end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wd123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

% this gives the equal weight for each omic data (used for comparison of clustering result)
V_tot=[]; for dd=1:K;  V_tot=[V_tot, 1*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wde123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

% this only uses the first leared target matrix for clustering analysis (used for comparison of clustering result)
V_tot=tresult{2}; 
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wdf123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations



%% this is our algorithm using equal weight to each omic data (no weight learning)
 K=length(data_set3);gg=ones(1,K);  c=0.1; rho=2;  lam=0.001; mu=1;  mu=1; eta=1;     
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update2(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n);   

%this use all the learned target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set]; 
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_ed123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

%this only use the first leared target matrix for clustering analysis
V_tot=tresult{2}; 
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_edf123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations



%% this is SIMLR 

% using rna for clustering
[S_simlr00, F1] = SIMLR_LARGE(double(data_set3{1}),CCC,30); %Y is already decreased matrix, %%S is the learned similarity, F is the latent embedding 
Clus_simlr1=litekmeans(F1,CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

% using mirna for clustering
[S_simlr00, F2] = SIMLR_LARGE(double(data_set3{2}),CCC,30); %Y is already decreased matrix, %%S is the learned similarity, F is the latent embedding 
Clus_simlr2=litekmeans(F2,CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

% using cna for clustering
[S_simlr00, F3] = SIMLR_LARGE(double(data_set3{3}),CCC,30); %Y is already decreased matrix, %%S is the learned similarity, F is the latent embedding 
Clus_simlr3=litekmeans(F3,CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

% using all omics data sets for clustering
Clus_simlr123=litekmeans([F1,F2,F3],CCC,'Replicates',50); %% run k-means on embeddings to get cell populations


%% Sparse spectral clustering

% using rna for clustering
[Ps,objs,errs,iters] = sparsesc(eye(n)-Wfc0s_euc_near_n_average{1}{1}{1},0.00001,CCC);
[V_tot1, temp, evs]=eig1(Ps, CCC); 
Clus_ind_sparsc1=litekmeans(double(V_tot1),CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

% using mirna for clustering
[Ps,objs,errs,iters] = sparsesc(eye(n)-Wfc0s_euc_near_n_average{2}{1}{1},0.00001,CCC);
[V_tot2, temp, evs]=eig1(Ps, CCC); 
Clus_ind_sparsc2=litekmeans(double(V_tot2),CCC,'Replicates',50); 

% using cna for clustering
[Ps,objs,errs,iters] = sparsesc(eye(n)-Wfc0s_euc_near_n_average{3}{1}{1},0.00001,CCC);
[V_tot3, temp, evs]=eig1(Ps, CCC); 
Clus_ind_sparsc3=litekmeans(double(V_tot3),CCC,'Replicates',50); 

% using all omics data sets for clustering
Clus_ind_sparsc123=litekmeans(double([V_tot1,V_tot2,V_tot3]),CCC,'Replicates',50); 

clus_ind_set12={ini_labs idxxk1 idxxk2 idxxk3 idxxk123 Clus_ind_spec1 Clus_ind_spec2 Clus_ind_spec3 Clus_ind_spec123, ...
    Clus_ind_sparsc1 Clus_ind_sparsc2 Clus_ind_sparsc3  Clus_ind_sparsc123, ...
    Clus_simlr1 Clus_simlr2 Clus_simlr3  Clus_simlr123, ...
    pi_index1 pi_index2 pi_index3  pi_index123, ...
    idx_cm123, idx_pm123, ...
    Clus_ind_wok123 Clus_ind_ed123 Clus_ind_edf123, ...
    Clus_ind_wd123, Clus_ind_wde123, Clus_ind_wdf123};
    

%% Name of clustering methods    
gp_title11={'Rand' 'kmean1' 'kmean2' 'kmean3' 'kmean123' 'spec1' 'spec2' 'spec3' 'spec123', ...
    's-spec1' 's-spec2' 's-spec3'  's-spec123', ...
    'SIM1' 'SIM2' 'SIM3'  'SIM123', ...
    'Cons1' 'Cons2' 'Cons3'  'Cons123', ...
    'Cen-M123'  'Pair-M123', ...
    'Ker-add123'  'Our-ed123' 'Our-edf123', ...
    'Our-wd123', 'Our-wde123', 'Our-wdf123'} 

%% weights of omics data sets assigned in the proposed clustering method
ck_set={ck123};   

%% save the file
save(sprintf('surv_analysis_mu1_lam0001_%d_norandom.mat', iii_ind),  'clus_ind_set12', 'gp_title11', 'surv_time', 'surv_stat', 'ck_set')
end



%%% Part 2: Perform survival analysis using the inferred clusters.

for iii=1:22; 

% load obtained files
load(sprintf('surv_analysis_mu1_lam0001_%d_norandom.mat', iii))
CCC=max(clus_ind_set12{1});

name_i=1; 
gp_title22={'Clus1','Clus2','Clus3','Clus4','Clus5', 'Clus6', 'ALL'}
colors_set={'-b','-g','-r','-k','-m', '--b','--g','--r','--k','--m'};
cl_title={'Blue', 'Green', 'Red', 'Black', 'Magenta','Blue--', 'Green--', 'Red--', 'Black--', 'Magenta--'  };


for jjj=1:length(clus_ind_set12)
    ind_set=cell(1,CCC);   
for ii=1:CCC
    ind_set{ii}=find(clus_ind_set12{jjj}==ii);
end
    ind_set_set{jjj}=ind_set;
end

cannot_use=zeros(1,length(clus_ind_set12));
for jjj=1:length(clus_ind_set12)
    bbm_set=cell(1,CCC); bb_surv=cell(1,CCC);
    for ijk=1:CCC
        bbm_set{ijk}= surv_time(ind_set_set{jjj}{ijk});
        bb_surv{ijk}= surv_stat(ind_set_set{jjj}{ijk}) ;
        if sum(bb_surv{ijk}==1)== length(bb_surv{ijk}); cannot_use(jjj)=1; end
    end
    bbm_set_set{jjj}=bbm_set; bb_surv_set{jjj}=bb_surv;     
end


for jjjj=1:length(clus_ind_set12); for kkkk=1:CCC
    bbm_set_set_m{jjjj}{kkkk}=bbm_set_set{jjjj}{kkkk}/30;
end; end

for jjjj=1:length(clus_ind_set12); for kkkk=1:CCC
    bbm_set_set_m{jjjj}{kkkk}=bbm_set_set_m{jjjj}{kkkk}.*(bbm_set_set_m{jjjj}{kkkk}>0);
end; end


diff_area_sum=zeros(1,length(gp_title11)); diff_area_min=zeros(1,length(gp_title11));
diff_avg_sum=zeros(1,length(gp_title11)); diff_avg_min=zeros(1,length(gp_title11));
diff_area_set=cell(1,length(gp_title11)); diff_avg_set=cell(1,length(gp_title11)); 
diff_area_set2=cell(1,length(gp_title11));

%% Fit Weibull survival model to each inferred cluster and compute the metrics that measure heterogeneity of survival times between inferred groups.

for jjjj=setdiff(1:length(clus_ind_set12),find(cannot_use))
bord=120; 
[pd1_A_set,pd1_B_set]=generate_surv_func_general(CCC, bord, bbm_set_set_m,bb_surv_set, jjjj, colors_set,gp_title11,gp_title22,cl_title, name_i)
[diff_area11]=diff_area_func(CCC, pd1_A_set,pd1_B_set, bord, 0.01);
close all
diff_area_set{jjjj}=diff_area11; diff_area_sum(jjjj)=sum(sum(diff_area11)); diff_area_min(jjjj)=min(min(diff_area11+10000*diag(ones(1,CCC))));
[diff_area22]=diff_area_func_ave(CCC, pd1_A_set,pd1_B_set, bord, 0.01);
diff_avg_set{jjjj}=diff_area22; diff_avg_sum(jjjj)=sum(sum(diff_area22))/CCC; diff_avg_min(jjjj)=min(diff_area22);
[diff_area33]=diff_area_func2(CCC, pd1_A_set,pd1_B_set, bord, 0.01);
diff_area_set2{jjjj}=diff_area33;
end

diff_area_set120=diff_area_set; diff_area_sum120=diff_area_sum;diff_area_min120=diff_area_min;
diff_avg_set120=diff_avg_set; diff_avg_sum120=diff_avg_sum;diff_avg_min120=diff_avg_min;
diff_area_set2_120=diff_area_set2;


area_prop=cell(1,length(clus_ind_set12)); area_prop_sum=cell(1,length(clus_ind_set12));  area_prop_min=cell(1,length(clus_ind_set12));
for jjjj=setdiff(1:length(clus_ind_set12),find(cannot_use))
area_p=diff_area_set2{jjjj}./diff_area_set{jjjj}; area_p(isnan(area_p))=0;
area_prop{jjjj}=area_p; area_prop_sum{jjjj}=sum(sum(area_p))/(CCC^2-CCC); area_prop_min{jjjj}=min(min(area_p+diag(ones(1,CCC)))); 
end

diff_area_sum120=diff_area_sum120/(CCC^2-CCC);


diff_area_sum120=[gp_title11;num2cell(diff_area_sum120)]
diff_area_min120=[gp_title11;num2cell(diff_area_min120)]
diff_avg_sum120=[gp_title11;num2cell(diff_avg_sum120)]
diff_avg_min120=[gp_title11;num2cell(diff_avg_min120)]
area_prop_sum=[gp_title11;(area_prop_sum)]
area_prop_min=[gp_title11;(area_prop_min)]

save(sprintf('surv_analysis_performance_mu1_lam0001_%d_norandom.mat',iii), ...
    'clus_ind_set12','CCC','bord','diff_area_set120','diff_area_sum120','diff_area_min120', ...
    'diff_avg_set120', 'diff_avg_sum120', 'diff_avg_min120','diff_area_set2_120','area_prop','area_prop_sum','area_prop_min')
end



%% this summarizes the results that will be used to generate figures.

iiii=setdiff(1:33, [10 11 12 16 20 24 25 27 29 31 33]); 

tot_num=29;
area_prop_min_tot=[]; area_prop_sum_tot=[]; diff_area_sum_tot=[]; diff_area_min_tot=[]; diff_avg_sum_tot=[];
for ii=1:length(iiii)
load(sprintf('surv_analysis_performance_mu1_lam0001_%d_norandom.mat', ii))
    for jj=1:tot_num
    if cell2mat(area_prop_sum(2,jj))>0.01;   a5=1; 
    else; area_prop_sum{2,jj}=quantile(cell2mat(area_prop_sum(2,:)),0.3);
    area_prop_min{2,jj}=quantile(cell2mat(area_prop_min(2,:)),0.3);
    diff_area_sum120{2,jj}=quantile(cell2mat(diff_area_sum120(2,:)),0.3);
    diff_area_min120{2,jj}=quantile(cell2mat(diff_area_min120(2,:)),0.3); 
    diff_avg_sum120{2,jj}=quantile(cell2mat(diff_avg_sum120(2,:)),0.3); 
    end
    end
    
area_prop_sum_tot=[area_prop_sum_tot;cell2mat(area_prop_sum(2,:))];   %24 by 59 matrix
area_prop_min_tot=[area_prop_min_tot;cell2mat(area_prop_min(2,:))];
diff_area_sum_tot=[diff_area_sum_tot;cell2mat(diff_area_sum120(2,:))];
diff_area_min_tot=[diff_area_min_tot;cell2mat(diff_area_min120(2,:))];
diff_avg_sum_tot=[diff_avg_sum_tot;cell2mat(diff_avg_sum120(2,:))];
end

comp_min=[]; comp_sum=[]; comp_area_min=[]; comp_area_sum=[]; comp_avg_sum=[];
for ii=1:29
    comp_min=[comp_min,area_prop_min_tot(:,ii)<=area_prop_min_tot(:,27)];
    comp_sum=[comp_sum,area_prop_sum_tot(:,ii)<=area_prop_sum_tot(:,27)];
    comp_area_min=[comp_area_min,diff_area_min_tot(:,ii)<=diff_area_min_tot(:,27)];
    comp_area_sum=[comp_area_sum,diff_area_sum_tot(:,ii)<=diff_area_sum_tot(:,27)];
    comp_avg_sum=[comp_avg_sum,diff_avg_sum_tot(:,ii)<=diff_avg_sum_tot(:,27)];
end

save('final_result_four_measures.mat', 'area_prop_sum_tot', 'area_prop_min_tot','diff_area_sum_tot','diff_area_min_tot', ...
 'diff_avg_sum_tot', 'iiii', 'gp_title11','comp_min','comp_sum', 'comp_avg_min', 'comp_avg_sum')


