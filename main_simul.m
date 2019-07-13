
%% This file runs a main simulation by identifying 22 different cancer types (Need for Fig1(B) type results)

addpath(genpath(pwd))

load('all_data')

ind_cancer=setdiff(1:33, [10 11 12 16 20 24 25 27 29 31 33]); 

ind_sel=[]; size_ind_sel=[]; size_ind_sel0=[]; size_cum=0;
for iii=1:length(ind_cancer)
    ind_sel=[ind_sel;find(all_pat==ind_cancer(iii))]; size_ind_sel=[size_ind_sel, sum(all_pat==ind_cancer(iii))];
    size_cum=size_cum+sum(all_pat==ind_cancer(iii));
    size_ind_sel0=[size_ind_sel0, size_cum];
end
size_ind_sel0=[0,size_ind_sel0];


%start from here
ksize=length(ind_cancer);  %%ksize is the number of cancer types.

nmi_cen1=[]; nmi_cen2=[]; nmi_cen3=[]; nmi_cen123=[];
nmi_cm123=[]; nmi_pm123=[];
nmik1=[]; nmik2=[]; nmik3=[]; nmik123=[];
nmis1=[]; nmis2=[]; nmis3=[]; nmis123=[];
nmi_keradd123=[];
nmi_wd123=[]; nmi_wde123=[]; nmi_wdf123=[];
nmi_ed123=[]; nmi_edf123=[];
nmi_sml1=[]; nmi_sml2=[]; nmi_sml3=[]; nmi_sml123=[];
nmi_ssc1=[]; nmi_ssc2=[]; nmi_ssc3=[]; nmi_ssc123=[];

pmi_cen1=[]; pmi_cen2=[]; pmi_cen3=[]; pmi_cen123=[];
pmi_cm123=[]; pmi_pm123=[];
pmik1=[]; pmik2=[]; pmik3=[]; pmik123=[];
pmis1=[]; pmis2=[]; pmis3=[]; pmis123=[];
pmi_keradd123=[];
pmi_wd123=[]; pmi_wde123=[]; pmi_wdf123=[];
pmi_ed123=[]; pmi_edf123=[];
pmi_sml1=[]; pmi_sml2=[]; pmi_sml3=[]; pmi_sml123=[];
pmi_ssc1=[]; pmi_ssc2=[]; pmi_ssc3=[]; pmi_ssc123=[];


rmi_cen1=[]; rmi_cen2=[]; rmi_cen3=[]; rmi_cen123=[];
rmi_cm123=[]; rmi_pm123=[];
rmik1=[]; rmik2=[]; rmik3=[]; rmik123=[];
rmis1=[]; rmis2=[]; rmis3=[]; rmis123=[];
rmi_keradd123=[];
rmi_wd123=[]; rmi_wde123=[]; rmi_wdf123=[];
rmi_ed123=[]; rmi_edf123=[];
rmi_sml1=[]; rmi_sml2=[]; rmi_sml3=[]; rmi_sml123=[];
rmi_ssc1=[]; rmi_ssc2=[]; rmi_ssc3=[]; rmi_ssc123=[];


parfor uuu=1:50    %% number of iterations
n4= 30*22;   %% consider total random 660 pateints; 660 should be multiple of ksize.
ttt10=[];
for tuu=1:ksize
ttt0=randperm(size_ind_sel(tuu)); ttt10=[ttt10, ttt0(1:(n4/ksize))+size_ind_sel0(tuu)];
end
ttt1=ind_sel(ttt10);

n=n4;

sel_exp=(all_exp(:,ttt1)); 
sel_cna=(all_cna(:,ttt1));   
sel_mirna=(all_mirna(:,ttt1));    
sel_clin=all_clin(ttt1,:); 
true_labs=all_pat(ttt1);   true_labs0=true_labs;

for jjj=1:ksize;  true_labs0(true_labs==ind_cancer(jjj))=jjj; end
true_labs=true_labs0;
C = max(true_labs); CCC=C; 

sel_exp2 = double(reg_pca(sel_exp',100));
sel_cna2 = double(reg_pca(sel_cna',100));
sel_mirna2 = double(reg_pca(sel_mirna',100));

data_set3={sel_exp2,sel_mirna2 sel_cna2};  
data_set33={double(sel_exp2),double(sel_mirna2), double(sel_cna2)}; 
in_Xi= [sel_exp2,sel_mirna2,sel_cna2]; 


%% consensus 
U = {'U_H','std',[]};   
r = 100;
w = ones(r,1); % the weight of each partitioning
rep = 10; % the number of ECC runs
maxIter = 40;
minThres = 1e-5;
utilFlag = 0;


IDX1 = BasicCluster_RPS(data_set33{1},r,CCC,'sqEuclidean',1);%Generate basic partitions
[pi_sumbest,pi_index1,pi_converge,pi_utility,cons1] = RunECC(IDX1,CCC,U,w,rep,maxIter,minThres,utilFlag); % run KCC for consensus clustering

IDX2 = BasicCluster_RPS(data_set33{2},r,CCC,'sqEuclidean',1);%Generate basic partitions
[pi_sumbest,pi_index2,pi_converge,pi_utility,cons2] = RunECC(IDX2,CCC,U,w,rep,maxIter,minThres,utilFlag); % run KCC for consensus clustering

IDX3 = BasicCluster_RPS(data_set33{3},r,CCC,'sqEuclidean',1);%Generate basic partitions
[pi_sumbest,pi_index3,pi_converge,pi_utility,cons2] = RunECC(IDX3,CCC,U,w,rep,maxIter,minThres,utilFlag); % run KCC for consensus clustering

IDXX123=[IDX1,IDX2,IDX3];
[pi_sumbest,pi_index123,pi_converge,pi_utility,cons3] = RunECC(IDXX123,CCC,U,ones(3*r,1),rep,maxIter,minThres,utilFlag); 

nmi_cen1=[nmi_cen1,Cal_NMI(pi_index1, true_labs)];pmi_cen1=[pmi_cen1,purity(CCC, pi_index1, true_labs)];
rmi_cen1=[rmi_cen1,RandIndex(true_labs, pi_index1)];
nmi_cen2=[nmi_cen2,Cal_NMI(pi_index2, true_labs)];pmi_cen2=[pmi_cen2,purity(CCC, pi_index2, true_labs)];
rmi_cen2=[rmi_cen2,RandIndex(true_labs, pi_index2)];
nmi_cen3=[nmi_cen3,Cal_NMI(pi_index3, true_labs)];pmi_cen3=[pmi_cen3,purity(CCC, pi_index3, true_labs)];
rmi_cen3=[rmi_cen3,RandIndex(true_labs, pi_index3)];

nmi_cen123=[nmi_cen123,Cal_NMI(pi_index123, true_labs)];pmi_cen123=[pmi_cen123,purity(CCC, pi_index123, true_labs)];
rmi_cen123=[rmi_cen123,RandIndex(true_labs, pi_index123)];


%% spectral_centroid_multiview
num_views=length(data_set33); sigma_c=10*ones(1,num_views); %true_labs=[ones(1,round(3/n)),1*ones(1,round(3/n))
[Fc Pc Rc nmi_cm2 avgent_cm AR_cm  idx_cm123] = spectral_centroid_multiview(data_set33,num_views,CCC,sigma_c, 0.01*ones(1,num_views),true_labs,30);
nmi_cm123=[nmi_cm123,Cal_NMI(idx_cm123, true_labs)];pmi_cm123=[pmi_cm123,purity(CCC, idx_cm123, true_labs)];
rmi_cm123=[rmi_cm123,RandIndex(true_labs, idx_cm123)];


%% spectral_pairwise_multview
num_views=length(data_set33); sigma_c=10*ones(1,num_views); %true_labs=[ones(1,round(3/n)),1*ones(1,round(3/n))
[Fc Pc Rc nmi_pm2 avgent_cm AR_cm  idx_pm123] = spectral_pairwise_multview(data_set33,num_views,CCC,sigma_c, 0.01,true_labs,30);
nmi_pm123=[nmi_pm123,Cal_NMI(idx_pm123, true_labs)];pmi_pm123=[pmi_pm123,purity(CCC, idx_pm123, true_labs)];
rmi_pm123=[rmi_pm123,RandIndex(true_labs, idx_pm123)];



%% this is for multiple cancer
K=length(data_set3);  gg=ones(1,K);
Wfc0s_euc_reg=cell(1,K); Wfc0s_euc_near=cell(1,K);  %Wfg_euc_reg_n=cell(1,K); 

sigma_set=1:0.25:2;   
k_set=10:2:30;

[Wfc0s_euc_near_n]=generate_sim_matrices2(K,data_set3,gg, true_labs,sigma_set,k_set);



Wfc0s_euc_near_n_average=cell(1,K);
for iii=1:K
    ini_average=0;
    for iiaverage=1:55
    Wfc0s_euc_near_n_average{iii}{1}{1}=ini_average+Wfc0s_euc_near_n{iii}{1}{iiaverage}/55;
    end
end


%% this is the k-menas using only first omic data
idxxk1=litekmeans(double(data_set3{1}),CCC,'Replicates',50); %% run k-means on embeddings to get cell populations
idxxk2=litekmeans(double(data_set3{2}),CCC,'Replicates',50); %% run k-means on embeddings to get cell populations
idxxk3=litekmeans(double(data_set3{3}),CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

%% this is the k-menas using all omic datasets
idxxk123=litekmeans(double(in_Xi),CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

nmik1=[nmik1,Cal_NMI(idxxk1, true_labs)];pmik1=[pmik1,purity(CCC, idxxk1, true_labs)];
rmik1=[rmik1,RandIndex(true_labs, idxxk1)];
nmik2=[nmik2,Cal_NMI(idxxk2, true_labs)];pmik2=[pmik2,purity(CCC, idxxk2, true_labs)];
rmik2=[rmik2,RandIndex(true_labs, idxxk2)];
nmik3=[nmik3,Cal_NMI(idxxk3, true_labs)];pmik3=[pmik3,purity(CCC, idxxk3, true_labs)];
rmik3=[rmik3,RandIndex(true_labs, idxxk3)];
nmik123=[nmik123,Cal_NMI(idxxk123, true_labs)];pmik123=[pmik123,purity(CCC, idxxk123, true_labs)];
rmik123=[rmik123,RandIndex(true_labs, idxxk123)];



%% this is the spectral clustering using the first omic data
[V_tot1, temp, evs]=eig1(double(Wfc0s_euc_near_n{1}{1}{22}), CCC); 
Clus_ind_spec1=litekmeans(V_tot1,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

[V_tot2, temp, evs]=eig1(double(Wfc0s_euc_near_n{2}{1}{22}), CCC); 
Clus_ind_spec2=litekmeans(V_tot2,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

[V_tot3, temp, evs]=eig1(double(Wfc0s_euc_near_n{3}{1}{22}), CCC); 
Clus_ind_spec3=litekmeans(V_tot3,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

Clus_ind_spec123=litekmeans([V_tot1,V_tot2,V_tot3],CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

nmis1=[nmis1,Cal_NMI(true_labs, Clus_ind_spec1)]; pmis1=[pmis1,purity(CCC, Clus_ind_spec1, true_labs)];
rmis1=[rmis1,RandIndex(true_labs, Clus_ind_spec1)];
nmis2=[nmis2,Cal_NMI(true_labs, Clus_ind_spec2)]; pmis2=[pmis2,purity(CCC, Clus_ind_spec2, true_labs)];
rmis2=[rmis2,RandIndex(true_labs, Clus_ind_spec2)];
nmis3=[nmis3,Cal_NMI(true_labs, Clus_ind_spec3)]; pmis3=[pmis3,purity(CCC, Clus_ind_spec3, true_labs)];
rmis3=[rmis3,RandIndex(true_labs, Clus_ind_spec3)];
nmis123=[nmis123,Cal_NMI(true_labs, Clus_ind_spec123)]; pmis123=[pmis123,purity(CCC, Clus_ind_spec123, true_labs)];
rmis123=[rmis123,RandIndex(true_labs, Clus_ind_spec123)];


%% this is our algorithm without kernel learning to each omic data
K=length(data_set3); gg=ones(1,K);  c=0.1; rho=2;  lam=0.001; mu=0.1; eta=1;     %K=1; gg=[1];
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update2_2HW(CCC, c,rho, n, K, 1,1, gg, lam, mu, eta, Wfc0s_euc_near_n_average);   %Wfc0_euc_reg); % Wfc0_euc_near_n) ;
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set]; 
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wok123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

nmi_keradd123=[nmi_keradd123,Cal_NMI(true_labs, Clus_ind_wok123)]; pmi_keradd123=[pmi_keradd123,purity(CCC, Clus_ind_wok123, true_labs)];
rmi_keradd123=[rmi_keradd123,RandIndex(true_labs, Clus_ind_wok123)];



%% this is our algorithm with learned weight to each omic data
K=length(data_set3); gg=ones(1,K); c=0.1; rho=2;  lam=0.001; mu=0.1; eta=1;     %K=1; gg=[1];
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update2_2HW(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n);   %Wfc0_euc_reg); % Wfc0_euc_near_n) ;
ck123=ck;
%this incorportates learned weight in the three target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wd123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

V_tot=[]; for dd=1:K;  V_tot=[V_tot, 1*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wde123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

%this uses the first leared target matrix for clustering analysis
V_tot=tresult{2}; 
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wdf123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

nmi_wd123=[nmi_wd123,Cal_NMI(true_labs, Clus_ind_wd123)]; pmi_wd123=[pmi_wd123,purity(CCC, Clus_ind_wd123, true_labs)];
rmi_wd123=[rmi_wd123,RandIndex(true_labs, Clus_ind_wd123)];
nmi_wde123=[nmi_wde123,Cal_NMI(true_labs, Clus_ind_wde123)]; pmi_wde123=[pmi_wde123,purity(CCC, Clus_ind_wde123, true_labs)];
rmi_wde123=[rmi_wde123,RandIndex(true_labs, Clus_ind_wde123)];
nmi_wdf123=[nmi_wdf123,Cal_NMI(true_labs, Clus_ind_wdf123)]; pmi_wdf123=[pmi_wdf123,purity(CCC, Clus_ind_wdf123, true_labs)];
rmi_wdf123=[rmi_wdf123,RandIndex(true_labs, Clus_ind_wdf123)];


%% this is our algorithm using equal weight to each omic data
 K=length(data_set3);gg=ones(1,K);  c=0.1; rho=2;  lam=0.001;  mu=1; eta=1;     
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update2_2HW0(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n);   %Wfc0_euc_reg); % Wfc0_euc_near_n) ;

%this use all the learned target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set]; 
%V_tot=[tresult{K+3}(1)*tresult{K+2}(:,1:CCC),tresult{K+3}(2)*tresult{K+2}(:,(CCC+1):(2*CCC)),tresult{K+3}(3)*tresult{K+2}(:,(2*CCC+1):(3*CCC))];
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_ed123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

%this only use the first leared target matrix for clustering analysis
V_tot=tresult{2}; 
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_edf123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

nmi_ed123=[nmi_ed123,Cal_NMI(true_labs, Clus_ind_ed123)]; pmi_ed123=[pmi_ed123,purity(CCC, Clus_ind_ed123, true_labs)];
rmi_ed123=[rmi_ed123,RandIndex(true_labs, Clus_ind_ed123)];
nmi_edf123=[nmi_edf123,Cal_NMI(true_labs, Clus_ind_edf123)]; pmi_edf123=[pmi_edf123,purity(CCC, Clus_ind_edf123, true_labs)];
rmi_edf123=[rmi_edf123,RandIndex(true_labs, Clus_ind_edf123)];

%% this is simlr1 using all the omic data
[S_simlr00, F1] = SIMLR_LARGE(double(data_set3{1}),CCC,30); %Y is already decreased matrix, %%S is the learned similarity, F is the latent embedding 
Clus_simlr1=litekmeans(F1,CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

[S_simlr00, F2] = SIMLR_LARGE(double(data_set3{2}),CCC,30); %Y is already decreased matrix, %%S is the learned similarity, F is the latent embedding 
Clus_simlr2=litekmeans(F2,CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

[S_simlr00, F3] = SIMLR_LARGE(double(data_set3{3}),CCC,30); %Y is already decreased matrix, %%S is the learned similarity, F is the latent embedding 
Clus_simlr3=litekmeans(F3,CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

Clus_simlr123=litekmeans([F1,F2,F3],CCC,'Replicates',50); %% run k-means on embeddings to get cell populations

nmi_sml1=[nmi_sml1,Cal_NMI(true_labs, Clus_simlr1)]; pmi_sml1=[pmi_sml1,purity(CCC, Clus_simlr1, true_labs)];
rmi_sml1=[rmi_sml1,RandIndex(true_labs, Clus_simlr1)];
nmi_sml2=[nmi_sml2,Cal_NMI(true_labs, Clus_simlr2)]; pmi_sml2=[pmi_sml2,purity(CCC, Clus_simlr2, true_labs)];
rmi_sml2=[rmi_sml2,RandIndex(true_labs, Clus_simlr2)];
nmi_sml3=[nmi_sml3,Cal_NMI(true_labs, Clus_simlr3)]; pmi_sml3=[pmi_sml3,purity(CCC, Clus_simlr3, true_labs)];
rmi_sml3=[rmi_sml3,RandIndex(true_labs, Clus_simlr3)];
nmi_sml123=[nmi_sml123,Cal_NMI(true_labs, Clus_simlr123)]; pmi_sml123=[pmi_sml123,purity(CCC, Clus_simlr123, true_labs)];
rmi_sml123=[rmi_sml123,RandIndex(true_labs, Clus_simlr123)];


%% sparse spect cluster
[Ps,objs,errs,iters] = sparsesc(eye(n)-Wfc0s_euc_near_n{1}{1}{22},0.00001,CCC);
[V_tot1, temp, evs]=eig1(Ps, CCC); 
Clus_ind_sparsc1=litekmeans(double(V_tot1),CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

[Ps,objs,errs,iters] = sparsesc(eye(n)-Wfc0s_euc_near_n{2}{1}{22},0.00001,CCC);
[V_tot2, temp, evs]=eig1(Ps, CCC); 
Clus_ind_sparsc2=litekmeans(double(V_tot2),CCC,'Replicates',50); 

[Ps,objs,errs,iters] = sparsesc(eye(n)-Wfc0s_euc_near_n{3}{1}{22},0.00001,CCC);
[V_tot3, temp, evs]=eig1(Ps, CCC); 
Clus_ind_sparsc3=litekmeans(double(V_tot3),CCC,'Replicates',50); 

Clus_ind_sparsc123=litekmeans(double([V_tot1,V_tot2,V_tot3]),CCC,'Replicates',50); 

nmi_ssc1=[nmi_ssc1,Cal_NMI(true_labs, Clus_ind_sparsc1)]; pmi_ssc1=[pmi_ssc1,purity(CCC, Clus_ind_sparsc1, true_labs)];
rmi_ssc1=[rmi_ssc1,RandIndex(true_labs, Clus_ind_sparsc1)];
nmi_ssc2=[nmi_ssc2,Cal_NMI(true_labs, Clus_ind_sparsc2)]; pmi_ssc2=[pmi_ssc2,purity(CCC, Clus_ind_sparsc2, true_labs)];
rmi_ssc2=[rmi_ssc2,RandIndex(true_labs, Clus_ind_sparsc2)];
nmi_ssc3=[nmi_ssc3,Cal_NMI(true_labs, Clus_ind_sparsc3)]; pmi_ssc3=[pmi_ssc3,purity(CCC, Clus_ind_sparsc3, true_labs)];
rmi_ssc3=[rmi_ssc3,RandIndex(true_labs, Clus_ind_sparsc3)];
nmi_ssc123=[nmi_ssc123,Cal_NMI(true_labs, Clus_ind_sparsc123)]; pmi_ssc123=[pmi_ssc123,purity(CCC, Clus_ind_sparsc123, true_labs)];
rmi_ssc123=[rmi_ssc123,RandIndex(true_labs, Clus_ind_sparsc123)];


end

save('similar_color_map.mat',   'Wfc0s_euc_near_n')


nmi_set{1}=nmi_cen1;nmi_set{2}=nmi_cen2;nmi_set{3}=nmi_cen3;nmi_set{4}=nmi_cen123; nmi_set{5}=nmi_cm123;
nmi_set{6}=nmi_pm123;nmi_set{7}=nmik1;nmi_set{8}=nmik2;nmi_set{9}=nmik3;nmi_set{10}=nmik123;
nmi_set{11}=nmis1;nmi_set{12}=nmis2;nmi_set{13}=nmis3;nmi_set{14}=nmis123;
nmi_set{15}=nmi_keradd123; nmi_set{16}=nmi_wd123;nmi_set{17}=nmi_wde123;nmi_set{18}=nmi_wdf123;
nmi_set{19}=nmi_ed123;nmi_set{20}=nmi_edf123;nmi_set{21}=nmi_sml1;nmi_set{22}=nmi_sml2;nmi_set{23}=nmi_sml3;
nmi_set{24}=nmi_sml123; nmi_set{25}=nmi_ssc1;nmi_set{26}=nmi_ssc2;nmi_set{27}=nmi_ssc3;
nmi_set{28}=nmi_ssc123;


pmi_set{1}=pmi_cen1;pmi_set{2}=pmi_cen2;pmi_set{3}=pmi_cen3;pmi_set{4}=pmi_cen123; pmi_set{5}=pmi_cm123;
pmi_set{6}=pmi_pm123;pmi_set{7}=pmik1;pmi_set{8}=pmik2;pmi_set{9}=pmik3;pmi_set{10}=pmik123;
pmi_set{11}=pmis1;pmi_set{12}=pmis2;pmi_set{13}=pmis3;pmi_set{14}=pmis123;
pmi_set{15}=pmi_keradd123; pmi_set{16}=pmi_wd123;pmi_set{17}=pmi_wde123;pmi_set{18}=pmi_wdf123;
pmi_set{19}=pmi_ed123;pmi_set{20}=pmi_edf123;pmi_set{21}=pmi_sml1;pmi_set{22}=pmi_sml2;pmi_set{23}=pmi_sml3;
pmi_set{24}=pmi_sml123;pmi_set{25}=pmi_ssc1;pmi_set{26}=pmi_ssc2;pmi_set{27}=pmi_ssc3;
pmi_set{28}=pmi_ssc123;


rmi_set{1}=rmi_cen1;rmi_set{2}=rmi_cen2;rmi_set{3}=rmi_cen3;rmi_set{4}=rmi_cen123; rmi_set{5}=rmi_cm123;
rmi_set{6}=rmi_pm123;rmi_set{7}=rmik1;rmi_set{8}=rmik2;rmi_set{9}=rmik3;rmi_set{10}=rmik123;
rmi_set{11}=rmis1;rmi_set{12}=rmis2;rmi_set{13}=rmis3;rmi_set{14}=rmis123;
rmi_set{15}=rmi_keradd123; rmi_set{16}=rmi_wd123;rmi_set{17}=rmi_wde123;rmi_set{18}=rmi_wdf123;
rmi_set{19}=rmi_ed123;rmi_set{20}=rmi_edf123;rmi_set{21}=rmi_sml1;rmi_set{22}=rmi_sml2;rmi_set{23}=rmi_sml3;
rmi_set{24}=rmi_sml123;rmi_set{25}=rmi_ssc1;rmi_set{26}=rmi_ssc2;rmi_set{27}=rmi_ssc3;
rmi_set{28}=rmi_ssc123;

save('simul_threedata_22.mat','rmi_set','nmi_set','pmi_set' )


