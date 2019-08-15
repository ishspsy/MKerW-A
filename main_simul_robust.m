%% This file perform robustness test of MKerW-A with respect to changes of regularization parameters (need for Figure S4)

addpath(genpath(pwd))

load('all_data')

rng(100)

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


nmi_wd123_set=[]; pmi_wd123_set=[]; rmi_wd123_set=[];
nmi_wd123_setm=[]; pmi_wd123_setm=[]; rmi_wd123_setm=[];
nmi_wd123_setc=[]; pmi_wd123_setc=[]; rmi_wd123_setc=[];
nmi_wd123_setr=[]; pmi_wd123_setr=[]; rmi_wd123_setr=[];


parfor uuu=1:50    %% number of iterations
n4= 30*22; 
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

sel_exp2 = reg_pca(sel_exp',100);
sel_cna2 = reg_pca(sel_cna',100);
sel_mirna2 = reg_pca(sel_mirna',100);


data_set3={sel_exp2,sel_mirna2 sel_cna2};  

data_set33={double(sel_exp2),double(sel_mirna2), double(sel_cna2)}; 
in_Xi= [sel_exp2,sel_mirna2,sel_cna2]; 



%% this is for multiple cancer
K=length(data_set3);  gg=ones(1,K);
Wfc0s_euc_reg=cell(1,K); Wfc0s_euc_near=cell(1,K);  

sigma_set=1:0.25:2;   
k_set=10:2:30;

[Wfc0s_euc_near_n]=generate_sim_matrices2(K,data_set3,gg, true_labs,sigma_set,k_set);


%% this is our algorithm with learned weight to each omic data

lam_set=[0.00001 0.0001 0.001 0.01 0.1];

nmi_wd123=[]; pmi_wd123=[]; rmi_wd123=[];

for ijk=1:length(lam_set);
lam=lam_set(ijk);    
K=length(data_set3); gg=ones(1,K); c=0.1; rho=2;  mu=0.1; eta=1;    
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n);  
ck123=ck;
%this incorportates learned weight in the three target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wd123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

nmi_wd123=[nmi_wd123,Cal_NMI(true_labs, Clus_ind_wd123)]; pmi_wd123=[pmi_wd123,purity(CCC, Clus_ind_wd123, true_labs)];
rmi_wd123=[rmi_wd123,RandIndex(true_labs, Clus_ind_wd123)];
end
nmi_wd123_set=[nmi_wd123_set;nmi_wd123];  pmi_wd123_set=[pmi_wd123_set;pmi_wd123]; rmi_wd123_set=[rmi_wd123_set;rmi_wd123];



mu_set=[0.001 0.01 0.1 1];
nmi_wd123=[]; pmi_wd123=[]; rmi_wd123=[];

for ijk=1:length(mu_set);  
K=length(data_set3); gg=ones(1,K); c=0.1; rho=2;  mu=mu_set(ijk);  lam=0.001;  eta=1;     
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n);  
ck123=ck;
%this incorportates learned weight in the three target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wd123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

nmi_wd123=[nmi_wd123,Cal_NMI(true_labs, Clus_ind_wd123)]; pmi_wd123=[pmi_wd123,purity(CCC, Clus_ind_wd123, true_labs)];
rmi_wd123=[rmi_wd123,RandIndex(true_labs, Clus_ind_wd123)];
end

nmi_wd123_setm=[nmi_wd123_setm;nmi_wd123];  pmi_wd123_setm=[pmi_wd123_setm;pmi_wd123]; rmi_wd123_setm=[rmi_wd123_setm;rmi_wd123];



c_set=[0.001 0.01 0.1 1];
nmi_wd123=[]; pmi_wd123=[]; rmi_wd123=[];

for ijk=1:length(c_set);
c=c_set(ijk);    
K=length(data_set3); gg=ones(1,K); rho=2;  mu=0.1; lam=0.001; eta=1;    
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n); 
ck123=ck;
%this incorportates learned weight in the three target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wd123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

nmi_wd123=[nmi_wd123,Cal_NMI(true_labs, Clus_ind_wd123)]; pmi_wd123=[pmi_wd123,purity(CCC, Clus_ind_wd123, true_labs)];
rmi_wd123=[rmi_wd123,RandIndex(true_labs, Clus_ind_wd123)];
end

nmi_wd123_setc=[nmi_wd123_setc;nmi_wd123];  pmi_wd123_setc=[pmi_wd123_setc;pmi_wd123]; rmi_wd123_setc=[rmi_wd123_setc;rmi_wd123];


rho_set=[0.1 0.5 1 2 5];
nmi_wd123=[]; pmi_wd123=[]; rmi_wd123=[];

for ijk=1:length(rho_set);
c=0.1;    
K=length(data_set3); gg=ones(1,K); rho=rho_set(ijk);  mu=0.1; lam=0.001; eta=1;    
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n); 
ck123=ck;
%this incorportates learned weight in the three target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wd123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

nmi_wd123=[nmi_wd123,Cal_NMI(true_labs, Clus_ind_wd123)]; pmi_wd123=[pmi_wd123,purity(CCC, Clus_ind_wd123, true_labs)];
rmi_wd123=[rmi_wd123,RandIndex(true_labs, Clus_ind_wd123)];
end

nmi_wd123_setr=[nmi_wd123_setr;nmi_wd123];  pmi_wd123_setr=[pmi_wd123_setr;pmi_wd123]; rmi_wd123_setr=[rmi_wd123_setr;rmi_wd123];

end

save('simul_robust.mat', 'nmi_wd123_set','pmi_wd123_set','rmi_wd123_set', ...
    'nmi_wd123_setm','pmi_wd123_setm','rmi_wd123_setm', ...
    'nmi_wd123_setc','pmi_wd123_setc','rmi_wd123_setc', ...
    'nmi_wd123_setr','pmi_wd123_setr','rmi_wd123_setr')
    

