%% This file infers the target cluster number  (need for Fig 1.(A))

addpath(genpath(pwd))

load('all_data')

ind_cancer=setdiff(1:33, [10 11 12 16 20 24 25 27 29 31 33]); 


pat_num=[];
for ii=1:33; pat_num=[pat_num,sum(all_pat==ii)];end;  find(pat_num>30)   %%this is the cancer index having at least 30 patients


ind_sel=[]; size_ind_sel=[]; size_ind_sel0=[]; size_cum=0;
for iii=1:length(ind_cancer)
    ind_sel=[ind_sel;find(all_pat==ind_cancer(iii))]; size_ind_sel=[size_ind_sel, sum(all_pat==ind_cancer(iii))];
    size_cum=size_cum+sum(all_pat==ind_cancer(iii));
    size_ind_sel0=[size_ind_sel0, size_cum];
end
size_ind_sel0=[0,size_ind_sel0];


%start from here
ksize=length(ind_cancer);  %%ksize is the number of cancer types.

nnmi_wd123=[]; ppmi_wd123=[]; rrmi_wd123=[]; 

parfor uuu=1:50    %% number of iterations
n4= ksize*30;  
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
for jjj=1:ksize;  true_labs0(true_labs==ind_cancer(jjj))=jjj; end; true_labs=true_labs0;
sel_exp2 = reg_pca(sel_exp',100);sel_cna2 = reg_pca(sel_cna',100);sel_mirna2 = reg_pca(sel_mirna',100);
data_set3={sel_exp2,sel_mirna2 sel_cna2};   

sigma_noise=1;
seln_exp=sel_exp+normrnd(0,sigma_noise,size(sel_exp,1),n);
seln_cna=sel_cna+normrnd(0,sigma_noise,size(sel_cna,1),n);
seln_mirna=sel_mirna+normrnd(0,sigma_noise,size(sel_mirna,1),n);

seln_exp2 = reg_pca(seln_exp',100);seln_cna2 = reg_pca(seln_cna',100);seln_mirna2 = reg_pca(seln_mirna',100);
datan_set3={seln_exp2,seln_mirna2 seln_cna2};  


%% this is for multiple cancer
K=length(data_set3);  gg=ones(1,K);
Wfc0s_euc_reg=cell(1,K); Wfc0s_euc_near=cell(1,K);   

sigma_set=1:0.25:2;  
k_set=10:2:30;

[Wfc0s_euc_near_n]= generate_sim_matrices2(K,data_set3,gg, true_labs,sigma_set,k_set);
[Wfc0s_euc_near_nn]=generate_sim_matrices2(K,datan_set3,gg, true_labs,sigma_set,k_set);


%% this is our algorithm with learned weight to each omic data
K=length(data_set3); gg=ones(1,K); c=0.1; rho=2;  lam=0.001; mu=1; eta=1;     %K=1; gg=[1];

nmi_wd123=[];  pmi_wd123=[]; rmi_wd123=[];
for CCC=2:30
[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update2_2HW(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_n);   %Wfc0_euc_reg); % Wfc0_euc_near_n) ;
ck123=ck;

%this incorportates learned weight in the three target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wd123=litekmeans(V_tot,CCC,'Replicates',50);  %% run k-means on embeddings to get cell populations

[P_set,V_set,V_tot,ck,W_set,Wg_set]=clus_sim_update2_2HW(CCC, c,rho, n, K, 5,11, gg, lam, mu, eta, Wfc0s_euc_near_nn);   %Wfc0_euc_reg); % Wfc0_euc_near_n) ;
ck123=ck;

%this incorportates learned weight in the three target matrices for clustering analysis
tresult=[P_set,V_set,V_tot,ck,W_set,Wg_set];  tresult_final=tresult;
V_tot=[]; for dd=1:K;  V_tot=[V_tot, tresult{K+1}(dd)*tresult{K+2}(:,(dd*CCC-CCC+1):(dd*CCC))]; end
V_tot=V_tot./ repmat(sqrt(sum(V_tot.^2,2)),1,size(V_tot,2));
Clus_ind_wd123n=litekmeans(V_tot,CCC,'Replicates',50); 

nmi_wd123=[nmi_wd123,Cal_NMI(Clus_ind_wd123n, Clus_ind_wd123)]; pmi_wd123=[pmi_wd123,purity(CCC, Clus_ind_wd123, Clus_ind_wd123n)];
rmi_wd123=[rmi_wd123,RandIndex(Clus_ind_wd123n, Clus_ind_wd123)];
end

nnmi_wd123=[nnmi_wd123;nmi_wd123]; ppmi_wd123=[ppmi_wd123;pmi_wd123];
rrmi_wd123=[rrmi_wd123;rmi_wd123];
end

[mean(nnmi_wd123);mean(ppmi_wd123);mean(rrmi_wd123)]
[std(nnmi_wd123);std(ppmi_wd123);std(rrmi_wd123)]

%save('sim_22cancers_cls2_number.mat', 'nnmi_wd123','ppmi_wd123','rrmi_wd123')



%% generating stability graph (target cluster number) 
aaa=importdata('sim_22cancers_cls2_number.mat')
load('norm_cls_number_const.mat')

normalized_stab=mean(aaa.ppmi_wd123)./norm_cls_purity_const(1,:);
standard_stab=std(aaa.ppmi_wd123)./norm_cls_purity_const(1,:);

normalized_nmi=[mean(aaa.nnmi_wd123)];   
standard_nmi=[std(aaa.nnmi_wd123)];     

normalized_ari=[mean(aaa.rrmi_wd123)];  
standard_ari=[std(aaa.rrmi_wd123)];     

mean_velocity = normalized_stab;   
std_velocity =standard_stab;    

mean_velocity2=normalized_nmi;    
std_velocity2=standard_nmi;

mean_velocity3=normalized_ari;    
std_velocity3=standard_ari;


clf
figure
hold on
ttt1=bar(2:30,mean_velocity2);  
ttt2=errorbar(2:30,mean_velocity2,std_velocity2,'.'); ttt2.Color='b'
ttt=xlabel('Target cluster number'); ttt.FontSize=15;
xlim([1 31]); ylim([0 1])
ylabel('NMI', 'FontSize', 15)
alpha 0.5
set(gca,'FontSize', 14);

print -depsc normalized_simul_22cancer_nmi
