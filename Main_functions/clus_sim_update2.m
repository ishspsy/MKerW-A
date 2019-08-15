function  [P_set,V_set,V_tot, ck, W_set,Wg_set]=clus_sim_update2(CCC, c,rho, n, K, p,q, gg, lam, mu, eta, W_euc_double)   % W_euc_nearest_double
%% This function solves our algorithm but using the equal weight to each omic data.


%% Input
% CCC is the target cluster number
% c is the small constant (represet \epsilon in the manuscript)
% rho is the penalty parameter controling multiple similarity matrices
% n is the number of samples (patients)
% K is the number of omics data sets considered (3 in our setting)
% p and q are number of parameters of \sigma and k used to construct similarity matrices (we use p=5 and q=11, so total 5*11=55 kernels)
% gg is the vector of ones (this is just for convenience when generalizing the algorithm to the case with feature selection) 
% lam and mu are main regularization parameters
% eta is the step-size in ADMM algorithm
% W_euc_double is the set of similarity matrices using multiple Gaussian kernels


%% Output
% P_set is the set of the obtained K target matrices
% V_set is the set of the first CCC number of eigenvectors of the obtained target matrices
% V_tot is the matrices constructed by merging components of V_set 
% ck is the learned weight for each omic data. It always have ck=ones(1,K) as the all data set is used equally.




%initial step
Wi_set=cell(1,K);
for k=1:K 
      Wi_set{k}=ones(p,q)/(p*q);
end
       
Wgi_set=cell(1,K);
for k=1:K
     Wgi_set{k}=ones(1,gg(k))/(gg(k));
end


kernel_ini_set=cell(K,max(gg));   kernel_ini_set2=cell(1,K);
for k=1:K
    ker_in=0;
    for jj=1:gg(k)
    ave_ker=0;
    for lll=1:(p*q);
    ave_ker=ave_ker+W_euc_double{k}{jj}{lll}/(p*q);
    end
    kernel_ini_set{k}{jj}=ave_ker;   %dimension of W_euc_double is K*gg*55 cells
    ker_in=ker_in+kernel_ini_set{k}{jj}/gg(k);
    end
    kernel_ini_set2{k}=ker_in;
end



Q_set=zeros(K,n,n); Gamma_set=cell(1,K);  
for k=1:K
[V, temp, evs]=eig1(eye(n)-kernel_ini_set2{k}, CCC, 0); 
LL= V(:,1:CCC) ;  Q=LL*LL';  Q_set(k,:,:)=Q;  Gamma_set{k}=zeros(n,n);
end

cki=ones(1,K);  %%

err=10;  Pi_set=Q_set; Qi_set=Q_set; Gammai_set=Gamma_set; W_set=Wi_set; Wg_set=Wgi_set; ck=cki;
rep=0; obj_fun=[];

while  (err>0.001) + (rep<10) >1
rep=rep+1;

%update P,Q,Gamma
[rep2, P_set, Q_set, Gamma_set, err0]=clus_sim_update0_2HW(CCC, ck, c, rho, lam, mu, eta, kernel_ini_set2);

err_P=squeeze(sum(sum(sum((P_set-Pi_set).^2))));


      
W_set=cell(1,K); err_W=0; 
for kk=1:K
  W0=zeros(p,q); W=zeros(p,q);
     for ii=1:(p*q)
         mid_trace=0;
         for jj=1:gg(k)
             mid_trace=mid_trace+ trace(ck(kk)*Wg_set{kk}(jj)*W_euc_double{kk}{jj}{ii}*squeeze(P_set(kk,:,:))/rho);
         end
                W0(ii)=exp(mid_trace);   
     end
     sw0=sum(sum(W0));
     for ii=1:(p*q)
       W(ii)=W0(ii)/ sw0;
     end
     sisnw=(sum(sum(isnan(W))));
       W(isnan(W))=  1/sisnw;
       sw2=sum(sum(W));
       for ii=1:(p*q)
       W(ii)=W(ii)/ sw2;
       end
       W_set{kk}=W;  err_W=err_W + norm(W_set{kk}- Wi_set{kk},'fro')^2;
end



 Wg_set=cell(1,K); err_Wg=0;
for kk=1:K
  W0=zeros(1,gg(kk)); W=zeros(1,gg(kk));
         for jj=1:gg(kk)
         mid_trace=0;
         for ii=1:(p*q)
             mid_trace=mid_trace+ trace(ck(kk)*W_set{kk}(ii)*W_euc_double{kk}{jj}{ii}*squeeze(P_set(kk,:,:))/rho);
         end
                W0(jj)=exp(mid_trace);   
         end
     sw0=sum(sum(W0));
     for jj=1:gg(kk)
       W(jj)=W0(jj)/ sw0;
     end
     sisw=(sum(sum(isnan(W))));
       W(isnan(W))= 1/sisw;
       swsw=sum(sum(W));
       for jj=1:gg(kk)
       W(jj)=W(jj)/ swsw;
       end
       Wg_set{kk}=W;  err_Wg=err_Wg + norm(Wg_set{kk}- Wgi_set{kk},'fro')^2;
end



ck=ones(1,K); cki=ones(1,K);  %%
err_ck=0;   %err_ck+norm(ck- cki,'fro')^2;



kernel_ini_set2=cell(1,K);
for kk=1:K
    ker_in=0;
    for jj=1:gg(kk)
    for lll=1:(p*q);
    ker_in=ker_in+W_euc_double{kk}{jj}{lll}*Wg_set{kk}(jj)*W_set{kk}(lll);
    end
    end
    kernel_ini_set2{kk}=ker_in;
end


err=err_P+err_W+err_Wg+err_ck; 
Wi_set=W_set; Wgi_set=Wg_set; Pi_set=P_set; Gammai_set=Gamma_set; Qi_set=Q_set; cki=ck;
err
end

%choose the first CCC eigenvectors
V_set=cell(1,K); V_tot=[];
for k=1:K
    [V, temp, evs]=eig1(squeeze(P_set(k,:,:)), CCC);   V=V./ repmat(sqrt(sum(V.^2,2)),1,CCC);   
    V_set{k}=V;  V_tot=[V_tot, V];
end
    


        
