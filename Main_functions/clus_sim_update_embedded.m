function  [rep, P_set, Q_set, Gamma_set ,err]=clus_sim_update_embedded(CCC,ck,c, rho, lam, mu, eta, kernel_ini_set)  
%% This function solves the embedded ADMM algorithm.

%% Input
% CCC is the target cluster number
% ck is the updated weight of the data
% c is the small constant (represet \epsilon in the manuscript)
% rho is the penalty parameter controling multiple similarity matrices
% lam and mu are main regularization parameters
% eta is the step-size in ADMM algorithm
% kernel_ini_set is the set of updated similarity matrices

%% Output
% rep is the number of runs
% P_set is the updated target matrices
% Q_set is the updated matrices corresponding to the surrogate P_set used in ADMM algorithm
% Gamma_set is the updated Lagrangian multipliers
% err is the final error of the algorithm



K=length(kernel_ini_set);    %1 by K
n=size(kernel_ini_set{1},1);  % n by n

Q_set=zeros(K,n,n); Gamma_set=cell(1,K);
for k=1:K
[V, temp, evs]=eig1(eye(n)-kernel_ini_set{k}, CCC, 0); 
LL= V(:,1:CCC) ;  Q=LL*LL';  Q_set(k,:,:)=Q;  Gamma_set{k}=zeros(n,n);
end

err=10;  Pi_set=Q_set; Qi_set=Q_set; Gammai_set=Gamma_set;
rep=0; obj_fun=[];


while  (err>0.001) + (rep<10) >1
rep=rep+1;
%update P
P_set=zeros(K,n,n);  err_P=0;
for k=1:K
    const0=1/(2*c*ck(k)+eta+2*mu*(K-1));
    T=ck(k)*kernel_ini_set{k}-Gamma_set{k}+eta* squeeze(Q_set(k,:,:))+2*mu* (squeeze(sum(Q_set,1))-squeeze(Q_set(k,:,:)));
    P=const0.*(T-lam).*(T>lam) + const0.*(T+lam).*(T<-lam);
    P_set(k,:,:) = P;  
    err_P= err_P+ norm(squeeze(P_set(k,:,:)-Pi_set(k,:,:)),'fro')^2;
end
    

%update Q
Q_set=zeros(K,n,n); const=1/(eta+2*mu*(K-1)); err_Q=0;
for k=1:K
    T=2*mu* (squeeze(sum(P_set,1))-squeeze(P_set(k,:,:))) + Gamma_set{k}+ eta* squeeze(P_set(k,:,:));
    B=const*T; C=(B+B')*0.5;   [V, D, U]=eig(C);  [X]=projection(diag(D),CCC);
    Q=V*diag(X)*V';
    Q_set(k,:,:) = Q;
    err_Q= err_Q+ norm(squeeze(Q_set(k,:,:)-Qi_set(k,:,:)),'fro')^2;
end


%update Gamma
Gamma_set=cell(1,K);   err_Gamma=0;
for k=1:K
 Gamma_set{k}=  Gammai_set{k} + eta*  squeeze(P_set(k,:,:)- Q_set(k,:,:));
 err_Gamma=err_Gamma+norm(Gamma_set{k} - Gammai_set{k}, 'fro')^2; 
end


err=err_Gamma+ err_Q+ err_P;
Pi_set=P_set; Gammai_set=Gamma_set; Qi_set=Q_set;
err;
end

        
