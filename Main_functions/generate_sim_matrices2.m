function[Wfc0s_euc_near_n]=generate_sim_matrices(K,data_set3,gg, true_labs,sigma_set,k_set);
%% This code generates multiple similarity matrices for each data set.

%% Input
% K is the number of data sets used in clustering analysis
% data_set3 is the set of data sets
% gg is the K-dimensional vector such that k-th element represents a number of subgroups of k-th data in clustering analysis. Since we do not split each data, we use gg=ones(1,K).
% true_labs is true labels.
% sigma_set is the set of \sigma values (a variance parameter used in Gaussian kernel) used in constructing similarity matrices.
% k_set is the set of k values (a nearest-neighbor parameter used in Gaussian kernel) used in constructing similarity matrices.

%% Output
% Wfc0s_euc_near_n is the set of similarity matrices for each data. e.g., Wfc0s_euc_near_n{1} shows simialrity matrices for the first data.


Wfc0s_euc_near_n=cell(1,K);
for iii=1:K
XX=data_set3{iii};
for gi=1:gg(iii)
X=XX;
[n p]=size(X); N=n;

Diff = sqrt(dist2(X));
[DDist,~]=sort(Diff,2);

W_euc_orig=pdist(X); W_euc0=squareform(W_euc_orig);
W_euc_nearest_set_n=cell(length(sigma_set), length(k_set));

for tts=1:length(sigma_set)
for ttk=1:length(k_set)
sigma=sigma_set(tts); cln=k_set(ttk);
muij=mean(DDist(:,2:cln),2); muij=max(muij,0.1);  logn=round(cln);
muij=repmat(muij,[1,length(muij)])+repmat(muij',[length(muij),1]);  muij=muij*0.5;
W_euc=exp(-(W_euc0.^2)./(2*(muij*sigma).^2 ))./ (muij*sigma*sqrt(2*pi))  ;
W_ind=zeros(N,N);
for ii=1:N
W_ind(ii,find(W_euc(ii,:)>quantile(W_euc(ii,:),(n-logn)/n)))=1;
end
W_ind=(W_ind+W_ind')*0.5; W_ind(W_ind>0)=1;
W_euc_nearest=single(W_euc.*W_ind);
WW2= W_euc_nearest; DWW2= (max(sum(WW2),eps)).^(-0.5); WW2=  single(diag(DWW2)*WW2*diag(DWW2));
W_euc_nearest_set_n{tts,ttk}=WW2;
end
end
Wfc0s_euc_near_n{iii}{gi}=W_euc_nearest_set_n; 
end
end
