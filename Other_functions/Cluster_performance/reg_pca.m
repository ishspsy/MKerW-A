function X = reg_pca(in_X, K)
%% Perform PCA analysis.

%% Input
% in_X is the omic data
% K is the dimension of the embedded space

%% Output
% X is the obtained score matrix with K dimension


in_X = in_X - repmat(mean(in_X),size(in_X,1),1);
[U, S, ~] = svd(in_X);
K = min(size(S,2),K);
X = U(:,1:K)*diag(sqrt(diag(S(1:K,1:K))));
X = X./repmat(sqrt(sum(X.*X,2)),1,K);
end
