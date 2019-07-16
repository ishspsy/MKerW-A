function labs=sc(X, ccc, varargin)
% simple spec cluster
[U,D,V]=svd(X);
spec=U(:,1:ccc);
labs=litekmeans(spec,ccc,'Replicates',20);