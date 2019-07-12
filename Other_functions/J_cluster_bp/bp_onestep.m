function [bp,labs]=bp_onestep(algo,X,ccc,varargin)
% first draw a bootstrap sample: bp, which is the bootstrap sample of (1:n)
% then return the cluster label of X(bp,:)
[n,p]=size(X);
bp=sort( randsample(n,n,true) );
X_bp=X(bp,:);
labs=algo(X_bp,ccc,varargin);
end
