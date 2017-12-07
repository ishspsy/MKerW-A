function P = gClusterDistribution(IDX,Ki,n)
    maxKi = max(Ki);
    counts = hist(IDX,0:maxKi);
    P = counts(2:maxKi+1,:)./repmat(n-counts(1,:),maxKi,1); 
end