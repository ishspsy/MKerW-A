function IDX = BasicCluster_RFS(Data,r,K,dist,nFeature)    
    [n,p] = size(Data);
    IDX = zeros(n,r);
    
    parfor i=1:r
        IDX(:,i) = kmeans(Data(:,randsample(p,nFeature)),K,...
        'distance',dist,'emptyaction','singleton','replicates',1);
    end
end