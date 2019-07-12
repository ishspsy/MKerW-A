function IDX = BasicCluster_RPS(Data,r,K,dist,randKi)
    [n,p] = size(Data);
    IDX = zeros(n,r);
    [n1,p1] = size(randKi);
    
    if n1>1
        Ki = randKi; % here randKi is the given Ki   
    elseif randKi==1&&sqrt(n)>K
        Ki = randsample(2:ceil(sqrt(n)),r,true); % here Ki is randomized
    else
        Ki = K*ones(r,1); % here Ki is equal to K
    end
    
    parfor i=1:r
        IDX(:,i) = kmeans(Data,Ki(i),'distance',dist,...
        'emptyaction','singleton','replicates',1);
    end  
end