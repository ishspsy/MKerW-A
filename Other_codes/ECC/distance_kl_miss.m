function D = distance_kl_miss(U,C,weight,n,r,K,sumKi,binIDX,M)
    D = zeros(n,K);
    
    for l=1:n
        idx = M(l,:);
        D(l,:) = -(log2(C(:,binIDX(l,idx))+eps)*weight(idx))'; % 前面有负号
    end
end