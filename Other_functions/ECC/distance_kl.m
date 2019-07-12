function D = distance_kl(U,C,weight,n,r,K,sumKi,binIDX)
    D = zeros(n,K);
    for l=1:n
        D(l,:) = -(log2(C(:,binIDX(l,:))+eps)*weight)';%ГЫвд -1
    end
end