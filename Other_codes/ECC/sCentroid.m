function C = sCentroid(idx,K,r,sumKi) 
    C = zeros(K,sumKi(r+1));
    for l = 1:K
        C(l,idx(l,:)+sumKi(1:r)') = 1;
    end  
end