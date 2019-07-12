function C = gCentroid(IDX,index,K,n,r,sumKi,Ki)
    C = zeros(K,sumKi(r+1));  
    num = zeros(K,1);         
    maxKi= max(Ki);           
    
    for k=1:K
        members = (index==k);
        
        if any(members) 
            num(k) = sum(members);
            idx = IDX(members,:); 
            counts = hist(idx,(1:maxKi)); 
            if size(counts,1) == 1
                C(k,idx+sumKi(1:r)') = 1;
            else
                for i = 1:r
                    C(k,sumKi(i)+1:sumKi(i+1)) = counts(1:Ki(i),i)'/num(k);
                end
            end
            
        else 
            C(k,:) = sCentroid(IDX(randsample(n,1),:),1,r,sumKi);
        end
    end
end