function C = sCentroid_miss(idx,K,r,Ki,sumKi)

%因为要频繁调用，因此把函数写成两个，一个用于处理缺失样本，另一个处理无缺失样本
%这个函数适用于缺失样本

    C = zeros(K,sumKi(r+1));
    for l = 1:K
        for i = 1:r
            if idx(l,i)>0 % 对于非缺失样本，对应位置为1
                C(l,idx(l,i)+sumKi(i)) = 1;
            else % 对于缺失样本，意味着中心退化，因此随机位置补1
                C(l,randsample(Ki(i),1)+sumKi(i)) = 1;
            end
        end
    end
    
end