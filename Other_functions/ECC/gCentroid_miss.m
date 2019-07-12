function C = gCentroid_miss(IDX,index,K,n,r,sumKi,Ki)

%返回的C不会退化

    C = zeros(K,sumKi(r+1));  % 存储中心
    num = zeros(K,1);         % 存储每个簇的样本量
    maxKi= max(Ki);
    
    for k=1:K
        members = (index==k); % 第K个簇样本的逻辑向量
        
        if any(members) % 如果没有出现退化，计算中心
            num(k) = sum(members); % 如果没有缺失样本的话，该簇样本量
            idx = IDX(members,:); % 取出该簇样本
            counts = hist(idx,(0:maxKi)); % counts为(maxKi+1)*r矩阵，存储样本在各基础聚类中的簇分布，第一行为在各簇缺失样本数
            if size(counts,1) == 1
                for i = 1:r
                    if idx(i)>0 % 对于非缺失样本，对应位置为1
                        C(k,idx(i)+sumKi(i)) = 1;
                    else % 对于缺失样本，意味着中心退化，因此随机位置补1
                        C(k,randsample(Ki(i),1)+sumKi(i)) = 1;
                    end
                end
            else
                for i = 1:r
                    if counts(1,i)==num(k) % 说明该簇所有样本在基础聚类i中都是缺失样本，意味着中心退化
                        C(k,randsample(Ki(i),1)+sumKi(i)) = 1; % 由于中心退化，因此随机位置补1
                    else % 说明该簇样本在基础聚类i中并非全是缺失样本
                        C(k,sumKi(i)+1:sumKi(i+1)) = counts(2:Ki(i)+1,i)'/...
                            (num(k)-counts(1,i));
                    end
                end
            end
            
        else % 如果出现退化空簇，找一个样本成为新的中心
            C(k,:) = sCentroid_miss(IDX(randsample(n,1),:),1,r,Ki,sumKi);
        end
        
    end
end