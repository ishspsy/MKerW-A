function C = gCentroid_miss(IDX,index,K,n,r,sumKi,Ki)

%���ص�C�����˻�

    C = zeros(K,sumKi(r+1));  % �洢����
    num = zeros(K,1);         % �洢ÿ���ص�������
    maxKi= max(Ki);
    
    for k=1:K
        members = (index==k); % ��K�����������߼�����
        
        if any(members) % ���û�г����˻�����������
            num(k) = sum(members); % ���û��ȱʧ�����Ļ����ô�������
            idx = IDX(members,:); % ȡ���ô�����
            counts = hist(idx,(0:maxKi)); % countsΪ(maxKi+1)*r���󣬴洢�����ڸ����������еĴطֲ�����һ��Ϊ�ڸ���ȱʧ������
            if size(counts,1) == 1
                for i = 1:r
                    if idx(i)>0 % ���ڷ�ȱʧ��������Ӧλ��Ϊ1
                        C(k,idx(i)+sumKi(i)) = 1;
                    else % ����ȱʧ��������ζ�������˻���������λ�ò�1
                        C(k,randsample(Ki(i),1)+sumKi(i)) = 1;
                    end
                end
            else
                for i = 1:r
                    if counts(1,i)==num(k) % ˵���ô����������ڻ�������i�ж���ȱʧ��������ζ�������˻�
                        C(k,randsample(Ki(i),1)+sumKi(i)) = 1; % ���������˻���������λ�ò�1
                    else % ˵���ô������ڻ�������i�в���ȫ��ȱʧ����
                        C(k,sumKi(i)+1:sumKi(i+1)) = counts(2:Ki(i)+1,i)'/...
                            (num(k)-counts(1,i));
                    end
                end
            end
            
        else % ��������˻��մأ���һ��������Ϊ�µ�����
            C(k,:) = sCentroid_miss(IDX(randsample(n,1),:),1,r,Ki,sumKi);
        end
        
    end
end