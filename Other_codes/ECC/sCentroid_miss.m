function C = sCentroid_miss(idx,K,r,Ki,sumKi)

%��ΪҪƵ�����ã���˰Ѻ���д��������һ�����ڴ���ȱʧ��������һ��������ȱʧ����
%�������������ȱʧ����

    C = zeros(K,sumKi(r+1));
    for l = 1:K
        for i = 1:r
            if idx(l,i)>0 % ���ڷ�ȱʧ��������Ӧλ��Ϊ1
                C(l,idx(l,i)+sumKi(i)) = 1;
            else % ����ȱʧ��������ζ�������˻���������λ�ò�1
                C(l,randsample(Ki(i),1)+sumKi(i)) = 1;
            end
        end
    end
    
end