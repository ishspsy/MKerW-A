function [Ki,sumKi,binIDX,missFlag,missMatrix,distance,Pvector,weight]=...
        Preprocess(IDX,U,n,r,w,utilFlag)
    Ki = max(IDX); 
    sumKi = zeros(r+1,1); 
    for i=1:r             
        sumKi(i+1) = sumKi(i)+Ki(i);
    end
    binIDX = IDX+repmat(sumKi(1:r)',n,1);
    
    if sum(any(IDX==0))>0  
        missFlag = 1; 
        missMatrix = IDX>0; 
    else 
        missFlag = 0;
        missMatrix = [];
    end
    
    if missFlag == 1
        switch lower(U{1,1})
            case 'u_c'
                distance = @distance_euc_miss;
            case 'u_h'
                distance = @distance_kl_miss;
            case 'u_cos'
                distance = @distance_cos_miss;
            case 'u_lp'
                distance = @distance_lp_miss;
            otherwise
                error('Preprocess:UnknowUtilityFunction',...
                    'Currently only support U_c,U_H,U_cos,U_Lp.');
        end
    else 
        switch lower(U{1,1})
            case 'u_c'
                distance = @distance_euc;
            case 'u_h'
                distance = @distance_kl;
            case 'u_cos'
                distance = @distance_cos;
            case 'u_lp'
                distance = @distance_lp;
            otherwise
                error('Preprocess:UnknowUtilityFunction',...
                    'Currently only support U_c,U_H,U_cos,U_Lp.');
        end
    end
    
    if (strcmpi(U{1,2},'norm')==1 || utilFlag==1) 
        P = gClusterDistribution(IDX,Ki,n); 
        switch lower(U{1,1}) 
            case 'u_c'
                Pvector = sum(P.^2);
            case 'u_h'
                Pvector = -sum(P.*log2(P+eps)); 
            case 'u_cos'
                Pvector = sqrt(sum(P.^2));
            case 'u_lp'
                p = U{1,3};
                Pvector = (sum(P.^p)).^(1/p);
            otherwise
                error('Preprocess:UnknowUtilityFunction',...
                    'Currently only support U_c,U_H,U_cos,U_Lp.');
        end
    else
        Pvector = [];
    end
    
    switch lower(U{1,2}) 
        case 'std'
            weight = w; 
        case 'norm'
            weight = w./Pvector'; 
        otherwise
             error('Preprocess:UnknownUtilityType',...
                'Only support two types of utility: std, norm');
    end
end