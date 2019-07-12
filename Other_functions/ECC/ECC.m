function [sumbest,index,converge,utility] =ECC(IDX,K,U,w,weight,...
    distance,maxIter,minThres,utilFlag,missFlag,missMatrix,...
    n,r,Ki,sumKi,binIDX,Pvector)

if (utilFlag==1 && missFlag==1) 
    C = sCentroid_miss(IDX(randsample(n,K),:),K,r,Ki,sumKi); 
    sumbest = inf;
    converge = zeros(100,1)-1;
    utility = zeros(100,2)-1;
            
    for i = 1:maxIter 
        D = feval(distance,U,C,weight,n,r,K,sumKi,binIDX,missMatrix);
        [d, idx] = min(D, [], 2); 
        totalsum = sum(d); 
        
        if abs(sumbest - totalsum) < minThres 
            break;           
        elseif totalsum < sumbest 
            index = idx; 
            C = gCentroid_miss(IDX,index,K,n,r,sumKi,Ki); 
            sumbest = totalsum; 
            converge(i) = sumbest; 
            utility(i,:) = (UCompute_miss(index,U,w,C,n,r,K,sumKi,Pvector,missMatrix))'; 
        else 
            warning('ECC:OptimizationException',...
                'The objective function value increases.');
            break;
        end 
    end
    
elseif  (utilFlag==1 && missFlag==0) 
    C = sCentroid(IDX(randsample(n,K),:),K,r,sumKi); 
    sumbest = inf;
    converge = zeros(100,1)-1;
    utility = zeros(100,2)-1;
            
    for i = 1:maxIter
        D = feval(distance,U,C,weight,n,r,K,sumKi,binIDX);
        [d, idx] = min(D, [], 2); 
        totalsum = sum(d); 
        
        if abs(sumbest - totalsum) < minThres
            break;           
        elseif totalsum < sumbest 
            index = idx; 
            C = gCentroid(IDX,index,K,n,r,sumKi,Ki);
            sumbest = totalsum; 
            converge(i) = sumbest; 
            utility(i,:) = (UCompute(index,U,w,C,n,r,K,sumKi,Pvector))'; 
        else 
            warning('ECC:OptimizationException',...
                'The objective function value increases.');
            break;
        end 
    end
    
elseif (utilFlag==0 && missFlag==1) 
    C = sCentroid_miss(IDX(randsample(n,K),:),K,r,Ki,sumKi); 
    sumbest = inf;
    converge = zeros(100,1)-1;
    utility = []; 
            
    for i = 1:maxIter
        D = feval(distance,U,C,weight,n,r,K,sumKi,binIDX,missMatrix);
        [d, idx] = min(D, [], 2); 
        totalsum = sum(d); 
        
        if abs(sumbest - totalsum) < minThres 
            break;           
        elseif totalsum < sumbest
            index = idx; 
            C = gCentroid_miss(IDX,index,K,n,r,sumKi,Ki); 
            sumbest = totalsum;
            converge(i) = sumbest; 
        else 
            warning('ECC:OptimizationException',...
                'The objective function value increases.');
            break;
        end 
    end
    
else 
    C = sCentroid(IDX(randsample(n,K),:),K,r,sumKi); 
    sumbest = inf;
    converge = zeros(100,1)-1;
    utility = []; 
            
    for i = 1:maxIter
        D = feval(distance,U,C,weight,n,r,K,sumKi,binIDX);
        [d, idx] = min(D, [], 2);
        totalsum = sum(d); 
        
        if abs(sumbest - totalsum) < minThres 
            break;           
        elseif totalsum < sumbest 
            index = idx; 
            C = gCentroid(IDX,index,K,n,r,sumKi,Ki);
            sumbest = totalsum; 
            converge(i) = sumbest; 
        else 
            warning('ECC:OptimizationException',...
                'The objective function value increases.');
            break;
        end 
    end
    
end

end