function [pi_sumbest,pi_index,pi_converge,pi_utility,t] = ...
    RunECC(IDX,K,U,w,rep,maxIter,minThres,utilFlag)
    tic; 
    rand('state', 0);    
    [n,r] = size(IDX);
    
    if nargin>8 
        error('RunECC:TooManyInputs',...
            'At most 8 input arguments: IDX,U,K,w,rep,maxIter,minThres,utilFlag.');
    elseif nargin<8
		utilFlag = 0;
    elseif nargin<7
        minThres = E-5;
    elseif nargin<6
        maxIter = 20;
    elseif nargin<5
        rep = 5;
    elseif nargin<4
        w = ones(r,1);
    elseif nargin<3
        U = {'u_h','std'};
    elseif nargin<2
        error('RunECC:TooFewInputs',...
            'At least 2 input arguments required: IDX,K.');
    end
    
    [Ki,sumKi,binIDX,missFlag,missMatrix,distance,Pvector,weight] = ...
        Preprocess(IDX,U,n,r,w,utilFlag); 

    l_sumbest = zeros(1,rep); 
    l_index = zeros(n,rep); 
    l_converge = zeros(100,rep);
        
    if utilFlag==1 
        l_utility = zeros(100,2*rep); 
   
        for p = 1:rep
            [sumbest,index,converge,utility] = ECC(IDX,K,U,w,weight,...
                distance,maxIter,minThres,utilFlag,missFlag,missMatrix,...
                n,r,Ki,sumKi,binIDX,Pvector); 
            l_sumbest(p) = sumbest;
            l_index(:,p) = index;
            l_converge(:,p) = converge;
            l_utility(:,(2*p-1):(2*p)) = utility;
        end
        
        [pi_sumbest,pos] = min(l_sumbest);
        pi_index = l_index(:,pos);
        pi_converge = l_converge(:,pos);
        pi_utility = l_utility(:,(2*pos-1):(2*pos));
        
    elseif utilFlag==0 
        for p = 1:rep
            [sumbest,index,converge] = ECC(IDX,K,U,w,weight,...
                distance,maxIter,minThres,utilFlag,missFlag,missMatrix,...
                n,r,Ki,sumKi,binIDX,Pvector);
            l_sumbest(p) = sumbest;
            l_index(:,p) = index;
            l_converge(:,p) = converge;
        end    
        [pi_sumbest,pos] = min(l_sumbest);
        pi_index = l_index(:,pos);
        pi_converge = l_converge(:,pos);
        pi_utility = [];
        
    else
        error('RunECC:UnknowUtilityFlag',...
                'utilFlag must be 1 or 0.');
        
    end

    t = toc;
end
