function simu_bp()
% fname=fopen('simu_bp.txt','w+');
%% paras
k_real=4; 
p=400; 
s=50; % number of informative features in p 
n=200;
sigma=6; % signal to noise ratio
repeat=2; % main repeat times
algo=@sc;

type=1; % type 1 or 0
bp_repeat=5; % bp repeat times
cmin=2;
cmax=8;

k_knn=10;

disp('bp_clr')

for run=1:repeat
    % Generate data
    [Bs,~,~]=svd( random('Normal',0,1,s,s) );
    Bs=Bs(1:k_real,:); % k_real*s
    B=[ sigma*Bs, zeros(k_real,p-s) ]; %k_real*p
    true_labs=random('unid',k_real,n,1);
    W=random('Normal',0,1,n,p);
    Z=zeros(n,k_real);
    for class=1:k_real
        Z(:,class)= (true_labs==class);
    end
    Y=Z*B+W;
    distance=bp_clr(type,bp_repeat,algo,Y,cmin,cmax);
    bp_score=mean(distance');
    % disp(bp_score)
    [~,imin]=min(bp_score);
    imin=imin+cmin-1;
    disp(imin)
end

disp('bp_clr_knn')

for run=1:repeat
    % Generate data
    [Bs,~,~]=svd( random('Normal',0,1,s,s) );
    Bs=Bs(1:k_real,:); % k_real*s
    B=[ sigma*Bs, zeros(k_real,p-s) ]; %k_real*p
    true_labs=random('unid',k_real,n,1);
    W=random('Normal',0,1,n,p);
    Z=zeros(n,k_real);
    for class=1:k_real
        Z(:,class)= (true_labs==class);
    end
    Y=Z*B+W;
    distance=bp_clr(type,bp_repeat,algo,Y,cmin,cmax,k_knn);
    bp_score=mean(distance');
    % disp(bp_score)
    [~,imin]=min(bp_score);
    imin=imin+cmin-1;
    disp(imin)
end

%% output
