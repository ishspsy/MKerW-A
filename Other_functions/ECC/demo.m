load('Alizadeh-2000-v1.mat');

U = {'U_H','std',[]};   
K = length(unique(gnd)); % number of clusters for consensus clustering
r = 100;
w = ones(r,1); % the weight of each partitioning
rep = 10; % the number of ECC runs
maxIter = 40;
minThres = 1e-5;
utilFlag = 0;


IDX = BasicCluster_RPS(fea,r,K,'sqEuclidean',1);%Generate basic partitions
[pi_sumbest,pi_index,pi_converge,pi_utility,t] = RunECC(IDX,K,U,w,rep,maxIter,minThres,utilFlag); % run KCC for consensus clustering

[~, ~, Rn, NMI] = exMeasure(pi_index,gnd)
