------------------------------------------------------------------------------------------
	               Readme for ECC package
	 		       version March 25, 2016
------------------------------------------------------------------------------------------

The package includes the MATLAB code of ECC package for gene expression data clustering.

This package was developed by Mr. Hongfu Liu (liu.hongf@husky.neu.edu). For any problem concerning the code, please feel free to contact Mr. Liu.

RunECC is the main function with the following inputs and outputs.
%  Input
%  IDX:       the set of basic partitions 
%  U:         the utility function (Here we use U = {'U_H','std',[]};) 
%  w:         the weight vector of basic partitions (Default setting w = ones(r,1), r is the number of basic partitions) 
%  rep:       the repetition time of K-means (Default setting rep = 20)
%  maxIter:   the maximum iteration number (Default setting maxIter = 50)
%  minThres:  the threshold for stopping criterion (Default setting minThres = 1e-5)
%  utilFlag:  the indicator to calculate the utility (Default setting utilFlag = 0)

%  Output
% pi_sumbest: the objective function value of K-means
% pi_index:   the partition for gene expression data
% pi_converge:the objective function value in each iteration
% pi_utility: the utility value (Here we do not calculate utility value)
% t:          the execution time         

We also provide two strategies to generate basic partitions, BasicCluster_RPS and BasicCluster_RFS. Besides, the evaluation function exMeasure including VIn, VDn, Rn and NMI, is also provided by this package

demo.m is an example.
