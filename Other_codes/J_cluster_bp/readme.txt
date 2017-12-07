Estimating Clustering instability using Bootstrap

"bp" always means bootstrap, it can refer to the bootstrap samples if it is the name of function output

"algo" in the code is the clustering algorithm to be evaluated
    algo should be a function like:
    labs=algo(X,ccc,varargin)

"bp_clr" is the main function of my method

"bp_clr" is the main function of the original method (using knn)

"bp_onestep" returns one bp sample and the labels of this sample

“simu_bp” is a simulation using the simplest gaussian mixture, print out the most stable cluster number. You can use it as a demo.