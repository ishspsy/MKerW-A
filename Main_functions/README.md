
# Overview

This directory includes main functions used in main clustering analysis

## Codes

- [generate_sim_matrices2.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/generate_sim_matrices2.m)
: this code constructs multiple similarity matrices using Gaussian kernels (Step 1).

- [clus_sim_update2_2HW.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/clus_sim_update2_2HW.m)
: this code is the main *MKerW-A* algorithm solving optimization problem (Step 2).

- [clus_sim_update0_2HW.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/clus_sim_update0_2HW.m)
: this code solves embedded ADMM algorithm in Step 2.

- [clus_sim_update2_2HW0.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/clus_sim_update2_2HW0.m)
: this code solves our algorithm but using the equal weight to each omic data. 

- [Contingency.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/Contingency.m)
: this code forms contingency matrix for two vectors. Used for computing ARI value.

- [diff_area_func_ave.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/diff_area_func_ave.m)
: this code computes the differences of the fitted survival curves of each cluster in terms of curves of all patients.

- [diff_area_func.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/diff_area_func.m)
: this code computes the absolute area between two fitted survival curves (used to compute 'L1-dif').

- [diff_area_func2.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/diff_area_func2.m)
: this code computes the relative area of the fitted survival curves (used to compute 'L1-dif').     

- [generate_surv_func_general.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/generate_surv_func_general.m)
: this function provides the fitted Weibull survival curve for each inferred group. 

- [purity.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/purity.m)
: this code computes Purity between two clustering labels.

- [RandIndex.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/RandIndex.m)
: this code calculates Rand Indices (ARI) to compare two partitions. This code is by David Corney (2000) 

- [reg_pca.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/reg_pca.m)
: this code performs PCA analysis. 

- [cal_normalized_purtiy.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/cal_normalized_purtiy.m)
: this code computes average values and standard deviations of the cluster performances (e.g. NMI, purity, ARI)
for the randomly assigned clusters given each target cluster numbers from 2 to 30.





	


