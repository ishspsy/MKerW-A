# Overview

This directory includes main functions that directly related to the proposed clustering algorithm.

## Codes

- [generate_sim_matrices.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/generate_sim_matrices.m)
: this code constructs multiple similarity matrices using Gaussian kernels (Step 1).

- [clus_sim_update.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/clus_sim_update.m)
: this code is the main *MKerW-A* algorithm solving optimization problem (Step 2).

- [clus_sim_update_embedded.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/clus_sim_update_embedded.m)
: this code solves the embedded ADMM algorithm in Step 2.

- [clus_sim_update2.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_functions/clus_sim_update2.m)
: this code solves the proposed algorithm but using the equal weight to each omic data. 






	


