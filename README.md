# MKerW-A: Integrating multidimensional data for clustering analysis



## Overview

*MKerW-A* is a novel multi-view spectral clustering framework to integrate different omics data types measured from the same subjects by treating each omic data type as a different informative representation between patients. *MKerW-A* learns the weight of each data type as well as a similarity measure between patients via a non-convex optimization framework. It solves the proposed non-convex problem iteratively using the ADMM algorithm.


### Main functions

- [generate_sim_matrices2.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_code/generate_sim_matrices2.m)
: Contruct multiple similarity matrices using Gaussian kernels (Step 1).

- [clus_sim_update2_2HW.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_code/clus_sim_update2_2HW.m)
: Main *MKerW-A* algorithm solving optimization problem (Step 2).

- [clus_sim_update0_2HW.m](https://github.com/ishspsy/MKerW-A/blob/master/Main_code/clus_sim_update0_2HW.m)
: Algorithm solving embedded ADMM algorithm in Step 2.

- [Example (BRCA)](https://github.com/ishspsy/MKerW-A/blob/master/example_BRCA.m)
: Simple example file of *MKerW-A* using BRCA cancer.


### Example files

Please follow the links to reproduce the clustering results of TCGA data sets

-  [Clustering analysis for each cancer type](https://github.com/ishspsy/MKerW-A/blob/master/main_real.m)
: Generate clustering results for each of the 22 cancer types.

-  [Clustering analysis for identifying 22 cancer types](https://github.com/ishspsy/MKerW-A/blob/master/main_simul.m)
: Generate clustering results related to identifying 22 cancer types.

-  [Analysis on target cluster number](https://github.com/ishspsy/MKerW-A/blob/master/main_simul_cls.m)
: Choose target cluster number.

-  [Sensitivity analysis](https://github.com/ishspsy/MKerW-A/blob/master/main_simul_robust.m)
: Sensitivity analysis with respect to changes of reglarization parameters.

-  [Sensitivity test w.r.t. additive noise](https://github.com/ishspsy/MKerW-A/blob/master/main_simul_pert.m)
: Sensitivity analysis with respect to additive noise.


**Note** Most of the simulations and TCGA data applications were implemented on an Apple MacBook Pro (2.7 GHz, 8 GB of memory) using the MATLAB 2016b. 







### Directory

- All the functions used in the proposed algorithm *MKerW-A* are located in the directory ["**Main_Code**"](https://github.com/ishspsy/MKerW-A/tree/master/Main_code).

- All the other supplementary files are located in the directory ["**Other_codes**"](https://github.com/ishspsy/MKerW-A/tree/master/Other_codes).

- All the codes generating figures presented in the manuscript are located in the directory ["**Generating_Figures**"](https://github.com/ishspsy/MKerW-A/tree/master/Generating_Figures).

- All the resulting files (e.g. .MAT and .eps) are located in the directory ["**Results_files**"](https://github.com/ishspsy/MKerW-A/tree/master/Resulting_files).




## Example data sets

The 22 TCGA cancer data sets saved in the matlab file can be obtained from the dropbox directory [**all_data.mat**](https://www.dropbox.com/s/v22fx0j2gnpeta6/all_data.mat?dl=0). 

Specifically, *all_clin* includes clinical information of the patients from the 22 cancer types.

- *all_exp* is the RNA data for all the patinets.

- *all_mirna* is the MicroRNA data for all the patinets.

- *all_cna* is the CNA data for all the patinets.

- *all_pat* is the index vector of patients indicating corresponding cancer types (See the first column of *all_clin* for the original cancer name).




## DOWNLOAD

We provide MATLAB implementations of *MKerW-A* in the MKerW-A branch. The 22 TCGA cancer data sets saved in the matlab file can be downloaded from the dropbox directory [**all_data.mat**](https://www.dropbox.com/s/v22fx0j2gnpeta6/all_data.mat?dl=0). 



## Authors

* [**Seyoung Park**](http://people.yale.edu/search/seyoung_park.profile), Hao Xu, and   [**Hongyu Zhao**](https://publichealth.yale.edu/biostat/people/hongyu_zhao.profile)

  Department of Biostatistics, School of Public Health, Yale University


## Contact

* seyoung.park@yale.edu

## License

This project is licensed under the MIT License.




