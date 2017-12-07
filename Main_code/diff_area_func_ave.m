function[diff_area]=diff_area_func_ave(CCC, pd1_A_set,pd1_B_set, bord, bordmin);
%% This code computes the differences of the fitted survival curves of each cluster in terms of curves of all patients

%% Input
% CCC is the number of cluster
% pd1_A_set and pd1_B_set are inferred parameters (scale and shape) of weibull distribution 
% bord and bordmin are the upper and lower limit time points used in the computation

%% Ouput
% diff_area is the CCC-dimensional vector, where ith component represents
% the area difference between curves corresponding to the i-th inferred clusters and that of all patients



diff_area=zeros(1,CCC);

for ii=1:CCC; 
    lam1=pd1_A_set(ii); kk1=pd1_B_set(ii); 
    for jj=CCC+1;
    lam2=pd1_A_set(jj); kk2=pd1_B_set(jj);   
    int_x=(lam2^(-kk2) *lam1^(kk1))^(1/(kk1-kk2));
    if int_x>bord
    area1=gammainc(lam1^(-kk1) *  bord^(kk1),   1/kk1)*lam1*gamma(1/kk1)/kk1 - gammainc(lam1^(-kk1) *  bordmin^(kk1),   1/kk1)*lam1*gamma(1/kk1)/kk1; area1=abs(area1);
    area2=gammainc(lam2^(-kk2) *  bord^(kk2),   1/kk2)*lam2*gamma(1/kk2)/kk2 - gammainc(lam2^(-kk2) *  bordmin^(kk2),   1/kk2)*lam2*gamma(1/kk2)/kk2; area2=abs(area2);
    diff_area(ii)=abs(area1-area2);
    else
    area1=gammainc(lam1^(-kk1) *  int_x^(kk1),   1/kk1)*lam1*gamma(1/kk1)/kk1 - gammainc(lam1^(-kk1) *  bordmin^(kk1),   1/kk1)*lam1*gamma(1/kk1)/kk1; area1=abs(area1);
    area2=gammainc(lam2^(-kk2) *  int_x^(kk2),   1/kk2)*lam2*gamma(1/kk2)/kk2 - gammainc(lam2^(-kk2) *  bordmin^(kk2),   1/kk2)*lam2*gamma(1/kk2)/kk2; area2=abs(area2);
    area_mid1=abs(area1-area2);
    area1=gammainc(lam1^(-kk1) *  int_x^(kk1),   1/kk1)*lam1*gamma(1/kk1)/kk1 - gammainc(lam1^(-kk1) *  bord^(kk1),   1/kk1)*lam1*gamma(1/kk1)/kk1; area1=abs(area1);
    area2=gammainc(lam2^(-kk2) *  int_x^(kk2),   1/kk2)*lam2*gamma(1/kk2)/kk2 - gammainc(lam2^(-kk2) *  bord^(kk2),   1/kk2)*lam2*gamma(1/kk2)/kk2; area2=abs(area2);
    area_mid2=abs(area1-area2);  
    diff_area(ii)=area_mid1+area_mid2;
    end
    end
end
        
        
   
   