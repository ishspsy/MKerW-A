function[pd1_A_set,pd1_B_set]=generate_surv_func_general(CCC, bord, bb_set1,surv_stat1, jjj, colors_set,gp_title11,gp_title22,cl_title, name_i)
%% This function provides the fitted Weibull survival curve for each inferred group.

%% Input
% CCC is the target cluster number
% bord is the end time (month) point for survival analysis.
% bb_set1 is the set of survival time of patients in each inferred group.
% surv_stat1 is the set of survival status (0 or 1, where 1 represents a censored one) of patients in each inferred group.
% jjj is the number of clustering methods used for comparisons.
% colors_set is the set of colors that will be used to drow survival curves.
% gp_title11 is the list of the names of considered clustering methods.
% gp_title22 is the set of the texts of legends that will be used to express fitted survival curves.

%% Output
% pd1_A_set is the fitted Weibull distribution parameters (scale)
% pd1_B_set is the fitted Weibull distribution parameters (shape)


bb_set=bb_set1{jjj}; for ii=1:CCC; bb_set{ii}(bb_set{ii}==0)=0.01; leng{ii}=length(bb_set{ii}); end
bb_surv=surv_stat1{jjj};

inc_size=[]; inc_size_sum=0;
for ii=1:CCC
    inc_size=[inc_size, length(bb_set{ii})]; inc_size_sum=inc_size_sum+length(bb_set{ii});
end
inc_size=[inc_size,inc_size_sum];

bb_set_all=[]; bb_surv_all=[];
for ii=1:CCC
    bb_set_all=[bb_set_all',bb_set{ii}']';  bb_surv_all=[bb_surv_all',bb_surv{ii}']'; 
end


ind_chose=1:CCC;  leng_chose=CCC;


pd1_A_set=[]; pd1_B_set=[];


figure
ax1 = gca;
[f1,x1]=ecdf(ax1,bb_set{ind_chose(1)},'Censoring',bb_surv{ind_chose(1)},'function','survivor');
stairs(x1,f1,colors_set{ind_chose(1)})
hold on
%pd1 = fitdist(bb_set{ind_chose(1)},'wbl','Censoring',bb_surv{ind_chose(1)}); pd1_A_set=[pd1_A_set,pd1.A];pd1_B_set=[pd1_B_set,pd1.B];
%plot(0:1:max(bb_set{ind_chose(1)}+10),1-cdf('wbl',0:1:max(bb_set{ind_chose(1)}+10),pd1.A,pd1.B), colors_set{ind_chose(1)}, 'LineWidth', 3);
%hold on
if leng_chose>1
    for ijk=2:leng_chose
[f2,x2] = ecdf(bb_set{ind_chose(ijk)},'Censoring',bb_surv{ind_chose(ijk)},'function','survivor');
stairs(x2,f2,colors_set{ind_chose(ijk)})
xlim([0 bord])
xlabel('Overall Survival')
ylabel('Survival Probability')
%pd1 = fitdist(bb_set{ind_chose(ijk)},'wbl','Censoring',bb_surv{ind_chose(ijk)});
%plot(0:1:max(bb_set{ind_chose(ijk)}+10),1-cdf('wbl',0:1:max(bb_set{ind_chose(ijk)}+10),pd1.A,pd1.B), colors_set{ind_chose(ijk)}, 'LineWidth', 3); hold on
    end
end
[f2,x2] = ecdf(bb_set_all,'Censoring',bb_surv_all,'function','survivor');
stairs(x2,f2,colors_set{CCC+1})
xlim([0 bord])
xlabel('Overall Survival')
ylabel('Survival Probability')
%pd1 = fitdist(bb_set_all,'wbl','Censoring',bb_surv_all);
%plot(0:1:max(bb_surv_all+10),1-cdf('wbl',0:1:max(bb_surv_all+10),pd1.A,pd1.B), colors_set{CCC+1}, 'LineWidth', 3);
hold on

[f1,x1]=ecdf(ax1,bb_set{ind_chose(1)},'Censoring',bb_surv{ind_chose(1)},'function','survivor');
%stairs(x1,f1,colors_set{ind_chose(1)})
%hold on
pd1 = fitdist(bb_set{ind_chose(1)},'wbl','Censoring',bb_surv{ind_chose(1)}); pd1_A_set=[pd1_A_set,pd1.A];pd1_B_set=[pd1_B_set,pd1.B];
plot(0:1:max(bb_set{ind_chose(1)}+10),1-cdf('wbl',0:1:max(bb_set{ind_chose(1)}+10),pd1.A,pd1.B), colors_set{ind_chose(1)}, 'LineWidth', 3);
hold on


if leng_chose>1
    for ijk=2:leng_chose
[f2,x2] = ecdf(bb_set{ind_chose(ijk)},'Censoring',bb_surv{ind_chose(ijk)},'function','survivor');
%stairs(x2,f2,colors_set{ind_chose(ijk)})
xlim([0 bord])
xlabel('Overall Survival')
ylabel('Survival Probability')
pd1 = fitdist(bb_set{ind_chose(ijk)},'wbl','Censoring',bb_surv{ind_chose(ijk)});
pd1_A_set=[pd1_A_set,pd1.A];pd1_B_set=[pd1_B_set,pd1.B];
plot(0:1:max(bb_set{ind_chose(ijk)}+10),1-cdf('wbl',0:1:max(bb_set{ind_chose(ijk)}+10),pd1.A,pd1.B), colors_set{ind_chose(ijk)}, 'LineWidth', 3); hold on
    end
end
[f2,x2] = ecdf(bb_set_all,'Censoring',bb_surv_all,'function','survivor');
%stairs(x2,f2,colors_set{CCC+1})
pd1 = fitdist(bb_set_all,'wbl','Censoring',bb_surv_all);
pd1_A_set=[pd1_A_set,pd1.A];pd1_B_set=[pd1_B_set,pd1.B];

plot(0:1:max(bb_set_all+10),1-cdf('wbl',0:1:max(bb_set_all+10),pd1.A,pd1.B), colors_set{CCC+1}, 'LineWidth', 3);




tt=title(sprintf('%s', gp_title11{jjj}));  %%0.12

string_vec=cell(1,length(ind_chose)+1);
for i4=1:(length(ind_chose)+1); string_vec{i4}='i4'; end


hlegend1=legend(ax1,  string_vec, 'Location','northeast');    
for i4=1:(length(ind_chose));  
hlegend1.String{i4}=[gp_title22{ind_chose(i4)} '-' num2str(inc_size(i4))];
end
hlegend1.String{CCC+1}=[gp_title22{CCC+1} '-' num2str(inc_size(CCC+1))];

%drawnow

print([sprintf('surv_curve_gp_%d_%s', name_i, gp_title11{jjj})],'-depsc')
image=2; 
%image=hlegend1;
