%% generating Fig 1 and S7
addpath(genpath(pwd))

aaa=importdata('sim_22cancers_cls2_number1.mat')
load('norm_cls_number_const.mat')

normalized_stab=mean(aaa.ppmi_wd123)./norm_cls_purity_const(1,:);
standard_stab=std(aaa.ppmi_wd123)./norm_cls_purity_const(1,:);
 

clf
figure
hold on
ttt1=bar(2:30, normalized_stab);   ttt1.FaceColor='k'
ttt2=errorbar(2:30, normalized_stab, standard_stab,'.');  ttt2.Color='k'
ttt=xlabel('Target cluster number'); ttt.FontSize=15;
xlim([1 31]); 
ylabel('Adjusted Purity', 'FontSize', 15)
alpha 0.2
set(gca,'FontSize', 14);
print -depsc normalized_simul_22cancer



%% generating stability graph (target cluster number) 
aaa=importdata('sim_22cancers_cls2_number2.mat')
load('norm_cls_number_const.mat')

normalized_stab=mean(aaa.ppmi_wd123)./norm_cls_purity_const(1,:);
standard_stab=std(aaa.ppmi_wd123)./norm_cls_purity_const(1,:);

clf
figure
hold on
ttt1=bar(2:30, normalized_stab);   ttt1.FaceColor='k'
ttt2=errorbar(2:30, normalized_stab, standard_stab,'.');  ttt2.Color='k'
ttt=xlabel('Target cluster number'); ttt.FontSize=15;
xlim([1 31]); 
ylabel('Adjusted Purity', 'FontSize', 15)
alpha 0.2
set(gca,'FontSize', 14);
print -depsc normalized_simul_22cancer2



