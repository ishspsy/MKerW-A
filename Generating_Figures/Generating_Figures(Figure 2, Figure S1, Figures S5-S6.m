%%% This file generates Figure 2, Figure S1, and Figures S5-S6.

addpath(genpath(pwd))

%% Genearing Figure S1
load('similar_color_map.mat')
close all
clf
subplot(1,3,1)
colormap(flipud(hot))
imagesc(Wfc0s_euc_near_n{1}{1}{1}); caxis([0 0.2])
title('RNA')
colorbar;  

subplot(1,3,2)
colormap(flipud(hot))
imagesc(Wfc0s_euc_near_n{2}{1}{1});caxis([0 0.2])
title('miRNA')
colorbar;  

subplot(1,3,3)
colormap(flipud(hot))
imagesc(Wfc0s_euc_near_n{3}{1}{1});caxis([0 0.2])
title('CNA')
colorbar;  
fig = gcf;
fgP = fig.Position;
fgP(4) = fgP(4)/2.5;
set(gcf,'PaperPositionMode','auto');

print -depsc multi_cancer_pn_gph_1118



%% Genearing Figure 2
load('simul_threedata_22.mat')
aaa=importdata('simiclu.mat');

meanmat_mut1=[]; stdmat_mut1=[];
for ii=1:28; meanmat_mut1=[meanmat_mut1,mean(nmi_set{ii})];  stdmat_mut1=[stdmat_mut1,std(nmi_set{ii})]; end
meanmat_mut2=[]; stdmat_mut2=[];
for ii=1:28; meanmat_mut2=[meanmat_mut2,mean(pmi_set{ii})];  stdmat_mut2=[stdmat_mut2,std(pmi_set{ii})]; end
meanmat_mut3=[]; stdmat_mut3=[];
for ii=1:28; meanmat_mut3=[meanmat_mut3,mean(rmi_set{ii})];  stdmat_mut3=[stdmat_mut3,std(rmi_set{ii})]; end


meanmat_mut=[meanmat_mut1;meanmat_mut2;meanmat_mut3];
stdmat_mut=[stdmat_mut1;stdmat_mut2;stdmat_mut3];
    
hindex={'Cons-R' 'Cons-M' 'Cons-C' 'Cons-A' 'C-A' 'P-A' 'K-R' 'K-M' 'K-C' 'K-A' 'S-R' 'S-M'  'S-C'  'S-A' 'Ker-A' 'MKerW-A' 'MKerequ-A' 'MKer-fst-ALL' 'MKer-A' 'MequKer-fst-A' 'SIM-R' 'SIM-M'  'SIM-C' 'SIM-A'  'SS-R' 'SS-M' 'SS-C' 'SS-A'}
hind2=[1:16, 19, 21:24];
hindex=hindex(hind2); meanmat_mut=meanmat_mut(:,hind2);  stdmat_mut=stdmat_mut(:,hind2);
meanmat_mut2=[hindex;num2cell(meanmat_mut)]

%nmi
meanmat_mut1=meanmat_mut(1,:);
meanmat_mut2=[1:length(hind2); (meanmat_mut1);meanmat_mut1];
meanmat_mut22=sortrows(meanmat_mut2',-2)
hind22=(meanmat_mut22(:,1));
hindex(hind22)
std_val1=stdmat_mut(1,hind22);

clear vars title

clf
hold on
val1=meanmat_mut1(hind22) 
ttt1=bar(1:21,val1);  
ttt2=errorbar(1:21,val1,std_val1,'.'); ttt2.Color='b'
set(gca,'XTick',1:1:21)
xlim([0 22]);  ylim([min(val1)-0.1  min(1,max(val1)+0.1)])    
ylabel('NMI', 'FontSize', 15)
ttt=xlabel('Method'); ttt.FontSize=15;
alpha 0.5
xticklabels(hindex(hind22));  
xtickangle(45)
htext=text(-0.03,1.03,'B','Units','normalized')
htext.FontSize=20;
set(gca,'FontSize', 14);
print -depsc multi_cancer_nmi_bar 


%% Genearing Figure S5
%purity
meanmat_mut1=meanmat_mut(2,:);
meanmat_mut2=[1:length(hind2); (meanmat_mut1);meanmat_mut1];
meanmat_mut22=sortrows(meanmat_mut2',-2)
hind22=(meanmat_mut22(:,1));
hindex(hind22)
std_val1=stdmat_mut(2,hind22);

clear vars title

clf
hold on
val1=meanmat_mut1(hind22) 
ttt1=bar(1:21,val1);  
ttt2=errorbar(1:21,val1,std_val1,'.'); ttt2.Color='k'
set(gca,'XTick',1:1:21)
xlim([0 22]);  ylim([min(val1)-0.1  min(1,max(val1)+0.1)])    
ylabel('Purity', 'FontSize', 15)
ttt=xlabel('Method'); ttt.FontSize=15;
alpha 0.5
xticklabels(hindex(hind22));   
xtickangle(45)
set(gca,'FontSize', 14);
print -depsc multi_cancer_purity_bar 


%% Genearing Figure S6
%ARI
meanmat_mut1=meanmat_mut(3,:);
meanmat_mut2=[1:length(hind2); (meanmat_mut1);meanmat_mut1];
meanmat_mut22=sortrows(meanmat_mut2',-2)
hind22=(meanmat_mut22(:,1));
hindex(hind22)
std_val1=stdmat_mut(3,hind22);

clear vars title

clf
hold on
val1=meanmat_mut1(hind22) 
ttt1=bar(1:21,val1);  
ttt2=errorbar(1:21,val1,std_val1,'.'); ttt2.Color='k'
set(gca,'XTick',1:1:21)
xlim([0 22]);  ylim([min(val1)-0.1  min(1,max(val1)+0.1)])   
ylabel('ARI', 'FontSize', 15)
ttt=xlabel('Method'); ttt.FontSize=15;
alpha 0.5
xticklabels(hindex(hind22));  
xtickangle(45)
set(gca,'FontSize', 14);
print -depsc multi_cancer_ari_bar 
