%% This file geneartes Fig 1(C), Fig S8, and Fig S9.

addpath(genpath(pwd))


load('simul_threedata_22_stability_mean.mat')


meanmat_mut1=[]; stdmat_mut1=[];
for ii=1:28; meanmat_mut1=[meanmat_mut1;mean(nnmi_set{ii}')];  stdmat_mut1=[stdmat_mut1;std(nnmi_set{ii}')]; end
meanmat_mut2=[]; stdmat_mut2=[];
for ii=1:28; meanmat_mut2=[meanmat_mut2;mean(ppmi_set{ii}')];  stdmat_mut2=[stdmat_mut2;std(ppmi_set{ii}')]; end
meanmat_mut3=[]; stdmat_mut3=[];
for ii=1:28; meanmat_mut3=[meanmat_mut3;mean(rrmi_set{ii}')];  stdmat_mut3=[stdmat_mut3;std(rrmi_set{ii}')]; end


aaa=importdata('simiclu.mat');aaa1=importdata('simiclu1.mat'); aaa2=importdata('simiclu2.mat'); aaa3=importdata('simiclu3.mat'); 
%icluster
nmi_icl1=[]; nmi_icl2=[]; nmi_icl3=[];
for jj=1:50; nmi_icl1=[nmi_icl1,Cal_NMI(aaa(:,jj),aaa1(:,jj))];
 nmi_icl2=[nmi_icl2,Cal_NMI(aaa(:,jj),aaa2(:,jj))];
nmi_icl3=[nmi_icl3,Cal_NMI(aaa(:,jj),aaa3(:,jj))];  end
meanmat_mut1=[meanmat_mut1;[mean(nmi_icl1),mean(nmi_icl2),mean(nmi_icl3)]];
stdmat_mut1=[stdmat_mut1;[std(nmi_icl1),std(nmi_icl2),std(nmi_icl3)]];
  
nmi_icl1=[]; nmi_icl2=[]; nmi_icl3=[];
for jj=1:50; nmi_icl1=[nmi_icl1,purity(22,aaa(:,jj),aaa1(:,jj))];
 nmi_icl2=[nmi_icl2,purity(22,aaa(:,jj),aaa2(:,jj))];
nmi_icl3=[nmi_icl3,purity(22,aaa(:,jj),aaa3(:,jj))];  end
meanmat_mut2=[meanmat_mut2;[mean(nmi_icl1),mean(nmi_icl2),mean(nmi_icl3)]];
stdmat_mut2=[stdmat_mut2;[std(nmi_icl1),std(nmi_icl2),std(nmi_icl3)]];

nmi_icl1=[]; nmi_icl2=[]; nmi_icl3=[];
for jj=1:50; nmi_icl1=[nmi_icl1,RandIndex(aaa(:,jj),aaa1(:,jj))];
 nmi_icl2=[nmi_icl2,RandIndex(aaa(:,jj),aaa2(:,jj))];
nmi_icl3=[nmi_icl3,RandIndex(aaa(:,jj),aaa3(:,jj))];  end
meanmat_mut3=[meanmat_mut3;[mean(nmi_icl1),mean(nmi_icl2),mean(nmi_icl3)]];
stdmat_mut3=[stdmat_mut3;[std(nmi_icl1),std(nmi_icl2),std(nmi_icl3)]];


hindex={'Cons-R' 'Cons-M' 'Cons-C' 'Cons-A' 'C-A' 'P-A' 'K-R' 'K-M' 'K-C' 'K-A' 'S-R' 'S-M'  'S-C'  'S-A' 'Ker-A' 'MKer-A' 'MKerequ-A' 'MKer-fst-ALL' 'MKerW-A' 'MequKer-fst-A' 'SIM-R' 'SIM-M'  'SIM-C' 'SIM-A' 'SS-R' 'SS-M' 'SS-C' 'SS-A' 'iCluster-A'}
hind2=[1:16, 19, 21:24,29];
hindex=hindex(hind2);

set_jj=[0.5 1 2]; set_ii={'NMI', 'Purity','ARI'};

for ii=1:3
    for jj=1:3

if ii==1
meanmat_mut11=meanmat_mut1(hind2,jj)';   stdmat_mut11=stdmat_mut1(hind2,jj)';
elseif ii==2
meanmat_mut11=meanmat_mut2(hind2,jj)';   stdmat_mut11=stdmat_mut2(hind2,jj)';
elseif ii==3
meanmat_mut11=meanmat_mut3(hind2,jj)';   stdmat_mut11=stdmat_mut3(hind2,jj)';
end    
meanmat_mut111=[1:length(hind2); (meanmat_mut11);meanmat_mut11];
meanmat_mut1111=sortrows(meanmat_mut111',-2)
hind22=(meanmat_mut1111(:,1))
xtiti=hindex(hind22);

close all

clf
hold on
val1=meanmat_mut11(hind22); std_val1=stdmat_mut11(hind22);
ttt1=bar(1:22,val1);  
ttt2=errorbar(1:22,val1,std_val1,'.');  
set(gca,'XTick',1:1:22)
xlim([0 23]);  ylim([min(val1)-0.1  min(1,max(val1)+0.1)])     %ylim([0.3 1])
if ii==1
ylabel('NMI', 'FontSize', 15)
elseif ii==2
ylabel('Purity', 'FontSize', 15)
elseif ii==3
ylabel('ARI', 'FontSize', 15)
end
ttt=xlabel('Method'); ttt.FontSize=15;
alpha 0.5
xticklabels(xtiti)
xtickangle(45)
if (ii==1)*(jj==2)==1
htext=text(-0.03,1.05,'C','Units','normalized')
htext.FontSize=20;
end
if jj==1
    htitle=title('\sigma=0.5'); htitle.FontSize=14;
elseif (jj==2)*(ii>1)==1
    htitle=title('\sigma=1'); htitle.FontSize=14;
elseif jj==3
    htitle=title('\sigma=2'); htitle.FontSize=14;   
end
set(gca,'FontSize', 14);

print([sprintf( 'multi_cancer_22_pert_%d_%d',ii,jj)],'-depsc')
    end
end
