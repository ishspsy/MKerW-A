%% This file geneartes Fig S4

addpath(genpath(pwd))

clear all
colors = {'b*-','r+-','kv-','gs-'};
load('simul_robust.mat')

colors = {'b*-','r+-','kv-','gs-'};
lam_set={'0.00001', '0.0001', '0.001', '0.01'};

clf
subplot(2,2,1)
tot_val=[mean(nmi_wd123_set);mean(pmi_wd123_set);mean(rmi_wd123_set)]; tot_val=tot_val(:,1:end-1);
for ii=1:3
pplo=plot(1:4,tot_val(ii,:),colors{ii});pplo.LineWidth=1.5; hold on
end
xlim([1 4])
set(gca,'XTick',1:4);
xticklabels(lam_set);  
hlegend=legend('NMI', 'Purity', 'ARI')
hlegend.FontSize=12; hlegend.Location='east';
htext=text(0.01,1.04, 'A','Units','normalized')
htext.FontSize=20;
xlabel('$\lambda$', 'Interpreter','latex');
%htitle=title(sprintf( '%s',dataset{i}))
set(gca, 'FontSize',15)

mu_set=[0.001 0.01 0.1 1];
subplot(2,2,2)
tot_val=[mean(nmi_wd123_setm);mean(pmi_wd123_setm);mean(rmi_wd123_setm)];  %tot_val=tot_val(:,1:end-1);
for ii=1:3
pplo=plot(1:4,tot_val(ii,:),colors{ii});pplo.LineWidth=1.5; hold on
end
xlim([1 4])
set(gca,'XTick',1:4);
xticklabels(mu_set);  
hlegend=legend('NMI', 'Purity', 'ARI')
hlegend.FontSize=12; hlegend.Location='east';
htext=text(0.01,1.04, 'B','Units','normalized')
htext.FontSize=20;
xlabel('$\mu$', 'Interpreter','latex');
set(gca, 'FontSize',15)

c_set=[0.001 0.01 0.1 1];
subplot(2,2,3)
tot_val=[mean(nmi_wd123_setm);mean(pmi_wd123_setm);mean(rmi_wd123_setm)];  %tot_val=tot_val(:,1:end-1);
for ii=1:3
pplo=plot(1:4,tot_val(ii,:),colors{ii});pplo.LineWidth=1.5; hold on
end
xlim([1 4])
set(gca,'XTick',1:4);
xticklabels(c_set);  
hlegend=legend('NMI', 'Purity', 'ARI')
hlegend.FontSize=12; hlegend.Location='east';
htext=text(0.01,1.04, 'C','Units','normalized')
htext.FontSize=20;
xlabel('$c$', 'Interpreter','latex');
set(gca, 'FontSize',15)

rho_set=[0.1 0.5 1 2 5];
subplot(2,2,4)
tot_val=[mean(nmi_wd123_setm);mean(pmi_wd123_setm);mean(rmi_wd123_setm)];  %tot_val=tot_val(:,1:end-1);
for ii=1:3
pplo=plot(1:4,tot_val(ii,:),colors{ii});pplo.LineWidth=1.5; hold on
end
xlim([1 4])
set(gca,'XTick',1:4);
xticklabels(rho_set);  
hlegend=legend('NMI', 'Purity', 'ARI')
hlegend.FontSize=12; hlegend.Location='east';
htext=text(0.01,1.04, 'D','Units','normalized')
htext.FontSize=20;
xlabel('$\rho$', 'Interpreter','latex');
set(gca, 'FontSize',15)

print -depsc cancer_robust_plot1

