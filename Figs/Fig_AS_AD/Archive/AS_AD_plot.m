%--------------------------------------------------------------------------
%File Name: AS_AD_plot.m
%Author: Philip Coyle
%Date Created: 10/17/2018
% /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/AS_AD_plot
% AS_AD_plot
%--------------------------------------------------------------------------

clear all 
close all
clc

load('AS_AD_demand_crisis.mat')

% Crisis State
fig(1) = figure;
subplot(2,2,1)
box on
grid on
hold on

plot(100*(as_c_demand1-1),400*(pi_grid1-1),'LineWidth',2,'Color','k','LineStyle','-');
plot(100*(ad_c_demand1-1),400*(pi_grid1-1),'LineWidth',2,'Color','b','LineStyle','-'); 
h = plot(100*(Sc(1,1)-1),400*(Sc(1,2)-1));
set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',10)

plot(100*(as_c_demand2-1),400*(pi_grid2-1),'LineWidth',2,'Color','k','LineStyle','--');
plot(100*(ad_c_demand2-1),400*(pi_grid2-1),'LineWidth',2,'Color','b','LineStyle','--');  
h = plot(100*(Sc(2,1)-1),400*(Sc(2,2)-1));
set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',10)

set(gca,'XLim',[-7.5 -5],'Ylim',[-2.5 0.75],'FontSize',15)
xlabel('Consumption (C)','FontSize',15)
ylabel('Inflation (\Pi)','FontSize',15)
title('Crisis State (Demand Shock Only)','FontSize',15,'FontWeight','normal')

% Deflationary Regime
load('AS_AD_sun_deflationary.mat')
subplot(2,2,2)
box on
grid on
hold on

h1(1) = plot(100*(as_c_sun1-1),400*(pi_grid1-1),'LineWidth',2,'Color','k','LineStyle','-');
h2(1) = plot(100*(ad_c_sun1-1),400*(pi_grid1-1),'LineWidth',2,'Color','b','LineStyle','-'); 
h = plot(100*(Sd(1,1)-1),400*(Sd(1,2)-1));
set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',10)

h1(2) = plot(100*(as_c_sun2-1),400*(pi_grid2-1),'LineWidth',2,'Color','k','LineStyle','--');
h2(2) = plot(100*(ad_c_sun2-1),400*(pi_grid2-1),'LineWidth',2,'Color','b','LineStyle','--');  
h = plot(100*(Sd(2,1)-1),400*(Sd(2,2)-1));
set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',10)

set(gca,'XLim',[-3 -0.5],'Ylim',[-1.5 -1],'FontSize',15)
xlabel('Consumption (C)','FontSize',15)
ylabel('Inflation (\Pi)','FontSize',15)
title('Deflationary Regime (Sunspot Shock Only)','FontSize',15,'FontWeight','normal')

L = legend([h1(1) h2(1) h1(2) h2(2)],'AS Curve \Pi^{targ} = 0%','AD Curve \Pi^{targ} = 0%','AS Curve \Pi^{targ} = 2%','AD Curve \Pi^{targ} = 2%');
set(L,'Position',[0.51 0.45 0.01 0.01],'Orientation','horizontal','FontSize',12)
legend('boxoff');
% set(L,'Location','Best','FontSize',15)

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','AS_AD_plot_cALPHAstar.eps');
