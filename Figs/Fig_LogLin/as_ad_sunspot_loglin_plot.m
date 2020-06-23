%--------------------------------------------------------------------------
% File Name: as_ad_sunspot_loglin_plot.m
% Author: Philip Coyle
% Date Created: 01/04/2019
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/SS_Opt_Inf/SS_Opt_Inf/AS_AD_Plot/InfTarg
% as_ad_sunspot_loglin_plot
%--------------------------------------------------------------------------

clear all
close all
clc

%% Parameters
cPHIpi      = 2;
cBET        = 1/(1.0025);
cSIGMA      = 1;
cCHIc       = 1;
cCHIn       = 1;
cTHETA      = 11;
cKAPPA      = 0.02;
cRstar      = 0.0025;

p_t = 1;
p_d = 0.975;

%% Main Code

style = {'-','--',':'};
cPIstar_grid = [0, 0.0050];

y_d = zeros(length(cPIstar_grid),1);
pi_d = zeros(length(cPIstar_grid),1);

pi_d_out1 = zeros(101,length(cPIstar_grid));
pi_d_out2 = zeros(101,length(cPIstar_grid));

for i = 1:length(cPIstar_grid)
    cPIstar= cPIstar_grid(i);
    A_t = [1-p_t, -cSIGMA*(p_t - cPHIpi), -(1-p_t), -cSIGMA*(1-p_t);
         -cKAPPA, 1-cBET*p_t, 0, -cBET*(1-p_t);
         -(1-p_d), -cSIGMA*(1-p_d), 1-p_d, -cSIGMA*p_d;
          0, -cBET*(1-p_d), -cKAPPA, 1-cBET*p_d];
    b_t = [-cSIGMA*cPIstar*(1-cPHIpi);cPIstar*(1-cBET);cSIGMA*cRstar;cPIstar*(1-cBET)];

    x_t = A_t\b_t;

    y_t = x_t(1);
    pi_t = x_t(2);
    y_d(i) = x_t(3);
    pi_d(i) = x_t(4);

    disp({'Inflation Target','Output (T)','Inflation (T)','Output (D)','Inflation (D)'});
    disp({num2str(400*(cPIstar)), num2str(100*y_t), num2str(400*pi_t),num2str(100*y_d(i)), num2str(400*pi_d(i))});

    y_d_grid = linspace(-0.05,0,101);
    pi_d_out1(:,i) = (y_d_grid*cKAPPA  + cBET*(1-p_d)*pi_t + cPIstar*(1-cBET))/(1-cBET*p_d);
    pi_d_out2(:,i) = (y_d_grid*(1-p_d) -cSIGMA*(pi_t*(1-p_d) + cRstar))/(cSIGMA*p_d);
end

%% Plotting
fig(1) = figure;
subplot(2,2,1)
box on 
grid on
hold on
for i = 1:length(cPIstar_grid)
    h1(i) = plot(100*y_d_grid,400*pi_d_out1(:,i),'LineWidth',2,'Color','k','LineStyle',style{i});
    h2(i) = plot(100*y_d_grid,400*pi_d_out2(:,i),'LineWidth',2,'Color','b','LineStyle',style{i});
    
    h = plot(100*y_d(i),400*pi_d(i));
    set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',6)
end
% L = legend([h1(1) h2(1) h1(2) h2(2)],'Phillips Curve(\pi* = 0%)','Euler Equation (\pi* = 0%)','Phillips Curve (\pi* = 2%)','Euler Equation (\pi* = 2%)');
% set(L,'Location','SouthEast','FontSize',12)

xlabel('Output (y)','FontSize',15)
ylabel('Inflation (\pi)','FontSize',15)
title('Deflationary Regime','FontSize',15,'FontWeight','normal')
set(gca,'xlim',[-2 -0],'Ylim',[-2 0],'FontSize',15)


savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end
set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','as_ad_DR_cPIstar_lin.eps');
print(fig(1),'-depsc',strcat(savedir,'as_ad_DR_cPIstar_lin.eps'));

