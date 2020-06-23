%--------------------------------------------------------------------------
% File Name: as_ad_demand_loglin_plot.m
% Author: Philip Coyle
% Date Created: 01/04/2019
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/SS_Opt_Inf/SS_Opt_Inf/AS_AD_Plot/InfTarg
% as_ad_demand_loglin_plot
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
cRz         = 0;
cRc         = 0.15833/10;

p_z = 1;
p_c = 0.75;

%% Main Code

style = {'-','--',':'};

cPIstar_grid = [0, 0.0050];

y_c = zeros(2,1);
pi_c = zeros(2,1);

pi_c_out1 = zeros(101,length(cPIstar_grid));
pi_c_out2 = zeros(101,length(cPIstar_grid));

for i = 1:length(cPIstar_grid)
    cPIstar= cPIstar_grid(i);    
    A_z = [1-p_z, -cSIGMA*(p_z - cPHIpi), -(1-p_z), -cSIGMA*(1-p_z);
         -cKAPPA, 1-cBET*p_z, 0, -cBET*(1-p_z);
         -(1-p_c), -cSIGMA*(1-p_c), 1-p_c, -cSIGMA*p_c;
          0, -cBET*(1-p_c), -cKAPPA, 1-cBET*p_c];
    b_z = [-cRz-cPIstar*cSIGMA*(1-cPHIpi);cPIstar*(1-cBET); cSIGMA*(cRstar-cRc);cPIstar*(1-cBET)];
    x_z = A_z\b_z;

    y_z = x_z(1);
    pi_z = x_z(2);
    y_c(i) = x_z(3);
    pi_c(i) = x_z(4);

    disp({'Inflation Target','Output (N)','Inflation (N)','Output (C)','Inflation (C)'});
    disp({num2str(400*(cPIstar)), num2str(100*y_z), num2str(400*pi_z),num2str(100*y_c(i)), num2str(400*pi_c(i))});
    
    y_c_grid = linspace(-0.1,0.02,101);
    pi_c_out1(:,i) = (y_c_grid*cKAPPA + cBET*(1-p_c)*pi_z + cPIstar*(1-cBET))/(1-cBET*p_c);
    pi_c_out2(:,i) = (y_c_grid*(1-p_c)-cSIGMA*(pi_z*(1-p_c) + cRstar - cRc))/(cSIGMA*p_c);
end

%% Plotting
fig(1) = figure;
subplot(2,2,1)
box on 
grid on
hold on
for i = 1:length(cPIstar_grid)
    h1(i) = plot(100*y_c_grid,400*pi_c_out1(:,i),'LineWidth',2,'Color','k','LineStyle',style{i});
    h2(i) = plot(100*y_c_grid,400*pi_c_out2(:,i),'LineWidth',2,'Color','b','LineStyle',style{i});
    
    h = plot(100*y_c(i),400*pi_c(i));
    set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',6)
end
L = legend([h1(1) h2(1) h1(2) h2(2)],'AS Curve(\pi* = 0%)','AD Curve (\pi* = 0%)','AS Curve (\pi* = 2%)','AD Curve (\pi* = 2%)');
% L = legend([h1(1) h2(1) h1(2) h2(2)],'Phillips Curve(\pi* = 0%)','Euler Equation (\pi* = 0%)','Phillips Curve (\pi* = 2%)','Euler Equation (\pi* = 2%)');
set(L,'Location','SouthEast','FontSize',10)

xlabel('Output (y)','FontSize',15)
ylabel('Inflation (\pi)','FontSize',15)
title('Crisis State','FontSize',15,'FontWeight','normal')
set(gca,'Xlim',[-8 -4],'Ylim',[-4 1],'FontSize',15)

savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end
set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','as_ad_CRISIS_cPIstar_lin.eps');
print(fig(1),'-depsc',strcat(savedir,'as_ad_CRISIS_cPIstar_lin.eps'));