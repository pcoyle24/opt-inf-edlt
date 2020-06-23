% --------------------------------------------------------------------------
% File Name: valuefunc_cheb.m
% Author: Philip Coyle
% Date Created: 11/05/2018
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/AR1/1d/sunspot
% valuefunc
% --------------------------------------------------------------------------

clear all
close all
clc

addpath('/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/AR1/1d/calib')

% Load parameters
P = parameters1;

pi_targ = 400*(P.pi_targ -1);

%% Load in Data (Demand shock only)
load 'exp_val_opt_inf_51del_cheb.mat'

%% Load in Data (Demand shock and sunpot shock)
load 'exp_val_opt_inf_sunspot_101del.mat'


%% Plotting
fig(1) = figure(1);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig(1),'defaultAxesColorOrder',[left_color; right_color]);
box on
grid on
hold on
for j = 1:2
    if j == 1
        yyaxis left
%         ylabel('Demand Shock Only','FontSize',15)
        h(j) = plot(pi_targ,value_d,'k','LineWidth',2);
    else
        yyaxis right
%         ylabel('Demand & Sunspot Shock','FontSize',15)
        h(j) = plot(pi_targ,value_ds,'b','LineWidth',2);
    end
    
end

xlabel('\Pi^{targ}','FontSize',16)
set(gca,'XLim',[pi_targ(1)-0.001 pi_targ(end)],'YLim',[min(value_ds)-0.2 max(value_ds)+ 0.2],'FontSize',16)

line([opt_inf_d opt_inf_d],get(gca,'YLim'),'Color','k','LineStyle','--','LineWidth',1)
line([opt_inf_ds opt_inf_ds],get(gca,'YLim'),'Color','b','LineStyle','--','LineWidth',1)

L = legend([h(1) h(2)],'Demand Shock Only (Left Axis)','Demand & Sunspot Shock (Right Axis)','Location','SouthEast');



set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','Valuefunc_sunspot_del51_cheb.eps');


