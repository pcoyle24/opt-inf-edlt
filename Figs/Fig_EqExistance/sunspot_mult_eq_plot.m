%--------------------------------------------------------------------------
% File Name: sunspot_mult_eq_plot.m
% Author: Philip Coyle
% Date Created: 03/04/2019
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/2State_Markov_Shock/AS_AD_Plot/MultipleZLBEq
% sunspot_mult_eq_plot.m
%--------------------------------------------------------------------------

clear all 
close all
clc

load('mult_equib_exist_sun.mat');
load('mult_equib_exist_sun_bisection.mat')
p_s_grid         = (0.9:0.001:1);
p_d_grid         = (0.9:0.001:1);
[ss,dd] = meshgrid(p_s_grid,p_d_grid);

fig(1) = figure;
% subplot(2,2,1)
box on
grid on
hold on
for i = 1:length(p_s_grid)
    disp(i);
    for j = 1:length(p_d_grid)
        if Rs(i,j) > 1 % Check that Steady State Nominal Interest Rate is above ZLB
            if flag_save(i,j) > 0 % Check that nonlinear equation solver found a solution
                if multeq_save(i,j) == 2
                    h1 = scatter(100*ss(i,j),100*dd(i,j),20,'b','fill');
                elseif multeq_save(i,j) == 3
                    h2 = scatter(100*ss(i,j),100*dd(i,j),20,'r','fill');
                end
            end
        end
    end 
end
set(gca,'XLim',[90 100],'YLim',[90 100],'FontSize',25)
xlabel('100p_D','FontSize',25)
ylabel('100p_T','FontSize',25)
L = legend([h1 h2],'Unique Sunspot Equilibrium','Two Sunspot Equilibria');
set(L,'Location','SouthWest','FontSize',15)


savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','mult_equib_exist_sun.eps')
print(fig(1),'-depsc',strcat(savedir,'mult_equib_exist_sun.eps'));