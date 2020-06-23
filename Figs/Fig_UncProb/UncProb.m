%--------------------------------------------------------------------------
% File Name: UncProb.m
% Author: Philip Coyle
% Date Created: 01/17/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/UncProb
% UncProb
%--------------------------------------------------------------------------

clear all 
close all
clc

%% Housekeeping
p_t         = linspace(0.99,1,101)';
p_d         = 0.975;
p_z         = 0.995;
p_c         = 0.75;

%% Unconditional Probability
unc_prob_crisis = zeros(length(p_t),1);
unc_prob_def = zeros(length(p_t),1);

for i = 1:length(p_t)
    unc_prob_crisis(i) = (1-p_z)/(2-p_z-p_c);
    unc_prob_def(i) = (1-p_t(i))./(2-p_t(i)-p_d);
end

unc_prob = [unc_prob_crisis, unc_prob_def, unc_prob_def./unc_prob_crisis];

%% Plotting
colors = {'k-.','k','k:'};

scale = max(unc_prob(:,2))/0.4;

fig(1) = figure(1);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig(1),'defaultAxesColorOrder',[left_color; right_color]);
subplot(2,2,1)
box on
grid on
hold on
for k = 1:size(unc_prob,2)
    if k == 3
        yyaxis right
        set(gca,'XLim',[min(p_t) max(p_t)],'XTick',(0.99:0.002:1),'Ylim',[0 max(unc_prob(:,3))/scale],'FontSize',12)
        ylabel({'[Relative to Uncond. Prob.', 'of Crisis State]'},'FontSize',12)
    else
        yyaxis left
        set(gca,'XLim',[min(p_t) max(p_t)],'XTick',(0.99:0.002:1),'Ylim',[0 max(unc_prob(:,2))/scale],'XTick',(0:0.02:0.4),'FontSize',12)          
        ylabel('Unconditional Probability','FontSize',12)
    end
    h1(k) = plot(p_t,unc_prob(:,k),colors{k},'LineWidth',2);
    xlabel('p_T (Persistence of Target Regime)','FontSize',12)
    line([0.995 0.995],get(gca,'YLim'),'Color','k','LineStyle','-','LineWidth',1);
    
end

L1 = legend([h1(1) h1(2)],'Crisis State','Deflationary Regime');
set(L1,'Location','NorthEast','FontSize',12)

savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','uncprob.eps');
print(fig(1),'-depsc',strcat(savedir,'uncprob.eps'));

