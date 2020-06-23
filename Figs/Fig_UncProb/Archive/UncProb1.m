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
p_t         = linspace(0.9,1,101)';
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

unc_prob = [unc_prob_crisis, unc_prob_def];%, unc_prob_def./unc_prob_crisis];

%% Plotting
colors = {'k-.','k'};

fig(1) = figure(1);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig(1),'defaultAxesColorOrder',[left_color; right_color]);
subplot(2,2,1)
box on
grid on
hold on
for k = 1:size(unc_prob,2)
    yyaxis left
    set(gca,'XLim',[min(p_t) max(p_t)],'XTick',(0.9:0.02:1),'Ylim',[0 1],'YTick',(0:0.2:1),'FontSize',15)      
    
    h1(k) = plot(p_t,unc_prob(:,k),colors{k},'LineWidth',2);
    xlabel('p_T','FontSize',15)
    ylabel('Unconditional Probability','FontSize',15)
    line([0.975 0.975],get(gca,'YLim'),'Color','k','LineStyle','-','LineWidth',1);
%     title('Unconditional Probability','FontWeight','normal','FontSize',15)
end
yyaxis right
set(gca,'XLim',[min(p_t) max(p_t)],'XTick',(0.9:0.02:1),'Ylim',[0 50],'YTick',(0:25:50),'FontSize',15) 
ylabel('DR:CS','FontSize',15)

L1 = legend([h1(1) h1(2)],'Crisis State','Deflationary Regime');
set(L1,'Location','SouthWest','FontSize',12)

savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','uncprob1.eps');
% print(fig(1),'-depsc',strcat(savedir,'uncprob1.eps'));

