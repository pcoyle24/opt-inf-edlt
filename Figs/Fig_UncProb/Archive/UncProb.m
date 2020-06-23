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
p_s         = [0.9, 0.995, 0.999]; 
p_d         = linspace(0.95,1,101)'; 
p_z         = [0.9, 0.9925, 0.995]; 
p_c         = linspace(0.7,0.85,101)'; 

%% Unconditional Probability
unc_prob_crisis = zeros(length(p_c),length(p_z));
unc_prob_def = zeros(length(p_d),length(p_s));

for i = 1:length(p_s)
    unc_prob_crisis(:,i) = (1-p_z(i))./(2-p_z(i)-p_c);
    unc_prob_def(:,i) = (1-p_s(i))./(2-p_s(i)-p_d);
end

unc_prob = zeros([size(unc_prob_def),2]);
unc_prob(:,:,1) = unc_prob_crisis;
unc_prob(:,:,2) = unc_prob_def;

%% Plotting
colors_sun = {'b--','k','r-.'};
colors_demand = {'b--','r-.','k'};
titles = {{'Unconditional Probability', 'Crisis State'},{'Unconditional Probability','Deflationary Regime'}};

fig(1) = figure(1);
for i = 1:size(unc_prob,3)
    subplot(2,2,i)
    box on
    grid on
    hold on
    
    for k = 1:size(unc_prob,2)
        if i == 1
            h1(k) = plot(p_c,unc_prob(:,k,i),colors_demand{k},'LineWidth',2);
            set(gca,'XLim',[min(p_c) max(p_c)],'Ylim',[0 1],'YTick',(0:0.2:1),'FontSize',15)  
            xlabel('p_C','FontSize',15)
            line([0.75 0.75],get(gca,'YLim'),'Color','k','LineStyle','-','LineWidth',1);        
        else
            h2(k) = plot(p_d,unc_prob(:,k,i),colors_sun{k},'LineWidth',2);
            set(gca,'XLim',[min(p_d) max(p_d)],'Ylim',[0 1],'YTick',(0:0.2:1),'FontSize',15)
            xlabel('p_D','FontSize',15)
            line([0.975 0.975],get(gca,'YLim'),'Color','k','LineStyle','-','LineWidth',1);        
        end
%         ylabel('Unconditional Probability','FontSize',15) 
    end
    title(titles{i},'FontSize',15,'FontWeight','Normal') 
end

L1 = legend([h1(1) h1(2) h1(3)],'p_N = 0.99','p_N = 0.9925','p_N = 0.995');
set(L1,'Location','NorthEast','FontSize',12)
L2 = legend([h2(1) h2(2) h2(3)],'p_T = 0.99','p_T = 0.995','p_T = 0.999');
set(L2,'Location','NorthWest','FontSize',12)
        

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
% print(fig(1),'-depsc',strcat(savedir,'uncprob.eps'));
