%--------------------------------------------------------------------------
%File Name: welfare_states.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_5
% welfare_states.m
%--------------------------------------------------------------------------

clear all 
close all
clc
% dbstop if error

if ispc
    addpath ..\..\Calib\
    addpath ..\..\min_inf_targ\
else
    addpath ../../Calib/
    addpath ../../min_inf_targ/
end

load('params.mat')
%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c
cBET        = 1/(1.0025);
cCHIc       = 1;
cCHIn       = 1;
cTHETA      = 11;
cTAU        = 1/cTHETA;
cVARPHI     = params(1);
cPHIpi      = 2;
cPHIy       = 0;
cRzlb       = 1;
cIOTA       = 1;
cALPHA      = params(3);

p_s         = 0.995; 
p_d         = 0.975; 
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
p_z         = 0.995;
p_c         = 0.75; 

%% Get Effecient Steady State Values
ess;

%% Demand Shock Only
cPItarg_min_ann = -2; 
cPItarg_max_ann = 4;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n);

% Get Standard Regime Welfare (Demand Shock Only)
welfare_demand;

% Unconditional Probablity Calculation: Fundamentals Shock only.]
unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

unc_welfare_demand = unc_prob_zero*Vsz + unc_prob_crisis*Vsc;

max_welf_demand = max(unc_welfare_demand);
inx_demand = find(max_welf_demand == unc_welfare_demand);
optinf = 400*(cPItarg(inx_demand)-1);
disp(char(strcat('The Optimal Inflation Target is',{' '},num2str(optinf),'% with a demand shock only.')))

%% Sunspot Shock Only
load eqm_exist_sun.mat
min_inf_targ_sun = min_inf_targ;

cPItarg_min_ann = min_inf_targ_sun; 
cPItarg_max_ann = 4;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n);

% Get Deflationary Regime Welfare (Sunspot Shock Only)
welfare_sun;

% Unconditional Probablity Calculation
% Standard and Deflationary Regime (No Sunpot Shock)
unc_prob_targ = (1-p_d)/(2-p_s-p_d);
unc_prob_def = (1-p_s)/(2-p_s-p_d);
unc_welfare_sun = unc_prob_targ*Vs + unc_prob_def*Vd;

max_welf_sun = max(unc_welfare_sun);
inx_sun = find(max_welf_sun == unc_welfare_sun);
optinf = 400*(cPItarg(inx_sun)-1);
disp(char(strcat('The Optimal Inflation Target is',{' '},num2str(optinf),'% with a sunspot shock only.')))

%% Sunspot and Demand Shock
load eqm_exist_sun_demand.mat
min_inf_targ_sun_demand = min_inf_targ;

cPItarg_min_ann = min_inf_targ_sun_demand; 
cPItarg_max_ann = 4;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n);

% Get Standard and Deflationary Regime Welfare (Demand Shock and Sunspot Shock)
welfare_sun_demand;

% Unconditional Probablity Calculation
% Standard and Deflationary Regime (With Sunpot Shock)
StateMat = [p_s, 1-p_s;
            1-p_d, p_d];
        
ShockMat = [p_z, 1-p_z;
            1-p_c, p_c];
        
TransMat = kron(StateMat,ShockMat);

unc_prob = limitdist(TransMat);

unc_prob_sz = unc_prob(1);
unc_prob_sc = unc_prob(2);
unc_prob_dz = unc_prob(3);
unc_prob_dc = unc_prob(4);

unc_welfare_sun_demand = unc_prob_sz*Vsz + unc_prob_sc*Vsc + unc_prob_dz*Vdz + unc_prob_dc*Vdc;

max_welf_sun_demand = max(unc_welfare_sun_demand);
inx_sun_demand = find(max_welf_sun_demand == unc_welfare_sun_demand);
optinf = 400*(cPItarg(inx_sun_demand)-1);
disp(char(strcat('The Optimal Inflation Target is',{' '},num2str(optinf),'% with a demand shock and sunspot shock.')))

%% Plotting
% [row, col] = size(unc_welfare_demand);
% X = zeros(row,col,3);
% X(:,:,1) = unc_welfare_demand;
% X(:,:,2) = unc_welfare_sun;
% X(:,:,3) = unc_welfare_sun_demand;
colors = {'k','b','r'};
% shock_string = {'Demand Shock Only (Right Axis)','Sunspot Shock Only (Left Axis)','Demand & Sunspot Shock (Left Axis)'};
% shock_string = {'Demand Shock Only','Sunspot Shock Only',{'Demand & Sunspot','Shock'}};
shock_string = {{'Model with Crisis','Shock Only'},{'Model with Sunspot','Shock Only'},{'Model with Crisis &','Sunspot Shocks'}};

inx = [inx_demand inx_sun inx_sun_demand];

fig(1) = figure(1);
for j = 1:3
    subplot(2,3,j)
    box on
    grid on
    hold on
%     plot(400*(cPItarg-1),X(:,:,j),colors{1},'LineWidth',2);
    if j == 1
        cPItarg_min_ann = -2; 
        cPItarg_max_ann = 4;
        cPItarg_min = cPItarg_min_ann/400 + 1;
        cPItarg_max = cPItarg_max_ann/400 + 1;
        cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
        cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n);
    
%         plot(400*(cPItarg-1),unc_welfare_demand,colors{1},'LineWidth',2);
        plot(400*(cPItarg-1),100*CT_demand,colors{1},'LineWidth',2);
        set(gca,'XLim',[cPItarg_min_ann cPItarg_max_ann],'XTick',[cPItarg_min_ann, 0, 2, 4],'FontSize',15)                
        
        line([400*(cPItarg(inx(j))-1) 400*(cPItarg(inx(j))-1)],get(gca,'YLim'),'Color','b','LineStyle','-','LineWidth',1);
    elseif j == 2
        cPItarg_min_ann = min_inf_targ_sun; 
        cPItarg_max_ann = 4;
        cPItarg_min = cPItarg_min_ann/400 + 1;
        cPItarg_max = cPItarg_max_ann/400 + 1;
        cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
        cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n);

%         plot(400*(cPItarg-1),unc_welfare_sun,colors{1},'LineWidth',2);
        plot(400*(cPItarg-1),100*CT_sun,colors{1},'LineWidth',2);
%         set(gca,'XLim',[cPItarg_min_ann cPItarg_max_ann],'XTick',[cPItarg_min_ann, 0, 2, 4],'FontSize',15)
        set(gca,'XLim',[cPItarg_min_ann cPItarg_max_ann],'XTick',[0, 2, 4],'FontSize',15)
                
        line([400*(cPItarg(inx(j))-1) 400*(cPItarg(inx(j))-1)],get(gca,'YLim'),'Color','b','LineStyle','-','LineWidth',1);
    elseif j ==3
        cPItarg_min_ann = min_inf_targ_sun_demand; 
        cPItarg_max_ann = 4;
        cPItarg_min = cPItarg_min_ann/400 + 1;
        cPItarg_max = cPItarg_max_ann/400 + 1;
        cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
        cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n);
        
%         plot(400*(cPItarg-1),unc_welfare_sun_demand,colors{1},'LineWidth',2);
        plot(400*(cPItarg-1),100*CT_demand_sun,colors{1},'LineWidth',2);
%         set(gca,'XLim',[cPItarg_min_ann cPItarg_max_ann],'XTick',[cPItarg_min_ann, 0, 2, 4],'FontSize',15)
        set(gca,'XLim',[cPItarg_min_ann cPItarg_max_ann],'XTick',[0, 2, 4],'FontSize',15)
                
        line([400*(cPItarg(inx(j))-1) 400*(cPItarg(inx(j))-1)],get(gca,'YLim'),'Color','b','LineStyle','-','LineWidth',1);
    end
%     set(gca,'XLim',[-2 4],'FontSize',15)
%     line([400*(cPItarg(inx(j))-1) 400*(cPItarg(inx(j))-1)],get(gca,'YLim'),'Color','k','LineStyle','-','LineWidth',1);
    xlabel('Inflation Target','FontSize',15)
    ylabel('Consumption Transfer (%)','FontSize',15) 
    title(shock_string{j},'FontSize',15,'FontWeight','normal')
    
end

% savedir = cd;
% savedir = fullfile(savedir, '..');
% if ispc 
%     savedir = strcat(savedir,'\Final\');
% else
%     savedir = strcat(savedir,'/Final/');
% end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','ct_states.eps');
% print(fig(1),'-depsc',strcat(savedir,'welfare_states.eps'));

