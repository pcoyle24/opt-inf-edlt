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
    addpath ..\Calib\
    addpath ..\min_inf_targ\
else
    addpath ../Calib/
    addpath ../min_inf_targ/
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

% Get Standard Regime Consumption Transfer (Demand Shock Only)
ct_demand;

max_ct_demand = max(CT_demand);
inx_demand = find(max_ct_demand == CT_demand);
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

% Get Deflationary Regime Consumption Transfer (Sunspot Shock Only)
ct_sun;

max_ct_sun = max(CT_sun);
inx_sun = find(max_ct_sun == CT_sun);
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

% Get Standard and Deflationary Regime Consumption Transfer (Demand Shock and Sunspot Shock)
ct_sun_demand;

max_ct_sun_demand = max(CT_demand_sun);
inx_sun_demand = find(max_ct_sun_demand == CT_demand_sun);
optinf = 400*(cPItarg(inx_sun_demand)-1);
disp(char(strcat('The Optimal Inflation Target is',{' '},num2str(optinf),'% with a demand shock and sunspot shock.')))

%% Plotting
colors = {'k','b','r'};
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
    
        plot(400*(cPItarg-1),100*CT_demand,colors{1},'LineWidth',2);
        set(gca,'XLim',[cPItarg_min_ann cPItarg_max_ann],'XTick',[cPItarg_min_ann, 0, 2, 4],'YLim',[-0.2 0],'FontSize',15)
        line([400*(cPItarg(inx(j))-1) 400*(cPItarg(inx(j))-1)],get(gca,'YLim'),'Color','b','LineStyle','-','LineWidth',1);        
    elseif j == 2
        cPItarg_min_ann = min_inf_targ_sun; 
        cPItarg_max_ann = 4;
        cPItarg_min = cPItarg_min_ann/400 + 1;
        cPItarg_max = cPItarg_max_ann/400 + 1;
        cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
        cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n);

        plot(400*(cPItarg-1),100*CT_sun,colors{1},'LineWidth',2);
        set(gca,'XLim',[cPItarg_min_ann cPItarg_max_ann],'XTick',[0, 2, 4],'YLim',[-2 0],'FontSize',15)
        line([400*(cPItarg(inx(j))-1) 400*(cPItarg(inx(j))-1)],get(gca,'YLim'),'Color','b','LineStyle','-','LineWidth',1);
    elseif j ==3
        cPItarg_min_ann = min_inf_targ_sun_demand; 
        cPItarg_max_ann = 4;
        cPItarg_min = cPItarg_min_ann/400 + 1;
        cPItarg_max = cPItarg_max_ann/400 + 1;
        cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
        cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n);
        
        plot(400*(cPItarg-1),100*CT_demand_sun,colors{1},'LineWidth',2);
        set(gca,'XLim',[cPItarg_min_ann cPItarg_max_ann],'XTick',[0, 2, 4],'YLim',[-2 0],'FontSize',15)                
        line([400*(cPItarg(inx(j))-1) 400*(cPItarg(inx(j))-1)],get(gca,'YLim'),'Color','b','LineStyle','-','LineWidth',1);        
    end

    xlabel('Inflation Target','FontSize',15)
    ylabel('Welfare','FontSize',15) 
    title(shock_string{j},'FontSize',15,'FontWeight','normal')
    
end

savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','ct_states.eps');
print(fig(1),'-depsc',strcat(savedir,'ct_states.eps'));

