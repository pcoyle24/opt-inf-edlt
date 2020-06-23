%--------------------------------------------------------------------------
%File Name: welfare_states.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd % /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/SS_Opt_Inf/Model_2\Ex_1'
% welfare_states.m
%--------------------------------------------------------------------------

clear all 
close all
clc
% dbstop if error

%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c
cBET        = 1/(1.0025);
cCHIc       = 1;
cCHIn       = 1;
cTHETA      = 11;
cTAU        = 1/cTHETA;
cVARPHI     = 2550;
cPHIpi      = 2;
cPHIy       = 0;
cRzlb       = 1;
cIOTA       = 1;
cALPHA      = 0.96;
cPItarg_min_ann = -2; 
cPItarg_max_ann = 4;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.05 + 1;
cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n); 
p_s         = 0.995; %(0.85:0.01:1);
p_d         = 0.98; %(0.85:0.01:1);
cDELz       = 1;
c           = 0.20657/10;%0.198/10;
cDELc       = 1 + c;
p_z         = 0.995;  % 0.99  
p_c         = 0.75;   % 0.75, 0.8, 5/6,  

%% Main Code
% Get Standard Regime Welfare (Demand Shock Only)
welfare_sss;

% Get Deflationary Regime Welfare (Demand Shock Only)
welfare_dss;

% Unconditional Probablity Calculation
% Standard and Deflationary Regime (No Sunpot Shock)
unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

unc_welfare_sss = unc_prob_zero*Vsz + unc_prob_crisis*Vsc;
unc_welfare_dss = unc_prob_zero*Vdz + unc_prob_crisis*Vdc;


% Get Standard and Deflationary Regime Welfare (Demand Shock and Sunspot Shock)
welfare_sun;

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

unc_welfare_sun = unc_prob_sz*Vsz + unc_prob_sc*Vsc + unc_prob_dz*Vdz + unc_prob_dc*Vdc;

%% Plotting
[row, col] = size(unc_welfare_sss);
X = zeros(row,col,3);
X(:,:,1) = unc_welfare_sss;
X(:,:,2) = unc_welfare_dss;
X(:,:,3) = unc_welfare_sun;
colors = {'k','b','r'};
shock_string = {'Target Regime (Demand Shock Only)','Deflationary Regime (Demand Shock Only)','Target & DeflationaryRegime (Demand & Sunspot Shock)'};

fig(1) = figure(1);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig(1),'defaultAxesColorOrder',[left_color; right_color]);
box on
grid on
hold on
for j = 1:3
    if j == 1
        yyaxis left
        ylabel('Target Regime','FontSize',15)
    else
        yyaxis right
        ylabel('Deflationary Regime & Sunspot Shock','FontSize',15)
    end
    h(j) = plot(400*(cPItarg-1),X(:,:,j),colors{j},'LineWidth',2);
    set(gca,'XLim',[-2 4],'FontSize',15)
    xlabel('Inflation Target','FontSize',15)
    
%     ylabel('Unconditional Welfare (TR)','FontSize',15)
end
L = legend([h(1) h(2) h(3)],shock_string{1},shock_string{2},shock_string{3},'Location', 'South');

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','welfare_states.eps');

