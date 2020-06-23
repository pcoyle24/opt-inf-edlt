% --------------------------------------------------------------------------
% File Name: accuracy_mex.m
% Author: Philip Coyle
% Date Created: 03/29/2018
% Date Updated: 04/02/2019
% Add Appropriate Paths %
% codetoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Codes/Code_Toolbox'; 
% mextoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Codes/mex_functions';
% addpath(genpath(codetoolbox));
% addpath(genpath(mextoolbox)); 
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_AR1/Accuracy
% matalb -nodesktop -nosplash -r accuracy_mex
% --------------------------------------------------------------------------
clear all
close all
clc

%% Housekeeping
addpath ../savedata101/
addpath ../common/

% Load in parameters
P = parameters2;

% Load in data
load(strcat('ChebPFs_PItarg',num2str(400*(P.pi_targ-1)),'_Ps',num2str(P.Ps),'_Pd',num2str(P.Pd),'.mat'))

O.e_pts = 50;

% Conditional Probability for Sun Spot Shock
TransMat = [P.Ps, 1-P.Ps;
            1-P.Pd, P.Pd]; 
        
% Load discrsetized state space
G = grids_cheb(O,P);

% Specify parmeters for simulation
sim = 101000;
burn = 1000;

% Generate Random demand and sunspot shocks
del = del_sim(sim + burn,P);
sun = sunspot_sim(sim + burn,P);

%% Compute Accuracy
% Allocate space for residuals
resid = nan(sim - burn,2);

for i = burn+1:sim+burn
    del_today = del(i);
    sun_today = sun(i);
    r_out = eqm_acc(del_today,sun_today,P,S,G,O,C,TransMat);
    
    % Store Residuals
    resid(i-burn,1) = r_out(1); % EE
    resid(i-burn,2) = r_out(2); % PC 
end

%% Output Accuracy Table
resid_10 = log10(abs(resid));

% Mean Error
ee_mu = round(mean(resid_10(:,1)),2);
pc_mu = round(mean(resid_10(:,2)),2);

% 95th Percentile Error
ee_95 = round(prctile(resid_10(:,1),95),2);
pc_95 = round(prctile(resid_10(:,2),95),2);



Mean_log_10_Residuals = [ee_mu;pc_mu];
Percentile_95_log_10_Residuals = [ee_95;pc_95];
PF = {'Comsumption Euler Equation'; 'Phillips Curve'};
T = table(PF,Mean_log_10_Residuals,Percentile_95_log_10_Residuals);

disp(T)

% Export Tale
writetable(T,'accuracy_50GH.txt');




