%--------------------------------------------------------------------------
%File Name: calib_cALPHA.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/Calib
% calib_cALPHA.m
%--------------------------------------------------------------------------

% clear all 
% close all
% clc

load 'params_temp.mat'

%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA  cPItarg cDELz cDELc p_z p_c cPItarg_init_ann cPItarg_init
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
if cPItarg_init_ann == 0.5
    cALPHA_guess  = 0.5;
elseif cPItarg_init_ann == 1
    cALPHA_guess  = 0.75;
elseif cPItarg_init_ann == 1.5
    cALPHA_guess  = 0.8;
else 
    cALPHA_guess  = 0.9;
end
cPItarg_min_ann = cPItarg_init_ann - 0.1; 
cPItarg_max_ann = cPItarg_init_ann + 0.1;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
if ~isinteger(cPItarg_n)
    cPItarg     = linspace(cPItarg_min,cPItarg_max,round(cPItarg_n)); 
else
    cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n); 
end
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
% p_z         = 0.995; 
% p_c         = 0.75;    
cPItarg_init = cPItarg_init_ann/400 + 1;

%% Main Code
% options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'iter');

x0 = cALPHA_guess;
func = @(x) param_solve_cALPHA(x);
x_out_param = func(x0);

cALPHA   = x_out_param(1);

params = [cVARPHI, c, cALPHA];

name = {'pc_07','pc_08','pc_085','pz_099','pz_09925','pz_09975','pitarg_05','pitarg_1','pitarg_15'};

save(strcat('params_',name{k},'.mat'),'params')
