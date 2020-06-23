%--------------------------------------------------------------------------
%File Name: calib_cALPHA.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/SA/OptInf05/SA_Calib
% calib_cALPHA.m
%--------------------------------------------------------------------------

clear all 
close all
% clc

load 'params.mat'

%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA  cPItarg cDELz cDELc p_z p_c cPItarg_init
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
cALPHA_guess  = 0.5;
cPItarg_min_ann = 0.4; 
cPItarg_max_ann = 0.6;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
if ~isinteger(cPItarg_n)
    cPItarg     = linspace(cPItarg_min,cPItarg_max,ceil(cPItarg_n)); 
else
    cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n); 
end
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
p_z         = 0.995; 
p_c         = 0.75;    

cPItarg_init = 0.5/400 + 1;

%% Main Code
% options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'iter');

x0 = cALPHA_guess;
func = @(x) param_solve_cALPHA(x);
x_out_param = func(x0);

cALPHA   = x_out_param(1);

params = [cVARPHI, c, cALPHA];

save('params.mat','params')
