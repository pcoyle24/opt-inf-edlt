%--------------------------------------------------------------------------
%File Name: calib_cVARPHI_cSHOCK.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/Calib
% calib_cVARPHI_cSHOCK.m
%--------------------------------------------------------------------------

clear all 
close all
clc
% dbstop if error

%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg p_s p_d cDELz p_z p_c c_init pi_init
cBET        = 1/(1.0025);
cCHIc       = 1;
cCHIn       = 1;
cTHETA      = 11;
cTAU        = 1/cTHETA;
cVARPHI_guess = 1250;
cPHIpi      = 2;
cPHIy       = 0;
cRzlb       = 1;
cIOTA       = 1;
cALPHA      = 1; %0.5; 
cPItarg_ann = 0; 
cPItarg     = cPItarg_ann/400 + 1;
p_s         = 0.995;
p_d         = 0.975;
cDELz       = 1;
c_guess     = 0.9/100;
p_z         = 1;  
p_c         = 0.75; 

c_init      = 1 - 7/100;
pi_init     = -2/400 + 1;

%% Main Code
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

x0 = [cVARPHI_guess;c_guess];
func = @(x) param_solve_cVARPHI_c(x);

x_out_param = fsolve(func,x0,options);

cVARPHI   = x_out_param(1);
c  = x_out_param(2);

params = [cVARPHI, c];

save('params.mat','params')