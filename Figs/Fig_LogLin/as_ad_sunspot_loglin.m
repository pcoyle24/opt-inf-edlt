%--------------------------------------------------------------------------
% File Name: as_ad_sunspot_loglin.m
% Author: Philip Coyle
% Date Created: 01/04/2019
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/SS_Opt_Inf/SS_Opt_Inf/AS_AD_Plot/InfTarg
% as_ad_sunspot_loglin
%--------------------------------------------------------------------------

clear all
close all
clc

syms p_t p_d cPHIpi cBET cKAPPA cRstar cPIstar cSIGMA

%% Simplified case when p_t is an absorbing state
p_t = 1;
A_t = [1-p_t, -cSIGMA*(p_t - cPHIpi), -(1-p_t), -cSIGMA*(1-p_t);
     -cKAPPA, 1-cBET*p_t, 0, -cBET*(1-p_t);
     -(1-p_d), -cSIGMA*(1-p_d), 1-p_d, -cSIGMA*p_d;
      0, -cBET*(1-p_d), -cKAPPA, 1-cBET*p_d];
b_t = [-cSIGMA*cPIstar*(1-cPHIpi);cPIstar*(1-cBET);cSIGMA*cRstar;cPIstar*(1-cBET)];
out = simplify(A_t\b_t,'steps',500);

disp(out);


