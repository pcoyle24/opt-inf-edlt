%--------------------------------------------------------------------------
% File Name: as_ad_demand_loglin.m
% Author: Philip Coyle
% Date Created: 01/04/2019
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/SS_Opt_Inf/SS_Opt_Inf/AS_AD_Plot/InfTarg
% as_ad_demand_loglin
%--------------------------------------------------------------------------

clear all
close all
clc

syms p_z p_c cPHIpi cBET cKAPPA cRstar cRz cRc cPIstar  cSIGMA

%% Let p_z be an absorbing state
p_z = 1;
cRz = 0;
A_z = [1-p_z, -cSIGMA*(p_z - cPHIpi), -(1-p_z), -cSIGMA*(1-p_z);
     -cKAPPA, 1-cBET*p_z, 0, -cBET*(1-p_z);
     -(1-p_c), -cSIGMA*(1-p_c), 1-p_c, -cSIGMA*p_c;
      0, -cBET*(1-p_c), -cKAPPA, 1-cBET*p_c];
b_z = [-cRz-cPIstar*cSIGMA*(1-cPHIpi);cPIstar*(1-cBET); cSIGMA*(cRstar-cRc);cPIstar*(1-cBET)];

out = simplify(A_z\b_z,'steps',500);

disp(out);





