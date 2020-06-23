%--------------------------------------------------------------------------
%File Name: calib.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/Calib
% calib.m
%--------------------------------------------------------------------------

clear all
close all
clc

warning off 

global p_z p_c cPItarg_init_ann

cases = [0.70, 0.995, 2;
% 0.75, 0.995, 2;
0.80, 0.995, 2;
0.85, 0.995, 2;
0.75, 0.99, 2;
0.75, 0.9925, 2;
% 0.75, 0.995, 2;
0.75, 0.9975, 2;
0.75, 0.995, 0.5;
0.75, 0.995, 1;
0.75, 0.995, 1.5];
% 0.75, 0.995, 2];

for k = 1:length(cases)
    msg  = sprintf('p_c = %.2f, p_z = %.4f,  pi_targ = %.1f.',cases(k,1),cases(k,2),cases(k,3));
    disp(msg);
    p_c = cases(k,1);
   
    calib_cVARPHI_cSHOCK;
    disp('Found correct cVARPHI and cSHOCK values')
    
    p_z = cases(k,2);
    cPItarg_init_ann = cases(k,3);
    calib_cALPHA;
    disp('Found correct cALPHA value')
    
    find_min_inftarg_sun_demand
    disp(' ')
end


