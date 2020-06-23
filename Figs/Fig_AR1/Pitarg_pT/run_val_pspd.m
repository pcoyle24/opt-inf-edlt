% --------------------------------------------------------------------------
% File Name: run_val_pspd.m
% Author: Philip Coyle
% Date Created: 02/06/2018
% Last Updated: 11/29/2018
% 
% Add Appropriate Paths %
% codetoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/Code_Toolbox'; 
% mextoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/mex_functions';
% addpath(genpath(codetoolbox));
% addpath(genpath(mextoolbox)); 
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_AR1/Pitarg_pT
% matlab -nodesktop -nosplash -r run_val_pspd
% --------------------------------------------------------------------------

clear all
close all
clc

run script_val_mex_pitarg_pspd_grid
run val_savedata_mex_pspd

exit