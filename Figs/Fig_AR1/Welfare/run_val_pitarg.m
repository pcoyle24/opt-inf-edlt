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
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_AR1/Welfare
% matlab -nodesktop -nosplash -r run_val_pitarg
% --------------------------------------------------------------------------

clear all
close all
clc

run script_val_mex_pitarg_grid
run val_savedata_mex

exit