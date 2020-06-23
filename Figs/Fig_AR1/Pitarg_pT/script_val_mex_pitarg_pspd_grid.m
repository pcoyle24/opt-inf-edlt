% --------------------------------------------------------------------------
% File Name: script_val_mex_pitarg_pspd_grid.m
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
% matlab -nodesktop -nosplash -r script_val_mex_pitarg_pspd_grid
% --------------------------------------------------------------------------

% Projection Method (Chebyshev basis):
% Canonical New Keynesian Model (Rotemberg Pricing) with AR1 Preference
% Shock and Sunspot shock. This script utilizes mex functions for speed 
% pickups and must be run in Linux.


clear all;
close all;
clc;

% Directory to save data
addpath ../common/
addpath ../mex_functions/
savedir = pwd;
savedir = strcat(pwd,'/savedata101/');

%% --------------------------------------------------------------------------
%  Initialize Policy Functions
%  --------------------------------------------------------------------------
% Load parameters, steady state and grids
% P = parameters2; % 10 bps grid interval
P = parameters3; % 5 bps grid interval

% Specify grid options
O.delbound = [1-4*P.bound 1+4*P.bound];
O.del_pts = 101;

%% ------------------------------------------------------------------------
%  Solving the Model
%  ------------------------------------------------------------------------
global pi_yesterday
pi_yesterday       = 1;
for d = 1:length(P.Pd_grid)
    P.Pd = P.Pd_grid(d);
    for s = 1:length(P.Ps_grid)
        P.Ps = P.Ps_grid(s);

        % Specify grid options (sunspot shock)
        O.sunbound = [P.Pd, P.Ps];
        O.sun_pts = 2;
        O.e_pts = 10;

        % Chebyshev polynomial order in each dimension (must <= pts in grid)
        O.n1 = 4;

        % Load discrsetized state space
        G = grids_cheb(O,P);

        % Initialize Chebyshev polynomial parameters
        C = chebpoly(G,O);

        % Conditional Probability for Sun Spot Shock
        TransMat = [P.Ps, 1-P.Ps;
                    1-P.Pd, P.Pd]; 

        for i = 1:length(P.pi_targ)
            S = steadystate(P,i);

            % Retrieve initial policy functions
            pf = guess_nzlb_val_mex(P,S,G,O,C); 

            % ----------------------------------------------------------------------
            % Initialize algorithm (Non-Linear Least Squares Minimization)
            % ----------------------------------------------------------------------
            % Calculate initial Chebyshev coefficients (Standard Regime, Non-ZLB)          
            C.Ac_s = Fchebweights11(O.n1,O.del_pts,pf.c_s,C.T,C.P,C.X);
            C.Ainf_s = Fchebweights11(O.n1,O.del_pts,pf.inf_s,C.T,C.P,C.X);

            % Calculate initial Chebyshev coefficients (Stanard Regime, ZLB)          
            C.Ac_zlb_s = Fchebweights11(O.n1,O.del_pts,pf.c_zlb_s,C.T,C.P,C.X);
            C.Ainf_zlb_s = Fchebweights11(O.n1,O.del_pts,pf.inf_zlb_s,C.T,C.P,C.X);

            % Calculate initial Chebyshev coefficients (Deflationary Regime, Non-ZLB)          
            C.Ac_d = Fchebweights11(O.n1,O.del_pts,pf.c_d,C.T,C.P,C.X);
            C.Ainf_d = Fchebweights11(O.n1,O.del_pts,pf.inf_d,C.T,C.P,C.X);

            % Calculate initial Chebyshev coefficients (Deflationary Regime, ZLB)          
            C.Ac_zlb_d = Fchebweights11(O.n1,O.del_pts,pf.c_zlb_d,C.T,C.P,C.X);
            C.Ainf_zlb_d = Fchebweights11(O.n1,O.del_pts,pf.inf_zlb_d,C.T,C.P,C.X);

            theta_old = [C.Ac_s,C.Ainf_s,C.Ac_zlb_s,C.Ainf_zlb_s,C.Ac_d,C.Ainf_d,C.Ac_zlb_d,C.Ainf_zlb_d];
            options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','display','iter','MaxFunEvals',300,'MaxIter',500,'TolFun',1e-12,'TolX',1e-10,'FunValCheck','On','SpecifyObjectiveGradient',true);
            func    = @(weights) eqm_mex(weights,P,S,G,O,C,TransMat);
            disp(char('Calculating the Policy Functions for Economy with Sunpot Shock'))

            [theta_new,norm,residual,exitflag,output,lambda,jacob] = lsqnonlin(func,theta_old,[],[],options);

            % Weights (Standard Regime)
            C.Ac_s        = theta_new(1:O.n1+1,1);
            C.Ainf_s      = theta_new(1:O.n1+1,2);
            C.Ac_zlb_s    = theta_new(1:O.n1+1,3);
            C.Ainf_zlb_s  = theta_new(1:O.n1+1,4);
            % Weights (Deflationary Regime)
            C.Ac_d        = theta_new(1:O.n1+1,5);
            C.Ainf_d      = theta_new(1:O.n1+1,6);
            C.Ac_zlb_d    = theta_new(1:O.n1+1,7);
            C.Ainf_zlb_d  = theta_new(1:O.n1+1,8);

            % Build out Policy functions
            % Standard Regime, NZLB
            pf.c_s = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C.Ac_s,C.max,C.T,C.P);    
            pf.inf_s = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C.Ainf_s,C.max,C.T,C.P); 
            pf.pitilde_s  = pf.inf_s/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
            pf.n_s        = (pf.c_s./(1-(P.varphi/2).*(pf.pitilde_s-1).^2));
            pf.y_s        = pf.n_s;
            pf.w_s        = pf.n_s.^P.chin.*pf.c_s.^P.chic;
            pf.r_s        = S.pi_targ./(P.beta*G.del_grid').*((pf.inf_s./S.pi_targ).^(P.phi_pi).*(pf.y_s./S.y).^(P.phi_y));

            % Standard Regime, ZLB
            pf.c_zlb_s = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C.Ac_zlb_s,C.max,C.T,C.P);    
            pf.inf_zlb_s = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C.Ainf_zlb_s,C.max,C.T,C.P); 
            pf.pitilde_zlb_s  = pf.inf_zlb_s/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
            pf.n_zlb_s        = (pf.c_zlb_s./(1-(P.varphi/2).*(pf.pitilde_zlb_s-1).^2));
            pf.y_zlb_s        = pf.n_zlb_s;
            pf.w_zlb_s        = pf.n_zlb_s.^P.chin.*pf.c_zlb_s.^P.chic;
            pf.r_zlb_s = ones(length(pf.r_s),1);

            % Deflationary Regime, NZLB
            pf.c_d = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C.Ac_d,C.max,C.T,C.P);    
            pf.inf_d = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C.Ainf_d,C.max,C.T,C.P); 
            pf.pitilde_d  = pf.inf_d/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
            pf.n_d        = (pf.c_d./(1-(P.varphi/2).*(pf.pitilde_d-1).^2));
            pf.y_d        = pf.n_d;
            pf.w_d        = pf.n_d.^P.chin.*pf.c_d.^P.chic;
            pf.r_d        = S.pi_targ./(P.beta*G.del_grid').*((pf.inf_d./S.pi_targ).^(P.phi_pi).*(pf.y_d./S.y).^(P.phi_y)); 

            % Deflationary Regime, ZLB
            pf.c_zlb_d = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C.Ac_zlb_d,C.max,C.T,C.P);    
            pf.inf_zlb_d = Fallcheb111(O.delbound,O.del_pts,G.del,O.n1,C.Ainf_zlb_d,C.max,C.T,C.P); 
            pf.pitilde_zlb_d  = pf.inf_zlb_d/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
            pf.n_zlb_d        = (pf.c_zlb_d./(1-(P.varphi/2).*(pf.pitilde_zlb_d-1).^2));
            pf.y_zlb_d        = pf.n_zlb_d;
            pf.w_zlb_d        = pf.n_zlb_d.^P.chin.*pf.c_zlb_d.^P.chic;
            pf.r_zlb_d = ones(length(pf.r_d),1);

            C.Ar_s = Fchebweights11(O.n1,O.del_pts,pf.r_s,C.T,C.P,C.X);
            C.Ar_d = Fchebweights11(O.n1,O.del_pts,pf.r_d,C.T,C.P,C.X);

            disp(char('Calculating Value Function'))
            %----------------------------------------------------------------------
            % Initialize algorithm (Value Function Iteration) 
            %----------------------------------------------------------------------
            % Calculate initial Chebyshev coefficients (Non-ZLB)          
            C.Av_s = Fchebweights11(O.n1,O.del_pts,pf.v_s,C.T,C.P,C.X);         
            C.Av_zlb_s = Fchebweights11(O.n1,O.del_pts,pf.v_zlb_s,C.T,C.P,C.X);
            C.Av_d = Fchebweights11(O.n1,O.del_pts,pf.v_d,C.T,C.P,C.X);         
            C.Av_zlb_d = Fchebweights11(O.n1,O.del_pts,pf.v_zlb_d,C.T,C.P,C.X);

            theta_old_val = [C.Av_s,C.Av_zlb_s,C.Av_d,C.Av_zlb_d];

            % Preallocate arrays to store policy function updates    
            func_val    = @(weights) eqm_val(weights,P,S,G,O,C,pf,TransMat);    
            options_val = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','display','iter','MaxFunEvals',300,'MaxIter',500,'TolFun',1e-10,'TolX',1e-10,'FunValCheck','On','SpecifyObjectiveGradient',false);

            [theta_new,norm,residual,exitflag,output,lambda,jacob] = lsqnonlin(func_val,theta_old_val,[],[],options_val);

            % Map Coefficients 
            C.Av_s     = theta_new(1:O.n1+1,1);
            C.Av_zlb_s = theta_new(1:O.n1+1,2);
            C.Av_d     = theta_new(1:O.n1+1,1);
            C.Av_zlb_d = theta_new(1:O.n1+1,2);

            pf.v_s      =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Av_s,C.max,C.T,C.P);    
            pf.v_zlb_s  =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Av_zlb_s,C.max,C.T,C.P); 
            pf.v_d      =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Av_d,C.max,C.T,C.P);    
            pf.v_zlb_d  =  Fallcheb111(O.delbound,O.del_pts,G.del_grid,O.n1,C.Av_zlb_d,C.max,C.T,C.P); 

            disp(char('Policy Functions Calculated'))

            % Save coefficients
            save(strcat(savedir,'ChebPFs_PItarg',num2str(400*(P.pi_targ(i)-1)),'_Ps',num2str(P.Ps),'_Pd',num2str(P.Pd),'.mat'),'C','O','pf','S');    
        end
    end
end