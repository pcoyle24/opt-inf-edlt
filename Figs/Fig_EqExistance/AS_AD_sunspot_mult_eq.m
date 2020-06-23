%--------------------------------------------------------------------------
% File Name: AS_AD_sunspot_mult_eq.m
% Author: Philip Coyle
% Date Created: 03/04/2019
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/2State_Markov_Shock/AS_AD_Plot/MultipleZLBEq
% AS_AD_sunspot_mult_eq.m
%--------------------------------------------------------------------------

clear all 
close all
clc

if ispc
    addpath ..\Calib\
else
    addpath ../Calib/
end
load params.mat

%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c
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
cALPHA      = params(3);
cPItarg_ann = 2; 
cPItarg     = cPItarg_ann/400 + 1;
p_s_grid         = (0.9:0.001:1);
p_d_grid         = (0.9:0.001:1);
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
c_check   = 0;

multeq_save   = nan(length(p_s_grid),length(p_d_grid));
it = 0;

%% Main Code
for s = 1:length(p_s_grid)
    p_s = p_s_grid(s);
    for d = 1:length(p_d_grid)
        p_d = p_d_grid(d);

        % Get Steady State Values
        [S,flag] = get_ss_sun;
        if flag > 0
            options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

            pi_grid = linspace(-15/400 + 1, 3/400 + 1,501);
            pi_grid = sort([pi_grid,cPItarg]);

            y_guess = ones(length(pi_grid),1)*S.Yd;
            as = zeros(length(pi_grid),1);
            ad = zeros(length(pi_grid),1);

            for i = 1:length(pi_grid)
                state = pi_grid(i);
                x_guess = y_guess(i);

                func_as = @(x_guess) get_as_sun(x_guess,state,S,c_check);
                as_out = fsolve(func_as,x_guess,options);
                func_ad = @(x_guess) get_ad_sun(x_guess,state,S,c_check);
                ad_out = fsolve(func_ad,x_guess,options);

                as(i) = as_out;
                ad(i) = ad_out;
            end

            sort_asad = sort(abs(as-ad));

%% Identify Equilbrium Points
            inx1 = find(abs(as-ad) == sort_asad(1));
            inx2 = find(abs(as-ad) == sort_asad(2));

            it = 2;
            while abs(as(inx1) - as(inx2)) < 1e-3 % Making sure there is no overlap
                inx2 = find(abs(as-ad) == sort_asad(it));
                it = it + 1;
            end

            inx3 = find(abs(as-ad) == sort_asad(3));
            it = 3;
            while abs(as(inx1) - as(inx3)) < 1e-1 || abs(as(inx2) - as(inx3)) < 1e-1 % Making sure there is no overlap
                inx3 = find(abs(as-ad) == sort_asad(it));

                it = it + 1;
            end
            inx = [inx1, inx2, inx3];

%% Record Multiple Equilbrium Existence
            if sum(abs(as(inx) - ad(inx)) < 1e-2) == length(inx)
                multeq_save(s,d) = length(inx);
            elseif sum(abs(as(inx) - ad(inx)) < 1e-2) == length(inx)-1
                multeq_save(s,d) = length(inx)-1;
            elseif sum(abs(as(inx) - ad(inx)) < 1e-2) == length(inx)-2
                multeq_save(s,d) = 0;
            end  
        end
%% Display Message
        it = it+1;
        if mod(it,10) == 0
            msg  = sprintf('p_s = %.4f, p_d = %.4f',p_s,p_d);
            disp(msg);
            disp(' ');
        end
    end
end

save('mult_equib_exist_sun_ASAD.mat','multeq_save');
