%--------------------------------------------------------------------------
% File Name: AS_AD_low_pitarg_plot.m
% Author: Philip Coyle
% Date Created: 03/04/2019
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/2State_Markov_Shock/AS_AD_Plot/MultipleZLBEq
% AS_AD_low_pitarg_plot.m
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
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg p_s p_d cDELz 
cBET            = 1/(1.0025);
cCHIc           = 1;
cCHIn           = 1;
cTHETA          = 11;
cTAU            = 1/cTHETA;
cVARPHI         = params(1);
cPHIpi          = 2;
cPHIy           = 0;
cRzlb           = 1;
cIOTA           = 1;
cALPHA          = params(3);
cPItarg         = -1/400 + 1; 
cDELz           = 1;
p_s             = 1;
p_d             = 0.975;
c_check_grid    = [0 1];

for cc = 1:length(c_check_grid)
    c_check = c_check_grid(cc);

    %% Main Code
    % Get Steady State Values
    S = get_ss_sun;
    options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
    
    pi_grid = linspace(-8/400 + 1, 8/400 + 1,501);
    pi_grid = sort([pi_grid,cPItarg]);

    if c_check == 0
        y_guess = ones(length(pi_grid),1)*S.Yd;
    elseif c_check == 1
        c_guess = ones(length(pi_grid),1)*S.Cd;
    end

    as = zeros(length(pi_grid),1);
    ad = zeros(length(pi_grid),1);

    for i = 1:length(pi_grid)
        state = pi_grid(i);
        if c_check == 0
            x_guess = y_guess(i);
        elseif c_check == 1
            x_guess = c_guess(i);
        end

        if c_check == 0
            func_as = @(x_guess) get_as_sun(x_guess,state,S,c_check);
        elseif c_check == 1
            func_as = @(x_guess) get_as_sun(x_guess,state,S,c_check);
        end
        as_out = fsolve(func_as,x_guess,options);

        state = pi_grid(i);

        if c_check == 0
            func_ad = @(x_guess) get_ad_sun(x_guess,state,S,c_check);
        elseif c_check == 1
            func_ad = @(x_guess) get_ad_sun(x_guess,state,S,c_check);
        end
        ad_out = fsolve(func_ad,x_guess,options);

        as(i) = as_out;
        ad(i) = ad_out;
    end

    sort_asad = sort(abs(as-ad));

    %% Plotting
    r = max(cPItarg/cBET*(pi_grid/cPItarg).^cPHIpi,1);
    zlb_inx = find(r > 1,1);
    
    fig(1) = figure;
    subplot(2,2,1)
    box on
    grid on
    hold on
    fill([-100,20,20,-100],[400*(pi_grid(zlb_inx-1)-1),400*(pi_grid(zlb_inx-1)-1),10,10],[0.85 0.85 0.85]);
    hh(1) = plot(100*(as-1),400*(pi_grid-1),'LineWidth',2,'Color','k');
    hh(2) = plot(100*(ad-1),400*(pi_grid-1),'LineWidth',2,'Color','b');
    set(gca,'Xlim',[-100 20],'Ylim',[-10 10],'FontSize',15)
    line([-2 -2],[-1.5 -0.5],'linewidth',1,'color','k')
    line([2 2],[-1.5 -0.5],'linewidth',1,'color','k')
    line([-2 2],[-1.5 -1.5],'linewidth',1,'color','k')
    line([-2 2],[-0.5 -0.5],'linewidth',1,'color','k')
    
    if c_check == 0
        xlabel('Output (Y)','FontSize',15)
    elseif c_check == 1
%         L = legend([hh(1) hh(2)],'Phillips Curve','Euler Equation');
        L = legend([hh(1) hh(2) ],'AS Curve','AD Curve');
        set(L,'Location','SouthEast','FontSize',15)
        xlabel('Consumption (C)','FontSize',15)
    end
    ylabel('Inflation (\Pi)','FontSize',15)
    
    axes('position',[.17 .63 .1 .1])
    box on % put box around new pair of axes
    hold on
    fill([-0.1,0.1,0.1,-0.1],[400*(pi_grid(zlb_inx-1)-1),400*(pi_grid(zlb_inx-1)-1),-0.99,-0.99],[0.85 0.85 0.85]);
    plot(100*(as-1),400*(pi_grid-1),'LineWidth',2,'Color','k');
    plot(100*(ad-1),400*(pi_grid-1),'LineWidth',2,'Color','b');
    set(gca,'Xlim',[-0.1 0.1],'Ylim',[-1.01 -0.99],'FontSize',10)
    
    
    set(fig(1),'PaperOrientation','Landscape');
    set(fig(1),'PaperPosition',[0 0 11 8.5]);
    
    savedir = cd;
    savedir = fullfile(savedir, '..');
    if ispc 
        savedir = strcat(savedir,'\Final\');
    else
        savedir = strcat(savedir,'/Final/');
    end
    
    if c_check == 0
        print(fig(1),'-depsc','AS_AD_low_pitarg_output_alt.eps');
        print(fig(1),'-depsc',strcat(savedir,'AS_AD_low_pitarg_output_alt.eps'));
    elseif c_check == 1
        print(fig(1),'-depsc','AS_AD_low_pitarg_consumption_alt.eps');
        print(fig(1),'-depsc',strcat(savedir,'AS_AD_low_pitarg_consumption_alt.eps'));
    end
end