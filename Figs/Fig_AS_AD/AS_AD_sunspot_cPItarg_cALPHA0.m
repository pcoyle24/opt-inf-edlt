%--------------------------------------------------------------------------
% File Name: AS_AD_sunspot_cPItarg_cALPHA0.m
% Author: Philip Coyle
% Date Created: 08/31/2018
% Date Updated: 10/17/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/AS_AD_plot
% AS_AD_sunspot_cPItarg_cALPHA0
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
cALPHA_grid      = 0;
cPItarg_grid     = [0/400 + 1, 2/400 + 1]; 
cDELz       = 1;
p_s         = 0.995;
p_d         = 0.975; 

%% Main Code: Generate Data
pi_grid = zeros(101,2);

Sd = zeros(length(cALPHA_grid)*length(cPItarg_grid),2);
St = zeros(length(cALPHA_grid)*length(cPItarg_grid),2);

it = 1;
for a = 1:length(cALPHA_grid)
    cALPHA = cALPHA_grid(a);
    for d = 1:length(cPItarg_grid)
        cPItarg = cPItarg_grid(d);
        % Get Steady State Values
        S = get_ss_sun;
        St(it,1) = S.Ct;
        St(it,2) = S.PIt;
        St(it,3) = S.Yt;
        Sd(it,1) = S.Cd;
        Sd(it,2) = S.PId;
        Sd(it,3) = S.Yd;

        options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

        pi_grid(:,it) = linspace(0.99*S.PId, S.PIt*1.01,101)'; 
        y_guess = ones(length(pi_grid),1)*S.Yd;
        c_guess = ones(length(pi_grid),1)*S.Cd;

        if it == 1
            as_c = zeros(length(pi_grid),2);
            ad_c = zeros(length(pi_grid),2);
            as_y = zeros(length(pi_grid),2);
            ad_y = zeros(length(pi_grid),2);
            alpha = zeros(length(pi_grid),2);
            pitarg = zeros(length(pi_grid),2);
        end

        for i = 1:length(pi_grid)
            state = pi_grid(i,it);
            x_guess_y = y_guess(i);
            x_guess_c = c_guess(i);

            func_as_c = @(x_guess) get_as_sun(x_guess,state,S,1);
            as_out_c = fsolve(func_as_c,x_guess_c,options);
            
            func_as_y = @(x_guess) get_as_sun(x_guess,state,S,0);
            as_out_y = fsolve(func_as_y,x_guess_y,options);

            func_ad_c = @(x_guess) get_ad_sun(x_guess,state,S,1);
            ad_out_c = fsolve(func_ad_c,x_guess_c,options);
            
            func_ad_y = @(x_guess) get_ad_sun(x_guess,state,S,0);
            ad_out_y = fsolve(func_ad_y,x_guess_y,options);

            as_c(i,it) = as_out_c;
            ad_c(i,it) = ad_out_c;
            
            as_y(i,it) = as_out_y;
            ad_y(i,it) = ad_out_y;
            
            alpha(i,it) = cALPHA;
            pitarg(i,it) = cPItarg;
        end
        it = it + 1;
    end
end

%% Prep Data for plotting
r = max(pitarg./cBET.*(pi_grid./pitarg).^cPHIpi,1);
r_inx(1) = find(r(:,1) > 1,1);
r_inx(2) = find(r(:,2) > 1,1);

as_c_sun1 = as_c(1:r_inx(1)-1,1,1);
as_c_sun2= as_c(1:r_inx(2)-1,2,1);
ad_c_sun1 = ad_c(1:r_inx(1)-1,1,1);
ad_c_sun2= ad_c(1:r_inx(2)-1,2,1);
pi_grid1 = pi_grid(1:r_inx(1)-1,1);
pi_grid2 = pi_grid(1:r_inx(2)-1,2);
save('AS_AD_sun_deflationary.mat','as_c_sun1','as_c_sun2','ad_c_sun1','ad_c_sun2','pi_grid1','pi_grid2','Sd');


%% Plotting

fig(1) = figure;
subplot(2,2,1)
box on
grid on
hold on

h1(1) = plot(100*(as_c_sun1-1),400*(pi_grid1-1),'LineWidth',2,'Color','k','LineStyle','-');
h2(1) = plot(100*(ad_c_sun1-1),400*(pi_grid1-1),'LineWidth',2,'Color','b','LineStyle','-'); 
h = plot(100*(Sd(1,1)-1),400*(Sd(1,2)-1));
set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',6)

h1(2) = plot(100*(as_c_sun2-1),400*(pi_grid2-1),'LineWidth',2,'Color','k','LineStyle','--');
h2(2) = plot(100*(ad_c_sun2-1),400*(pi_grid2-1),'LineWidth',2,'Color','b','LineStyle','--');  
h = plot(100*(Sd(2,1)-1),400*(Sd(2,2)-1));
set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',6)

% set(gca,'XLim',[-3 -0.5],'Ylim',[-1.5 -1],'YTick',(-1.5:0.1:-1),'FontSize',25)
set(gca,'XLim',[-3.25 -0.25],'Ylim',[-5 -0],'FontSize',15)
xlabel('Consumption (C)','FontSize',15)
ylabel('Inflation (\Pi)','FontSize',15)
title('Deflationary Regime','FontSize',15,'FontWeight','normal')

L = legend([h1(1) h2(1) h1(2) h2(2)],'AS Curve (\Pi^{targ} = 0%)','AD Curve (\Pi^{targ} = 0%)','AS Curve (\Pi^{targ} = 2%)','AD Curve (\Pi^{targ} = 2%)');
% L = legend([h1(1) h2(1) h1(2) h2(2)],'Phillips Curve (\Pi^{targ} = 0%)','Euler Equation (\Pi^{targ} = 0%)','Phillips Curve (\Pi^{targ} = 2%)','Euler Equation (\Pi^{targ} = 2%)');
set(L,'Location','SouthEast','FontSize',10)

savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','AS_AD_plot_cALPHAstar_sunspot_cALPHA0.eps');
print(fig(1),'-depsc',strcat(savedir,'AS_AD_plot_cALPHAstar_sunspot_cALPHA0.eps'));
