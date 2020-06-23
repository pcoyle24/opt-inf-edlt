%--------------------------------------------------------------------------
% File Name: AS_AD_demand_cPItarg_cALPHA1.m
% Author: Philip Coyle
% Date Created: 08/31/2018
% Date Updated: 10/17/2018
% /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/AS_AD_plot
% AS_AD_demand_cPItarg_cALPHA1
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
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c
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
cALPHA_grid      = 1;
cPItarg_grid     = [0/400 + 1, 2/400 + 1]; 
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
p_z         = 0.995;
p_c         = 0.75;   
%% Main Code: Generate Data
z_grid = [0,1];

Sz = zeros(length(cALPHA_grid)*length(cPItarg_grid),2);
Sc = zeros(length(cALPHA_grid)*length(cPItarg_grid),2);

for z = 0:1
    it = 1;
    for a = 1:length(cALPHA_grid)
        cALPHA = cALPHA_grid(a);
        for d = 1:length(cPItarg_grid)
            cPItarg = cPItarg_grid(d);
            % Get Steady State Values
            S = get_ss_demand;
            if z == 1
                Sz(it,1) = S.Ctz;
                Sz(it,2) = S.PItz;
                Sz(it,3) = S.Ytz;
                Sc(it,1) = S.Ctc;
                Sc(it,2) = S.PItc;
                Sc(it,3) = S.Ytc;
            end

            options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

    %          pi_grid(:,it) = linspace(0.99*S.PItc,S.PItz*1.01,101)'; 
            pi_grid = linspace(-5/400+1,5/400+1,101)'; 
            y_guess = ones(length(pi_grid),1)*S.Ytz;
            c_guess = ones(length(pi_grid),1)*S.Ctz;

            if z == 0 && it == 1
                as_c = zeros(length(pi_grid),4,2);
                ad_c = zeros(length(pi_grid),4,2);
                as_y = zeros(length(pi_grid),4,2);
                ad_y = zeros(length(pi_grid),4,2);
                alpha = zeros(length(pi_grid),4,2);
                pitarg = zeros(length(pi_grid),4,2);
            end

            for i = 1:length(pi_grid)
                state = pi_grid(i);
                x_guess_y = y_guess(i);
                x_guess_c = c_guess(i);

                func_as_c = @(x_guess) get_as_demand(x_guess,state,S,z,1);
                as_out_c = fsolve(func_as_c,x_guess_c,options);
                
                func_as_y = @(x_guess) get_as_demand(x_guess,state,S,z,0);
                as_out_y = fsolve(func_as_y,x_guess_y,options);

                func_ad_c = @(x_guess) get_ad_demand(x_guess,state,S,z,1);
                ad_out_c = fsolve(func_ad_c,x_guess_c,options);
                
                func_ad_y = @(x_guess) get_ad_demand(x_guess,state,S,z,0);
                ad_out_y = fsolve(func_ad_y,x_guess_y,options);

                as_c(i,it,z+1) = as_out_c;
                ad_c(i,it,z+1) = ad_out_c;
                
                as_y(i,it,z+1) = as_out_y;
                ad_y(i,it,z+1) = ad_out_y;


                alpha(i,it,z+1) = cALPHA;
                pitarg(i,it,z+1) = cPItarg;

            end
            it = it + 1;
        end
    end
end

%% Prep Data for Plotting
r_inx = zeros(length(cALPHA_grid)*length(cPItarg_grid),2);

for z = 0:1
    it = 1;
    for a = 1:length(cALPHA_grid)
        cALPHA = cALPHA_grid(a);
        for p = 1:length(cPItarg_grid)
            cPItarg = cPItarg_grid(p);
            if z == 1
                r = max(cPItarg./(cDELz*cBET).*(pi_grid./cPItarg).^cPHIpi,1);
                r_inx(it,2) = find(r > 1,1);
            else
                r = max(cPItarg./(cDELc*cBET).*(pi_grid./cPItarg).^cPHIpi,1);
                r_inx(it,1) = find(r > 1,1);
            end           
            it = it + 1;
        end
    end
end

as_c_demand1 = as_c(1:r_inx(1,1)-1,1,1);
as_c_demand2= as_c(1:r_inx(2,1)-1,2,1);
ad_c_demand1 = ad_c(1:r_inx(1,1)-1,1,1);
ad_c_demand2= ad_c(1:r_inx(2,1)-1,2,1);
pi_grid1 = pi_grid(1:r_inx(1,1)-1);
pi_grid2 = pi_grid(1:r_inx(2,1)-1);

save('AS_AD_demand_crisis.mat','as_c_demand1','as_c_demand2','ad_c_demand1','ad_c_demand2','pi_grid1','pi_grid2','Sc');

%% Plotting
fig(1) = figure;
subplot(2,2,1)
box on
grid on
hold on

h1(1) = plot(100*(as_c_demand1-1),400*(pi_grid1-1),'LineWidth',2,'Color','k','LineStyle','-');
h2(1) = plot(100*(ad_c_demand1-1),400*(pi_grid1-1),'LineWidth',2,'Color','b','LineStyle','-'); 
h = plot(100*(Sc(1,1)-1),400*(Sc(1,2)-1));
set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',6)

h1(2) = plot(100*(as_c_demand2-1),400*(pi_grid2-1),'LineWidth',2,'Color','k','LineStyle','--');
h2(2) = plot(100*(ad_c_demand2-1),400*(pi_grid2-1),'LineWidth',2,'Color','b','LineStyle','--');  
h = plot(100*(Sc(2,1)-1),400*(Sc(2,2)-1));
set(h,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',6)

% set(gca,'XLim',[-8 -4],'Ylim',[-3 1],'FontSize',25)
set(gca,'XLim',[-10 5],'Ylim',[-5 2],'FontSize',15)
xlabel('Consumption (C)','FontSize',15)
ylabel('Inflation (\Pi)','FontSize',15)
title('Crisis State','FontSize',15,'FontWeight','normal')

% L = legend([h1(1) h2(1) h1(2) h2(2)],'AS Curve \Pi^{targ} = 0%','AD Curve \Pi^{targ} = 0%','AS Curve \Pi^{targ} = 2%','AD Curve \Pi^{targ} = 2%');
% L = legend([h1(1) h2(1) h1(2) h2(2)],'Phillips Curve \Pi^{targ} = 0%','Euler Equation \Pi^{targ} = 0%','Phillips Curve \Pi^{targ} = 2%','Euler Equation \Pi^{targ} = 2%');
% set(L,'Position',[0.51 0.45 0.01 0.01],'Orientation','horizontal','FontSize',12)
% legend('boxoff');
% set(L,'Location','SouthEast','FontSize',20)

savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','AS_AD_plot_cALPHAstar_demand_cALPHA1.eps');
print(fig(1),'-depsc',strcat(savedir,'AS_AD_plot_cALPHAstar_demand_cALPHA1.eps'));