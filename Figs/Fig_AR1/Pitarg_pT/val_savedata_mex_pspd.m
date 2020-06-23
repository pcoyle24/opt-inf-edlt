% --------------------------------------------------------------------------
% File Name: val_savedata_mex_pspd.m
% Author: Philip Coyle
% Date Created: 02/06/2018
% Last Updated: 07/19/2018
% 
% Add Appropriate Paths %
% codetoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/Code_Toolbox'; 
% mextoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/mex_functions';
% addpath(genpath(codetoolbox));
% addpath(genpath(mextoolbox)); 
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_AR1/Pitarg_pT
% matlab -nodesktop -nosplash -r val_savedata_mex_pspd
% --------------------------------------------------------------------------


clear all
close all 
clc

% Add path to saved mat files
addpath ../common/
addpath ../mex_functions/
addpath ../savedata101/

% Load parameters
% P = parameters2; % 10 bps grid interval
P = parameters3; % 5 bps grid interval

run_first_time = 0;

if run_first_time == 1
    % Start Parallel Pool
    gcp;
    
    P.bound = (P.sigma^2/(1-P.rho^2))^0.5;
    O.delbound = [1-4*P.bound 1+4*P.bound];
    O.del_pts = 101;
    O.e_pts = 10;

    % Chebyshev polynomial order in each dimension (must <= pts in grid)
    O.n1 = 4;
    value = zeros(length(P.pi_targ),length(P.Ps_grid),length(P.Pd_grid));
    value_inx = zeros(length(P.Ps_grid),length(P.Pd_grid));
    for d = 1:length(P.Pd_grid)
        P.Pd = P.Pd_grid(d);
        disp(strcat('P_d =',num2str(P.Pd)));  

        for s = length(P.Ps_grid):-1:1
            P.Ps = P.Ps_grid(s);
            disp(strcat('P_s =',num2str(P.Ps))); 

            for i = 1:length(P.pi_targ)
                load(strcat('ChebPFs_PItarg',num2str(400*(P.pi_targ(i)-1)),'_Ps',num2str(P.Ps),'_Pd',num2str(P.Pd),'.mat'))

                E_val = sim_val_zlb(P,O,C);
                value(i,s,d) = E_val;

                disp(strcat('Value calculated for PItarg =',num2str(400*(P.pi_targ(i)-1))))
            end
            max_value = max(value(:,s,d));
            value_inx(s,d) = find(max_value == value(:,s,d));
    %         if value_inx(s,d) == 1
    %             value_inx(1:s-1,d) = value_inx(s,d);
    %             break;
    %         end       
        end
    end

    opt_inf = 400*(P.pi_targ(value_inx)-1);
    save('Ev_OptInf_sun_5bps.mat','opt_inf');
else
    load 'Ev_OptInf_sun_5bps.mat'
end

%% Plotting
colors = {'-k'};

fig(1) = figure(1);
subplot(2,2,1)
box on
hold on
grid on
for i = 1:length(colors)
    h(i) = plot(P.Ps_grid(1:1:end),opt_inf(1:1:end),colors{i},'LineWidth',2);
end
line([0.995 0.995],get(gca,'YLim'),'Color','k','LineStyle','-','LineWidth',1)
line([0.999 0.999],get(gca,'YLim'),'Color','k','LineStyle','-.','LineWidth',1)

xlabel('p_T (Persistance of Target Regime)','FontSize',12)
ylabel('Optimal Inflation Target','FontSize',12)
set(gca,'XLim',[P.Ps_grid(1) P.Ps_grid(end)],'XTick',(P.Ps_grid(1):0.002:P.Ps_grid(end)),'Ylim',[1 2],'FontSize',12)

savedir = cd;
savedir = fullfile(savedir, '../..');
savedir = strcat(savedir,'/Final/');

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','OptInf_pTpD_sun_5bps.eps');
print(fig(1),'-depsc',strcat(savedir,'OptInf_pTpD_sun_5bps.eps'));
