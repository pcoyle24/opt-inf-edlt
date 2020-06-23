% --------------------------------------------------------------------------
% File Name: val_savedata_mex.m
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
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_AR1/Welfare
% matlab -nodesktop -nosplash -r val_savedata_mex
% --------------------------------------------------------------------------


clear all
close all 
clc

% Directory to save data
addpath ../common/
addpath ../mex_functions/
addpath ../savedata101/

% Load parameters
P = parameters1;
P.bound = (P.sigma^2/(1-P.rho^2))^0.5;
O.delbound = [1-4*P.bound 1+4*P.bound];
O.del_pts = 101;
O.e_pts = 10;

run_first_time = 0;

if run_first_time == 1
    % Start Parallel Pool
    gcp;
    
    % Chebyshev polynomial order in each dimension (must <= pts in grid)
    O.n1 = 4;
    value = zeros(length(P.pi_targ),1);

    for i = 1:length(P.pi_targ)
        load(strcat('ChebPFs_PItarg',num2str(400*(P.pi_targ(i)-1)),'_Ps',num2str(P.Ps),'_Pd',num2str(P.Pd),'.mat'))

        E_val = sim_val_zlb(P,O,C);
        value(i,1) = E_val;

        disp(strcat('Value calculated for PItarg =',num2str(400*(P.pi_targ(i)-1))))
    end

    max_value = max(value(:,1));
    value_inx = find(max_value == value(:,1));
    opt_inf = 400*(P.pi_targ(value_inx)-1);
    disp(opt_inf);


    %% Save results
    value_ds = value;
    opt_inf_ds = opt_inf;

    save('Ev_sun_10bps.mat','value_ds','opt_inf_ds')
else
    load('Ev_sun_10bps.mat');
    opt_inf = opt_inf_ds;  
    value = value_ds;
end

%% Plotting
fig(1) = figure(1);
subplot(2,2,1)
box on
hold on
grid on
plot(400*(P.pi_targ(1:1:end)-1),value,'k','LineWidth',2)
xlabel('\Pi^{targ}','FontSize',16)
set(gca,'XLim',[400*(P.pi_targ(1)-1) 400*(P.pi_targ(end)-1)],'YLim',[min(value) - 0.2 max(value) + 0.2],'FontSize',16)
line([opt_inf opt_inf],get(gca,'YLim'),'Color','b','LineStyle','-','LineWidth',1.5)




savedir = cd;
savedir = fullfile(savedir, '../..');
savedir = strcat(savedir,'/Final/');

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','welfare_sun.eps');
print(fig(1),'-depsc',strcat(savedir,'welfare_sun.eps'));