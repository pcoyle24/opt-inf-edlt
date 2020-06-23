% --------------------------------------------------------------------------
% File Name: plot_irfs.m
% Author: Philip Coyle
% Date Created: 11/09/2018
% 
% Add Appropriate Paths 
% codetoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/Code_Toolbox'; 
% mextoolbox = '/mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/Codes/mex_functions';
% addpath(genpath(codetoolbox));
% addpath(genpath(mextoolbox)); 
% 
% Run Code %
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_AR1/PFs_and_IRFs
% plot_irfs
% --------------------------------------------------------------------------


clear all
close all 
clc

% Directory to save data
addpath ../common/
addpath ../mex_functions/
addpath ../savedata101/

%% ------------------------------------------------------------------------
%  Initialize Policy Functions
%  ------------------------------------------------------------------------
% Load parameters
P = parameters;

%% ------------------------------------------------------------------------
%  IRF Analysis
%  ------------------------------------------------------------------------
time = 1:40;
per = 40;
c_s_irf = zeros(per,2);
pi_s_irf = zeros(per,2);
r_s_irf = zeros(per,2);

c_d_irf = zeros(per,2);
pi_d_irf = zeros(per,2);
r_d_irf = zeros(per,2);

X = zeros(per,2,6);

shock = 3.75*P.bound;

for i = 1:length(P.pi_targ)
    load(strcat('ChebPFs_PItarg',num2str(400*(P.pi_targ(i)-1)),'.mat'));
    
    c_s_init = S.c_s;
    pi_s_init = S.inf_s;
    r_s_init = S.r_s;
    
    c_d_init = S.c_d;
    pi_d_init = S.inf_d;
    r_d_init = S.r_d;
    
    del_yesterday = 1;
    
    for k = 1:per
        if k == 1
            del_today = P.rho*(del_yesterday - 1) + 1 + shock;
        else
            del_today = P.rho*(del_yesterday - 1) + 1;
        end
        % Target Regime
        r_s_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ar_s,C.max,C.T,C.P);
        if r_s_today >=1
            c_s_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_s,C.max,C.T,C.P);
            pi_s_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_s,C.max,C.T,C.P);
        else
            r_s_today = 1;
            c_s_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_zlb_s,C.max,C.T,C.P);
            pi_s_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_zlb_s,C.max,C.T,C.P);
        end
        
        c_s_irf(k,i) = 100*(c_s_today-1);
        pi_s_irf(k,i) = 400*(pi_s_today-1);
        r_s_irf(k,i) = 400*(r_s_today-1);
        
        % Deflationary Regime
        r_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ar_d,C.max,C.T,C.P);
        if r_d_today >=1
            c_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_d,C.max,C.T,C.P);
            pi_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_d,C.max,C.T,C.P);
        else
            r_d_today = 1;
            c_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_zlb_d,C.max,C.T,C.P);
            pi_d_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_zlb_d,C.max,C.T,C.P);
        end
        
        c_d_irf(k,i) = 100*(c_d_today-1);
        pi_d_irf(k,i) = 400*(pi_d_today-1);
        r_d_irf(k,i) = 400*(r_d_today-1);
        
        del_yesterday = del_today;
    end
end 
X(:,:,1) = c_s_irf;
X(:,:,2) = pi_s_irf;
X(:,:,3) = r_s_irf;
X(:,:,4) = c_d_irf;
X(:,:,5) = pi_d_irf;
X(:,:,6) = r_d_irf;


%% ------------------------------------------------------------------------
%  Plot IRFs
%  ------------------------------------------------------------------------

title_fig = {'Consumption (T)', 'Inflation Rate (T)', 'Interest Rate (T)','Consumption (D)', 'Inflation Rate (D)', 'Interest Rate (D)'};
colors = {'k', 'b'};

fig(1) = figure(1);
for i = 1:length(title_fig)
    for k = 1:length(colors)
        subplot(2,3,i)
        box on
        hold on
        grid on
        if i == 6
            h(k) = plot(time,X(:,k,i),colors{k},'LineWidth',2);
        else
            plot(time,X(:,k,i),colors{k},'LineWidth',2)
        end
        xlabel('Periods','FontSize',16)
        title(title_fig{i},'FontSize',16,'FontWeight','Normal')
        set(gca,'XLim',[time(1) time(end)],'FontSize',16)
    end
end
L = legend([h(1) h(2)],'\Pi^{targ} = 1.5%', '\Pi^{targ} = 2%');
set(L,'Location','SouthEast','FontSize',15)

savedir = cd;
savedir = fullfile(savedir, '../..');
savedir = strcat(savedir,'/Final/');

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','IRFs_sunspot.eps');
print(fig(1),'-depsc',strcat(savedir,'IRFs_sunspot.eps'));