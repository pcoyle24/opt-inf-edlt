%--------------------------------------------------------------------------
%File Name: IRFs_demand_shock.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/SA/LongerCrisis0_85/SA_Fig_3
% IRFs_demand_shock.m
%--------------------------------------------------------------------------

clear all 
close all
clc
% dbstop if error

if ispc
    addpath ..\SA_Calib\
else
    addpath ../SA_Calib/
end

load('params.mat')

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
cPItarg_min_ann = 0; 
cPItarg_max_ann = 2;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg     = [cPItarg_min cPItarg_max]; 
p_s         = 0.995;
p_d         = 0.975; 
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
p_z         = 0.995; 
p_c         = 0.85;%0.75; 

%% Main Code
%Allocate Space for functions
Csz   = zeros(length(cALPHA),length(cPItarg));
PIsz  = zeros(length(cALPHA),length(cPItarg));
Ysz   = zeros(length(cALPHA),length(cPItarg));
Nsz   = zeros(length(cALPHA),length(cPItarg));
Wsz   = zeros(length(cALPHA),length(cPItarg));
Rsz   = zeros(length(cALPHA),length(cPItarg));
Vsz   = zeros(length(cALPHA),length(cPItarg));

Csc   = zeros(length(cALPHA),length(cPItarg));
PIsc  = zeros(length(cALPHA),length(cPItarg));
Ysc   = zeros(length(cALPHA),length(cPItarg));
Nsc   = zeros(length(cALPHA),length(cPItarg));
Wsc   = zeros(length(cALPHA),length(cPItarg));
Rsc   = zeros(length(cALPHA),length(cPItarg));
Rsc_check   = zeros(length(cALPHA),length(cPItarg));
Vsc   = zeros(length(cALPHA),length(cPItarg));

it_count = numel(Vsc);

options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
tic;

it = 0;

for i = 1:length(cALPHA)
    for j = 1:length(cPItarg)
        x0 = ones(10,1);
        func = @(x) demand_solve(x,i,j);

        x_out_markov_sss = fsolve(func,x0,options);

        Csz_out   = x_out_markov_sss(1);
        PIsz_out  = x_out_markov_sss(2);
        Ysz_out   = x_out_markov_sss(3);
        Nsz_out   = x_out_markov_sss(4);
        Vsz_out   = x_out_markov_sss(5);

        Wsz_out = Nsz_out^cCHIn*Csz_out^cCHIc;
        Rsz_out = (cPItarg(j)/(cDELz*cBET))*((PIsz_out/cPItarg(j))^cPHIpi)*((Ysz_out/Ysz_out)^cPHIy);
        if Rsz_out < 1
            Rsz_out = cRzlb;
        end


        Csc_out   = x_out_markov_sss(6);
        PIsc_out  = x_out_markov_sss(7);
        Ysc_out   = x_out_markov_sss(8);
        Nsc_out   = x_out_markov_sss(9);
        Vsc_out   = x_out_markov_sss(10);

        Wsc_out = Nsc_out^cCHIn*Csc_out^cCHIc;
        Rsc_out_check = (cPItarg(j)/(cDELc*cBET))*((PIsc_out/cPItarg(j))^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);

        Rsc_out = (cPItarg(j)/(cDELc*cBET))*((PIsc_out/cPItarg(j))^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);
        if Rsc_out < 1
            Rsc_out = cRzlb;
        end


        Csz(i,j)   = Csz_out;
        PIsz(i,j)  = PIsz_out;
        Ysz(i,j)   = Ysz_out;
        Nsz(i,j)   = Nsz_out;
        Wsz(i,j)   = Wsz_out;
        Rsz(i,j)   = Rsz_out;
        Vsz(i,j)   = Vsz_out;

        Csc(i,j)   = Csc_out;
        PIsc(i,j)  = PIsc_out;
        Ysc(i,j)   = Ysc_out;
        Nsc(i,j)   = Nsc_out;
        Wsc(i,j)   = Wsc_out;
        Rsc(i,j)   = Rsc_out;
        Rsc_check(i,j)   = Rsc_out_check;
        Vsc(i,j)   = Vsc_out;

        it = it+1;

    end
end

%% Plotting
Xsz = zeros(5,2,3);
Xsc = zeros(5,2,3);
X = zeros(15,2,3);

for i = 1:2
    Xsz(:,i,1) = 400*(PIsz(i)-1);
%     Xsz(:,i,2) = 100*(Ysz(i)-1);
    Xsz(:,i,2) = 100*(Csz(i)-1);
    Xsz(:,i,3) = 400*(Rsz(i)-1);


    Xsc(:,i,1) = 400*(PIsc(i)-1);
%     Xsc(:,i,2) = 100*(Ysc(i)-1);
    Xsc(:,i,2) = 100*(Csc(i)-1);
    Xsc(:,i,3) = 400*(Rsc(i)-1);
    
    X(1:5,i,1) =  Xsz(:,i,1);
    X(6:10,i,1) =  Xsc(:,i,1);
    X(11:15,i,1) =  Xsz(:,i,1);
    
    X(1:5,i,2) =  Xsz(:,i,2);
    X(6:10,i,2) =  Xsc(:,i,2);
    X(11:15,i,2) =  Xsz(:,i,2);
    
    
    X(1:5,i,3) =  Xsz(:,i,3);
    X(6:10,i,3) =  Xsc(:,i,3);
    X(11:15,i,3) =  Xsz(:,i,3); 
end
 


header = {'Inflation','Consumption','Policy Rate'};
pi_targ_string = {'0% Inflation Target', '2% Inflation Target'};
colors = {'k','b'};
% y_label = { 'Annualized Percent',{'Percent Deveation from','Efficient Steady State'}, {'Percent Deveation from','Efficient Steady State'}, '', '', ''};

period = 1:1:15;

fig(1) = figure(1);
for i = 1:length(cPItarg)
    for j = 1:length(header)
        subplot(2,3,j)
        box on
        grid on
        hold on
        
        if j == 2
            h(i) = plot(period,X(:,i,j),colors{i},'LineWidth',2);
            set(gca,'YLim',[-10 5],'FontSize',15)
        else
            plot(period,X(:,i,j),colors{i},'LineWidth',2);
        end
            
        set(gca,'XLim',[1 15],'FontSize',15)
        if i == length(cALPHA)  
            xlabel('Periods','FontSize',15)
%             ylabel(y_label{j},'FontSize',15)
            title(header{j},'FontSize',15,'FontWeight','normal')
        end        
    end
end
    
% L = legend([h(1) h(2) h(3) h(4) ],alpha_string{1},alpha_string{2},alpha_string{3},alpha_string{4},'Location', 'NorthWest');
% L = legend([h(1) h(2) h(3) ],alpha_string{1},alpha_string{2},alpha_string{3},'Location', 'NorthWest');
L = legend([h(1) h(2)],pi_targ_string{1},pi_targ_string{2},'Location', 'North');
% L = legend([h(1) h(2)],pi_targ_string{1},pi_targ_string{2},'Position', [0.9 0.9 0.01 0.01]);


set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','demand_only_IRFs.eps');

