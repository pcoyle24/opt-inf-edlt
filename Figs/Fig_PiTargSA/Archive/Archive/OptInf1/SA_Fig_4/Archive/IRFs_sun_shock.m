%--------------------------------------------------------------------------
%File Name: IRFs_sun_shock.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/SS_Opt_Inf/Model_2\Ex_1'
% IRFs_sun_shock.m
%--------------------------------------------------------------------------

clear all 
close all
clc
dbstop if error

%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c
cBET        = 1/(1.0025);
cCHIc       = 1;
cCHIn       = 1;
cTHETA      = 11;
cTAU        = 1/cTHETA;
cVARPHI     = 1250; % For an approx 200 bps decline in inflation 
% cVARPHI     = 2550; % For an approx 100 bps decline in inflation 
cPHIpi      = 2;
cPHIy       = 0;
cRzlb       = 1;
cIOTA       = 1;
cALPHA      = 0.876; % For an approx 200 bps decline in inflation 
% cALPHA      = 0.915; % For an approx 100 bps decline in inflation 
cPItarg_min_ann = 0; 
cPItarg_max_ann = 2;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg     = [cPItarg_min cPItarg_max]; 
p_s         = 0.995; %(0.85:0.01:1);
p_d         = 0.975; %(0.85:0.01:1);
cDELz       = 1;
c           = 0.20657/10;%0.198/10;
cDELc       = 1 + c;
p_z         = 0.995;  % 0.99  
p_c         = 0.75;   % 0.75, 0.8, 5/6,  

%Allocate Space for functions
Cs   = zeros(length(cALPHA),length(cPItarg));
PIs  = zeros(length(cALPHA),length(cPItarg));
Ys   = zeros(length(cALPHA),length(cPItarg));
Ns   = zeros(length(cALPHA),length(cPItarg));
Ws   = zeros(length(cALPHA),length(cPItarg));
Rs   = zeros(length(cALPHA),length(cPItarg));
Vs   = zeros(length(cALPHA),length(cPItarg));

Cd   = zeros(length(cALPHA),length(cPItarg));
PId  = zeros(length(cALPHA),length(cPItarg));
Yd   = zeros(length(cALPHA),length(cPItarg));
Nd   = zeros(length(cALPHA),length(cPItarg));
Wd   = zeros(length(cALPHA),length(cPItarg));
Rd   = zeros(length(cALPHA),length(cPItarg));
Vd   = zeros(length(cALPHA),length(cPItarg));

it_count = numel(Vd);

options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
tic;

it = 0;

for i = 1:length(cALPHA)
    for j = 1:length(cPItarg)
        x0 = ones(10,1);
        func = @(x) sun_solve(x,i,j);

        x_out_sun = fsolve(func,x0,options);

        Cs_out   = x_out_sun(1);
        PIs_out  = x_out_sun(2);
        Ys_out   = x_out_sun(3);
        Ns_out   = x_out_sun(4);
        Vs_out   = x_out_sun(5);

        Ws_out = Ns_out^cCHIn*Cs_out^cCHIc;
        Rs_out = (cPItarg(j)/cBET)*((PIs_out/cPItarg(j))^cPHIpi)*((Ys_out/Ys_out)^cPHIy);

        Cd_out   = x_out_sun(6);
        PId_out  = x_out_sun(7);
        Yd_out   = x_out_sun(8);
        Nd_out   = x_out_sun(9);
        Vd_out   = x_out_sun(10);

        Wd_out = Nd_out^cCHIn*Cd_out^cCHIc;
        Rd_out = cRzlb;              


        Cs(i,j)   = Cs_out;
        PIs(i,j)  = PIs_out;
        Ys(i,j)   = Ys_out;
        Ns(i,j)   = Ns_out;
        Ws(i,j)   = Ws_out;
        Rs(i,j)   = Rs_out;
        Vs(i,j)   = Vs_out;

        Cd(i,j)   = Cd_out;
        PId(i,j)  = PId_out;
        Yd(i,j)   = Yd_out;
        Nd(i,j)   = Nd_out;
        Wd(i,j)   = Wd_out;
        Rd(i,j)   = Rd_out;
        Vd(i,j)   = Vd_out;

        it = it+1;
    end
end

%% Plotting
Xs = zeros(5,2,3);
Xd = zeros(5,2,3);
X = zeros(15,2,3);

for i = 1:2
    Xs(:,i,1) = 400*(PIs(i)-1);
%     Xs(:,i,2) = 100*(Ys(i)-1);
    Xs(:,i,2) = 100*(Cs(i)-1);
    Xs(:,i,3) = 400*(Rs(i)-1);


    Xd(:,i,1) = 400*(PId(i)-1);
%     Xd(:,i,2) = 100*(Yd(i)-1);
    Xd(:,i,2) = 100*(Cd(i)-1);
    Xd(:,i,3) = 400*(Rd(i)-1);
    
    X(1:5,i,1) =  Xs(:,i,1);
    X(6:10,i,1) =  Xd(:,i,1);
    X(11:15,i,1) =  Xs(:,i,1);
    
    X(1:5,i,2) =  Xs(:,i,2);
    X(6:10,i,2) =  Xd(:,i,2);
    X(11:15,i,2) =  Xs(:,i,2);
    
    
    X(1:5,i,3) =  Xs(:,i,3);
    X(6:10,i,3) =  Xd(:,i,3);
    X(11:15,i,3) =  Xs(:,i,3); 
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
            set(gca,'YLim',[-5 10],'FontSize',15)
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

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
if cVARPHI == 1250 && cALPHA == 0.876
    print(fig(1),'-depsc','sunspot_only_IRFs_200bps.eps');
elseif cVARPHI == 2550 && cALPHA == 0.915
    print(fig(1),'-depsc','sunspot_only_IRFs_100bps.eps');
end
