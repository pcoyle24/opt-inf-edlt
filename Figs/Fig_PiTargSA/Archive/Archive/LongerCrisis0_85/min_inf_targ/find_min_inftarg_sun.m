%--------------------------------------------------------------------------
%File Name: find_min_inftarg_sun.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/SA/LongerCrisis0_85/min_inf_targ
% find_min_inftarg_sun.m
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
cPItarg_min_ann = -2; 
cPItarg_max_ann = 0;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n); 
p_s         = 0.995; 
p_d         = 0.975; 
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
p_z         = 0.995; 
p_c         = 0.85;%0.75; 

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

eqm_exist = find(Rs >= 1,1);
min_inf_targ = round(400*(cPItarg(eqm_exist)-1),2);
disp(char(strcat('The minimum inflation target consistent with equilibrium existence is ',num2str(min_inf_targ),'%')));

%% Plotting
header = {'Inflation','Consumption','Policy Rate'};

X = [PIs;Cs;Rs];

fig(1) = figure(1);
for j = 1:length(header)
    subplot(2,3,j)
    box on
    grid on
    hold on

    if j == 1 && 3
        h(i) = plot(400*(cPItarg-1),400*(X(j,:)-1),'k','LineWidth',2);
        ax = get(gca,'YLim');
        set(gca,'XTick',[-2,min_inf_targ,0],'YLim',[ax(1) ax(2)],'FontSize',15)
        line([400*(cPItarg(eqm_exist)-1) 400*(cPItarg(eqm_exist)-1)],[ax(1) ax(2)],'Color','k','LineStyle','-','LineWidth',1);
    else
        plot(400*(cPItarg-1),100*(X(j,:)-1),'k','LineWidth',2);
        ax = get(gca,'YLim');
        set(gca,'XTick',[-2,min_inf_targ,0],'YLim',[ax(1) ax(2)],'FontSize',15)
        line([400*(cPItarg(eqm_exist)-1) 400*(cPItarg(eqm_exist)-1)],[ax(1) ax(2)],'Color','k','LineStyle','-','LineWidth',1);
    end

    set(gca,'XLim',[-2 0],'FontSize',15)
    if i == length(cALPHA)  
        xlabel('Inflation Target','FontSize',15)
        title(header{j},'FontSize',15,'FontWeight','normal')
    end        
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','min_inftarg_sun.eps');

save('eqm_exist_sun','min_inf_targ');