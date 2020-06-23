%--------------------------------------------------------------------------
%File Name: equib_exist_demand.m
%Author: Philip Coyle
%Date Created: 10/17/2018
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/Eq_Existance
% equib_exist_demand
%--------------------------------------------------------------------------

clear all 
close all
clc
% dbstop if error

if ispc
    addpath ..\Calib\
else
    addpath ../Calib/
end

load('params.mat')

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
cALPHA      = params(3);
cPItarg_ann = 2; 
cPItarg     = cPItarg_ann/400 + 1;
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
p_z_grid         = (0:0.01:1);
p_c_grid         = (0:0.01:1);

%% Main Code
%Allocate Space for functions
Csz   = zeros(length(p_z_grid),length(p_c_grid));
PIsz  = zeros(length(p_z_grid),length(p_c_grid));
Ysz   = zeros(length(p_z_grid),length(p_c_grid));
Nsz   = zeros(length(p_z_grid),length(p_c_grid));
Wsz   = zeros(length(p_z_grid),length(p_c_grid));
Rsz   = zeros(length(p_z_grid),length(p_c_grid));

Csc   = zeros(length(p_z_grid),length(p_c_grid));
PIsc  = zeros(length(p_z_grid),length(p_c_grid));
Ysc   = zeros(length(p_z_grid),length(p_c_grid));
Nsc   = zeros(length(p_z_grid),length(p_c_grid));
Wsc   = zeros(length(p_z_grid),length(p_c_grid));
Rsc   = zeros(length(p_z_grid),length(p_c_grid));
flag_save   = nan(length(p_z_grid),length(p_c_grid));

it_count = numel(Csc);

options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000,'TolFun', 1e-32, 'Display', 'none');
tic;

it = 0;

for i = 1:length(p_z_grid)
    p_z = p_z_grid(i);
    for j = 1:length(p_c_grid)
        p_c = p_c_grid(j);
%         x0 = ones(10,1);
        x0 = ones(8,1);
        func = @(x) demand_solve(x);

        [x_out_markov_sss,~,flag] = fsolve(func,x0,options);
        flag_save(i,j) = flag;

        Csz_out   = x_out_markov_sss(1);
        PIsz_out  = x_out_markov_sss(2);
        Ysz_out   = x_out_markov_sss(3);
        Nsz_out   = x_out_markov_sss(4);
%         Vsz_out   = x_out_markov_sss(5);

        Wsz_out = Nsz_out^cCHIn*Csz_out^cCHIc;
        Rsz_out = (cPItarg/(cDELz*cBET))*((PIsz_out/cPItarg)^cPHIpi)*((Ysz_out/Ysz_out)^cPHIy);
        if Rsz_out < 1
            Rsz_out = cRzlb;
        end


        Csc_out   = x_out_markov_sss(5);
        PIsc_out  = x_out_markov_sss(6);
        Ysc_out   = x_out_markov_sss(7);
        Nsc_out   = x_out_markov_sss(8);
%         Vsc_out   = x_out_markov_sss(10);

        Wsc_out = Nsc_out^cCHIn*Csc_out^cCHIc;
        
        Rsc_out = (cPItarg/(cDELc*cBET))*((PIsc_out/cPItarg)^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);
        if Rsc_out < 1
            Rsc_out = cRzlb;
        end


        Csz(i,j)   = Csz_out;
        PIsz(i,j)  = PIsz_out;
        Ysz(i,j)   = Ysz_out;
        Nsz(i,j)   = Nsz_out;
        Wsz(i,j)   = Wsz_out;
        Rsz(i,j)   = Rsz_out;
%         Vsz(i,j)   = Vsz_out;

        Csc(i,j)   = Csc_out;
        PIsc(i,j)  = PIsc_out;
        Ysc(i,j)   = Ysc_out;
        Nsc(i,j)   = Nsc_out;
        Wsc(i,j)   = Wsc_out;
        Rsc(i,j)   = Rsc_out;
%         Rsc_check(i,j)   = Rsc_out_check;
%         Vsc(i,j)   = Vsc_out;

        it = it+1;
        if mod(it,50) == 0
            msg  = sprintf('p_z = %.2f, p_c = %.2f',p_z,p_c);
            disp(msg);
            disp(' ');
        end     

    end
end

save('equib_exist_demand.mat','flag_save','Rsz');
%% Plotting
[zz,cc] = meshgrid(p_z_grid,p_c_grid);

fig(1) = figure;
box on
grid on
hold on
for i = 1:length(p_z_grid)
    disp(i)
    for j = 1:length(p_c_grid)
%         disp(j)
        if flag_save(i,j) > 0
            scatter(zz(i,j),cc(i,j),20,'b','fill');
%         else
%             scatter(cc(i,j),zz(i,j),10,'r','fill');
        end
    end 
end
set(gca,'XLim',[0 1],'YLim',[0 1],'FontSize',15)
xlabel('p_C','FontSize',15)
ylabel('p_N','FontSize',15)



set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','equib_exist_demand.eps');

