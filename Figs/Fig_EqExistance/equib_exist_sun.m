%--------------------------------------------------------------------------
%File Name: equib_exist_sun.m
%Author: Philip Coyle
%Date Created: 10/17/2018
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/draft/Figs/Eq_Existance
% equib_exist_sun
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
cPItarg_ann = 2; 
cPItarg     = cPItarg_ann/400 + 1;
p_s_grid         = (0:0.01:1);
p_d_grid         = (0:0.01:1);
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;


%% Main Code
%Allocate Space for functions
Cs   = zeros(length(p_s_grid),length(p_d_grid));
PIs  = zeros(length(p_s_grid),length(p_d_grid));
Ys   = zeros(length(p_s_grid),length(p_d_grid));
Ns   = zeros(length(p_s_grid),length(p_d_grid));
Ws   = zeros(length(p_s_grid),length(p_d_grid));
Rs   = zeros(length(p_s_grid),length(p_d_grid));

Cd   = zeros(length(p_s_grid),length(p_d_grid));
PId  = zeros(length(p_s_grid),length(p_d_grid));
Yd   = zeros(length(p_s_grid),length(p_d_grid));
Nd   = zeros(length(p_s_grid),length(p_d_grid));
Wd   = zeros(length(p_s_grid),length(p_d_grid));
Rd   = zeros(length(p_s_grid),length(p_d_grid));
flag_save   = nan(length(p_s_grid),length(p_d_grid));

it_count = numel(Cd);

options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000,'TolFun', 1e-32, 'Display', 'none');
tic;

it = 0;

for i = 1:length(p_s_grid)
    p_s = p_s_grid(i);
    for j = 1:length(p_d_grid)
        p_d = p_d_grid(j);
        
%         x0 = ones(10,1);
        x0 = ones(8,1);
        func = @(x) sun_solve(x);

        [x_out_sun,~,flag] = fsolve(func,x0,options);
        flag_save(i,j) = flag;
        
        Cs_out   = x_out_sun(1);
        PIs_out  = x_out_sun(2);
        Ys_out   = x_out_sun(3);
        Ns_out   = x_out_sun(4);
%         Vs_out   = x_out_sun(5);

        Ws_out = Ns_out^cCHIn*Cs_out^cCHIc;
        Rs_out = (cPItarg/cBET)*((PIs_out/cPItarg)^cPHIpi)*((Ys_out/Ys_out)^cPHIy);

        Cd_out   = x_out_sun(5);
        PId_out  = x_out_sun(6);
        Yd_out   = x_out_sun(7);
        Nd_out   = x_out_sun(8);
%         Vd_out   = x_out_sun(10);

        Wd_out = Nd_out^cCHIn*Cd_out^cCHIc;
        Rd_out = cRzlb;              


        Cs(i,j)   = Cs_out;
        PIs(i,j)  = PIs_out;
        Ys(i,j)   = Ys_out;
        Ns(i,j)   = Ns_out;
        Ws(i,j)   = Ws_out;
        Rs(i,j)   = Rs_out;

        Cd(i,j)   = Cd_out;
        PId(i,j)  = PId_out;
        Yd(i,j)   = Yd_out;
        Nd(i,j)   = Nd_out;
        Wd(i,j)   = Wd_out;
        Rd(i,j)   = Rd_out;

        it = it+1;
        if mod(it,50) == 0
            msg  = sprintf('p_s = %.2f, p_d = %.2f',p_s,p_d);
            disp(msg);
            disp(' ');
        end
    end
end
save('equib_exist_sun.mat','flag_save','Rs');
%% Plotting
[ss,dd] = meshgrid(p_s_grid,p_d_grid);

flag_plot = nan(size(flag_save));

fig(1) = figure;
box on
grid on
hold on
for i = 1:length(p_s_grid)
    disp(i);
    for j = 1:length(p_d_grid)
        if flag_save(i,j) > 0
            if Rs(i,j) > 1
                if i > 80 && j > 80
                    scatter(ss(i,j),dd(i,j),20,'b','fill');
                end
%                 flag_plot(i,j) = 1;
%             else
%                 flag_plot(i,j) = 0;
            end
%         else
%             flag_plot(i,j) = 0;
        end
    end 
end
set(gca,'XLim',[0 1],'YLim',[0 1],'FontSize',15)
xlabel('p_D','FontSize',15)
ylabel('p_S','FontSize',15)

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','equib_exist_sun.eps')