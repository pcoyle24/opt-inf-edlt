%--------------------------------------------------------------------------
% File Name: mult_equib_exist_sun.m
% Author: Philip Coyle
% Date Created: 10/17/2018
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Eq_Existance
% mult_equib_exist_sun
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
p_s_grid         = (0.9:0.001:1);
p_d_grid         = (0.9:0.001:1);
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;


%% Main Code
%Allocate Space for functions
Rs   = zeros(length(p_s_grid),length(p_d_grid));
flag_save   = nan(length(p_s_grid),length(p_d_grid));
PId_mult_out = zeros(2,1);
it_count = numel(flag_save);

options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000,'TolFun', 1e-32, 'Display', 'none');
it = 0;
for i = 1:length(p_s_grid)
    p_s = p_s_grid(i);
    for j = 1:length(p_d_grid)
        p_d = p_d_grid(j);
        
        x0 = ones(8,1);
        func = @(x) sun_solve(x);
        [x_out_sun,~,flag] = fsolve(func,x0,options);
        flag_save(i,j) = flag;
        
        PIs_out  = x_out_sun(2);
        Ys_out   = x_out_sun(3);
        Rs_out = (cPItarg/cBET)*((PIs_out/cPItarg)^cPHIpi)*((Ys_out/Ys_out)^cPHIy);
        Rs(i,j)   = Rs_out;
        PId_out  = x_out_sun(6);
        
%         % Check if multiple equilbriums exist
%         if flag_save(i,j) > 0 
%             for c = 1:2
%                 if c == 1
%                     x1 = ones(8,1);
%                 else
% %                     if j == 1
%                         x1 = [-70/100 + 1; -20/400 + 1; 175/100 + 1 ;175/100 + 1;-75/100 + 1; -15/400 + 1; 150/100 + 1 ;150/100 + 1];
% %                     else
% %                         x1 = x_out_mult_sun2;
% %                     end
%                 end
%                 func = @(x) multeq_sun_solve(x);
%                 [x_out_mult_sun,~,flag] = fsolve(func,x1,options);
%                 if c == 2
%                     x_out_mult_sun2 = x_out_mult_sun;
%                 end
% 
%                 PId_mult_out(c)  = x_out_mult_sun(6);
%                 
%                 if c == 2
%                     if abs(PId_mult_out(1) - PId_mult_out(2)) < 1e-10
%                          flag_save(i,j) = 10;
%                     end
%                 end
%             end
%         end

        it = it+1;
        if mod(it,50) == 0
            msg  = sprintf('p_s = %.3f, p_d = %.3f',p_s,p_d);
            disp(msg);
            disp(' ');
        end
    end
end
save('mult_equib_exist_sun.mat','flag_save','Rs');
% %% Plotting
% [ss,dd] = meshgrid(p_s_grid,p_d_grid);
% 
% fig(1) = figure;
% box on
% grid on
% hold on
% for i = 1:length(p_s_grid)
%     disp(i);
%     for j = 1:length(p_d_grid)
%         if Rs(i,j) > 1
%             if flag_save(i,j) > 0
%                 scatter(100*ss(i,j),100*dd(i,j),20,'r','fill');
% %                 scatter(100*ss(j),100*dd(j),20,'r','fill');
%                 if flag_save(i,j) == 10
%                     scatter(100*ss(i,j),100*dd(i,j),20,'b','fill');
% %                     scatter(100*ss(j),100*dd(j),20,'b','fill');
%                 end
%             end
%         end
%     end 
% end
% set(gca,'XLim',[90 100],'YLim',[90 100],'FontSize',15)
% xlabel('100p_D','FontSize',15)
% ylabel('100p_S','FontSize',15)
% 
% set(fig(1),'PaperOrientation','Landscape');
% set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','mult_equib_exist_sun.eps')