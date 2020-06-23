%--------------------------------------------------------------------------
%File Name: equib_exist_plot.m
%Author: Philip Coyle
%Date Created: 10/17/2018
% cd /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Eq_Existance
% equib_exist_plot
%--------------------------------------------------------------------------

clear all 
close all
clc

p_z_grid         = 100*(0:0.01:1);
p_c_grid         = 100*(0:0.01:1);
p_s_grid         = 100*(0:0.01:1);
p_d_grid         = 100*(0:0.01:1);

load('equib_exist_demand.mat');
[zz,cc] = meshgrid(p_z_grid,p_c_grid);

fig(1) = figure(1);
% subplot(2,2,1)
subplot(1,2,1)
box on
grid on
hold on
axis square
for i = 1:length(p_z_grid)
    disp(i)
    for j = 1:length(p_c_grid)
        if flag_save(i,j) > 0
            if Rsz(i,j) >= 1
                scatter(zz(i,j),cc(i,j),5,'b','fill');
            end
        end
    end 
end
% set(gca,'XLim',[0 1],'XTick',(0:0.2:1),'YLim',[0 1],'YTick',(0:0.2:1),'FontSize',15)
set(gca,'XLim',[0 100],'XTick',(0:20:100),'YLim',[0 100],'YTick',(0:20:100),'FontSize',15)
xlabel('100p_C','FontSize',15)
ylabel('100p_N','FontSize',15)
title('Existence under Demand Shock Only','FontSize',15,'FontWeight','normal');

load('equib_exist_sun.mat','flag_save','Rs');
[ss,dd] = meshgrid(p_s_grid,p_d_grid);

% subplot(2,2,2)
subplot(1,2,2)
box on
grid on
hold on
axis square
for i = 1:length(p_s_grid)
    disp(i);
    for j = 1:length(p_d_grid)
        if flag_save(i,j) > 0
            if Rs(i,j) >= 1
                if i > 80 && j > 80
                    scatter(ss(i,j),dd(i,j),5,'b','fill');
                end
            end
        end
    end 
end
% set(gca,'XLim',[0 1],'XTick',(0:0.2:1),'YLim',[0 1],'YTick',(0:0.2:1),'FontSize',15)
set(gca,'XLim',[0 100],'XTick',(0:20:100),'YLim',[0 100],'YTick',(0:20:100),'FontSize',15)
xlabel('100p_D','FontSize',15)
ylabel('100p_T','FontSize',15)
title('Existence under Sunspot Shock Only','FontSize',15,'FontWeight','normal');

savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','equib_exist_sun_demand.eps');
print(fig(1),'-dpdf','equib_exist_sun_demand.pdf');
print(fig(1),'-dpdf',strcat(savedir,'equib_exist_sun_demand.pdf'));



% The code struggles to crop out the white-sapce. So, to manaully do this
% go to the appropriate path in linux and run the code: 
% pdfcrop --margins '-45 -140 -60 -135'  equib_exist_sun_demand.pdf equib_exist_sun_demand_crop.pdf
% to get a cropped pdf. Run:
% pdftops equib_exist_sun_demand_crop.pdf equib_exist_sun_demand.eps

