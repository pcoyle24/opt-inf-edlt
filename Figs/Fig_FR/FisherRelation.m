%% About File
%--------------------------------------------------------------------------
% filename     : FisherRelation.m
% author       : Philip Coyle
% date created : 08/28/2018
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Codes/StylizedModel/SS_Opt_Inf/FisherRelation
% FisherRelation.m
%--------------------------------------------------------------------------

clear all
close all
clc

%% Paramaters
cPHIpi = 2;
cRstar = 1/400;

%% Define Fisher Relation
% Set Grid Intervals
pi_m_low = -1.5/400;
pi_m_high = 0.5/400;
pts = 1001;
pi_m = linspace(pi_m_low,pi_m_high, pts)';

% Taylor Rule
tr = max(0,cRstar + cPHIpi*pi_m);
% 'Riskless' Fisher Relation
fr = pi_m + cRstar;

%% Get DSS
[~,inx] = sort(abs(tr - fr));
DSS = inx(1:2);
it = 2;
while abs(DSS(1) - DSS(2)) <= 10
    DSS(2) = inx(it + 1);
    it = it + 1;
end

%% Plotting
fig(1) = figure(1);
subplot(2,2,1)
grid off
box off
hold on

% Plot Lines
h(1) = plot(400*pi_m, 400*tr,'Color','k','LineWidth',2);
h(2) = plot(400*pi_m, 400*fr,'Color','r','LineWidth',2);

% Plot DSS
pdss1 = plot(400*pi_m(DSS(1)),400*tr(DSS(1)));
pdss2 = plot(400*pi_m(DSS(2)),400*tr(DSS(2))); 
set(pdss1,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8)
set(pdss2,'Marker','d','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8)

% Clean up x and y axis
set(gca,'XLim',[400*pi_m(1), 400*pi_m(end)],'XTick',[-400*cRstar 0],'XTickLabel',{'\beta','\Pi^{targ}' },'YLim',[-1 3],'YTick',[0],'YTickLabel',{'1'},'FontSize',20)
inf = text(400*pi_m(end),-1.2,'\bf \Pi','Interpreter','tex');
R = text(-1.6,2.9,'R','Interpreter','tex');
set(inf,'FontWeight','bold','HorizontalAlignment','center','FontSize',20)
set(R,'FontWeight','normal ','HorizontalAlignment','center','FontSize',20)
% annotation('arrow',[0.13 .92],[0.11 0.11])
% annotation('arrow',[0.13 0.13],[0.11 0.94])
annotation('arrow',[0.13 .49],[0.585 0.585])
annotation('arrow',[0.13 0.13],[0.585 0.94])

% Make legend
L = legend([h(1) h(2) pdss1 pdss2],'Taylor Rule','Fisher Relation','Target Steady State','Deflationary Steady State');
set(L,'Location','NorthWest','Fontsize',12)
set([pdss1 pdss2],'linestyle','none')
legend('boxoff')

%% Save Figure
savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','FisherRelation.eps');
print(fig(1),'-depsc',strcat(savedir,'FisherRelation.eps'));
