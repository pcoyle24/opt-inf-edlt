%% About File
%--------------------------------------------------------------------------
% filename     : FisherRelation.m
% author       : Philip Coyle
% date created : 08/28/2018
% cd
% /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/Codes/StylizedModel/SS_Opt_Inf/FisherRelation
% FisherRelation.m
%--------------------------------------------------------------------------

%% Main Code

clear all
close all
clc

x = 0:1:40;
y1 = 0.5.*x - 3;
% y2 = 0.6.*x - 1.5;
y2 = 0.7.*x - 6;
y3(1,(1:length(x))) = 0;
y4(1,(1:length(x))) = 0;

for i = 1:length(x)
    if i <= 20
        y3(1,i) = 0;
        y4(1,i) = 0;
    elseif i <= 20
        y3(1,i) = 1.5*x(1,i) - 20;
        y4(1,i) = 0;
    else
        y3(1,i) = 1.5*(x(1,i) - 20);
        y4(1,i) = x(1,i)*15.125/14.25 - 22*15.125/14.25;       
    end		
end

fig1 = figure(1);
box off
hold on
grid off
h1 = plot(x,y2,'k',x,y3,'r','LineWidth',2);
set(gca,'XLim',[0 50],'XTick',[8.5 30],'XTickLabel',{'\beta','\Pi^{targ}' },'YLim',[-10 35],'YTick',[0],'YTickLabel',{'1'},'FontSize',22)
plot([30 30],[15 15],'.k','MarkerSize',25);
plot([8.5 8.5],[0 0],'.k','MarkerSize',25);
tr = text(25.2,18.1,{'Target', 'Steady State'});
dr = text(12,-3.5,{'Deflationary','Steady Steady'});
inf = text(52.5,-10,'\bf \Pi','Interpreter','tex');
R = text(-2,35,'\bf R','Interpreter','tex');
set(inf,'FontWeight','bold','HorizontalAlignment','center','FontSize',30)
set(R,'FontWeight','bold','HorizontalAlignment','center','FontSize',30)
set(tr,'FontWeight','bold','HorizontalAlignment','center','FontSize',22)
set(dr,'FontWeight','bold','HorizontalAlignment','center','FontSize',22)
c = [0.37,0.45];
d = [0.55,0.47];
arrow_2 = annotation('textarrow',c,d,'String',({'Fisher Relation'}));
set(arrow_2,'FontWeight','bold','HorizontalAlignment','center','FontSize',22)
e = [0.57,0.45];
f = [0.20,0.30];
arrow_3 = annotation('textarrow',e,f,'String',{'Taylor Rule'});
set(arrow_3,'Color','r','FontWeight','bold','HorizontalAlignment','center','FontSize',22)
annotation('arrow',[0.13 .92],[0.11 0.11])
annotation('arrow',[0.13 0.13],[0.11 0.94])


savedir = cd;
savedir = fullfile(savedir, '..');
if ispc 
    savedir = strcat(savedir,'\Final\');
else
    savedir = strcat(savedir,'/Final/');
end

set(fig1,'PaperOrientation','Landscape');
set(fig1,'PaperPosition',[0 0 11 8.5]);
print(fig1,'-depsc','FisherRelation.eps');
print(fig1,'-depsc',strcat(savedir,'FisherRelation.eps'));
