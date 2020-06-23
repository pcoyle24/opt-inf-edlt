%--------------------------------------------------------------------------
%File Name: Markov_SSS_Inf_Targ.m
%Author: Philip Coyle
%Date Created: 02/02/2018
% cd % /mq/philiprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/SS_Opt_Inf/Model_2\Ex_1'
% Markov_SSS_Inf_Targ.m
%--------------------------------------------------------------------------

clear all 
close all
clc
dbstop if error

%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c
cBET        = 1/(1.0025);
cCHIc       = 1;
cCHIn       = 1;
cTHETA      = 11;%6;
cTAU        = 1/cTHETA;
cVARPHI     = 2550;%1850;
cPHIpi      = 4;%1.5;
cPHIy       = 0;
cRzlb       = 1;
cIOTA       = 1;
% cALPHA      = [0.765, 0.77, 0.775];
cALPHA      = 0.96;%0.835;
cPItarg_min_ann = -2; 
cPItarg_max_ann = 4;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.05 + 1;
cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n); 
cDELz       = 1;
c           = 0.20657/10;%0.198/10;
cDELc       = 1 + c;
p_z         = 0.995;  % 0.99  
p_c         = 0.75;   % 0.75, 0.8, 5/6, 
run_first_time = 1;

%% Main Code
if run_first_time == 1
    %Allocate Space for functions
    Cdz   = zeros(length(cALPHA),length(cPItarg));
    PIdz  = zeros(length(cALPHA),length(cPItarg));
    Ydz   = zeros(length(cALPHA),length(cPItarg));
    Ndz   = zeros(length(cALPHA),length(cPItarg));
    Wdz   = zeros(length(cALPHA),length(cPItarg));
    Rdz   = zeros(length(cALPHA),length(cPItarg));
    Vdz   = zeros(length(cALPHA),length(cPItarg));

    Cdc   = zeros(length(cALPHA),length(cPItarg));
    PIdc  = zeros(length(cALPHA),length(cPItarg));
    Ydc   = zeros(length(cALPHA),length(cPItarg));
    Ndc   = zeros(length(cALPHA),length(cPItarg));
    Wdc   = zeros(length(cALPHA),length(cPItarg));
    Rdc   = zeros(length(cALPHA),length(cPItarg));
    Rdc_check   = zeros(length(cALPHA),length(cPItarg));
    Vdc   = zeros(length(cALPHA),length(cPItarg));

    it_count = numel(Vdc);

    options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
    tic;

    it = 0;

    for i = 1:length(cALPHA)
        for j = 1:length(cPItarg)
            x0 = ones(10,1);
            func = @(x) markov_dss_solve(x,i,j);

            x_out_markov_dss = fsolve(func,x0,options);

            Cdz_out   = x_out_markov_dss(1);
            PIdz_out  = x_out_markov_dss(2);
            Ydz_out   = x_out_markov_dss(3);
            Ndz_out   = x_out_markov_dss(4);
            Vdz_out   = x_out_markov_dss(5);

            Wdz_out = Ndz_out^cCHIn*Cdz_out^cCHIc;
            Rdz_out = cRzlb;
%             Rdz_out = (cPItarg(j)/cBET)*((PIdz_out/cPItarg(j))^cPHIpi)*((Ydz_out/Ydz_out)^cPHIy);
%             if Rdz_out < 1
%                 Rdz_out = cRzlb;
%             end

            
            Cdc_out   = x_out_markov_dss(6);
            PIdc_out  = x_out_markov_dss(7);
            Ydc_out   = x_out_markov_dss(8);
            Ndc_out   = x_out_markov_dss(9);
            Vdc_out   = x_out_markov_dss(10);

            Wdc_out = Ndc_out^cCHIn*Cdc_out^cCHIc;
            Rdc_out_check = (cPItarg(j)/cBET)*((PIdc_out/cPItarg(j))^cPHIpi)*((Ydc_out/Ydc_out)^cPHIy);
            
%             Rdc_out = (cPItarg(j)/cBET)*((PIdc_out/cPItarg(j))^cPHIpi)*((Ydc_out/Ydc_out)^cPHIy);
%             if Rdc_out < 1
                Rdc_out = cRzlb;
%             end


            Cdz(i,j)   = Cdz_out;
            PIdz(i,j)  = PIdz_out;
            Ydz(i,j)   = Ydz_out;
            Ndz(i,j)   = Ndz_out;
            Wdz(i,j)   = Wdz_out;
            Rdz(i,j)   = Rdz_out;
            Vdz(i,j)   = Vdz_out;

            Cdc(i,j)   = Cdc_out;
            PIdc(i,j)  = PIdc_out;
            Ydc(i,j)   = Ydc_out;
            Ndc(i,j)   = Ndc_out;
            Wdc(i,j)   = Wdc_out;
            Rdc(i,j)   = Rdc_out;
            Rdc_check(i,j)   = Rdc_out_check;
            Vdc(i,j)   = Vdc_out;

            it = it+1;

            if mod(it,50) == 0
                display(char(strcat('At cALPHA Grid Iteration = ',num2str(i), ' of ', num2str(length(cALPHA)))));
                display(char(strcat('At cPItarg Grid Iteration = ',num2str(j),' of ', num2str(length(cPItarg)))));
                display((' '));
            end
        end
    end
    toc

    save('markov_dss_data.mat','Cdz','PIdz','Ydz','Ndz','Wdz','Rdz','Vdz','Cdc','PIdc','Ydc','Ndc','Wdc','Rdc','Rdc_check','Vdc');
else
    load markov_dss_data.mat
end

%% Unconditional Probablity

unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

unc_welfare = unc_prob_zero*Vdz + unc_prob_crisis*Vdc;
%% Plotting DZ Case
[row, col] = size(Cdz);
Xdz = zeros(row,col,5);
Xdz(:,:,1) = 100*(Cdz-1);
Xdz(:,:,2) = 400*(PIdz-1);
Xdz(:,:,3) = 100*(Ydz-1);
Xdz(:,:,4) = 400*(Rdz-1);
Xdz(:,:,5) = Vdz;
Xdz(:,:,6) = unc_welfare;


v_max_dz = zeros(length(cALPHA),1);
v_max_dz_index = zeros(length(cALPHA),1);

for i = 1:length(cALPHA)
    v_max_dz(i) = max(Vdz(i,:));
    v_max_dz_index(i) = find(v_max_dz(i) == Vdz(i,:));
end

v_max_unc = zeros(length(cALPHA),1);
v_max_unc_index = zeros(length(cALPHA),1);

for i = 1:length(cALPHA)
    v_max_unc(i) = max(unc_welfare(i,:));
    v_max_unc_index(i) = find(v_max_unc(i) == unc_welfare(i,:));
end



header = {'Consumption (DZ)','Inflation (DZ)','Output (DZ)','Nominal Interest Rate (DZ)','Value (DZ)','Unconditonal Welfare (DZ)'};
colors = {'k','b','r','g'};
% alpha_string = {'\alpha = 0','\alpha = 0.33', '\alpha = 0.67','\alpha = 1'};
alpha_string = {'\alpha = 0','\alpha = 0.5', '\alpha = 1'};
y_label = {{'Percent Deveation from','Efficient Steady State'}, '', {'Percent Deveation from','Efficient Steady State'}, '', '', ''};

fig(1) = figure(1);
for i = 1:length(cALPHA)
    for j = 1:length(header)
        subplot(2,3,j)
        box on
        grid on
        hold on
        if j == 4
            h(i) = plot(400*(cPItarg-1),Xdz(i,:,j),colors{i},'LineWidth',2);
        else 
            plot(400*(cPItarg-1),Xdz(i,:,j),colors{i},'LineWidth',2)
        end
        set(gca,'XLim',[-2 4],'FontSize',15)
        if i == length(cALPHA)  
            xlabel('Infaltion Target','FontSize',15)
            ylabel(y_label{j},'FontSize',15)
            title(header{j},'FontSize',15)
        end
        if j == 5
            line([400*(cPItarg(v_max_dz_index(i))-1) 400*(cPItarg(v_max_dz_index(i))-1)],get(gca,'YLim'),'Color',colors{i},'LineStyle','--','LineWidth',1);
        end
        if j == 6
            line([400*(cPItarg(v_max_unc_index(i))-1) 400*(cPItarg(v_max_unc_index(i))-1)],get(gca,'YLim'),'Color',colors{i},'LineStyle','--','LineWidth',1);
        end  
    end
end
    
% L = legend([h(1) h(2) h(3) h(4) ],alpha_string{1},alpha_string{2},alpha_string{3},alpha_string{4},'Location', 'NorthWest');
% L = legend([h(1) h(2) h(3) ],alpha_string{1},alpha_string{2},alpha_string{3},'Location', 'NorthWest');


set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','Markov_DSS_Inf_Targ_p_z.eps');

%% Plotting DC Case
[row, col] = size(Cdc);
Xdc = zeros(row,col,5);
Xdc(:,:,1) = 100*(Cdc-1);
Xdc(:,:,2) = 400*(PIdc-1);
Xdc(:,:,3) = 100*(Ydc-1);
Xdc(:,:,4) = 400*(Rdc-1);
Xdc(:,:,5) = Vdc;

v_max_dc = zeros(length(cALPHA),1);
v_max_dc_index = zeros(length(cALPHA),1);

for i = 1:length(cALPHA)
    v_max_dc(i) = max(Vdc(i,:));
    v_max_dc_index(i) = find(v_max_dc(i) == Vdc(i,:));
end

header = {'Consumption (DC)','Inflation (DC)','Output (DC)','Nominal Interest Rate (DC)','Value (DC)'};
colors = {'k','b','r','g'};
% alpha_string = {'\alpha = 0','\alpha = 0.33', '\alpha = 0.67','\alpha = 1'};
alpha_string = {'\alpha = 0','\alpha = 0.5', '\alpha = 1'};
y_label = {{'Percent Deveation from','Efficient Steady State'}, '', {'Percent Deveation from','Efficient Steady State'}, '', ''};


fig(2) = figure(2);
for i = 1:length(cALPHA)
    for j = 1:length(header)
        subplot(2,3,j)
        box on
        grid on
        hold on
        if j == 4
            h(i) = plot(400*(cPItarg-1),Xdc(i,:,j),colors{i},'LineWidth',2);
        else 
            plot(400*(cPItarg-1),Xdc(i,:,j),colors{i},'LineWidth',2)
        end
        set(gca,'XLim',[-2 4],'FontSize',15)
        if i == length(cALPHA)  
            xlabel('Infaltion Target','FontSize',15)
            ylabel(y_label{j},'FontSize',15)
            title(header{j},'FontSize',15)
        end
        if j == 5
            line([400*(cPItarg(v_max_dc_index(i))-1) 400*(cPItarg(v_max_dc_index(i))-1)],get(gca,'YLim'),'Color',colors{i},'LineStyle','--','LineWidth',1);
        end
    end
end
    
% L = legend([h(1) h(2) h(3) h(4) ],alpha_string{1},alpha_string{2},alpha_string{3},alpha_string{4},'Location', 'NorthWest');
% L = legend([h(1) h(2) h(3) ],alpha_string{1},alpha_string{2},alpha_string{3},'Location', 'NorthWest');


set(fig(2),'PaperOrientation','Landscape');
set(fig(2),'PaperPosition',[0 0 11 8.5]);
% print(fig(2),'-depsc','Markov_DSS_Inf_Targ_p_c.eps');