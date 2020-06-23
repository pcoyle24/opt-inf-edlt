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
% dbstop if error

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
            func = @(x) markov_sss_solve(x,i,j);

            x_out_markov_sss = fsolve(func,x0,options);

            Csz_out   = x_out_markov_sss(1);
            PIsz_out  = x_out_markov_sss(2);
            Ysz_out   = x_out_markov_sss(3);
            Nsz_out   = x_out_markov_sss(4);
            Vsz_out   = x_out_markov_sss(5);

            Wsz_out = Nsz_out^cCHIn*Csz_out^cCHIc;
            Rsz_out = (cPItarg(j)/cBET)*((PIsz_out/cPItarg(j))^cPHIpi)*((Ysz_out/Ysz_out)^cPHIy);
            if Rsz_out < 1
                Rsz_out = cRzlb;
            end

            
            Csc_out   = x_out_markov_sss(6);
            PIsc_out  = x_out_markov_sss(7);
            Ysc_out   = x_out_markov_sss(8);
            Nsc_out   = x_out_markov_sss(9);
            Vsc_out   = x_out_markov_sss(10);

            Wsc_out = Nsc_out^cCHIn*Csc_out^cCHIc;
            Rsc_out_check = (cPItarg(j)/cBET)*((PIsc_out/cPItarg(j))^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);
            
            Rsc_out = (cPItarg(j)/cBET)*((PIsc_out/cPItarg(j))^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);
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

            if mod(it,50) == 0
                display(char(strcat('At cALPHA Grid Iteration = ',num2str(i), ' of ', num2str(length(cALPHA)))));
                display(char(strcat('At cPItarg Grid Iteration = ',num2str(j),' of ', num2str(length(cPItarg)))));
                display((' '));
            end
        end
    end
    toc

    save('markov_sss_data.mat','Csz','PIsz','Ysz','Nsz','Wsz','Rsz','Vsz','Csc','PIsc','Ysc','Nsc','Wsc','Rsc','Rsc_check','Vsc');
else
    load markov_sss_data.mat
end

%% Unconditional Probablity

unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

unc_welfare = unc_prob_zero*Vsz + unc_prob_crisis*Vsc;

%% Plotting SZ Case
[row, col] = size(Csz);
Xsz = zeros(row,col,5);
Xsz(:,:,1) = 100*(Csz-1);
Xsz(:,:,2) = 400*(PIsz-1);
Xsz(:,:,3) = 100*(Ysz-1);
Xsz(:,:,4) = 400*(Rsz-1);
Xsz(:,:,5) = Vsz;
Xsz(:,:,6) = unc_welfare;

v_max_sz = zeros(length(cALPHA),1);
v_max_sz_index = zeros(length(cALPHA),1);

for i = 1:length(cALPHA)
    v_max_sz(i) = max(Vsz(i,:));
    v_max_sz_index(i) = find(v_max_sz(i) == Vsz(i,:));
end

v_max_unc = zeros(length(cALPHA),1);
v_max_unc_index = zeros(length(cALPHA),1);

for i = 1:length(cALPHA)
    v_max_unc(i) = max(unc_welfare(i,:));
    v_max_unc_index(i) = find(v_max_unc(i) == unc_welfare(i,:));
end



header = {'Consumption','Inflation','Output','Nominal Interest Rate','Value','Unconditonal Welfare'};
colors = {'k','b','r','g'};
% alpha_string = {'\alpha = 0','\alpha = 0.33', '\alpha = 0.67','\alpha = 1'};
alpha_string = {'\alpha = 0.77','\alpha = 0.85', '\alpha = 0.9'};
y_label = {{'Percent Deveation from','Efficient Steady State'}, '', {'Percent Deveation from','Efficient Steady State'}, '', '', ''};


fig(1) = figure(1);
for i = 1:length(cALPHA)
    for j = 1:length(header)
        subplot(2,3,j)
        box on
        grid on
        hold on
        if j == 4
            h(i) = plot(400*(cPItarg-1),Xsz(i,:,j),colors{i},'LineWidth',2);
        else 
            plot(400*(cPItarg-1),Xsz(i,:,j),colors{i},'LineWidth',2)
        end
        set(gca,'XLim',[-2 4],'FontSize',15)
        if i == length(cALPHA)  
            xlabel('Infaltion Target','FontSize',15)
            ylabel(y_label{j},'FontSize',15)
            title(header{j},'FontSize',15)
        end
        if j == 5
            line([400*(cPItarg(v_max_sz_index(i))-1) 400*(cPItarg(v_max_sz_index(i))-1)],get(gca,'YLim'),'Color',colors{i},'LineStyle','--','LineWidth',1);
        end
        if j == 6
            line([400*(cPItarg(v_max_unc_index(i))-1) 400*(cPItarg(v_max_unc_index(i))-1)],get(gca,'YLim'),'Color',colors{i},'LineStyle','--','LineWidth',1);
        end        
    end
end
    
% L = legend([h(1) h(2) h(3) h(4) ],alpha_string{1},alpha_string{2},alpha_string{3},alpha_string{4},'Location', 'NorthWest');
% L = legend([h(1) h(2) h(3) ],alpha_string{1},alpha_string{2},alpha_string{3},'Location', 'NorthWest');
L = legend(h(1),alpha_string{1},'Location', 'NorthWest');


set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
% print(fig(1),'-depsc','Markov_SSS_Inf_Targ_p_z.eps');

%% Plotting SC Case
[row, col] = size(Csc);
Xsc = zeros(row,col,5);
Xsc(:,:,1) = 100*(Csc-1);
Xsc(:,:,2) = 400*(PIsc-1);
Xsc(:,:,3) = 100*(Ysc-1);
Xsc(:,:,4) = 400*(Rsc-1);
Xsc(:,:,5) = Vsc;

v_max_sc = zeros(length(cALPHA),1);
v_max_sc_index = zeros(length(cALPHA),1);

for i = 1:length(cALPHA)
    v_max_sc(i) = max(Vsc(i,:));
    v_max_sc_index(i) = find(v_max_sc(i) == Vsc(i,:));
end

header = {'Consumption','Inflation','Output','Nominal Interest Rate','Value'};
colors = {'k','b','r','g'};
% alpha_string = {'\alpha = 0','\alpha = 0.33', '\alpha = 0.67','\alpha = 1'};
alpha_string = {'\alpha = 0.77','\alpha = 0.85', '\alpha = 0.9'};
y_label = {{'Percent Deveation from','Efficient Steady State'}, '', {'Percent Deveation from','Efficient Steady State'}, '', ''};


fig(2) = figure(2);
for i = 1:length(cALPHA)
    for j = 1:length(header)
        subplot(2,3,j)
        box on
        grid on
        hold on
        if j == 4
            h(i) = plot(400*(cPItarg-1),Xsc(i,:,j),colors{i},'LineWidth',2);
        else 
            plot(400*(cPItarg-1),Xsc(i,:,j),colors{i},'LineWidth',2)
        end
        if j == 2
            line([-2 4],[-1 -1],'Color',colors{i},'LineStyle','--','LineWidth',1);
        end
        if j == 3
            line([-2 4],[-7 -7],'Color',colors{i},'LineStyle','--','LineWidth',1);
        end
        set(gca,'XLim',[-2 4],'FontSize',15)
        if i == length(cALPHA)  
            xlabel('Infaltion Target','FontSize',15)
            ylabel(y_label{j},'FontSize',15)
            title(header{j},'FontSize',15)
        end
        if j == 5
            line([400*(cPItarg(v_max_sc_index(i))-1) 400*(cPItarg(v_max_sc_index(i))-1)],get(gca,'YLim'),'Color',colors{i},'LineStyle','--','LineWidth',1);
        end
    end
end
    
% L = legend([h(1) h(2) h(3) h(4) ],alpha_string{1},alpha_string{2},alpha_string{3},alpha_string{4},'Location', 'NorthWest');
% L = legend([h(1) h(2) h(3) ],alpha_string{1},alpha_string{2},alpha_string{3},'Location', 'NorthWest');
L = legend(h(1),alpha_string{1},'Location', 'NorthWest');


set(fig(2),'PaperOrientation','Landscape');
set(fig(2),'PaperPosition',[0 0 11 8.5]);
% print(fig(2),'-depsc','Markov_SSS_Inf_Targ_p_c.eps');