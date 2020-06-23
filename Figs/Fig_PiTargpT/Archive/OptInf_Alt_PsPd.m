%--------------------------------------------------------------------------
%File Name: OptInf_Alt_PsPd.m
%Author: Philip Coyle
%Date Created: 01/17/2018
% cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/OptInf/draft/Figs/Fig_6/Archive
% matlab -nodesktop -nosplash -r OptInf_Alt_PsPd
%--------------------------------------------------------------------------

clear all 
close all
clc
% dbstop if error

if ispc
    addpath ..\Calib\
    addpath ..\min_inf_targ\
    
    load eqm_exist_sun_demand.mat
    min_inf_targ_sun_demand = min_inf_targ;
else
    addpath ../../Calib/
    addpath ../../min_inf_targ/
    load eqm_exist_sun_demand.mat
    min_inf_targ_sun_demand = min_inf_targ;
end

load('params.mat')
%% Parameters
global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg p_s p_d cDELz cDELc p_z p_c
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
cPItarg_min_ann = min_inf_targ_sun_demand; 
cPItarg_max_ann = 2;
cPItarg_min = cPItarg_min_ann/400 + 1;
cPItarg_max = cPItarg_max_ann/400 + 1;
cPItarg_n   = (cPItarg_max_ann - cPItarg_min_ann)/0.01 + 1;
cPItarg     = linspace(cPItarg_min,cPItarg_max,cPItarg_n); 
p_s         = (0.99:0.0001:1); 
p_d         = [0.935,0.95,0.965,0.975]; 
cDELz       = 1;
c           = params(2);
cDELc       = 1 + c;
p_z         = 0.995;  
p_c         = 0.75;

%% Generate or Load Data
run_first_time = 1;

if run_first_time == 1
    %Allocate Space for functions
    Csz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    PIsz  = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Ysz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Nsz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Wsz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Rsz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Vsz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));

    Csc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    PIsc  = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Ysc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Nsc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Wsc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Rsc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Rsc_check   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Vsc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));

    Cdz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    PIdz  = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Ydz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Ndz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Wdz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Rdz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Vdz   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));

    Cdc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    PIdc  = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Ydc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Ndc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Wdc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Rdc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Rdc_check   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));
    Vdc   = zeros(length(cALPHA),length(cPItarg),length(p_s), length(p_d));

    it_count = numel(Vsc);

    options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000,'TolFun', 1e-32, 'Display', 'none');
    tic;

    it = 0;

    for i = 1:length(cALPHA)
        for j = 1:length(cPItarg)
            for k = 1:length(p_s)
                for q = 1:length(p_d)

                    if j == 1
                        x0 = ones(20,1);
                    else 
                        x0 = x_out_markov_sun;
                    end
                    func = @(x) sun_demand_solve(x,i,j,k,q);

                    x_out_markov_sun = fsolve(func,x0,options);

                    Csz_out   = x_out_markov_sun(1);
                    PIsz_out  = x_out_markov_sun(2);
                    Ysz_out   = x_out_markov_sun(3);
                    Nsz_out   = x_out_markov_sun(4);
                    Vsz_out   = x_out_markov_sun(5);

                    Wsz_out = Nsz_out^cCHIn*Csz_out^cCHIc;
                    Rsz_out = (cPItarg(j)/cBET)*((PIsz_out/cPItarg(j))^cPHIpi)*((Ysz_out/Ysz_out)^cPHIy);
                    if Rsz_out < 1
                        Rsz_out = cRzlb;
                    end


                    Csc_out   = x_out_markov_sun(6);
                    PIsc_out  = x_out_markov_sun(7);
                    Ysc_out   = x_out_markov_sun(8);
                    Nsc_out   = x_out_markov_sun(9);
                    Vsc_out   = x_out_markov_sun(10);

                    Wsc_out = Nsc_out^cCHIn*Csc_out^cCHIc;
                    Rsc_out_check = (cPItarg(j)/(cDELc*cBET))*((PIsc_out/cPItarg(j))^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);

                    Rsc_out = (cPItarg(j)/(cDELc*cBET))*((PIsc_out/cPItarg(j))^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);
                    if Rsc_out < 1
                        Rsc_out = cRzlb;
                    end

                    Cdz_out   = x_out_markov_sun(11);
                    PIdz_out  = x_out_markov_sun(12);
                    Ydz_out   = x_out_markov_sun(13);
                    Ndz_out   = x_out_markov_sun(14);
                    Vdz_out   = x_out_markov_sun(15);

                    Wdz_out = Ndz_out^cCHIn*Cdz_out^cCHIc;
                    Rdz_out = cRzlb;
        %             Rdz_out = (cPItarg(j)/cBET)*((PIdz_out/cPItarg(j))^cPHIpi)*((Ydz_out/Ydz_out)^cPHIy);
        %             if Rdz_out < 1
        %                 Rdz_out = cRzlb;
        %             end


                    Cdc_out   = x_out_markov_sun(16);
                    PIdc_out  = x_out_markov_sun(17);
                    Ydc_out   = x_out_markov_sun(18);
                    Ndc_out   = x_out_markov_sun(19);
                    Vdc_out   = x_out_markov_sun(20);

                    Wdc_out = Ndc_out^cCHIn*Cdc_out^cCHIc;
                    Rdc_out_check = (cPItarg(j)/(cDELc*cBET))*((PIdc_out/cPItarg(j))^cPHIpi)*((Ydc_out/Ydc_out)^cPHIy);

                    Rdc_out = (cPItarg(j)/(cDELc*cBET))*((PIdc_out/cPItarg(j))^cPHIpi)*((Ydc_out/Ydc_out)^cPHIy);
                    if Rdc_out < 1
                        Rdc_out = cRzlb;
                    end

                    Csz(i,j,k,q)   = Csz_out;
                    PIsz(i,j,k,q)  = PIsz_out;
                    Ysz(i,j,k,q)   = Ysz_out;
                    Nsz(i,j,k,q)   = Nsz_out;
                    Wsz(i,j,k,q)   = Wsz_out;
                    Rsz(i,j,k,q)   = Rsz_out;
                    Vsz(i,j,k,q)   = Vsz_out;

                    Csc(i,j,k,q)   = Csc_out;
                    PIsc(i,j,k,q)  = PIsc_out;
                    Ysc(i,j,k,q)   = Ysc_out;
                    Nsc(i,j,k,q)   = Nsc_out;
                    Wsc(i,j,k,q)   = Wsc_out;
                    Rsc(i,j,k,q)   = Rsc_out;
                    Rsc_check(i,j,k,q)   = Rsc_out_check;
                    Vsc(i,j,k,q)   = Vsc_out;



                    Cdz(i,j,k,q)   = Cdz_out;
                    PIdz(i,j,k,q)  = PIdz_out;
                    Ydz(i,j,k,q)   = Ydz_out;
                    Ndz(i,j,k,q)   = Ndz_out;
                    Wdz(i,j,k,q)   = Wdz_out;
                    Rdz(i,j,k,q)   = Rdz_out;
                    Vdz(i,j,k,q)   = Vdz_out;

                    Cdc(i,j,k,q)   = Cdc_out;
                    PIdc(i,j,k,q)  = PIdc_out;
                    Ydc(i,j,k,q)   = Ydc_out;
                    Ndc(i,j,k,q)   = Ndc_out;
                    Wdc(i,j,k,q)   = Wdc_out;
                    Rdc(i,j,k,q)   = Rdc_out;
                    Rdc_check(i,j,k,q)   = Rdc_out_check;
                    Vdc(i,j,k,q)   = Vdc_out;

                    it = it+1;

                    if mod(it,50) == 0
                        display(char(strcat('At cALPHA Grid Iteration = ',num2str(i), ' of ', num2str(length(cALPHA)))));
                        display(char(strcat('At cPItarg Grid Iteration = ',num2str(j),' of ', num2str(length(cPItarg)))));
                        display(char(strcat('At p_s Grid Iteration = ',num2str(k),' of ', num2str(length(p_s)))));
                        display(char(strcat('At p_d Grid Iteration = ',num2str(q),' of ', num2str(length(p_d)))));
                        display((' '));
                    end
                end
            end
        end
    end
    save('alt_PsPd_data.mat','Csz','PIsz','Ysz','Nsz','Wsz','Rsz','Vsz','Csc','PIsc','Ysc','Nsc','Wsc','Rsc','Rsc_check','Vsc','Cdz','PIdz','Ydz','Ndz','Wdz','Rdz','Vdz','Cdc','PIdc','Ydc','Ndc','Wdc','Rdc','Rdc_check','Vdc');
else
    display('Loading in already solved for allocations')
    load('alt_PsPd_data.mat')
end
%% Main Code
gamma = nan(4,16^2);
unc_prob = nan(4,16^2);

opt_inf_sss = get_opt_inf_baseline;
opt_inf_sss = ones(length(p_s),length(p_d))*opt_inf_sss;

%Preallocate space for optimal inflation target
opt_inf_unc = nan(length(p_s),length(p_d));
opt_inf_unc_1 = nan(length(p_s),length(p_d));
it = 1;
for i = 1:length(p_s)
    for j = 1:length(p_d)
        % We first calculate the unconditional, or steady state,
        % probabilities of being in both regimes.
        if (p_s(i) == 1 && p_d(j) == 1) || p_d(j) == 1
            % p_s = 1 and p_d = 1 is an edge case we cannot consider, 
            % because they are both absorbing states and unconditional 
            % probability is undefined
%             opt_inf_unc(i,j) = nan;
            opt_inf_unc_1(i,j) = nan;
        else
%             % V1 (Tai's method)
%             gamma_z = (1-p_c)/(2-p_z-p_c);
%             gamma_c = (1-p_z)/(2-p_z-p_c);
%             gamma_s = (1-p_d(j))/(2-p_s(i)-p_d(j));
%             gamma_d = (1-p_s(i))/(2-p_s(i)-p_d(j));
%             
%             gamma_sz = gamma_s*gamma_z;
%             gamma_sc = gamma_s*gamma_c;
%             gamma_dz = gamma_d*gamma_z;
%             gamma_dc = gamma_d*gamma_c;
%             
%             gamma_vec = [gamma_sz; gamma_sc; gamma_dz; gamma_dc];
%             gamma(:,it) = gamma_vec;
                      
            
            % V2 (code method)
            StateMat = [p_s(i), 1-p_s(i);
                        1-p_d(j), p_d(j)];

            ShockMat = [p_z, 1-p_z;
                        1-p_c, p_c];

            TransMat = kron(StateMat,ShockMat);

            unc_prob_vec = limitdist(TransMat);

            unc_prob_sz = unc_prob_vec(1);
            unc_prob_sc = unc_prob_vec(2);
            unc_prob_dz = unc_prob_vec(3);
            unc_prob_dc = unc_prob_vec(4);
            
            unc_prob(:,it) = unc_prob_vec;
            
            V_sz = Vsz(1,:,i,j);
            V_sc = Vsc(1,:,i,j);
            V_dz = Vdz(1,:,i,j);
            V_dc = Vdc(1,:,i,j);

            unc_welf = unc_prob_sz*V_sz + unc_prob_sc*V_sc + unc_prob_dz*V_dz + unc_prob_dc*V_dc;
%             unc_welf_1 = gamma_sz*V_sz + gamma_sc*V_sc + gamma_dz*V_dz + gamma_dc*V_dc;
            
            

            %Determine Optimal Inflation Target
            unc_welf_max = max(unc_welf);
            unc_welf_max_index = find(unc_welf_max == unc_welf);
            opt_inf_unc(i,j) = 400*(cPItarg(unc_welf_max_index)-1); 

%             %Determine Optimal Inflation Target
%             unc_welf_max_1 = max(unc_welf_1);
%             unc_welf_max_index_1 = find(unc_welf_max_1 == unc_welf_1);
%             opt_inf_unc_1(i,j) = 400*(cPItarg(unc_welf_max_index_1)-1);
            
            it = it + 1;
        end
    end
end

% Difference between the optimal inflation target in model with crisis
% shock only and model with both crisis shock and sunspot shock.
opt_inf_net = opt_inf_unc - opt_inf_sss;
% toc

%% Plotting
line_vals = [0.995 0.998 0.999];

colors = {'k:','b--','r-.','k'};
pd_string = {'p_D = 0.935','p_D = 0.95','p_D = 0.965','p_D = 0.975'};

fig(1) = figure(1);
box on
grid on
hold on
for i = 1:length(p_d)
    h(i) = plot(p_s,opt_inf_unc(:,i),colors{i},'LineWidth',2);
    xlabel('p_T','FontSize',20)
    ylabel('Optimal Inflation','FontSize',20)
    line([line_vals(i) line_vals(i)],get(gca,'YLim'),'Color','k','LineWidth',1)
    set(gca,'XLim',[0.99 1],'Ylim',[min_inf_targ_sun_demand cPItarg_max_ann],'FontSize',20)
    
end
L = legend([h(1) h(2) h(3) h(4)],pd_string{1},pd_string{2},pd_string{3},pd_string{4},'Location', 'NorthWest');
set(L,'FontSize',20)



% savedir = cd;
% savedir = fullfile(savedir, '..');
% if ispc 
%     savedir = strcat(savedir,'\Final\');
% else
%     savedir = strcat(savedir,'/Final/');
% end


set(fig(1),'PaperOrientation','Landscape');
set(fig(1),'PaperPosition',[0 0 11 8.5]);
print(fig(1),'-depsc','opt_inf_vary_ps_pd.eps');
% print(fig(1),'-depsc',strcat(savedir,'opt_inf_vary_ps_pd.eps'));

