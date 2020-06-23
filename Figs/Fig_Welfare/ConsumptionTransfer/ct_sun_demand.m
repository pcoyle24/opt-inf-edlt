global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c

display((' '));
disp('Solving for Allocations with Sunspot Shock and Demand Shock')

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

CT_demand_sun   = zeros(length(cALPHA),length(cPItarg));

it_count = numel(CT_demand_sun);

options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
tic;

it = 0;

for i = 1:length(cALPHA)
    for j = 1:length(cPItarg)
        x0 = ones(20,1);
        func = @(x) sun_demand_solve(x,i,j);

        x_out_markov_sun = fsolve(func,x0,options);

        Csz_out   = x_out_markov_sun(1);
        PIsz_out  = x_out_markov_sun(2);
        Ysz_out   = x_out_markov_sun(3);
        Nsz_out   = x_out_markov_sun(4);
        Vsz_out   = x_out_markov_sun(5);

        Wsz_out = Nsz_out^cCHIn*Csz_out^cCHIc;
        Rsz_out = (cPItarg(j)/cBET)*((PIsz_out/cPItarg(j))^cPHIpi)*((Ysz_out/Ysz_out)^cPHIy);
%         if Rsz_out < 1
%             Rsz_out = cRzlb;
%         end


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

        % Consumption Transfer
        params = [Csz_out,Nsz_out,Vsz_out,Csc_out,Nsc_out,Vsc_out,Cdz_out,Ndz_out,Vdz_out,Cdc_out,Ndc_out,Vdc_out];
        params_ess = [Cess_ds, Vess_ds];
        x1 = 0;
        func = @(x) ct_solve_sun_demand(x,params,params_ess);
        x_out_markov_sun_demand = fsolve(func,x1,options);
        CT_out = x_out_markov_sun_demand;
        
        
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
        
        CT_demand_sun(i,j)   = CT_out;

        it = it+1;

        if mod(it,100) == 0
            display('Solving...');
        end
    end
end