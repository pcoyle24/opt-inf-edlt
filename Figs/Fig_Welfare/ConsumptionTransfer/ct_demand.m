global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c

display((' '));
disp('Solving for Allocations with Demand Shock Only')

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

CT_demand   = zeros(length(cALPHA),length(cPItarg));

it_count = numel(CT_demand);

options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
tic;

it = 0;

for i = 1:length(cALPHA)
    for j = 1:length(cPItarg)
        x0 = ones(10,1);
        func = @(x) demaned_solve(x,i,j);

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
        Rsc_out_check = (cPItarg(j)/(cDELc*cBET))*((PIsc_out/cPItarg(j))^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);

        Rsc_out = (cPItarg(j)/(cDELc*cBET))*((PIsc_out/cPItarg(j))^cPHIpi)*((Ysc_out/Ysc_out)^cPHIy);
        if Rsc_out < 1
            Rsc_out = cRzlb;
        end
        
        % Consumption Transfer
        params_sz = [Csz_out,Nsz_out,Vsz_out,Csc_out,Nsc_out,Vsc_out];
        params_ess = [Cess_d, Vess_d];
        x1 = 0;
        func = @(x) ct_solve_demand(x,params_sz,params_ess);
        x_out_markov_sss = fsolve(func,x1,options);
        CT_out = x_out_markov_sss;


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
        
        CT_demand(i,j)   = CT_out;

        it = it+1;

        if mod(it,100) == 0
            display('Solving...');
        end
    end
end


Usz = log(Csz) - Nsz.^(1+cCHIn)./(1+cCHIn);
Usc = log(Csc) - Nsc.^(1+cCHIn)/(1+cCHIn);
