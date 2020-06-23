global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c

display((' '));
disp('Solving for Deflationary Regime Values')

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
            display('Solving...');
        end
    end
end