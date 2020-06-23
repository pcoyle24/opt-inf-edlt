global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz p_s p_d
display((' '));
disp('Solving for Allocations with Sunspot Shock Only')

%Allocate Space for functions
Cs   = zeros(length(cALPHA),length(cPItarg));
PIs  = zeros(length(cALPHA),length(cPItarg));
Ys   = zeros(length(cALPHA),length(cPItarg));
Ns   = zeros(length(cALPHA),length(cPItarg));
Ws   = zeros(length(cALPHA),length(cPItarg));
Rs   = zeros(length(cALPHA),length(cPItarg));
Vs   = zeros(length(cALPHA),length(cPItarg));

Cd   = zeros(length(cALPHA),length(cPItarg));
PId  = zeros(length(cALPHA),length(cPItarg));
Yd   = zeros(length(cALPHA),length(cPItarg));
Nd   = zeros(length(cALPHA),length(cPItarg));
Wd   = zeros(length(cALPHA),length(cPItarg));
Rd   = zeros(length(cALPHA),length(cPItarg));
Vd   = zeros(length(cALPHA),length(cPItarg));

CT_sun   = zeros(length(cALPHA),length(cPItarg));

it_count = numel(CT_sun);

options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
tic;

it = 0;

for i = 1:length(cALPHA)
    for j = 1:length(cPItarg)
        x0 = ones(10,1)*0.9;
        func = @(x) sun_solve(x,i,j);

        x_out_sun = fsolve(func,x0,options);

        Cs_out   = x_out_sun(1);
        PIs_out  = x_out_sun(2);
        Ys_out   = x_out_sun(3);
        Ns_out   = x_out_sun(4);
        Vs_out   = x_out_sun(5);

        Ws_out = Ns_out^cCHIn*Cs_out^cCHIc;
        Rs_out = (cPItarg(j)/cBET)*((PIs_out/cPItarg(j))^cPHIpi)*((Ys_out/Ys_out)^cPHIy);

        Cd_out   = x_out_sun(6);
        PId_out  = x_out_sun(7);
        Yd_out   = x_out_sun(8);
        Nd_out   = x_out_sun(9);
        Vd_out   = x_out_sun(10);

        Wd_out = Nd_out^cCHIn*Cd_out^cCHIc;
        Rd_out = cRzlb;      
        
        % Consumption Transfer
        params_s = [Vs_out,Vd_out];
        params_ess = [Cesst,Nesst,Cessd,Nessd];
        x1 = 0;
        func = @(x) ct_solve_sun(x,params_s,params_ess);
        x_out_markov_sss = fsolve(func,x1,options);
        CT_out = x_out_markov_sss;


        Cs(i,j)   = Cs_out;
        PIs(i,j)  = PIs_out;
        Ys(i,j)   = Ys_out;
        Ns(i,j)   = Ns_out;
        Ws(i,j)   = Ws_out;
        Rs(i,j)   = Rs_out;
        Vs(i,j)   = Vs_out;

        Cd(i,j)   = Cd_out;
        PId(i,j)  = PId_out;
        Yd(i,j)   = Yd_out;
        Nd(i,j)   = Nd_out;
        Wd(i,j)   = Wd_out;
        Rd(i,j)   = Rd_out;
        Vd(i,j)   = Vd_out;
        
        CT_sun(i,j)   = CT_out;

        it = it+1;

        if mod(it,100) == 0
            display('Solving...');
        end
    end
end

Us = log(Cs) - Ns.^(1+cCHIn)./(1+cCHIn);
Ud = log(Cd) - Nd.^(1+cCHIn)/(1+cCHIn);
