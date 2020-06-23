global cBET cCHIc cCHIn cPHIpi cPHIy cRzlb cDELz cDELc p_z p_c p_s p_d

cPItarg_ess = 0/400 + 1;
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

%% Demand Shock Only
disp('Solving for Efficient Steady State Allocations Demand Shock Only')

x0 = ones(10,1);
func = @(x) ess_demand_solve(x);
x_out_ess = fsolve(func,x0,options);

Cessz   = x_out_ess(1);
PIessz  = x_out_ess(2);
Yessz   = x_out_ess(3);
Nessz   = x_out_ess(4);
Vessz   = x_out_ess(5);
Wessz = Nessz^cCHIn*Cessz^cCHIc;
Ressz = (cPItarg_ess/(cDELz*cBET))*((PIessz/cPItarg_ess)^cPHIpi)*((Yessz/Yessz)^cPHIy);

Cessc   = x_out_ess(6);
PIessc  = x_out_ess(7);
Yessc   = x_out_ess(8);
Nessc   = x_out_ess(9);
Vessc   = x_out_ess(10);
Wessc = Nessc^cCHIn*Cessc^cCHIc;
Ressc = (cPItarg_ess/cDELc*cBET)*((PIessc/cPItarg_ess)^cPHIpi)*((Yessc/Yessc)^cPHIy);


%% Sunspot Shock Only
disp('Solving for Efficient Steady State Allocations Sunspot Shock Only')

x1 = ones(10,1);
func = @(x) ess_sun_solve(x);
x_out_ess = fsolve(func,x1,options);

Cesst   = x_out_ess(1);
PIesst  = x_out_ess(2);
Yesst   = x_out_ess(3);
Nesst   = x_out_ess(4);
Vesst   = x_out_ess(5);
Wesst = Nesst^cCHIn*Cesst^cCHIc;
Resst = (cPItarg_ess/(cDELz*cBET))*((PIesst/cPItarg_ess)^cPHIpi)*((Yesst/Yesst)^cPHIy);

Cessd   = x_out_ess(6);
PIessd  = x_out_ess(7);
Yessd   = x_out_ess(8);
Nessd   = x_out_ess(9);
Vessd   = x_out_ess(10);
Wessd = Nessd^cCHIn*Cessd^cCHIc;
Ressd = (cPItarg_ess/(cDELz*cBET))*((PIessd/cPItarg_ess)^cPHIpi)*((Yessd/Yessd)^cPHIy);
%% Demand and Sunspot Shock
disp('Solving for Efficient Steady State Allocations Demand & Sunspot Shock')

x2 = ones(20,1);
func = @(x) ess_sun_demand_solve(x);
x_out_ess = fsolve(func,x2,options);

Cesstz   = x_out_ess(1);
PIesstz  = x_out_ess(2);
Yesstz   = x_out_ess(3);
Nesstz   = x_out_ess(4);
Vesstz   = x_out_ess(5);
Wesstz = Nesstz^cCHIn*Cesstz^cCHIc;
Resstz = (cPItarg_ess/(cDELz*cBET))*((PIesstz/cPItarg_ess)^cPHIpi)*((Yesstz/Yesstz)^cPHIy);

Cesstc   = x_out_ess(6);
PIesstc  = x_out_ess(7);
Yesstc   = x_out_ess(8);
Nesstc   = x_out_ess(9);
Vesstc   = x_out_ess(10);
Wesstc = Nesstc^cCHIn*Cesstc^cCHIc;
Resstc = (cPItarg_ess/(cDELc*cBET))*((PIesstc/cPItarg_ess)^cPHIpi)*((Yesstc/Yesstc)^cPHIy);

Cessdz   = x_out_ess(11);
PIessdz  = x_out_ess(12);
Yessdz   = x_out_ess(13);
Nessdz   = x_out_ess(14);
Vessdz   = x_out_ess(15);
Wessdz = Nessdz^cCHIn*Cessdz^cCHIc;
Ressdz = (cPItarg_ess/(cDELz*cBET))*((PIessdz/cPItarg_ess)^cPHIpi)*((Yessdz/Yessdz)^cPHIy);

Cessdc   = x_out_ess(16);
PIessdc  = x_out_ess(17);
Yessdc   = x_out_ess(18);
Nessdc   = x_out_ess(19);
Vessdc   = x_out_ess(20);
Wessdc = Nessdc^cCHIn*Cessdc^cCHIc;
Ressdc = (cPItarg_ess/(cDELc*cBET))*((PIessdc/cPItarg_ess)^cPHIpi)*((Yessdc/Yessdc)^cPHIy);

%% Unconditional Probablity Calculation: Demand Shock
unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

Cess_d = unc_prob_zero.*Cessz + unc_prob_crisis.*Cessc;
Vess_d = unc_prob_zero.*Vessz + unc_prob_crisis.*Vessc;

%% Unconditional Probablity Calculation: Sunspot Shock
unc_prob_targ = (1-p_d)/(2-p_s-p_d);
unc_prob_def = (1-p_s)/(2-p_s-p_d);

Cess_s = unc_prob_targ.*Cesst + unc_prob_def.*Cessd;
Vess_s = unc_prob_targ.*Vesst + unc_prob_def.*Vessd;

%% Unconditional Probablity Calculation: Demand & Sunspot Shock
StateMat = [p_s, 1-p_s;
            1-p_d, p_d];
        
ShockMat = [p_z, 1-p_z;
            1-p_c, p_c];
        
TransMat = kron(StateMat,ShockMat);

unc_prob = limitdist(TransMat);

unc_prob_sz = unc_prob(1);
unc_prob_sc = unc_prob(2);
unc_prob_dz = unc_prob(3);
unc_prob_dc = unc_prob(4);

Cess_ds = unc_prob_sz*Cesstz + unc_prob_sc*Cesstc + unc_prob_dz*Cessdz + unc_prob_dc*Cessdc;
Vess_ds = unc_prob_sz*Vesstz + unc_prob_sc*Vesstc + unc_prob_dz*Vessdz + unc_prob_dc*Vessdc;

