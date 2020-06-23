function residuals_ct_demand = ct_solve_demand(x_guess_ct_demand,params_sz,params_ess)

global cBET cCHIc cCHIn cDELz cDELc p_z p_c

CT   = x_guess_ct_demand(1);

Vz = params_sz(1);
Vc = params_sz(2);

Cessz = params_ess(1);
Nessz = params_ess(2);
Cessc = params_ess(3);
Nessc = params_ess(4);

% Unconditional Probablity Calculation: Fundamentals Shock only.
unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

Cess = unc_prob_zero*Cessz + unc_prob_crisis*Cessc;

if cCHIc == 1
    Usz = log(Cessz + CT*Cess) - Nessz^(1+cCHIn)/(1+cCHIn);
    Usc = log(Cessc + CT*Cess) - Nessc^(1+cCHIn)/(1+cCHIn);
else
    Usz = (Cessz + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nessz^(1+cCHIn)/(1+cCHIn);
    Usc = (Cessc + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nessc^(1+cCHIn)/(1+cCHIn);
end



A = [1-cBET*cDELz*p_z, -cBET*cDELz*(1-p_z);...
    -cBET*cDELc*(1-p_c), 1-cBET*cDELc*p_c];

b = [Usz;Usc];
x = A\b;

Vessz = x(1);
Vessc = x(2);

Vess = unc_prob_zero*Vessz + unc_prob_crisis*Vessc;
V = unc_prob_zero*Vz + unc_prob_crisis*Vc;

residuals_ct_demand = V - Vess;

