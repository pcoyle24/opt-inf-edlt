function residuals_ct_demand = ct_solve_demand(x_guess_ct_demand,params_sz,params_ess)

global cBET cCHIc cCHIn cDELz cDELc p_z p_c

CT   = x_guess_ct_demand(1);

Csz = params_sz(1);
Nsz = params_sz(2);
Vsz = params_sz(3);
Csc = params_sz(4);
Nsc = params_sz(5);
Vsc = params_sz(6);

Cess = params_ess(1);
Vess = params_ess(2);

% Unconditional Probablity Calculation: Fundamentals Shock only.
unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

if cCHIc == 1
    Usz = log(Csz + CT*Cess) - Nsz^(1+cCHIn)/(1+cCHIn);
    Usc = log(Csc + CT*Cess) - Nsc^(1+cCHIn)/(1+cCHIn);
else
    Usz = (Csz + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nsz^(1+cCHIn)/(1+cCHIn);
    Usc = (Csc + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nsc^(1+cCHIn)/(1+cCHIn);
end



A = [1-cBET*cDELz*p_z, -cBET*cDELz*(1-p_z);...
    -cBET*cDELc*(1-p_c), 1-cBET*cDELc*p_c];

b = [Usz;Usc];
x = A\b;

Vsz = x(1);
Vsc = x(2);


V = unc_prob_zero*Vsz + unc_prob_crisis*Vsc;

residuals_ct_demand = V - Vess;

