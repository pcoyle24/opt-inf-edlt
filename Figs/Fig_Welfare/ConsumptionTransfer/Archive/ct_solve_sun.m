function residuals_ct_sun = ct_solve_sun(x_guess_ct_sun,params_s,params_ess)

global cBET cCHIc cCHIn cDELz cDELc p_s p_d

CT   = x_guess_ct_sun(1);

Cs = params_s(1);
Ns = params_s(2);
Vs = params_s(3);
Cd = params_s(4);
Nd = params_s(5);
Vd = params_s(6);

Cess = params_ess(1);
Vess = params_ess(2);

% Unconditional Probablity Calculation: Sunspot Shock only.
unc_prob_targ = (1-p_d)/(2-p_s-p_d);
unc_prob_def = (1-p_s)/(2-p_s-p_d);

if cCHIc == 1
    Us = log(Cs + CT*Cess) - Ns^(1+cCHIn)/(1+cCHIn);
    Ud = log(Cd + CT*Cess) - Nd^(1+cCHIn)/(1+cCHIn);
else
    Us = (Cs + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Ns^(1+cCHIn)/(1+cCHIn);
    Ud = (Cd + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nd^(1+cCHIn)/(1+cCHIn);
end



A = [1-cBET*cDELz*p_s, -cBET*cDELz*(1-p_s);...
    -cBET*cDELz*(1-p_d), 1-cBET*cDELz*p_d];

b = [Us;Ud];
x = A\b;

Vs = x(1);
Vd = x(2);


V = unc_prob_targ*Vs + unc_prob_def*Vd;

residuals_ct_sun = V - Vess;

