function residuals_ct_sun = ct_solve_sun(x_guess_ct_sun,params_s,params_ess)

global cBET cCHIc cCHIn cDELz cDELc p_s p_d

CT   = x_guess_ct_sun(1);

Vt = params_s(1);
Vd = params_s(2);

Cesst = params_ess(1);
Nesst = params_ess(2);
Cessd = params_ess(3);
Nessd = params_ess(4);

% Unconditional Probablity Calculation: Sunspot Shock only.
unc_prob_targ = (1-p_d)/(2-p_s-p_d);
unc_prob_def = (1-p_s)/(2-p_s-p_d);

Cess = unc_prob_targ*Cesst + unc_prob_def*Nesst;

if cCHIc == 1
    Us = log(Cesst + CT*Cess) - Nesst^(1+cCHIn)/(1+cCHIn);
    Ud = log(Cessd + CT*Cess) - Nessd^(1+cCHIn)/(1+cCHIn);
else
    Us = (Cesst + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nesst^(1+cCHIn)/(1+cCHIn);
    Ud = (Cessd + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nessd^(1+cCHIn)/(1+cCHIn);
end



A = [1-cBET*cDELz*p_s, -cBET*cDELz*(1-p_s);...
    -cBET*cDELz*(1-p_d), 1-cBET*cDELz*p_d];

b = [Us;Ud];
x = A\b;

Vesst = x(1);
Vessd = x(2);


Vess = unc_prob_targ*Vesst + unc_prob_def*Vessd;
V = unc_prob_targ*Vt + unc_prob_def*Vd;

residuals_ct_sun = V - Vess;

