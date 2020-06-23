function residuals_ct_sun_demand = ct_solve_sun_demand(x_guess_ct_sun_demand,params,params_ess)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c

CT   = x_guess_ct_sun_demand(1);

Csz = params(1);
Nsz = params(2);
Vsz = params(3);
Csc = params(4);
Nsc = params(5);
Vsc = params(6);

Cdz = params(7);
Ndz = params(8);
Vdz = params(9);
Cdc = params(10);
Ndc = params(11);
Vdc = params(12);

Cess = params_ess(1);
Vess = params_ess(2);

% Unconditional Probablity Calculation: Demand and Sunspot Shock.
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

if cCHIc == 1
    Usz = log(Csz + CT*Cess) - Nsz^(1+cCHIn)/(1+cCHIn);
    Usc = log(Csc + CT*Cess) - Nsc^(1+cCHIn)/(1+cCHIn);
    Udz = log(Cdz + CT*Cess) - Ndz^(1+cCHIn)/(1+cCHIn);
    Udc = log(Cdc + CT*Cess) - Ndc^(1+cCHIn)/(1+cCHIn);
else
    Usz = (Csz + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nsz^(1+cCHIn)/(1+cCHIn);
    Usc = (Csc + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nsc^(1+cCHIn)/(1+cCHIn);
    Udz = (Cdz + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Ndz^(1+cCHIn)/(1+cCHIn);
    Udc = (Cdc + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Ndc^(1+cCHIn)/(1+cCHIn);
end


A = [1-cBET*cDELz*p_z*p_s,        -cBET*cDELz*(1-p_z)*p_s,      -cBET*cDELz*p_z*(1-p_s),       -cBET*cDELz*(1-p_z)*(1-p_s);...
    -cBET*cDELc*(1-p_c)*p_s,      1-cBET*cDELc*p_c*p_s,         -cBET*cDELc*(1-p_c)*(1-p_s),   -cBET*cDELc*p_c*(1-p_s);...
    -cBET*cDELz*p_z*(1-p_d),      -cBET*cDELz*(1-p_z)*(1-p_d),  1-cBET*cDELz*p_z*p_d,          -cBET*cDELz*(1-p_z)*p_d;...
    -cBET*cDELc*(1-p_c)*(1-p_d),  -cBET*cDELc*p_c*(1-p_d),      -cBET*cDELc*(1-p_c)*p_d,       1-cBET*cDELc*p_c*p_d];

b = [Usz;Usc;Udz;Udc];
x = A\b;

Vsz = x(1);
Vsc = x(2);
Vdz = x(3);
Vdc = x(4);

V = unc_prob_sz*Vsz + unc_prob_sc*Vsc + unc_prob_dz*Vdz + unc_prob_dc*Vdc;

residuals_ct_sun_demand = V - Vess;

