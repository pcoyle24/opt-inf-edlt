function residuals_ct_sun_demand = ct_solve_sun_demand(x_guess_ct_sun_demand,params,params_ess)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c

CT   = x_guess_ct_sun_demand(1);

Vtz = params(1);
Vsc = params(2);
Vdz = params(3);
Vdc = params(4);

Cesstz = params_ess(1);
Nesstz = params_ess(2);
Cesstc = params_ess(3);
Nesstc = params_ess(4);
Cessdz = params_ess(5);
Nessdz = params_ess(6);
Cessdc = params_ess(7);
Nessdc = params_ess(8);


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

Cess = unc_prob_sz*Cesstz + unc_prob_sc*Cesstc + unc_prob_dz*Cessdz + unc_prob_dc*Cessdc;

if cCHIc == 1
    Usz = log(Cesstz + CT*Cess) - Nesstz^(1+cCHIn)/(1+cCHIn);
    Usc = log(Cesstc + CT*Cess) - Nesstc^(1+cCHIn)/(1+cCHIn);
    Udz = log(Cessdz + CT*Cess) - Nessdz^(1+cCHIn)/(1+cCHIn);
    Udc = log(Cessdc + CT*Cess) - Nessdc^(1+cCHIn)/(1+cCHIn);
else
    Usz = (Cesstz + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nsz^(1+cCHIn)/(1+cCHIn);
    Usc = (Cesstc + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Nsc^(1+cCHIn)/(1+cCHIn);
    Udz = (Cessdz + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Ndz^(1+cCHIn)/(1+cCHIn);
    Udc = (Cessdc + CT*Cess)^(1-cCHIc)/(1-cCHIc) - Ndc^(1+cCHIn)/(1+cCHIn);
end


A = [1-cBET*cDELz*p_z*p_s,        -cBET*cDELz*(1-p_z)*p_s,      -cBET*cDELz*p_z*(1-p_s),       -cBET*cDELz*(1-p_z)*(1-p_s);...
    -cBET*cDELc*(1-p_c)*p_s,      1-cBET*cDELc*p_c*p_s,         -cBET*cDELc*(1-p_c)*(1-p_s),   -cBET*cDELc*p_c*(1-p_s);...
    -cBET*cDELz*p_z*(1-p_d),      -cBET*cDELz*(1-p_z)*(1-p_d),  1-cBET*cDELz*p_z*p_d,          -cBET*cDELz*(1-p_z)*p_d;...
    -cBET*cDELc*(1-p_c)*(1-p_d),  -cBET*cDELc*p_c*(1-p_d),      -cBET*cDELc*(1-p_c)*p_d,       1-cBET*cDELc*p_c*p_d];

b = [Usz;Usc;Udz;Udc];
x = A\b;

Vesstz = x(1);
Vesstc = x(2);
Vessdz = x(3);
Vessdc = x(4);

Vess = unc_prob_sz*Vesstz + unc_prob_sc*Vesstc + unc_prob_dz*Vessdz + unc_prob_dc*Vessdc;
V = unc_prob_sz*Vtz + unc_prob_sc*Vsc + unc_prob_dz*Vdz + unc_prob_dc*Vdc;


residuals_ct_sun_demand = V - Vess;

