function residuals_ess = ess_sun_demand_solve(x_guess_ess)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c

cPItarg_ess = 0/400 + 1;

Cesstz   = x_guess_ess(1);
PIesstz  = x_guess_ess(2);
Yesstz   = x_guess_ess(3);
Nesstz   = x_guess_ess(4);
Vesstz   = x_guess_ess(5);
Wesstz  = Nesstz^cCHIn*Cesstz^cCHIc;
Resstz  = (cPItarg_ess/(cDELz*cBET))*((PIesstz/cPItarg_ess)^cPHIpi)*((Yesstz/Yesstz)^cPHIy);

Cesstc   = x_guess_ess(6);
PIesstc  = x_guess_ess(7);
Yesstc   = x_guess_ess(8);
Nesstc   = x_guess_ess(9);
Vesstc   = x_guess_ess(10);
Wesstc  = Nesstc^cCHIn*Cesstc^cCHIc;
Resstc  = (cPItarg_ess/(cDELc*cBET))*((PIesstc/cPItarg_ess)^cPHIpi)*((Yesstc/Yesstc)^cPHIy);


Cessdz   = x_guess_ess(11);
PIessdz  = x_guess_ess(12);
Yessdz   = x_guess_ess(13);
Nessdz   = x_guess_ess(14);
Vessdz   = x_guess_ess(15);
Wessdz  = Nessdz^cCHIn*Cessdz^cCHIc;
Ressdz  = (cPItarg_ess/(cDELz*cBET))*((Cessdz/cPItarg_ess)^cPHIpi)*((Yessdz/Yessdz)^cPHIy);

Cessdc   = x_guess_ess(16);
PIessdc  = x_guess_ess(17);
Yessdc   = x_guess_ess(18);
Nessdc   = x_guess_ess(19);
Vessdc   = x_guess_ess(20);
Wessdc  = Nessdc^cCHIn*Cessdc^cCHIc;
Ressdc  = (cPItarg_ess/(cDELc*cBET))*((PIessdc/cPItarg_ess)^cPHIpi)*((Yessdc/Yessdc)^cPHIy);


%% TZ
residuals_ess(1) = Cesstz^(-cCHIc) - cBET*cDELz*Resstz*(p_s*(p_z*(Cesstz^(-cCHIc)*PIesstz^(-1)) + (1-p_z)*(Cesstc^(-cCHIc)*PIesstc^(-1))) +...
                    (1-p_s)*(p_z*(Cessdz^(-cCHIc)*PIessdz^(-1)) + (1-p_z)*(Cessdc^(-cCHIc)*PIessdc^(-1))));
residuals_ess(2) = Yesstz/Cesstz^(cCHIc)*(cVARPHI*(PIesstz/(cPItarg_ess^cALPHA) - 1)*PIesstz/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wesstz) - ...
                (cBET*cDELz*cVARPHI*(p_s*(p_z*((Yesstz/Cesstz^(cCHIc))*(PIesstz/(cPItarg_ess^cALPHA) - 1)*PIesstz/(cPItarg_ess^cALPHA)) + ...
                    (1-p_z)*((Yesstc/Cesstc^(cCHIc))*(PIesstc/(cPItarg_ess^cALPHA) - 1)*PIesstc/(cPItarg_ess^cALPHA))) + ...
                (1-p_s)*(p_z*((Yessdz/Cessdz^(cCHIc))*(PIessdz/(cPItarg_ess^cALPHA) - 1)*PIessdz/(cPItarg_ess^cALPHA)) + ...
                    (1-p_z)*((Yessdc/Cessdc^(cCHIc))*(PIessdc/(cPItarg_ess^cALPHA) - 1)*PIessdc/(cPItarg_ess^cALPHA)))));
residuals_ess(3) = Yesstz - Cesstz - cVARPHI/2*(PIesstz/(cPItarg_ess^cALPHA) - 1)^2*Yesstz;
residuals_ess(4) = Yesstz - Nesstz;
if cCHIc == 1
    residuals_ess(5) = Vesstz - log(Cesstz) + Nesstz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_s*(p_z*Vesstz + (1-p_z)*Vesstc) + ...
                        (1-p_s)*(p_z*Vessdz + (1-p_z)*Vessdc));
else
    residuals_ess(5) = Vesstz - Cesstz^(1-cCHIc)/(1-cCHIc) + Nesstz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_s*(p_z*Vesstz + (1-p_z)*Vesstc) + ...
                        (1-p_s)*(p_z*Vessdz + (1-p_z)*Vessdc));
end

%% TC
residuals_ess(6) = Cesstc^(-cCHIc) - cBET*cDELc*Resstc*(p_s*((1-p_c)*(Cesstz^(-cCHIc)*PIesstz^(-1)) + p_c*(Cesstc^(-cCHIc)*PIesstc^(-1))) + ...
                    (1-p_s)*((1-p_c)*(Cessdz^(-cCHIc)*PIessdz^(-1)) + p_c*(Cessdc^(-cCHIc)*PIessdc^(-1))));
residuals_ess(7) = Yesstc/Cesstc^(cCHIc)*(cVARPHI*(PIesstc/(cPItarg_ess^cALPHA) - 1)*PIesstc/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wesstc)  - ...
                (cBET*cDELc*cVARPHI*(p_s*((1-p_c)*((Yesstz/Cesstz^(cCHIc))*(PIesstz/(cPItarg_ess^cALPHA) - 1)*PIesstz/(cPItarg_ess^cALPHA)) + ...
                    p_c*((Yesstc/Cesstc^(cCHIc))*(PIesstc/(cPItarg_ess^cALPHA) - 1)*PIesstc/(cPItarg_ess^cALPHA))) + ...
                (1-p_s)*((1-p_c)*((Yessdz/Cessdz^(cCHIc))*(PIessdz/(cPItarg_ess^cALPHA) - 1)*PIessdz/(cPItarg_ess^cALPHA)) + ...
                    p_c*((Yessdc/Cessdc^(cCHIc))*(PIessdc/(cPItarg_ess^cALPHA) - 1)*PIessdc/(cPItarg_ess^cALPHA)))));
residuals_ess(8) = Yesstc - Cesstc - cVARPHI/2*(PIesstc/(cPItarg_ess^cALPHA) - 1)^2*Yesstc;
residuals_ess(9) = Yesstc - Nesstc;
if cCHIc == 1
    residuals_ess(10) = Vesstc - log(Cesstc) + Nesstc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_s*((1-p_c)*Vesstz + p_c*Vesstc) + ...
                        (1-p_s)*((1-p_c)*Vessdz + p_c*Vessdc));
else
    residuals_ess(10) = Vesstc - Cesstc^(1-cCHIc)/(1-cCHIc) + Nesstc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_s*((1-p_c)*Vesstz + p_c*Vesstc) + ...
                        (1-p_s)*((1-p_c)*Vessdz + p_c*Vessdc));
end

%% DZ
residuals_ess(11) = Cessdz^(-cCHIc) - cBET*cDELz*Ressdz*(p_d*(p_z*(Cessdz^(-cCHIc)*PIessdz^(-1)) + (1-p_z)*(Cessdc^(-cCHIc)*PIessdc^(-1))) + ...
                    (1-p_d)*(p_z*(Cesstz^(-cCHIc)*PIesstz^(-1)) + (1-p_z)*(Cesstc^(-cCHIc)*PIesstc^(-1))));
residuals_ess(12) = Yessdz/Cessdz^(cCHIc)*(cVARPHI*(PIessdz/(cPItarg_ess^cALPHA) - 1)*PIessdz/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wessdz) - ...
                (cBET*cDELz*cVARPHI*(p_d*(p_z*((Yessdz/Cessdz^(cCHIc))*(PIessdz/(cPItarg_ess^cALPHA) - 1)*PIessdz/(cPItarg_ess^cALPHA)) + ...
                    (1-p_z)*((Yessdc/Cessdc^(cCHIc))*(PIessdc/(cPItarg_ess^cALPHA) - 1)*PIessdc/(cPItarg_ess^cALPHA))) + ...
                (1-p_d)*(p_z*((Yesstz/Cesstz^(cCHIc))*(PIesstz/(cPItarg_ess^cALPHA) - 1)*PIesstz/(cPItarg_ess^cALPHA)) + ...
                    (1-p_z)*((Yesstc/Cesstc^(cCHIc))*(PIesstc/(cPItarg_ess^cALPHA) - 1)*PIesstc/(cPItarg_ess^cALPHA)))));
residuals_ess(13) = Yessdz - Cessdz - cVARPHI/2*(PIessdz/(cPItarg_ess^cALPHA) - 1)^2*Yessdz;
residuals_ess(14) = Yessdz - Nessdz;
if cCHIc == 1
    residuals_ess(15) = Vessdz - log(Cessdz) + Nessdz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_d*(p_z*Vessdz + (1-p_z)*Vessdc) + ...
                        (1-p_d)*(p_z*Vesstz + (1-p_z)*Vesstc));
else
    residuals_ess(15) = Vessdz - Cessdz^(1-cCHIc)/(1-cCHIc) + Nessdz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_d*(p_z*Vessdz + (1-p_z)*Vessdc) + ...
                        (1-p_d)*(p_z*Vesstz + (1-p_z)*Vesstc));
end

%% DC
residuals_ess(16) = Cessdc^(-cCHIc) - cBET*cDELc*Ressdc*(p_d*((1-p_c)*(Cessdz^(-cCHIc)*PIessdz^(-1)) + p_c*(Cessdc^(-cCHIc)*PIessdc^(-1))) + ...
                    (1-p_d)*((1-p_c)*(Cesstz^(-cCHIc)*PIesstz^(-1)) + p_c*(Cesstc^(-cCHIc)*PIesstc^(-1))));
residuals_ess(17) = Yessdc/Cessdc^(cCHIc)*(cVARPHI*(PIessdc/(cPItarg_ess^cALPHA) - 1)*PIessdc/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wessdc)  - ...
                (cBET*cDELc*cVARPHI*(p_d*((1-p_c)*((Yessdz/Cessdz^(cCHIc))*(PIessdz/(cPItarg_ess^cALPHA) - 1)*PIessdz/(cPItarg_ess^cALPHA)) + ...
                    p_c*((Yessdc/Cessdc^(cCHIc))*(PIessdc/(cPItarg_ess^cALPHA) - 1)*PIessdc/(cPItarg_ess^cALPHA))) + ...
                (1-p_d)*((1-p_c)*((Yesstz/Cesstz^(cCHIc))*(PIesstz/(cPItarg_ess^cALPHA) - 1)*PIesstz/(cPItarg_ess^cALPHA)) + ...
                    p_c*((Yesstc/Cesstc^(cCHIc))*(PIesstc/(cPItarg_ess^cALPHA) - 1)*PIesstc/(cPItarg_ess^cALPHA)))));
residuals_ess(18) = Yessdc - Cessdc - cVARPHI/2*(PIessdc/(cPItarg_ess^cALPHA) - 1)^2*Yessdc;
residuals_ess(19) = Yessdc - Nessdc;
if cCHIc == 1
    residuals_ess(20) = Vessdc - log(Cessdc) + Nessdc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_d*((1-p_c)*Vessdz + p_c*Vessdc) + ...
                        (1-p_d)*((1-p_c)*Vesstz + p_c*Vesstc));
else
    residuals_ess(20) = Vessdc - Cessdc^(1-cCHIc)/(1-cCHIc) + Nessdc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_d*((1-p_c)*Vessdz + p_c*Vessdc) + ...
                        (1-p_d)*((1-p_c)*Vesstz + p_c*Vesstc));
end

