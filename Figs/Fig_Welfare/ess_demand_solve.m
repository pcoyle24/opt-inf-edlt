function residuals_ess = ess_demand_solve(x_guess_ess)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy  cALPHA  cDELz cDELc p_z p_c

cPItarg_ess = 0/400 + 1;

Cessz   = x_guess_ess(1);
PIessz  = x_guess_ess(2);
Yessz   = x_guess_ess(3);
Nessz   = x_guess_ess(4);
Vessz   = x_guess_ess(5);
Wessz  = Nessz^cCHIn*Cessz^cCHIc;
Ressz  = (cPItarg_ess/(cDELz*cBET))*((PIessz/cPItarg_ess)^cPHIpi)*((Yessz/Yessz)^cPHIy);

Cessc   = x_guess_ess(6);
PIessc  = x_guess_ess(7);
Yessc   = x_guess_ess(8);
Nessc   = x_guess_ess(9);
Vessc   = x_guess_ess(10);
Wessc  = Nessc^cCHIn*Cessc^cCHIc;
Ressc  = (cPItarg_ess/(cDELc*cBET))*((PIessc/cPItarg_ess)^cPHIpi)*((Yessc/Yessc)^cPHIy);


residuals_ess(1) = Cessz^(-cCHIc) - cBET*cDELz*Ressz*(p_z*(Cessz^(-cCHIc)*PIessz^(-1)) + (1-p_z)*(Cessc^(-cCHIc)*PIessc^(-1)));
residuals_ess(2) = Yessz/Cessz^(cCHIc)*(cVARPHI*(PIessz/(cPItarg_ess^cALPHA) - 1)*PIessz/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wessz) - ...
                        (cBET*cDELz*cVARPHI*(p_z*((Yessz/Cessz^(cCHIc))*(PIessz/(cPItarg_ess^cALPHA) - 1)*PIessz/(cPItarg_ess^cALPHA)) + ...
                        (1-p_z)*((Yessc/Cessc^(cCHIc))*(PIessc/(cPItarg_ess^cALPHA) - 1)*PIessc/(cPItarg_ess^cALPHA))));
residuals_ess(3) = Yessz - Cessz - cVARPHI/2*(PIessz/(cPItarg_ess^cALPHA) - 1)^2*Yessz;
residuals_ess(4) = Yessz - Nessz;
if cCHIc == 1
    residuals_ess(5) = Vessz - log(Cessz) + Nessz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_z*Vessz + (1-p_z)*Vessc);
else
    residuals_ess(5) = Vessz - Cessz^(1-cCHIc)/(1-cCHIc) + Nessz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_z*Vessz + (1-p_z)*Vessc);
end

residuals_ess(6) = Cessc^(-cCHIc) - cBET*cDELc*Ressc*(p_c*(Cessc^(-cCHIc)*PIessc^(-1)) + (1-p_c)*(Cessz^(-cCHIc)*PIessz^(-1)));
residuals_ess(7) = Yessc/Cessc^(cCHIc)*(cVARPHI*(PIessc/(cPItarg_ess^cALPHA) - 1)*PIessc/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wessc) - ...
                        (cBET*cDELc*cVARPHI*(p_c*((Yessc/Cessc^(cCHIc))*(PIessc/(cPItarg_ess^cALPHA) - 1)*PIessc/(cPItarg_ess^cALPHA)) + ...
                        (1-p_c)*((Yessz/Cessz^(cCHIc))*(PIessz/(cPItarg_ess^cALPHA) - 1)*PIessz/(cPItarg_ess^cALPHA))));
residuals_ess(8) = Yessc - Cessc - cVARPHI/2*(PIessc/(cPItarg_ess^cALPHA) - 1)^2*Yessc;
residuals_ess(9) = Yessc - Nessc;
if cCHIc == 1
    residuals_ess(10) = Vessc - log(Cessc) + Nessc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_c*Vessc + (1-p_c)*Vessz);
else
    residuals_ess(10) = Vessc - Cessc^(1-cCHIc)/(1-cCHIc) + Nessc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_c*Vessc + (1-p_c)*Vessz);
end
