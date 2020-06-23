function residuals_ess = ess_solve(x_guess_ess)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cALPHA  

cPItarg_ess = 0/400 + 1;

Cess   = x_guess_ess(1);
PIess  = x_guess_ess(2);
Yess   = x_guess_ess(3);
Ness   = x_guess_ess(4);
Vess   = x_guess_ess(5);

Wess  = Ness^cCHIn*Cess^cCHIc;
Ress  = (cPItarg_ess/cBET)*((PIess/cPItarg_ess)^cPHIpi)*((Yess/Yess)^cPHIy);

residuals_ess(1) = Cess^(-cCHIc) - cBET*Ress*(Cess^(-cCHIc)*PIess^(-1));
residuals_ess(2) = Yess/Cess^(cCHIc)*(cVARPHI*(PIess/(cPItarg_ess^cALPHA) - 1)*PIess/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wess) - ...
                    (cBET*cVARPHI*((Yess/Cess^(cCHIc))*(PIess/(cPItarg_ess^cALPHA) - 1)*PIess/(cPItarg_ess^cALPHA)));
residuals_ess(3) = Yess - Cess - cVARPHI/2*(PIess/(cPItarg_ess^cALPHA) - 1)^2*Yess;
residuals_ess(4) = Yess - Ness;

if cCHIc == 1
    residuals_ess(5) = Vess - log(Cess) + Ness^(1+cCHIn)/(1+cCHIn) - cBET*Vess;
else
    residuals_ess(5) = Vess - Cess^(1-cCHIc)/(1-cCHIc) + Ness^(1+cCHIn)/(1+cCHIn) - cBET*Vess;
end

