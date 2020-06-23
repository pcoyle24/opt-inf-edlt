function residuals_ess = ess_sun_solve(x_guess_ess)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz p_s p_d

cPItarg_ess = 0/400 + 1;

Cesst   = x_guess_ess(1);
PIesst  = x_guess_ess(2);
Yesst   = x_guess_ess(3);
Nesst   = x_guess_ess(4);
Vesst   = x_guess_ess(5);
Wesst  = Nesst^cCHIn*Cesst^cCHIc;
Resst  = (cPItarg_ess/(cDELz*cBET))*((PIesst/cPItarg_ess)^cPHIpi)*((Yesst/Yesst)^cPHIy);

Cessd   = x_guess_ess(6);
PIessd  = x_guess_ess(7);
Yessd   = x_guess_ess(8);
Nessd   = x_guess_ess(9);
Vessd   = x_guess_ess(10);
Wessd  = Nessd^cCHIn*Cessd^cCHIc;
Ressd  = (cPItarg_ess/(cDELz*cBET))*((PIessd/cPItarg_ess)^cPHIpi)*((Yessd/Yessd)^cPHIy);


residuals_ess(1) = Cesst^(-cCHIc) - cBET*cDELz*Resst*(p_s*(Cesst^(-cCHIc)*PIesst^(-1)) + (1-p_s)*(Cessd^(-cCHIc)*PIessd^(-1)));
residuals_ess(2) = Yesst/Cesst^(cCHIc)*(cVARPHI*(PIesst/(cPItarg_ess^cALPHA) - 1)*PIesst/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wesst) - ...
                (cBET*cDELz*cVARPHI*(p_s*((Yesst/Cesst^(cCHIc))*(PIesst/(cPItarg_ess^cALPHA) - 1)*PIesst/(cPItarg_ess^cALPHA)) + ...
                (1-p_s)*((Yessd/Cessd^(cCHIc))*(PIessd/(cPItarg_ess^cALPHA) - 1)*PIessd/(cPItarg_ess^cALPHA))));
residuals_ess(3) = Yesst - Cesst - cVARPHI/2*(PIesst/(cPItarg_ess^cALPHA) - 1)^2*Yesst;
residuals_ess(4) = Yesst - Nesst;
if cCHIc == 1
    residuals_ess(5) = Vesst - log(Cesst) + Nesst^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_s*Vesst + (1-p_s)*Vessd);
else
    residuals_ess(5) = Vesst - Cesst^(1-cCHIc)/(1-cCHIc) + Nesst^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_s*Vesst + (1-p_s)*Vessd);
end

residuals_ess(6) = Cessd^(-cCHIc) - cBET*cDELz*Ressd*(p_d*(Cessd^(-cCHIc)*PIessd^(-1)) + (1-p_d)*(Cesst^(-cCHIc)*PIesst^(-1)));
residuals_ess(7) = Yessd/Cessd^(cCHIc)*(cVARPHI*(PIessd/(cPItarg_ess^cALPHA) - 1)*PIessd/(cPItarg_ess^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wessd) - ...
                (cBET*cDELz*cVARPHI*(p_d*((Yessd/Cessd^(cCHIc))*(PIessd/(cPItarg_ess^cALPHA) - 1)*PIessd/(cPItarg_ess^cALPHA)) + ...
                (1-p_d)*((Yesst/Cesst^(cCHIc))*(PIesst/(cPItarg_ess^cALPHA) - 1)*PIesst/(cPItarg_ess^cALPHA))));
residuals_ess(8) = Yessd - Cessd - cVARPHI/2*(PIessd/(cPItarg_ess^cALPHA) - 1)^2*Yessd;
residuals_ess(9) = Yessd - Nessd;
if cCHIc == 1
    residuals_ess(10) = Vessd - log(Cessd) + Nessd^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_d*Vessd + (1-p_d)*Vesst);
else
    residuals_ess(10) = Vessd - Cessd^(1-cCHIc)/(1-cCHIc) + Nessd^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_d*Vessd + (1-p_d)*Vesst);
end
