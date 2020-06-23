function residuals_sun = sun_solve(x_guess_sun,i,j)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz p_s p_d

Cs   = x_guess_sun(1);
PIs  = x_guess_sun(2);
Ys   = x_guess_sun(3);
Ns   = x_guess_sun(4);
Vs   = x_guess_sun(5);

Cd   = x_guess_sun(6);
PId  = x_guess_sun(7);
Yd   = x_guess_sun(8);
Nd   = x_guess_sun(9);
Vd   = x_guess_sun(10);

Ws  = Ns^cCHIn*Cs^cCHIc;
Rs  = (cPItarg(j)/cBET)*((PIs/cPItarg(j))^cPHIpi)*((Ys/Ys)^cPHIy);

Wd  = Nd^cCHIn*Cd^cCHIc;
Rd  = cRzlb;


residuals_sun(1) = Cs^(-cCHIc) - cBET*cDELz*Rs*(p_s*(Cs^(-cCHIc)*PIs^(-1)) + (1-p_s)*(Cd^(-cCHIc)*PId^(-1)));
residuals_sun(2) = Ys/Cs^(cCHIc)*(cVARPHI*(PIs/((cPItarg(j)^cIOTA*PIs^(1-cIOTA))^cALPHA(i)) - 1)*PIs/((cPItarg(j)^cIOTA*PIs^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Ws) - (cBET*cDELz*cVARPHI*(p_s*((Ys/Cs^(cCHIc))*(PIs/((cPItarg(j)^cIOTA*PIs^(1-cIOTA))^cALPHA(i)) - 1)*PIs/((cPItarg(j)^cIOTA*PIs^(1-cIOTA))^cALPHA(i))) + (1-p_s)*((Yd/Cd^(cCHIc))*(PId/((cPItarg(j)^cIOTA*PId^(1-cIOTA))^cALPHA(i)) - 1)*PId/((cPItarg(j)^cIOTA*PId^(1-cIOTA))^cALPHA(i)))));
residuals_sun(3) = Ys - Cs - cVARPHI/2*(PIs/((cPItarg(j)^cIOTA*PIs^(1-cIOTA))^cALPHA(i)) - 1)^2*Ys;
residuals_sun(4) = Ys - Ns;
if cCHIc == 1
    residuals_sun(5) = Vs - log(Cs) + Ns^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_s*Vs + (1-p_s)*Vd);
else
    residuals_sun(5) = Vs - Cs^(1-cCHIc)/(1-cCHIc) + Ns^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_s*Vs + (1-p_s)*Vd);
end

residuals_sun(6) = Cd^(-cCHIc) - cBET*cDELz*Rd*(p_d*(Cd^(-cCHIc)*PId^(-1)) + (1-p_d)*(Cs^(-cCHIc)*PIs^(-1)));
residuals_sun(7) = Yd/Cd^(cCHIc)*(cVARPHI*(PId/((cPItarg(j)^cIOTA*PId^(1-cIOTA))^cALPHA(i)) - 1)*PId/((cPItarg(j)^cIOTA*PId^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wd) - (cBET*cDELz*cVARPHI*(p_d*((Yd/Cd^(cCHIc))*(PId/((cPItarg(j)^cIOTA*PId^(1-cIOTA))^cALPHA(i)) - 1)*PId/((cPItarg(j)^cIOTA*PId^(1-cIOTA))^cALPHA(i))) + (1-p_d)*((Ys/Cs^(cCHIc))*(PIs/((cPItarg(j)^cIOTA*PIs^(1-cIOTA))^cALPHA(i)) - 1)*PIs/((cPItarg(j)^cIOTA*PIs^(1-cIOTA))^cALPHA(i)))));
residuals_sun(8) = Yd - Cd - cVARPHI/2*(PId/((cPItarg(j)^cIOTA*PId^(1-cIOTA))^cALPHA(i)) - 1)^2*Yd;
residuals_sun(9) = Yd - Nd;
if cCHIc == 1
    residuals_sun(10) = Vd - log(Cd) + Nd^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_d*Vd + (1-p_d)*Vs);
else
    residuals_sun(10) = Vd - Cd^(1-cCHIc)/(1-cCHIc) + Nd^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_d*Vd + (1-p_d)*Vs);
end
