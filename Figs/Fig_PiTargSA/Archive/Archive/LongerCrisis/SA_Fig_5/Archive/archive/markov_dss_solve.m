function residuals_markov_dss = markov_dss_solve(x_guess_markov_dss,i,j)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c

Cdz   = x_guess_markov_dss(1);
PIdz  = x_guess_markov_dss(2);
Ydz   = x_guess_markov_dss(3);
Ndz   = x_guess_markov_dss(4);
Vdz   = x_guess_markov_dss(5);

Cdc   = x_guess_markov_dss(6);
PIdc  = x_guess_markov_dss(7);
Ydc   = x_guess_markov_dss(8);
Ndc   = x_guess_markov_dss(9);
Vdc   = x_guess_markov_dss(10);

Wdz  = Ndz^cCHIn*Cdz^cCHIc;
Rdz  = cRzlb;

Wdc  = Ndc^cCHIn*Cdc^cCHIc;

% Rdc  = (cPItarg(j)/cBET)*((PIdc/cPItarg(j))^cPHIpi)*((Ydc/Ydc)^cPHIy);
% if Rdc < 1
    Rdc = cRzlb;
% end


residuals_markov_dss(1) = Cdz^(-cCHIc) - cBET*cDELz*Rdz*(p_z*(Cdz^(-cCHIc)*PIdz^(-1)) + (1-p_z)*(Cdc^(-cCHIc)*PIdc^(-1)));
residuals_markov_dss(2) = Ydz/Cdz^(cCHIc)*(cVARPHI*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)*PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wdz) - (cBET*cDELz*cVARPHI*(p_z*((Ydz/Cdz^(cCHIc))*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)*PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i))) + (1-p_z)*((Ydc/Cdc^(cCHIc))*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)*PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)))));
residuals_markov_dss(3) = Ydz - Cdz - cVARPHI/2*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)^2*Ydz;
residuals_markov_dss(4) = Ydz - Ndz;
if cCHIc == 1
    residuals_markov_dss(5) = Vdz - log(Cdz) + Ndz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_z*Vdz + (1-p_z)*Vdc);
else
    residuals_markov_dss(5) = Vdz - Cdz^(1-cCHIc)/(1-cCHIc) + Ndz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_z*Vdz + (1-p_z)*Vdc);
end

residuals_markov_dss(6) = Cdc^(-cCHIc) - cBET*cDELc*Rdc*(p_c*(Cdc^(-cCHIc)*PIdc^(-1)) + (1-p_c)*(Cdz^(-cCHIc)*PIdz^(-1)));
residuals_markov_dss(7) = Ydc/Cdc^(cCHIc)*(cVARPHI*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)*PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wdc) - (cBET*cDELc*cVARPHI*(p_c*((Ydc/Cdc^(cCHIc))*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)*PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i))) + (1-p_c)*((Ydz/Cdz^(cCHIc))*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)*PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)))));
residuals_markov_dss(8) = Ydc - Cdc - cVARPHI/2*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)^2*Ydc;
residuals_markov_dss(9) = Ydc - Ndc;
if cCHIc == 1
    residuals_markov_dss(10) = Vdc - log(Cdc) + Ndc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_c*Vdc + (1-p_c)*Vdz);
else
    residuals_markov_dss(10) = Vdc - Cdc^(1-cCHIc)/(1-cCHIc) + Ndc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_c*Vdc + (1-p_c)*Vdz);
end
