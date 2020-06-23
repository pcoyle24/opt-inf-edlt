function residuals_markov_sss = markov_sss_solve(x_guess_markov_sss)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c

Cz   = x_guess_markov_sss(1);
PIz  = x_guess_markov_sss(2);
Yz   = x_guess_markov_sss(3);
Nz   = x_guess_markov_sss(4);
Vz   = x_guess_markov_sss(5);

Cc   = x_guess_markov_sss(6);
PIc  = x_guess_markov_sss(7);
Yc   = x_guess_markov_sss(8);
Nc   = x_guess_markov_sss(9);
Vc   = x_guess_markov_sss(10);

Wz  = Nz^cCHIn*Cz^cCHIc;
Rz  = (cPItarg/cBET)*((PIz/cPItarg)^cPHIpi)*((Yz/Yz)^cPHIy);
% if Rz < 1
%     Rz = cRzlb ;
% end

Wc  = Nc^cCHIn*Cc^cCHIc;

Rc  = (cPItarg/(cDELc*cBET))*((PIc/cPItarg)^cPHIpi)*((Yc/Yc)^cPHIy);
if Rc < 1
    Rc = cRzlb;
end


residuals_markov_sss(1) = Cz^(-cCHIc) - cBET*cDELz*Rz*(p_z*(Cz^(-cCHIc)*PIz^(-1)) + (1-p_z)*(Cc^(-cCHIc)*PIc^(-1)));
residuals_markov_sss(2) = Yz/Cz^(cCHIc)*(cVARPHI*(PIz/((cPItarg^cIOTA*PIz^(1-cIOTA))^cALPHA) - 1)*PIz/((cPItarg^cIOTA*PIz^(1-cIOTA))^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wz) - (cBET*cDELz*cVARPHI*(p_z*((Yz/Cz^(cCHIc))*(PIz/((cPItarg^cIOTA*PIz^(1-cIOTA))^cALPHA) - 1)*PIz/((cPItarg^cIOTA*PIz^(1-cIOTA))^cALPHA)) + (1-p_z)*((Yc/Cc^(cCHIc))*(PIc/((cPItarg^cIOTA*PIc^(1-cIOTA))^cALPHA) - 1)*PIc/((cPItarg^cIOTA*PIc^(1-cIOTA))^cALPHA))));
residuals_markov_sss(3) = Yz - Cz - cVARPHI/2*(PIz/((cPItarg^cIOTA*PIz^(1-cIOTA))^cALPHA) - 1)^2*Yz;
residuals_markov_sss(4) = Yz - Nz;
if cCHIc == 1
    residuals_markov_sss(5) = Vz - log(Cz) + Nz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_z*Vz + (1-p_z)*Vc);
else
    residuals_markov_sss(5) = Vz - Cz^(1-cCHIc)/(1-cCHIc) + Nz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_z*Vz + (1-p_z)*Vc);
end

residuals_markov_sss(6) = Cc^(-cCHIc) - cBET*cDELc*Rc*(p_c*(Cc^(-cCHIc)*PIc^(-1)) + (1-p_c)*(Cz^(-cCHIc)*PIz^(-1)));
residuals_markov_sss(7) = Yc/Cc^(cCHIc)*(cVARPHI*(PIc/((cPItarg^cIOTA*PIc^(1-cIOTA))^cALPHA) - 1)*PIc/((cPItarg^cIOTA*PIc^(1-cIOTA))^cALPHA) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wc) - (cBET*cDELc*cVARPHI*(p_c*((Yc/Cc^(cCHIc))*(PIc/((cPItarg^cIOTA*PIc^(1-cIOTA))^cALPHA) - 1)*PIc/((cPItarg^cIOTA*PIc^(1-cIOTA))^cALPHA)) + (1-p_c)*((Yz/Cz^(cCHIc))*(PIz/((cPItarg^cIOTA*PIz^(1-cIOTA))^cALPHA) - 1)*PIz/((cPItarg^cIOTA*PIz^(1-cIOTA))^cALPHA))));
residuals_markov_sss(8) = Yc - Cc - cVARPHI/2*(PIc/((cPItarg^cIOTA*PIc^(1-cIOTA))^cALPHA) - 1)^2*Yc;
residuals_markov_sss(9) = Yc - Nc;
if cCHIc == 1
    residuals_markov_sss(10) = Vc - log(Cc) + Nc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_c*Vc + (1-p_c)*Vz);
else
    residuals_markov_sss(10) = Vc - Cc^(1-cCHIc)/(1-cCHIc) + Nc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_c*Vc + (1-p_c)*Vz);
end
