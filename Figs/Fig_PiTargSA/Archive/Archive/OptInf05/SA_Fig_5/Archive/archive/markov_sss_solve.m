function residuals_markov_sss = markov_sss_solve(x_guess_markov_sss,i,j)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c

Csz   = x_guess_markov_sss(1);
PIsz  = x_guess_markov_sss(2);
Ysz   = x_guess_markov_sss(3);
Nsz   = x_guess_markov_sss(4);
Vsz   = x_guess_markov_sss(5);

Csc   = x_guess_markov_sss(6);
PIsc  = x_guess_markov_sss(7);
Ysc   = x_guess_markov_sss(8);
Nsc   = x_guess_markov_sss(9);
Vsc   = x_guess_markov_sss(10);

Wsz  = Nsz^cCHIn*Csz^cCHIc;
Rsz  = (cPItarg(j)/cBET)*((PIsz/cPItarg(j))^cPHIpi)*((Ysz/Ysz)^cPHIy);
% if Rsz < 1
%     Rsz = cRzlb ;
% end

Wsc  = Nsc^cCHIn*Csc^cCHIc;

Rsc  = (cPItarg(j)/(cDELc*cBET))*((PIsc/cPItarg(j))^cPHIpi)*((Ysc/Ysc)^cPHIy);
if Rsc < 1
    Rsc = cRzlb;
end


residuals_markov_sss(1) = Csz^(-cCHIc) - cBET*cDELz*Rsz*(p_z*(Csz^(-cCHIc)*PIsz^(-1)) + (1-p_z)*(Csc^(-cCHIc)*PIsc^(-1)));
residuals_markov_sss(2) = Ysz/Csz^(cCHIc)*(cVARPHI*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)*PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wsz) - (cBET*cDELz*cVARPHI*(p_z*((Ysz/Csz^(cCHIc))*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)*PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i))) + (1-p_z)*((Ysc/Csc^(cCHIc))*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)*PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)))));
residuals_markov_sss(3) = Ysz - Csz - cVARPHI/2*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)^2*Ysz;
residuals_markov_sss(4) = Ysz - Nsz;
if cCHIc == 1
    residuals_markov_sss(5) = Vsz - log(Csz) + Nsz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_z*Vsz + (1-p_z)*Vsc);
else
    residuals_markov_sss(5) = Vsz - Csz^(1-cCHIc)/(1-cCHIc) + Nsz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_z*Vsz + (1-p_z)*Vsc);
end

residuals_markov_sss(6) = Csc^(-cCHIc) - cBET*cDELc*Rsc*(p_c*(Csc^(-cCHIc)*PIsc^(-1)) + (1-p_c)*(Csz^(-cCHIc)*PIsz^(-1)));
residuals_markov_sss(7) = Ysc/Csc^(cCHIc)*(cVARPHI*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)*PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wsc) - (cBET*cDELc*cVARPHI*(p_c*((Ysc/Csc^(cCHIc))*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)*PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i))) + (1-p_c)*((Ysz/Csz^(cCHIc))*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)*PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)))));
residuals_markov_sss(8) = Ysc - Csc - cVARPHI/2*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)^2*Ysc;
residuals_markov_sss(9) = Ysc - Nsc;
if cCHIc == 1
    residuals_markov_sss(10) = Vsc - log(Csc) + Nsc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_c*Vsc + (1-p_c)*Vsz);
else
    residuals_markov_sss(10) = Vsc - Csc^(1-cCHIc)/(1-cCHIc) + Nsc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_c*Vsc + (1-p_c)*Vsz);
end
