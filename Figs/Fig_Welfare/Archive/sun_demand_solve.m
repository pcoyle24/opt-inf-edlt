function residuals_markov_sss = sun_demand_solve(x_guess_markov_sun,i,j)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg DELbar p_s p_d cDELz cDELc p_z p_c

%% Load In Guesses
Csz   = x_guess_markov_sun(1);
PIsz  = x_guess_markov_sun(2);
Ysz   = x_guess_markov_sun(3);
Nsz   = x_guess_markov_sun(4);
Vsz   = x_guess_markov_sun(5);

Csc   = x_guess_markov_sun(6);
PIsc  = x_guess_markov_sun(7);
Ysc   = x_guess_markov_sun(8);
Nsc   = x_guess_markov_sun(9);
Vsc   = x_guess_markov_sun(10);

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

Cdz   = x_guess_markov_sun(11);
PIdz  = x_guess_markov_sun(12);
Ydz   = x_guess_markov_sun(13);
Ndz   = x_guess_markov_sun(14);
Vdz   = x_guess_markov_sun(15);

Cdc   = x_guess_markov_sun(16);
PIdc  = x_guess_markov_sun(17);
Ydc   = x_guess_markov_sun(18);
Ndc   = x_guess_markov_sun(19);
Vdc   = x_guess_markov_sun(20);

Wdz  = Ndz^cCHIn*Cdz^cCHIc;
Rdz  = cRzlb;

Wdc  = Ndc^cCHIn*Cdc^cCHIc;

Rdc  = (cPItarg(j)/(cDELc*cBET))*((PIdc/cPItarg(j))^cPHIpi)*((Ydc/Ydc)^cPHIy);
% if Rdc < 1
    Rdc = cRzlb;
% end

%% SZ
residuals_markov_sss(1) = Csz^(-cCHIc) - cBET*cDELz*Rsz*(p_s*(p_z*(Csz^(-cCHIc)*PIsz^(-1)) + (1-p_z)*(Csc^(-cCHIc)*PIsc^(-1))) + (1-p_s)*(p_z*(Cdz^(-cCHIc)*PIdz^(-1)) + (1-p_z)*(Cdc^(-cCHIc)*PIdc^(-1))));
residuals_markov_sss(2) = Ysz/Csz^(cCHIc)*(cVARPHI*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)*PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wsz) - (cBET*cDELz*cVARPHI*(p_s*(p_z*((Ysz/Csz^(cCHIc))*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)*PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i))) + (1-p_z)*((Ysc/Csc^(cCHIc))*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)*PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)))) + (1-p_s)*(p_z*((Ydz/Cdz^(cCHIc))*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)*PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i))) + (1-p_z)*((Ydc/Cdc^(cCHIc))*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)*PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i))))));
residuals_markov_sss(3) = Ysz - Csz - cVARPHI/2*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)^2*Ysz;
residuals_markov_sss(4) = Ysz - Nsz;
if cCHIc == 1
    residuals_markov_sss(5) = Vsz - log(Csz) + Nsz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_s*(p_z*Vsz + (1-p_z)*Vsc) + (1-p_s)*(p_z*Vdz + (1-p_z)*Vdc));
else
    residuals_markov_sss(5) = Vsz - Csz^(1-cCHIc)/(1-cCHIc) + Nsz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_s*(p_z*Vsz + (1-p_z)*Vsc) + (1-p_s)*(p_z*Vdz + (1-p_z)*Vdc));
end

%% SC
residuals_markov_sss(6) = Csc^(-cCHIc) - cBET*cDELc*Rsc*(p_s*((1-p_c)*(Csz^(-cCHIc)*PIsz^(-1)) + p_c*(Csc^(-cCHIc)*PIsc^(-1))) + (1-p_s)*((1-p_c)*(Cdz^(-cCHIc)*PIdz^(-1)) + p_c*(Cdc^(-cCHIc)*PIdc^(-1))));
residuals_markov_sss(7) = Ysc/Csc^(cCHIc)*(cVARPHI*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)*PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wsc)  - (cBET*cDELc*cVARPHI*(p_s*((1-p_c)*((Ysz/Csz^(cCHIc))*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)*PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i))) + p_c*((Ysc/Csc^(cCHIc))*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)*PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)))) + (1-p_s)*((1-p_c)*((Ydz/Cdz^(cCHIc))*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)*PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i))) + p_c*((Ydc/Cdc^(cCHIc))*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)*PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i))))));
residuals_markov_sss(8) = Ysc - Csc - cVARPHI/2*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)^2*Ysc;
residuals_markov_sss(9) = Ysc - Nsc;
if cCHIc == 1
    residuals_markov_sss(10) = Vsc - log(Csc) + Nsc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_s*((1-p_c)*Vsz + p_c*Vsc) + (1-p_s)*((1-p_c)*Vdz + p_c*Vdc));
else
    residuals_markov_sss(10) = Vsc - Csc^(1-cCHIc)/(1-cCHIc) + Nsc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_s*((1-p_c)*Vsz + p_c*Vsc) + (1-p_s)*((1-p_c)*Vdz + p_c*Vdc));
end

%% DZ
residuals_markov_sss(11) = Cdz^(-cCHIc) - cBET*cDELz*Rdz*(p_d*(p_z*(Cdz^(-cCHIc)*PIdz^(-1)) + (1-p_z)*(Cdc^(-cCHIc)*PIdc^(-1))) + (1-p_d)*(p_z*(Csz^(-cCHIc)*PIsz^(-1)) + (1-p_z)*(Csc^(-cCHIc)*PIsc^(-1))));
residuals_markov_sss(12) = Ydz/Cdz^(cCHIc)*(cVARPHI*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)*PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wdz) - (cBET*cDELz*cVARPHI*(p_d*(p_z*((Ydz/Cdz^(cCHIc))*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)*PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i))) + (1-p_z)*((Ydc/Cdc^(cCHIc))*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)*PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)))) + (1-p_d)*(p_z*((Ysz/Csz^(cCHIc))*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)*PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i))) + (1-p_z)*((Ysc/Csc^(cCHIc))*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)*PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i))))));
residuals_markov_sss(13) = Ydz - Cdz - cVARPHI/2*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)^2*Ydz;
residuals_markov_sss(14) = Ydz - Ndz;
if cCHIc == 1
    residuals_markov_sss(15) = Vdz - log(Cdz) + Ndz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_d*(p_z*Vdz + (1-p_z)*Vdc) + (1-p_d)*(p_z*Vsz + (1-p_z)*Vsc));
else
    residuals_markov_sss(15) = Vdz - Cdz^(1-cCHIc)/(1-cCHIc) + Ndz^(1+cCHIn)/(1+cCHIn) - cBET*cDELz*(p_d*(p_z*Vdz + (1-p_z)*Vdc) + (1-p_d)*(p_z*Vsz + (1-p_z)*Vsc));
end

%% DC
residuals_markov_sss(16) = Cdc^(-cCHIc) - cBET*cDELc*Rdc*(p_d*((1-p_c)*(Cdz^(-cCHIc)*PIdz^(-1)) + p_c*(Cdc^(-cCHIc)*PIdc^(-1))) + (1-p_d)*((1-p_c)*(Csz^(-cCHIc)*PIsz^(-1)) + p_c*(Csc^(-cCHIc)*PIsc^(-1))));
residuals_markov_sss(17) = Ydc/Cdc^(cCHIc)*(cVARPHI*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)*PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - (1 - cTHETA) - cTHETA*(1-cTAU)*Wdc)  - (cBET*cDELc*cVARPHI*(p_d*((1-p_c)*((Ydz/Cdz^(cCHIc))*(PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i)) - 1)*PIdz/((cPItarg(j)^cIOTA*PIdz^(1-cIOTA))^cALPHA(i))) + p_c*((Ydc/Cdc^(cCHIc))*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)*PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)))) + (1-p_d)*((1-p_c)*((Ysz/Csz^(cCHIc))*(PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i)) - 1)*PIsz/((cPItarg(j)^cIOTA*PIsz^(1-cIOTA))^cALPHA(i))) + p_c*((Ysc/Csc^(cCHIc))*(PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i)) - 1)*PIsc/((cPItarg(j)^cIOTA*PIsc^(1-cIOTA))^cALPHA(i))))));
residuals_markov_sss(18) = Ydc - Cdc - cVARPHI/2*(PIdc/((cPItarg(j)^cIOTA*PIdc^(1-cIOTA))^cALPHA(i)) - 1)^2*Ydc;
residuals_markov_sss(19) = Ydc - Ndc;
if cCHIc == 1
    residuals_markov_sss(20) = Vdc - log(Cdc) + Ndc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_d*((1-p_c)*Vdz + p_c*Vdc) + (1-p_d)*((1-p_c)*Vsz + p_c*Vsc));
else
    residuals_markov_sss(20) = Vdc - Cdc^(1-cCHIc)/(1-cCHIc) + Ndc^(1+cCHIn)/(1+cCHIn) - cBET*cDELc*(p_d*((1-p_c)*Vdz + p_c*Vdc) + (1-p_d)*((1-p_c)*Vsz + p_c*Vsc));
end

