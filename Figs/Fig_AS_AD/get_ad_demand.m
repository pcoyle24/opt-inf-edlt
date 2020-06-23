function R = get_ad_demand(x_guess,state,S,z,c_check)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c p_c_star

pi = state;

if c_check == 0;
    y = x_guess;
    c = y*(1-cVARPHI/2*(pi/(cPItarg^cALPHA)-1)^2);
else
    c = x_guess;
    y = c/(1-cVARPHI/2*(pi/(cPItarg^cALPHA)-1)^2);
end

% r = max(cPItarg/cBET*(pi/cPItarg)^cPHIpi,1);
if z == 1
    r = max(cPItarg/(cDELz*cBET)*(pi/cPItarg)^cPHIpi,1);
else
    r = max(cPItarg/(cDELc*cBET)*(pi/cPItarg)^cPHIpi,1);
end


if p_z == 1
    exp_ad = p_z*(c^(-cCHIc)*pi^(-1)) + (1-p_z)*(S.Ctc^(-cCHIc)*S.PItc^(-1));
else
    exp_ad = p_c*(c^(-cCHIc)*pi^(-1)) + (1-p_c)*(S.Ctz^(-cCHIc)*S.PItz^(-1));
end

if z == 1
    R = c^(-cCHIc) - cBET*cDELz*r*exp_ad;
else
    R = c^(-cCHIc) - cBET*cDELc*r*exp_ad;
end