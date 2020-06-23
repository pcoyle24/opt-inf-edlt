function R = get_ad_sun(x_guess,state,S,c_check)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg p_s p_d cDELz 

pi = state;

if c_check == 0
    y = x_guess;
    c = y*(1-cVARPHI/2*(pi/(cPItarg^cALPHA)-1)^2);
else
    c = x_guess;
end

r = max(cPItarg/cBET*(pi/cPItarg)^cPHIpi,1);

if p_s == 1
    exp_ad = p_d*(c^(-cCHIc)*pi^(-1)) + (1-p_d)*(S.Ct^(-cCHIc)*S.PIt^(-1));
elseif p_d == 1
    exp_ad = p_s*(c^(-cCHIc)*pi^(-1)) + (1-p_s)*(S.Cd^(-cCHIc)*S.PId^(-1));
else 
    exp_ad = p_d*(c^(-cCHIc)*pi^(-1)) + (1-p_d)*(S.Ct^(-cCHIc)*S.PIt^(-1));
end

R = c^(-cCHIc) - cBET*cDELz*r*exp_ad;