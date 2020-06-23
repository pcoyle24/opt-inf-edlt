function R = get_as_sun(x_guess,state,S,c_check)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg p_s p_d cDELz 

pi = state;
y = x_guess;
c = y*(1-cVARPHI/2*(pi/(cPItarg^cALPHA)-1)^2);


n = y;
w = n^cCHIn*c^cCHIc;

if p_s == 1
    exp_as = p_d*((y/(c^cCHIc))*cVARPHI*(pi/(cPItarg^cALPHA)-1)*pi/(cPItarg^cALPHA)) + (1-p_d)*((S.Yt/(S.Ct^cCHIc))*cVARPHI*(S.PIt/(cPItarg^cALPHA)-1)*S.PIt/(cPItarg^cALPHA));
elseif p_d == 1
    exp_as = p_s*((y/(c^cCHIc))*cVARPHI*(pi/(cPItarg^cALPHA)-1)*pi/(cPItarg^cALPHA)) + (1-p_s)*((S.Yd/(S.Cd^cCHIc))*cVARPHI*(S.PId/(cPItarg^cALPHA)-1)*S.PId/(cPItarg^cALPHA));
else
    exp_as = p_d*((y/(c^cCHIc))*cVARPHI*(pi/(cPItarg^cALPHA)-1)*pi/(cPItarg^cALPHA)) + (1-p_d)*((S.Yt/(S.Ct^cCHIc))*cVARPHI*(S.PIt/(cPItarg^cALPHA)-1)*S.PIt/(cPItarg^cALPHA));
end

R = y/(c^cCHIc)*(cVARPHI*(pi/(cPItarg^cALPHA)-1)*pi/(cPItarg^cALPHA) - (1-cTHETA)-(1-cTAU)*(cTHETA)*w) - cBET*cDELz*exp_as;