function R = get_as_demand(x_guess,state,S,z,c_check)

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c p_c_star

pi = state;
if c_check == 0
    y = x_guess;
    c = y*(1-cVARPHI/2*(pi/(cPItarg^cALPHA)-1)^2);
else
    c = x_guess;
    y = c/(1-cVARPHI/2*(pi/(cPItarg^cALPHA)-1)^2);
end

n = y;
w = n^cCHIn*c^cCHIc;

if p_z == 1%p_c == p_c_star
    exp_as = p_z*((y/(c^cCHIc))*cVARPHI*(pi/(cPItarg^cALPHA)-1)*pi/(cPItarg^cALPHA)) + (1-p_z)*((S.Ytc/(S.Ctc^cCHIc))*cVARPHI*(S.PItc/(cPItarg^cALPHA)-1)*S.PItc/(cPItarg^cALPHA));
else
    exp_as = p_c*((y/(c^cCHIc))*cVARPHI*(pi/(cPItarg^cALPHA)-1)*pi/(cPItarg^cALPHA)) + (1-p_c)*((S.Ytz/(S.Ctz^cCHIc))*cVARPHI*(S.PItz/(cPItarg^cALPHA)-1)*S.PItz/(cPItarg^cALPHA));
end

if z == 1
    R = y/(c^cCHIc)*(cVARPHI*(pi/(cPItarg^cALPHA)-1)*pi/(cPItarg^cALPHA) - (1-cTHETA)-(1-cTAU)*(cTHETA)*w) - cBET*cDELz*exp_as;
else
    R = y/(c^cCHIc)*(cVARPHI*(pi/(cPItarg^cALPHA)-1)*pi/(cPItarg^cALPHA) - (1-cTHETA)-(1-cTAU)*(cTHETA)*w) - cBET*cDELc*exp_as;
end