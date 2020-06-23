function S = get_ss_demand 

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg cDELz cDELc p_z p_c

options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
x0 = ones(10,1);
func = @(x) markov_sss_solve(x);

x_out_markov_sss = fsolve(func,x0,options);

% x0 = x_out_markov_sss;
% x_out_markov_sss = fsolve(func,x0,options);


Ctz_out   = x_out_markov_sss(1);
PItz_out  = x_out_markov_sss(2);
Ytz_out   = x_out_markov_sss(3);
Ntz_out   = x_out_markov_sss(4);
Vtz_out   = x_out_markov_sss(5);


Wtz_out = Ntz_out^cCHIn*Ctz_out^cCHIc;
Rtz_out = (cPItarg/cBET)*((PItz_out/cPItarg)^cPHIpi)*((Ytz_out/Ytz_out)^cPHIy);
if Rtz_out < 1
    Rtz_out = cRzlb;
end


Ctc_out   = x_out_markov_sss(6);
PItc_out  = x_out_markov_sss(7);
Ytc_out   = x_out_markov_sss(8);
Ntc_out   = x_out_markov_sss(9);
Vtc_out   = x_out_markov_sss(10);

Wtc_out = Ntc_out^cCHIn*Ctc_out^cCHIc;
Rtc_out_check = (cPItarg/(cDELc*cBET))*((PItc_out/cPItarg)^cPHIpi)*((Ytc_out/Ytc_out)^cPHIy);

Rtc_out = (cPItarg/(cDELc*cBET))*((PItc_out/cPItarg)^cPHIpi)*((Ytc_out/Ytc_out)^cPHIy);
if Rtc_out < 1
    Rtc_out = cRzlb;
end


S.Ctz   = Ctz_out;
S.PItz  = PItz_out;
S.Ytz   = Ytz_out;
Ntz   = Ntz_out;
Wtz   = Wtz_out;
Rtz   = Rtz_out;
Vtz   = Vtz_out;

S.Ctc   = Ctc_out;
S.PItc  = PItc_out;
S.Ytc   = Ytc_out;
Ntc   = Ntc_out;
Wtc   = Wtc_out;
Rtc   = Rtc_out;
Rtc_check   = Rtc_out_check;
Vtc   = Vtc_out;