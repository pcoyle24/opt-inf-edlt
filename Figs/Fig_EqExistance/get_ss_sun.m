function [S,flag] = get_ss_sun

global cBET cCHIc cCHIn cTHETA cTAU cVARPHI cPHIpi cPHIy cRzlb cIOTA cALPHA cPItarg p_s p_d cDELz 


options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

x0 = ones(8,1)*1;
func = @(x) sun_solve(x);

[x_out_sun,~,flag] = fsolve(func,x0,options);

Ct_out   = x_out_sun(1);
PIt_out  = x_out_sun(2);
Yt_out   = x_out_sun(3);
Nt_out   = x_out_sun(4);
% Vt_out   = x_out_sun(5);

Wt_out = Nt_out^cCHIn*Ct_out^cCHIc;
Rt_out = (cPItarg/cBET)*((PIt_out/cPItarg)^cPHIpi)*((Yt_out/Yt_out)^cPHIy);

Cd_out   = x_out_sun(5);
PId_out  = x_out_sun(6);
Yd_out   = x_out_sun(7);
Nd_out   = x_out_sun(8);
% Vd_out   = x_out_sun(10);

Wd_out = Nd_out^cCHIn*Cd_out^cCHIc;
Rd_out = cRzlb;              


S.Ct   = Ct_out;
S.PIt  = PIt_out;
S.Yt   = Yt_out;
Nt   = Nt_out;
Wt   = Wt_out;
Rt   = Rt_out;

S.Cd   = Cd_out;
S.PId  = PId_out;
S.Yd   = Yd_out;
Nd   = Nd_out;
Wd   = Wd_out;
Rd   = Rd_out;
