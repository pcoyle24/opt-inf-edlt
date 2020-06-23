global cBET cCHIc cCHIn cPHIpi cPHIy cRzlb cDELz cDELc p_z p_c

cPItarg_ess = 0/400 + 1;

disp('Solving for Efficient Steady State Allocations')

options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');
x0 = ones(4,1);
x0 = [x0;-200];
func = @(x) ess_solve(x);
x_out_ess = fsolve(func,x0,options);

Cess_out   = x_out_ess(1);
PIess_out  = x_out_ess(2);
Yess_out   = x_out_ess(3);
Ness_out   = x_out_ess(4);
% Vess_out   = x_out_ess(5);

Wess_out = Ness_out^cCHIn*Cess_out^cCHIc;
Ress_out = (cPItarg_ess/cBET)*((PIess_out/cPItarg_ess)^cPHIpi)*((Yess_out/Yess_out)^cPHIy);
if Ress_out < 1
    Ress_out = cRzlb;
end

Cess   = Cess_out;
PIess  = PIess_out;
Yess   = Yess_out;
Ness   = Ness_out;
Wess   = Wess_out;
Ress   = Ress_out;
% Vess   = Vess_out;

Uess = log(Cess) - Ness.^(1+cCHIn)./(1+cCHIn);

% Unconditional Probablity Calculation: Fundamentals Shock only.
unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

A = [1-cBET*cDELz*p_z, -cBET*cDELz*(1-p_z);...
    -cBET*cDELc*(1-p_c), 1-cBET*cDELc*p_c];

b = [Uess;Uess];
x = A\b;

Vz = x(1);
Vc = x(2);


Vess = unc_prob_zero.*Vz + unc_prob_crisis.*Vc;
