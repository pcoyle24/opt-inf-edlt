function residuals_param = param_solve_cVARPHI_c(x_guess_param)

global c_init pi_init

cVARPHI   = x_guess_param(1);
c         = x_guess_param(2);

cDELc       = 1 + c;

x0 = ones(10,1);
func = @(x) demand_solve(x,cVARPHI,cDELc);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

x_out_markov_sss = fsolve(func,x0,options);
Csc   = x_out_markov_sss(6);
PIsc  = x_out_markov_sss(7);


residuals_param(1) = c_init - Csc;
residuals_param(2) = pi_init- PIsc ;
