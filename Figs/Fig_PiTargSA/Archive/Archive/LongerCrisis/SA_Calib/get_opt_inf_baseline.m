function opt_inf_sss = get_opt_inf_baseline(cALPHA)
%% Parameters
global cPItarg p_z p_c

%% Main Code
%Allocate Space for functions
Vsz   = zeros(1,length(cPItarg));
Vsc   = zeros(1,length(cPItarg));
options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000,'TolFun', 1e-32, 'Display', 'none');

for j = 1:length(cPItarg)
    x0 = ones(10,1);
    func = @(x) demand_solve_alpha(x,j,cALPHA);

    x_out_markov_sss = fsolve(func,x0,options);
    Vsz_out   = x_out_markov_sss(5);
    Vsc_out   = x_out_markov_sss(10);

    Vsz(j)   = Vsz_out;
    Vsc(j)   = Vsc_out;
end

%% Unconditional Probablity

unc_prob_zero = (1-p_c)/(2-p_z-p_c);
unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

unc_welfare = unc_prob_zero*Vsz + unc_prob_crisis*Vsc;

v_max_unc = max(unc_welfare(1,:));
v_max_unc_index = find(v_max_unc == unc_welfare(1,:));


opt_inf_sss = cPItarg(v_max_unc_index);
end