function cALPHAout = param_solve_cALPHA(x_guess_param)

global cPItarg_init cPItarg p_z p_c

converged = 0;
cALPHA   = x_guess_param(1);
step = 0.02;
optinf_past = cPItarg_init;
while converged == 0
    % Main Script
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

    % Unconditional Probablity
    unc_prob_zero = (1-p_c)/(2-p_z-p_c);
    unc_prob_crisis = (1-p_z)/(2-p_z-p_c);

    unc_welfare = unc_prob_zero*Vsz + unc_prob_crisis*Vsc;

    v_max_unc = max(unc_welfare);
    v_max_unc_index = find(v_max_unc == unc_welfare(1,:));

    opt_inf = cPItarg(v_max_unc_index);
%     disp(400*(opt_inf-1))
    
    if opt_inf > cPItarg_init
        if optinf_past < cPItarg_init
            step = abs(step)/2;
        end
        cALPHA = cALPHA - step;        
    elseif opt_inf < cPItarg_init
        if optinf_past >cPItarg_init
            step = abs(step)/2;
        end   
        cALPHA = cALPHA + step;
    elseif cPItarg_init == opt_inf
        converged = 1;
        cALPHAout = cALPHA;
    end
    
    optinf_past = opt_inf;
end


