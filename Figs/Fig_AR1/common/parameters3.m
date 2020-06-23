function P = parameters3

% P = parameters
%   Defines model parameters
% Output:
%   P : Structure of model parameters

% Household (Annual Calibration)
P.beta     = (1/1.0025);   % Discount factor
P.chic     = 1;            % Constant of relative risk aversion
P.chin     = 1;            % Constant of relative risk aversion
P.theta    = 11;           % CES technology production


% Firm
P.varphi   = 1038;       % Rotemberg adjustment cost coefficient
P.tau      = 1/P.theta; 
P.alpha    = 0.9425;
P.iota     = 1;

pi_targ_min_annual = 1.2;
pi_targ_max_annual = 2;
pi_targ_min = pi_targ_min_annual/400 + 1;
pi_targ_max = pi_targ_max_annual/400 + 1;
pi_targ_n = round((pi_targ_max_annual - pi_targ_min_annual)/0.05 + 1);
P.pi_targ  = linspace(pi_targ_min, pi_targ_max, pi_targ_n);


% Monetary Policy 
P.phi_pi     = 2;       % Inflation coefficient: active interest rate rule
P.phi_y      = 0;       % Inflation coefficient: active interest rate rule

% Fiscal Policy
P.rho     = 0.8;       % Debt coefficient: passive tax rule
P.sigma   = 0.285/100;   % Standard deviation of demand shock
P.bound = (P.sigma^2/(1-P.rho^2))^0.5;



% Transition Probabilities 
% P.Ps      = 0.995; % Probability of transitioning from SSS to SSS 
% P.Ps_grid = (0.99:0.001:1); 
P.Ps_grid = (0.99:0.0005:1); 
% P.Pd      = 0.975; % Probability of transitioning from DSS to DSS
P.Pd_grid = 0.975; %[0.95 0.975 0.99]; 

% For Value Function Iteration
P.tol = 1e-10;
P.maxiter = 10000;