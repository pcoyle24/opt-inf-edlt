function R = eqm_acc(del_today,sun_today,P,S,G,O,C,CondProb)

% R= eqm_acc(weights,P,S,G,O,C,CondProb)
%   Outputs residuals and Analytical Jacobian of the equilibrium system of 
%   equations for least squares minimization/Chebyshev interpolation method.
% Inputs:
%   weights :   Parameter Weights for Cheb Poly (Matrix)
%   P       :   Structure of parameters
%   S       :   Structure of steady state values
%   G       :   Structure of grids
%   O       :   Structure of options
%   C       :   Structure of Chebyshev polynomial parameters
%   CondProb:   Sun Spot Transition Matrix
% Output:
%   R       :   Residuals
%   J       :   Jacobian


if sun_today == 1
    r_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ar_s,C.max,C.T,C.P);
    if r_today > 1
        c_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_s,C.max,C.T,C.P);
        inf_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_s,C.max,C.T,C.P);
    else
        r_today = 1;
        c_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_zlb_s,C.max,C.T,C.P);
        inf_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_zlb_s,C.max,C.T,C.P);
    end    
else
    r_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ar_d,C.max,C.T,C.P);
    if r_today > 1
        c_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_d,C.max,C.T,C.P);
        inf_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_d,C.max,C.T,C.P);
    else
        r_today = 1;
        c_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ac_zlb_d,C.max,C.T,C.P);
        inf_today = Fallcheb111(O.delbound,1,del_today,O.n1,C.Ainf_zlb_d,C.max,C.T,C.P);
    end   
end

pi_tilde_today     = inf_today/(S.pi_targ)^P.alpha;
n_today            = (c_today/(1-(P.varphi/2)*(pi_tilde_today-1)^2));
y_today            = n_today;
w_today            = n_today^P.chin*c_today^P.chic;

% Get Conditional Probabability 
cond_prob_sun  = CondProb(sun_today,:);

% Get Tomorrow's Shock 
del_tomorrow   = P.rho*(del_today - 1) + 1 + G.e_nodes;

%For GH Integration
exp_ee_s_int_temp = 0;
exp_pc_s_int_temp = 0;
exp_ee_d_int_temp = 0;
exp_pc_d_int_temp = 0;
exp_ee_int = 0;
exp_pc_int = 0;


for j = 1:length(G.e_weight)
     for s = 1:G.sun_nodes
         if s == 1
             r_s_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ar_s,C.max,C.T,C.P);
             if r_s_tomorrow >= 1
                 % Get Interpolated Values
                 %Standard Regime
                 c_s_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ac_s,C.max,C.T,C.P);
                 inf_s_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ainf_s,C.max,C.T,C.P);

                 % Build out other PF Values
                 pi_tilde_s_tomorrow     = inf_s_tomorrow/(S.pi_targ)^P.alpha;
                 n_s_tomorrow            = (c_s_tomorrow/(1-(P.varphi/2)*(pi_tilde_s_tomorrow-1)^2));
                 y_s_tomorrow            = n_s_tomorrow;

                 % Compute GH Integration
                 exp_ee_s             = c_s_tomorrow^(-P.chic)*inf_s_tomorrow^(-1);
                 exp_pc_s             = (y_s_tomorrow/c_s_tomorrow^(P.chic))*P.varphi*(pi_tilde_s_tomorrow - 1)*pi_tilde_s_tomorrow;
                 exp_ee_s_int_temp         = exp_ee_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_s);
                 exp_pc_s_int_temp         = exp_pc_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_s);               
             else 
                 % Get Interpolated Values
                 %Standard Regime
                 c_s_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ac_zlb_s,C.max,C.T,C.P);
                 inf_s_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ainf_zlb_s,C.max,C.T,C.P);

                 % Build out other PF Values
                 pi_tilde_s_tomorrow_zlb     = inf_s_tomorrow_zlb/(S.pi_targ)^P.alpha;
                 n_s_tomorrow_zlb            = (c_s_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_s_tomorrow_zlb-1)^2));
                 y_s_tomorrow_zlb            = n_s_tomorrow_zlb;

                 % Compute GH Integration
                 exp_ee_zlb_s             = c_s_tomorrow_zlb^(-P.chic)*inf_s_tomorrow_zlb^(-1);
                 exp_pc_zlb_s             = (y_s_tomorrow_zlb/c_s_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_s_tomorrow_zlb - 1)*pi_tilde_s_tomorrow_zlb;

                 exp_ee_s_int_temp         = exp_ee_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_s);
                 exp_pc_s_int_temp         = exp_pc_s_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_s);

             end 
         elseif s == 2
             r_d_tomorrow = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ar_d,C.max,C.T,C.P);
             if r_d_tomorrow >= 1
                 % Get Interpolated Values
                 %Deflationary Regime
                 c_d_tomorrow         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ac_d,C.max,C.T,C.P);
                 inf_d_tomorrow       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ainf_d,C.max,C.T,C.P);

                 % Build out other PF Values
                 pi_tilde_d_tomorrow     = inf_d_tomorrow/(S.pi_targ)^P.alpha;
                 n_d_tomorrow            = (c_d_tomorrow/(1-(P.varphi/2)*(pi_tilde_d_tomorrow-1)^2));
                 y_d_tomorrow            = n_d_tomorrow;

                 % Compute GH Integration
                 exp_ee_d             = c_d_tomorrow^(-P.chic)*inf_d_tomorrow^(-1);
                 exp_pc_d             = (y_d_tomorrow/c_d_tomorrow^(P.chic))*P.varphi*(pi_tilde_d_tomorrow - 1)*pi_tilde_d_tomorrow;

                 exp_ee_d_int_temp         = exp_ee_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_d);
                 exp_pc_d_int_temp         = exp_pc_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_d);
             else
                 %Deflationary Regime
                 c_d_tomorrow_zlb         = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ac_zlb_d,C.max,C.T,C.P);
                 inf_d_tomorrow_zlb       = Fallcheb111(O.delbound,1,del_tomorrow(j),O.n1,C.Ainf_zlb_d,C.max,C.T,C.P);

                 % Build out other PF Values
                 pi_tilde_d_tomorrow_zlb     = inf_d_tomorrow_zlb/(S.pi_targ)^P.alpha;
                 n_d_tomorrow_zlb            = (c_d_tomorrow_zlb/(1-(P.varphi/2)*(pi_tilde_d_tomorrow_zlb-1)^2));
                 y_d_tomorrow_zlb            = n_d_tomorrow_zlb;

                 % Compute GH Integration
                 exp_ee_zlb_d             = c_d_tomorrow_zlb^(-P.chic)*inf_d_tomorrow_zlb^(-1);
                 exp_pc_zlb_d             = (y_d_tomorrow_zlb/c_d_tomorrow_zlb^(P.chic))*P.varphi*(pi_tilde_d_tomorrow_zlb - 1)*pi_tilde_d_tomorrow_zlb;

                 exp_ee_d_int_temp         = exp_ee_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_ee_zlb_d);
                 exp_pc_d_int_temp         = exp_pc_d_int_temp + pi^(-0.5)*G.e_weight(j)*(exp_pc_zlb_d);
             end              
         end
     end
     exp_ee_int = cond_prob_sun*[exp_ee_s_int_temp;exp_ee_d_int_temp];
     exp_pc_int = cond_prob_sun*[exp_pc_s_int_temp;exp_pc_d_int_temp];
end
 
 % Residuals
 res1 = c_today^(-P.chic) - P.beta*del_today*r_today*exp_ee_int;
 res2 = (y_today/c_today^(P.chic))*(P.varphi*(pi_tilde_today-1)*pi_tilde_today-(1-P.theta)-P.theta*(1-P.tau)*w_today) - P.beta*del_today*exp_pc_int;             
 
 R = [res1, res2];
