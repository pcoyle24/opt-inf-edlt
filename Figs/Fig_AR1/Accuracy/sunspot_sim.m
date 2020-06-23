function sunspot = sunspot_sim(sim,P)

% Allocate space for output and assume starting off target regime
sunspot = nan(sim,1);
sunspot(1) = 1;
ss = rand(sim,1);

% 1 = Target Regime
% 2 = Deflationary Regime

% Determine Regime
for i = 2:sim
    if sunspot(i-1) == 1
        tol = P.Ps;
        if ss(i) <= tol
            sunspot(i) = 1;
        else
            sunspot(i) = 2;
        end
    elseif sunspot(i-1) == 2
        tol = P.Pd;
        if ss(i) <= tol
            sunspot(i) = 2;
        else
            sunspot(i) = 1;
        end
    end
end
