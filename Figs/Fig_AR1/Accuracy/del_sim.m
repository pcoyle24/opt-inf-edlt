function del = del_sim(sim,P)

del = zeros(sim,1);
del(1) = 1;
rn = rand(sim,2);

for i = 2:sim
    u1 = rn(i,1);
    u2 = rn(i,2);
    
    ed_today = P.sigma*(-2*log(u1))^(1/2)*cos(2*pi*u2);
    del(i) = P.rho*(del(i-1) - 1) + 1 + ed_today;     
end
