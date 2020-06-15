F{1} = @(e, y) 0.6*y(3);
F{2} = @(e, y) -0.6*y(3)+0.001407*y(4)*y(3);
F{3} = @(e, y) -0.001407*y(4)*y(3);

steps = 24;

sol1 = coupled_ode(0,F,steps,0,24,[0 5 995],"euler");
sol2 = coupled_ode(0,F,steps,0,24,[0 5 995],"heun");
sol3 = coupled_ode(0,F,steps,0,24,[0 5 995],"rk4");

plot(sol1(:,1),sol1(:,4),sol2(:,1),sol2(:,4),sol3(:,1),sol3(:,4));
legend("Euler", "Heun", "RK4");