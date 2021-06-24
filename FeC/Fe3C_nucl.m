% Calculate nucleation & growth of metastable iron carbide (Fe3C) in Fe-C alloy
clear
% add path to pcpsim functions
addpath("../m");

% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]
Rgas = 8.314; % gas constant [J/K*mol]

% annealing conditions
Ta = 473; % ageing T (K)
dt = 1e5; % ageing time in s

% Model parameters
afe = 0.286e-9; % lattice parameter [m]
na = 2; % 2 atoms per cubic cell
gam = 0.174; % surface tension [J/m^2]

% C diffusion in Fe
D0 = 6.2e-7; %Pre-exp. diffusion [m^2/s]
Qd = 80000; % Activation energy for diffusion [J/mol]
D = D0*exp(-Qd/(Rgas*Ta)); %diffusion coefficient [m^2/s]

% alloy data
X0 = 7e-4; % Initial C concentration
Xeq = 0.01*exp(-28400/(Rgas*Ta)); % equilibrium concentration 
Xp = 1/4; % precipitate carbon concentration

% Options
incub=1; % Calc. incubation time for nucleation

% log time grid - 10 pts per decade
nt = (log10(dt)+3)*10 + 1;
t = logspace(log10(0.001),log10(dt),nt);

% get calculation parameters
[rat,tau,b0,dG0,R0] = ngparam(afe,na,gam,D,Ta);

% initial condition
S = Xp*log(X0/Xeq) +(1-Xp).*log((1-X0)/(1-Xeq));
x0 = [0; 1.05*R0/S; X0];

% call dae solver
solver = 'daspk'; % 'ode15i', 'daspk' (octave only)
[x,F,S] = ngdae(t/tau,x0,...
  Xp,Xeq,b0,dG0,R0,incub,solver);

clf
subplot(2,2,1)
semilogx(t,x(:,1),'.-')
title('Clusters per atom');

subplot(2,2,2)
loglog(t,x(:,2).*rat,'.-',t,R0./S*rat,'.-')
ylabel('nm');
title('radius R and critical radius R^*')
legend('R','R^*','location','northwest') 

subplot(2,2,3)
loglog(t,x(:,3)/Xeq - 1,'.-')
title('X / X_{eq} - 1 ');
xlabel('t (s) ');

subplot(2,2,4)
semilogx(t,F,'.-')
title('Transformed volume fraction ');
xlabel('t (s) ');

print -dpdfcrop Fe3C_nucl
