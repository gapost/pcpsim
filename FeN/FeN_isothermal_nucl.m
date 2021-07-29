% Calculate Fe - 880 ppm N nucleation & growth at 373 K
% Exp: Abiko & Imai 1977
% model by mean R equations
clear


% annealing conditions
Ta = 373; % ageing T (K)
dt = 1e6; % ageing time in s

% Options
incub=1; % Calc. incubation time for nucleation

% alloy data
gam = 0.062; % surface tension [J/m^2]
X0 = 8.8e-4; % Initial N concentration


% log time grid - 10 pts per decade
% t in s
nt = (log10(dt)+1)*10 + 1;
t = logspace(-1,log10(dt),nt);

[x,F] = FeN_model(t,Ta,X0,gam,incub);

clf
subplot(2,2,1)
semilogx(t,x(:,1),'.-')
title('Clusters per atom');

subplot(2,2,2)
loglog(t,x(:,2),'.-')
ylabel('nm');
title('radius R')

subplot(2,2,3)
loglog(t,x(:,3),'.-')
title('matrix X_N');
xlabel('t (s)');

subplot(2,2,4)
semilogx(t,F,'.-')
title('Transformed volume fraction ');
xlabel('t (s)');

% print -dpng FeN_isothermal_nucl


