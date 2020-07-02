clear
clf

Ta = 473; %ageing temperature [K]
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]

Xq0s = 0.01;%Pre-exp. solubility limit 
Qp = 28400; %energy for the solubility limit [J/mol]
P.Xeqs = Xq0s* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix

D0 = 6.2e-7; %Pre-exp. diffusion [m^2/s]
Qd = 80000; %energy for diffusion [J/mol]
D =D0*exp(-Qd/(Rgas*Ta))*1e+18; %diffusion coefficient [nm^2/s]

P.a = 2.86e-10*1e9; %latrice parameter [nm]
Z = 1/20; %Zeldovich factor
gs = 0.174*1e-18; %surfuce tension [J/m^2]

P.Xc0 = 0.0007; %initial value of carbon 
P.Xp = 0.25;  %mole fraction of precipitation

Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
N0 = 1/Vat; %Number of sites per unit volume [1/nm^3]

P.R0s = (2*gs*Vat)/(Kb*Ta); %[nm] 

P.b0 = (4*pi*P.R0s^2*P.Xc0*Z)./P.a^2;
P.dG0 = (4/3)*pi*P.R0s^2*gs/Kb/Ta;
P.h = (4/3)*pi*N0*P.R0s^3;
P.S0 = P.Xp*log(P.Xc0./P.Xeqs) +(1-P.Xp).*log((1-P.Xc0)./(1-P.Xeqs));

% u = D*t'./P.a^2;

u = logspace(0, 5, 100);

% find lowest u so that t/tau ~ 0.1
u0 = 0.2*P.S0^2/2/P.b0;

i=find(u>=u0);
u=u(i);
u0=u(1)


ifunc = @(x,u) ng(x,u,P);

B = P.S0.^2 / 2 / P.b0;
A = (P.b0./P.S0.^2) .*exp(-P.dG0./P.S0.^2);
Na = -A*B*expint(B/u(1))+ A* exp(-B/u(1))*u(1);
logNa0 = log(A) + log(B) + 2*log(u(1)/B) - B/u(1);

xdot0 = ng([logNa0 1.05/P.S0 P.Xc0]', u(1), P);


tic
x = lsode (ifunc, [logNa0 1.05/P.S0+xdot0(2)*u(1) P.Xc0], u); % [0 0.88 0.9209537139] 
toc

[xdot,S,F] = func(x',u,P);

%u = D*t'./P.a^2;

u = u*P.a^2/D/60;

subplot(3,2,1)
loglog(u,x(:,2)*P.R0s,'.-')
hold on
loglog(u,1./S*P.R0s,'.-')
hold off
xlabel('u ');
ylabel('R ');

subplot(3,2,2)
semilogx(u,xdot(2,:),'.-')
xlabel('u ');
ylabel('dR/du ');

subplot(3,2,3)
loglog(u,x(:,3),'.-')
xlabel('u');
ylabel('Solute mole fraction');

subplot(3,2,4)
semilogx(u,xdot(3,:),'.-')
xlabel('u ');
ylabel('dC/du ');

subplot(3,2,5)
loglog(u,exp(x(:,1))*N0,'.-')
%semilogx(u,F,'.-')
xlabel('u');
ylabel('Density ');

subplot(3,2,6)
semilogx(u,F,'.-')
xlabel('u ');
ylabel('\Phi ');

u1 = u(1);
