clear
clf

Ta = 373; %ageing temperature [K]
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]
Kb1 = 8.617e-5; %boltzmann constant [eV/K]

% equilibrium solute mole fraction of m in the matrix
P.Xeq1 = 10.^(3.12-2160./Ta) *1e-2;
P.Xeq2 = 10.^(2.48-1770./Ta) *1e-2;
P.Xeqs = 10.^(2.43 - 1840./Ta) *1e-2;


D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
D =D0*exp(-Qd/(Kb1*Ta))*1e+18; %diffusion coefficient [nm^2/s]

P.a = 2.86e-10*1e9; %latrice parameter [nm]
Z = 1/20; %Zeldovich factor
gs = 0.062*1e-18; %surfuce tension [J/m^2] 0.174

%P.Xc0 = 466e-6; %initial value of carbon 
P.Xc0 = 0.00088; %initial value of carbon, abiko
P.Xp = 1/9;  %mole fraction of precipitation

Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
N0 = 1/Vat; %Number of sites per unit volume [1/nm^3]

P.R0s = (2*gs*Vat)/(Kb*Ta); %[nm] 

P.b0 = (4*pi*P.R0s^2*P.Xc0*Z)./P.a^2;
P.dG0 = (4/3)*pi*P.R0s^2*gs/Kb/Ta;
P.h = (4/3)*pi*N0*P.R0s^3;
P.S0 = P.Xp*log(P.Xc0./P.Xeqs) +(1-P.Xp).*log((1-P.Xc0)./(1-P.Xeqs));

%u = D*t'./P.a^2;
u = logspace(0, 10, 51);


%ifunc = @(x,u) nuclea_coars(x,u,P);
ifunc = @(x,u) nuclea(x,u,P);
%ifunc = @(x,u) nuclea_coars_calder(x,u,P);

B = P.S0.^2 / 2 / P.b0;
A = (P.b0./P.S0.^2) .*exp(-P.dG0./P.S0.^2);

x1 = u(1)/B;
if x1<0.002,
  u1=0.002*B;
  i = find(u>u1);
  u=[u1 u(i)];
endif
x1 = u(1)/B;

Na = -A*B*expint(B/u(1))+ A* exp(-B/u(1))*u(1);

tic
x = lsode (ifunc, [Na 1.05/P.S0 P.Xc0], u); 
toc

%[xdot, fcoars, S] = nuclea_coars(x',u,P);
[xdot, F, S] = nuclea(x',u,P);
%[xdot, fcoars, S] = nuclea_coars_calder(x',u,P);

t = u*P.a^2/D/60;

figure 1
subplot(3,2,1)
loglog(t,x(:,2).*P.R0s,'.-')
hold on
loglog(t,P.R0s./S,'.-')
hold off
xlabel('t (min) ');
ylabel('R (nm) ');

subplot(3,2,2)
semilogx(t,xdot(2,:),'.-')
xlabel('t (min) ');
ylabel('dR/du ');

subplot(3,2,3)
loglog(t,x(:,3),'.-')
xlabel('t (min)');
ylabel('Solute mole fraction');

subplot(3,2,4)
semilogx(t,xdot(3,:),'.-')
xlabel('t (min)');
ylabel('dC/du ');

subplot(3,2,5)
semilogx(t,x(:,1)*N0*1e9,'.-')
xlabel('t (min)');
ylabel('Density (\mu m^3)');

subplot(3,2,6)
%semilogx(u,xdot(1,:),'.-')
semilogx(t,F,'.-')
xlabel('t (min)');
ylabel('Transformed volume fraction ');

