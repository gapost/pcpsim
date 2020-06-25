clear
clf

t =logspace(-2,2,50);

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

u = D*t'./P.a^2;
u = linspace(1e0,1e4,1e5);
u = logspace(0, 9, 50);

function [xdot, Xc] = func(x,u,P)

Xc = (P.Xc0 - P.Xp*P.h*x(1).*x(2).^3) ./ (1 - P.h*x(1).*x(2).^3) ;
S = P.Xp*log(Xc./P.Xeqs) +(1-P.Xp).*log((1-Xc)./(1-P.Xeqs));

xdot(1) = (P.b0/S.^2) .*exp(-P.dG0/S.^2) .*exp(-(S.^2)./(2*P.b0*u)); 

if (x(1)== 0)
  y = 0;
else
  y = xdot(1)./ x(1);
endif

xdot(2) = (P.a^2/P.R0s^2/x(2)) .*((Xc - P.Xeqs*exp(1/P.Xp/x(2))) ./ (P.Xp - P.Xeqs*exp(1/P.Xp/x(2)))) - y.*(1.05/S - x(2));
endfunction

ifunc = @(x,u) func(x,u,P);

x = lsode (ifunc, [1e-12 0.88], u); % [1e-5 0.7] 

AXc= (P.Xc0 - P.Xp*P.h*x(:,1).*x(:,2).^3) ./ (1 - P.h*x(:,1).*x(:,2).^3) ;
AS = P.Xp*log(AXc./P.Xeqs) +(1-P.Xp).*log((1-AXc)./(1-P.Xeqs));
ARsS = P.R0s./AS;
ARsS1 = ARsS/P.R0s;

subplot(3,1,1)
loglog(u,x(:,2),'.-')
hold on
loglog(u,ARsS1,'.-')
xlabel('t (sec)');
ylabel('R (nm)');

subplot(3,1,2)
loglog(u,AXc,'.-')
xlabel('t (sec)');
ylabel('Solute mole fraction');

subplot(3,1,3)
semilogx(u,x(:,1),'.-')
xlabel('t (sec)');
ylabel('Density ');
