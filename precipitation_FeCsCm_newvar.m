clear
clf

t =logspace(-2,2,50);

Ta = 473; %ageing temperature [K]
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]

Xq0s = 0.01;%Pre-exp. solubility limit 
Xq0m = 0.026;%Pre-exp. solubility limit 
Xq0 = [0.01 0.026]; %Pre-exp. solubility limit

Qp = 28400; %energy for the solubility limit [J/mol]
P.Xeq = Xq0* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix

D0 = 6.2e-7; %Pre-exp. diffusion [m^2/s]
Qd = 80000; %energy for diffusion [J/mol]
D =D0*exp(-Qd/(Rgas*Ta))*1e+18; %diffusion coefficient [nm^2/s]

P.a = 2.86e-10*1e9; %latrice parameter [nm]
Z = 1/20; %Zeldovich factor

gm = 0.147*1e-18; %surfuce tension [J/m^2]
gs = 0.174*1e-18; %surfuce tension [J/m^2]
gamma = [0.174*1e-18 0.1475*1e-18]; %surfuce tension [J/m^2]

P.Xc0 = 0.0007; %initial value of carbon 
P.Xp = 0.25;  %mole fraction of precipitation

Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
N0 = 1/Vat; %Number of sites per unit volume [1/nm^3]

P.R0 = (2*gamma*Vat)/(Kb*Ta); %[nm] 

P.b0 = (4*pi*P.R0.^2*P.Xc0*Z)./P.a^2;
P.dG0 = (4/3)*pi*(P.R0.^2).*gamma/Kb/Ta;
P.h = (4/3)*pi*N0*P.R0.^3;

u = D*t'./P.a^2;
u = linspace(1e0,1e4,1e5);
u = logspace(0, 9.27, 100)';

function [xdot, Xc] = func(x,u,P)

Xc = (P.Xc0 - P.Xp*P.h(1)*x(1).*x(2).^3 - P.Xp*P.h(2)*x(3).*x(4).^3) ./ (1 - P.h(1)*(x(1).*x(2).^3 - P.h(2)*x(3).*x(4).^3)) ;

S = P.Xp*log(Xc./P.Xeq) +(1-P.Xp).*log((1-Xc)./(1-P.Xeq));

xdot(1) = (P.b0(1)/S(1).^2) .*exp(-P.dG0(1)/S(1).^2) .*exp(-(S(1).^2)./(2*P.b0(1)*u)); 
xdot(2) = (P.a^2/P.R0(1).^2./x(2)) .*((Xc - P.Xeq(1)*exp(1/P.Xp/x(2))) ./ (P.Xp - P.Xeq(1)*exp(1/P.Xp/x(2)))) - xdot(1)./ x(1).*(1.05/S(1) - x(2));

xdot(3) = (P.b0(2)/S(2).^2) .*exp(-P.dG0(2)/S(2).^2) .*exp(-(S(2).^2)./(2*P.b0(2)*u)); 
xdot(4) = (P.a^2/P.R0(2).^2/x(4)) .*((Xc - P.Xeq(2)*exp(1/P.Xp/x(4))) ./ (P.Xp - P.Xeq(2)*exp(1/P.Xp/x(4)))) - xdot(3)./ x(3).*(1.05/S(2) - x(4));

endfunction

ifunc = @(x,u) func(x,u,P);

x = lsode (ifunc, [1e-12 0.88 1e-12 1.12], u); % [1e-5 0.7] 
%x = lsode (ifunc, [1e-4 0.88 1e-4 1.12], u); % [1e-5 0.7] 

%AXc= (P.Xc0 - P.Xp*P.h*x(:,1).*x(:,2).^3) ./ (1 - P.h*x(:,1).*x(:,2).^3) ;
AXc = (P.Xc0 - P.Xp*P.h(1)*x(:,1).*x(:,2).^3 - P.Xp*P.h(2)*x(:,3).*x(:,4).^3) ./ (1 - P.h(1)*(x(:,1).*x(:,2).^3 - P.h(2)*x(:,3).*x(:,4).^3)) ;
AS = P.Xp*log(AXc./P.Xeq) +(1-P.Xp).*log((1-AXc)./(1-P.Xeq));
ARs = P.R0./AS;
ARs1 = ARs./P.R0;

subplot(3,1,1)
loglog(u,x(:,2),'.-')
hold on
loglog(u,x(:,4),'.-')
loglog(u,ARs1(:,1),'.-')
loglog(u,ARs1(:,2),'.-')
xlabel('t (sec)');
ylabel('R (nm)');

subplot(3,1,2)
loglog(u,AXc,'.-')
xlabel('t (sec)');
ylabel('Solute mole fraction');

subplot(3,1,3)
semilogx(u,x(:,1),'.-')
hold on
semilogx(u,x(:,3),'.-')
hold off
xlabel('t (sec)');
ylabel('Density ');
