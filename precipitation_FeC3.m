clear
clf

t =logspace(-4,-2,1000);

Ta = 473; %ageing temperature [k]
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]

Xq0m = 0.026;%Pre-exp. solubility limit 
Qp = 24800; %energy for the solubility limit [J/mol]

P.Xeqm = Xq0m* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix

D0 = 6.2e-7; %Pre-exp. diffusion [m^2/s]
Qd = 80000; %energy for diffusion [J/mol]
P.D =D0*exp(-Qd/(Rgas*Ta))*1e+18; %diffusion coefficient [nm^2/s]

a = 2.86e-10*1e9; %latrice parameter [nm]
P.Z = 1/20; %Zeldovich factor
gm = 0.147*1e-18; %surfuce tension [J/m^2]

P.Xc0 = 0.0007; %initial value of carbon 
P.Xp = 0.25;  %mole fraction of precipitation

Vat = a^3/2;
P.N0 = 1/Vat;
P.R0m = (2*gm*Vat)/(Kb*Ta); %[nm] 

P.b0 = (4*pi*P.R0m^2*P.D*P.Xc0)./a^4;
P.dG0 = (4/3)*pi*P.R0m^2*gm/Kb/Ta;


##Sm = Xp*log(Xc/Xeqm)+(1-Xp)*log((1-Xc)/(1-Xeqm)); %thermodynamical function giving the driving force for nucleation 
%Dgm = - P.ga*Sm; %driving force [J/m^3]
##DGmS = (16/3)*pi*(gm^3/Dgm^2); % Gibbs energy[J]
##RmS = R0m./Sm; %critical radius [m]
##bm = (4*pi*RmS^2*D*Xc0)./a^4; %[1/s] 
##tm = 1/(2*bm*Z); %[t]
function [xdot, Xc] = func(x,t,P)
Xc= (P.Xc0 - (4/3)*pi*(P.Xp*x(1).*x(2).^3))./(1 - (4/3)*pi*(P.Xp*x(1).*x(2).^3));
Sm = P.Xp*log(Xc./P.Xeqm) +(1-P.Xp)*log((1-Xc)/(1-P.Xeqm));
xdot(1) = P.Z*P.b0*exp(-P.dG0/Sm.^2)*exp(-(P.b0/2/P.Z)/t);
xdot(2) = (P.D/x(2))*((Xc - P.Xeqm*exp(P.R0m/P.Xp/x(2)))/(1 - P.Xeqm*exp(P.R0m/P.Xp/x(2)))) + xdot(1)./x(1)*(1.05*P.R0m/Sm-x(2));
%xdot(2) = (P.D/x(2))*((Xc - P.Xeqm*exp(P.R0m/P.Xp/x(2)))/(1 - P.Xeqm*exp(P.R0m/P.Xp/x(2))));

endfunction

ifunc = @(x,t) func(x,t,P);

x = lsode (ifunc, [5e-4; 0.55], t);
loglog(t,x(:,2))

AXc= (P.Xc0 - (4/3)*pi*(P.Xp*x(:,1).*x(:,2).^3))./(1 - (4/3)*pi*(P.Xp*x(:,1).*x(:,2).^3));
ASm = P.Xp*log(AXc./P.Xeqm) +(1-P.Xp)*log((1-AXc)/(1-P.Xeqm));