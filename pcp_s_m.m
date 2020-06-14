clear
clf

t =logspace(-4,10,20)';

Ta = 473; %ageing temperature [K]
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]

Xq0m = 0.026;%Pre-exp. solubility limit 
Xq0s = 0.01;%Pre-exp. solubility limit 
Xq0 = [0.01 0.026]; %Pre-exp. solubility limit

Qp = 24800; %energy for the solubility limit [J/mol]

P.Xeqm = Xq0m* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix
P.Xeqs = Xq0s* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix
P.Xeq = Xq0.* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix

D0 = 6.2e-7; %Pre-exp. diffusion [m^2/s]
Qd = 80000; %energy for diffusion [J/mol]
P.D =D0*exp(-Qd/(Rgas*Ta))*1e+18; %diffusion coefficient [nm^2/s]

P.a = 2.86e-10*1e9; %latrice parameter [nm]
P.Z = 1/20; %Zeldovich factor
gm = 0.147*1e-18; %surfuce tension [J/m^2]
gs = 0.174*1e-18; %surfuce tension [J/m^2]
gamma = [0.174*1e-18 0.147*1e-18]; %surfuce tension [J/m^2]

P.Xc0 = 0.0007; %initial value of carbon 
P.Xp = 0.25;  %mole fraction of precipitation

Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
P.N0 = 1/Vat; %Number of sites per unit volume [1/nm^3]
P.R0m = (2*gm*Vat)/(Kb*Ta); %[nm] 
P.R0s = (2*gs*Vat)/(Kb*Ta); %[nm]
P.R0 = (2*gamma*Vat)/(Kb*Ta); %[nm]

P.b0 = (4*pi*P.R0.^2*P.Xc0)./P.a^2;
P.dG0 = (4/3)*pi*gamma.*P.R0.^2/Kb/Ta;

function [xdot, Xc] = func(x,t,P)
u = sqrt(P.D*t)./P.a;
%Xc= (P.Xc0 - (4/3)*pi*((P.Xp*x(1).*x(2).^3)+(P.Xp*x(3).*x(4).^3)))./(1 - (4/3)*pi*((P.Xp*x(1).*x(2).^3)+P.Xp*x(3).*x(4).^3));
Xc= (P.Xc0 - (4/3)*pi*(P.Xp*x(1).*x(2).^3))./(1 - (4/3)*pi*(P.Xp*x(1).*x(2).^3));
S = P.Xp*log(Xc./P.Xeq) +(1-P.Xp).*log((1-Xc)./(1-P.Xeq));
##xdot(1) = (P.Z*P.b0(1).*(u.^2)./t./S(1).^2).*exp(-P.dG0(1)./S(1).^2).*exp(-(S(1).^2/P.b0(1)/2/P.Z./u.^2));
##xdot(2) = ((u.^2)*P.a^2/x(2)./t).*((Xc - P.Xeq(1)*exp(P.R0(1)/P.Xp/x(2)))/(1 - P.Xeq(1)*exp(P.R0(1)/P.Xp/x(2)))) + xdot(1)./x(1).*(1.05*P.R0(1)/S(1)-x(2));
##xdot(3) = (P.Z*P.b0(2).*(u.^2)./t./S(2).^2).*exp(-P.dG0(2)./S(2).^2).*exp(-(S(2).^2/P.b0(2)/2/P.Z./u.^2));
##xdot(4) = ((u.^2)*P.a^2/x(4)./t).*((Xc - P.Xeq(2)*exp(P.R0(2)/P.Xp/x(4)))/(1 - P.Xeq(2)*exp(P.R0(2)/P.Xp/x(4)))) + xdot(3)./x(3).*(1.05*P.R0(2)/S(2)-x(4));
xdot(1) = (P.Z*P.b0(2).*(u.^2)./t./S(2).^2).*exp(-P.dG0(2)./S(2).^2).*exp(-(S(2).^2/P.b0(2)/2/P.Z./u.^2));
xdot(2) = ((u.^2)*P.a^2/x(2)./t).*((Xc - P.Xeq(2)*exp(P.R0(2)/P.Xp/x(2)))/(1 - P.Xeq(2)*exp(P.R0(2)/P.Xp/x(2)))) + xdot(1)./x(1).*(1.05*P.R0(1)/S(1)-x(2));

endfunction

ifunc = @(x,t) func(x,t,P);

%x = lsode (ifunc, [1e-12 0.7 1e-12 0.85 ], t); % [1e-5 0.7] 
x = lsode (ifunc, [1e-12 0.85], t); % [1e-5 0.7] 
##
AXc= (P.Xc0 - (4/3)*pi*(P.Xp*x(:,1).*x(:,2).^3))./(1 - (4/3)*pi*(P.Xp*x(:,1).*x(:,2).^3));
##ASs = P.Xp*log(AXc./P.Xeq(1)) +(1-P.Xp)*log((1-AXc)/(1-P.Xeq(1)));
##AN = P.Z*(P.b0(1)./ASs.^2).*exp(-P.dG0(1)./ASs.^2).*exp(-(ASs.^2/P.b0(1)/2/P.Z)./t');
##ARsS = P.R0(1)./ASs;
ASs = P.Xp*log(AXc./P.Xeq(2)) +(1-P.Xp)*log((1-AXc)/(1-P.Xeq(2)));
AN = P.Z*(P.b0(2)./ASs.^2).*exp(-P.dG0(2)./ASs.^2).*exp(-(ASs.^2/P.b0(2)/2/P.Z)./t');
ARsS = P.R0(2)./ASs;
##
##AXc= (P.Xc0 - (4/3)*pi*((P.Xp*x(:,1).*x(:,2).^3)+(P.Xp*x(:,3).*x(:,4).^3)))./(1 - (4/3)*pi*((P.Xp*x(:,1).*x(:,2).^3)+P.Xp*x(:,3).*x(:,4).^3));
##AS = P.Xp*log(AXc./P.Xeq) +(1-P.Xp)*log((1-AXc)./(1-P.Xeq));
##AN = P.Z*(P.b0./AS.^2).*exp(-P.dG0./AS.^2).*exp(-(AS.^2./P.b0/2/P.Z)./t);
##ARsS = P.R0./AS;
##
subplot(3,1,1)
loglog(t,x(:,2),'.-')
hold on
##loglog(t,x(:,4),'.-')
loglog(t,ARsS,'.-')
xlabel('t (sec)');
ylabel('R (nm)');

subplot(3,1,2)
loglog(t,AXc,'.-')
xlabel('t (sec)');
ylabel('Solute mole fraction');

subplot(3,1,3)
loglog(t,x(:,1),'.-')
##hold on
##loglog(t,x(:,3),'.-')
##hold off
xlabel('t (sec)');
ylabel('Density ');
