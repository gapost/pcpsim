clear
clf

t =logspace(-4,10,50);

Ta = 473; %ageing temperature [K]
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]

Xq0m = 0.026;%Pre-exp. solubility limit 
Xq0s = 0.01;%Pre-exp. solubility limit 
Qp = 24800; %energy for the solubility limit [J/mol]

P.Xeqm = Xq0m* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix
P.Xeqs = Xq0s* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix

D0 = 6.2e-7; %Pre-exp. diffusion [m^2/s]
Qd = 80000; %energy for diffusion [J/mol]
P.D =D0*exp(-Qd/(Rgas*Ta))*1e+18; %diffusion coefficient [nm^2/s]

a = 2.86e-10*1e9; %latrice parameter [nm]
P.Z = 1/20; %Zeldovich factor
gm = 0.147*1e-18; %surfuce tension [J/m^2]
gs = 0.174*1e-18; %surfuce tension [J/m^2]

P.Xc0 = 0.0007; %initial value of carbon 
P.Xp = 0.25;  %mole fraction of precipitation

Vat = a^3/2;
P.N0 = 1/Vat;
P.R0m = (2*gm*Vat)/(Kb*Ta); %[nm] 
P.R0s = (2*gs*Vat)/(Kb*Ta); %[nm] 

P.b0 = (4*pi*P.R0s^2*P.D*P.Xc0)./a^4;
P.dG0 = (4/3)*pi*P.R0s^2*gs/Kb/Ta;


##Sm = Xp*log(Xc/Xeqm)+(1-Xp)*log((1-Xc)/(1-Xeqm)); %thermodynamical function giving the driving force for nucleation 
%Dgm = - P.ga*Sm; %driving force [J/m^3]
##DGmS = (16/3)*pi*(gm^3/Dgm^2); % Gibbs energy[J]
##RmS = R0m./Sm; %critical radius [m]
##bm = (4*pi*RmS^2*D*Xc0)./a^4; %[1/s] 
##tm = 1/(2*bm*Z); %[t]
function [xdot, Xc] = func(x,t,P)
Xc= (P.Xc0 - (4/3)*pi*(P.Xp*x(1).*x(2).^3))./(1 - (4/3)*pi*(P.Xp*x(1).*x(2).^3));
Ss = P.Xp*log(Xc./P.Xeqs) +(1-P.Xp).*log((1-Xc)./(1-P.Xeqs));
xdot(1) = P.Z*(P.b0./Ss.^2).*exp(-P.dG0./Ss.^2).*exp(-(Ss.^2/P.b0/2/P.Z)./t);
xdot(2) = (P.D/x(2)).*((Xc - P.Xeqs*exp(P.R0s/P.Xp/x(2)))/(1 - P.Xeqs*exp(P.R0s/P.Xp/x(2)))) + xdot(1)./x(1).*(1.05*P.R0s/Ss-x(2));

endfunction

ifunc = @(x,t) func(x,t,P);

x = lsode (ifunc, [1e-16; 0.7], t); % [1e-5 0.7] 

AXc= (P.Xc0 - (4/3)*pi*(P.Xp*x(:,1).*x(:,2).^3))./(1 - (4/3)*pi*(P.Xp*x(:,1).*x(:,2).^3));
ASs = P.Xp*log(AXc./P.Xeqs) +(1-P.Xp)*log((1-AXc)/(1-P.Xeqs));
AN = P.Z*(P.b0./ASs.^2).*exp(-P.dG0./ASs.^2).*exp(-(ASs.^2/P.b0/2/P.Z)./t');
ARsS = P.R0s./ASs;

subplot(3,1,1)
loglog(t,x(:,2),'.-')
hold on
loglog(t,ARsS,'.-')
xlabel('t (sec)');
ylabel('R (nm)');

subplot(3,1,2)
loglog(t,AXc,'.-')
xlabel('t (sec)');
ylabel('Solute mole fraction');

subplot(3,1,3)
loglog(t,x(:,1),'.-')
xlabel('t (sec)');
ylabel('Density ');
