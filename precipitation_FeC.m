clear
clf

t =logspace(-4,4,1000);

function xdot = f (x,t)
 
Ta = 473; %ageing temperature [k]
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]

Xq0m = 0.026;%Pre-exp. solubility limit 
Xq0s = 0.01; %Pre-exp. solubility limit 
Qp = 24800; %energy for the solubility limit [J/mol]

Xeqm = Xq0m* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix
Xeqs = Xq0s* exp(-Qp/(Rgas*Ta)); % equilibrium solute mole fraction of m in the matrix

D0 = 6.2e-7; %Pre-exp. diffusion [m^2/s]
Qd = 80000; %energy for diffusion [J/mol]
D =D0*exp(-Qd/(Rgas*Ta)); %diffusion coefficient [m^2/s]

a = 2.86e-10; %latrice parameter [m]
Z = 1/20; %Zeldovich factor
gm = 0.147; %surfuce tension [J/m^2]
gs = 0.174; %surfuce tension [J/m^2]

Xc0 = 0.00007; %initial value of carbon 
Xp = 0.25; % mole fraction of precipitation

Vat = a^3/2;
N0 = 1/Vat;
R0m = (2*gm*Vat)/(Kb*Ta); %[m] 
R0s = (2*gm*Vat)/(Kb*Ta); %[m] 

Xc = (Xc0 - (4/3)*pi*(Xp*x(1).*x(2).^3))./(1 - (4/3)*pi*(Xp*x(1).*x(2).^3));
Sm = Xp*log(Xc/Xeqm)+(1-Xp)*log((1-Xc)/(1-Xeqm)); %thermodynamical function giving the driving force for nucleation 
Dgm = - (Kb*Ta*Sm)/Vat; %driving force [J/m^3]
DGmS = (16/3)*pi*(gm^3/Dgm^2); % Gibbs energy[J]
RmS = R0m./Sm; %critical radius [m]
bm = (4*pi*RmS^2*D*Xc0)./a^4; %[1/s] 
tm = 1/(2*bm*Z); %[t]

xdot(1) = N0*Z*bm*exp(-DGmS/Kb/Ta)*exp(-tm/t);
xdot(2) = (D/x(2))*((Xc - Xeqm*exp(R0m/Xp/x(2)))/(Xp - Xeqm*exp(R0m/Xp/x(2)))) + xdot(1)./x(1)*(1.05*RmS-x(2));

endfunction

x = lsode ("f", [5e-8; 5e-10], t);

semilogx(t,x(:,2))