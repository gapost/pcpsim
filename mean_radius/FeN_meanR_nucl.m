% Calculate Fe - 880 ppm N nucleation & growth at 373 K
% Exp: Abiko & Imai 1977
% model by mean R equations
clear

% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]

% annealing conditions
Ta = 373; % ageing T (K)
dt = 1e6; % ageing time in s

% lattice data
afe = 0.286; % lattice parameter [nm]
Vat = afe^3/2; %Atomic volume of bcc iron [nm^3]
rat = (3*Vat/4/pi)^(1/3); % atomic radius [nm]

% N diffusion
D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
D =D0*exp(-Qd./(kb*Ta))*1e+18; %diffusion coefficient [nm^2/s]
gam = D/rat^2; 

% alloy data
gs = 0.062; %surface tension [J/m^2]
X0 = 8.8e-4; % Initial N concentration, Abiko
Xeq = 10.^(2.43 - 1840./Ta) *1e-2; % solubility
Xp = 1/9; % precipitate

% Options
incub=1; % Calc. incubation time for nucleation
dbg=0; % debug level 

% derived quantities
gs *= 6.24150913; % convert to eV/nm2
R0 = (2*gs*Vat)./(kb*Ta)/rat; 
dG0 = 0.5*R0^3;
Z = 1/20;
S0 = Xp*log(X0./Xeq)+(1-Xp)*log((1-X0)./(1-Xeq));
Rc = R0./S0;
b0 = 4*pi*R0.^2*Z*(rat/afe)^4;

% log time grid - 10 pts per decade
% t in s
nt = (log10(dt)+1)*10 + 1;
t = logspace(-1,log10(dt),nt);

% clear figure
figure 1
clf

ifunc = @(t,x) mean_radius_ng(t,x,Xp,Xeq,b0,dG0,R0,incub,dbg);

odeopt = odeset('InitialStep', 0.001*gam,...
  'AbsTol',[1e-32, 1e-6, 1e-6]',...
  'NonNegative',[1 1 1]');

tic
[dd,x] = ode23(ifunc,t*gam,[0 1.05*R0/S0 X0],odeopt);
toc  

F = x(:,2).^3.*x(:,1);
S = Xp*log(x(:,3)./Xeq)+(1-Xp)*log((1-x(:,3))./(1-Xeq));
xdot = zeros(size(x));
for i=1:size(x,1)
  xdot(i,:) = mean_radius_ng(t(i)*gam,x(i,:)',Xp,Xeq,b0,dG0,R0,incub,0);
end


subplot(3,2,1)
loglog(t/60,x(:,2).*rat,'.-',t/60,R0./S*rat,'.-')
ylabel('R, R^* (nm) ');

subplot(3,2,2)
semilogx(t/60,xdot(:,2)*gam*rat,'.-')
ylabel('dR/dt (nm/s) ');

subplot(3,2,3)
semilogx(t/60,x(:,3),'.-')
ylabel('X');

subplot(3,2,4)
semilogx(t/60,F,'.-')
ylabel('Transformed volume fraction ');

subplot(3,2,5)
semilogx(t/60,x(:,1),'.-')
xlabel('t (min) ');
ylabel('Clusters per atom');

subplot(3,2,6)
semilogx(t/60,xdot(:,1)*gam,'.-')
xlabel('t (min) ');
ylabel('Jn (s^{-1})');

print2pdf(gcf,[20 24],'FeN_meanR_nucl')
