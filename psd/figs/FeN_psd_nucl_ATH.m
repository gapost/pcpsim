% Calculate FeN nucleation & growth at one temperature
clear

% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]

% annealing conditions
Ta = 373; % ageing T (K)
dt = 1e9; % ageing time in s

% lattice data
afe = 0.286; % lattice parameter [nm]
Vat = afe^3/2; %Atomic volume of bcc iron [nm^3]
rat = (3*Vat/4/pi)^(1/3); % atomic radius [nm]



for w=1:3,
  
% N diffusion
D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd1 = [0.95*0.76 0.76 1.05*0.76]; %energy for diffusion [ev]
Qd = Qd1(2);
D =D0*exp(-Qd./(kb*Ta))*1e+18; %diffusion coefficient [nm^2/s]
gam = D/rat^2;   
  
% alloy data
gs1 = [0.060 0.061 0.062]; %surface tension [J/m^2]
gs = gs1(3);
X0 = 8.8e-4; % Initial N concentration, Abiko
Xeq1 = 10.^(2.43 - 1840./Ta) *1e-2; % solubility
Xeq = (0.8+2*0.1)*Xeq1;
Xp = 1/9; % precipitate

% Options
incub=1; % Calc. incubation time for nucleation
dbg=1; % debug level (0: off, 1: messages, 2: msg+plots, 3: msg+more plots)

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

% R grid - R in units of rat
% find maximum radius
zfunc = @(r) time2R(abs(r),X0,Xp,Xeq,R0)-dt*gam;
Rmax = sqrt(Rc^2 + 2*(X0-Xeq)/(Xp-Xeq)*dt*gam); % guess
[Rmax, fval, info] = fsolve(zfunc,Rmax); % refine
% log R grid
R = logspace(log10(0.9*Rc),log10(Rmax/5),41);

% initial values
cutoff=1;
f0 = zeros(size(R));

% integrate the PDF
tic
[f, X, dfdt, cutoff] = ...
  psd_integ(f0,X0,t*gam,Xp,Xeq,R,R0,dG0,b0,incub,cutoff,dbg);
toc

% post-process
S = Xp*log(X./Xeq)+(1-Xp)*log((1-X)./(1-Xeq));
R1 = R0/Xp./log(X/Xeq);
tau = S.^2./X/2/b0;
Js = b0*X./S.^2 .* exp(-dG0./S.^2);
if incub,
  Js .*= exp(-tau./t/gam);
end

% sums 
dR = diff([0 R]);
dR(1) = dR(2);
dR2 = 0.5*diff([0 R.^2]);
dR4 = 0.25*diff([0 R.^4]);
 
Nt = (f*dR')';
Rm = (f*dR2')'./Nt;
F = (f*dR4')';

% plot

##figure 1
##clf
##
##subplot(3,2,1)
##loglog(t,Rm*rat,'.-',t,R1*rat,'.-')
##ylabel('R (nm) ');
##
##subplot(3,2,2)
##semilogx(t,[0 diff(Rm)./diff(t)*rat],'.-')
##ylabel('dR/dt (nm/s) ');
##
##subplot(3,2,3)
##loglog(t,X/Xeq-1,'.-')
##ylabel('X / X_{eq} - 1');
##
##subplot(3,2,4)
##semilogx(t,F,'.-')
##ylabel('Transformed volume fraction ');
##
##subplot(3,2,5)
##semilogx(t,Nt,'.-')
##xlabel('t (s)) ');
##ylabel('Clusters per atom');
##
##subplot(3,2,6)
##semilogx(t,Js*gam,'.-')
##xlabel('t (s) ');
##ylabel('Jn (s^{-1})');
##
##figure 2
##clf
##
##idx=41:10:101;
##lbls = {};
##for k=1:length(idx),
##  lbls = [lbls num2str(t(idx(k)),2)];
##endfor
##
##plot_psd

figure 3
semilogx(t/60,(X-X0)./X0*100)
hold on

endfor

xlim([0.1 1e4])
grid
hold off
legend('0.9*Xeq', 'Xeq', '1.1*Xeq');
title('Annealing T = 100 C');
##A = [t' X'];
##save -ascii tvsX.dat A