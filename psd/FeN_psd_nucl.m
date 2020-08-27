% Calculate FeN nucleation & growth at one temperature
clear

% options
Ta = 373; % ageing T (K)
dt = 10*60; % ageing time in s
gs = 0.062; %surface tension [J/m^2] 
X0 = 8.8e-4; % Initial N concentration, Abiko
incub=1; % Calc. incubation time for nucleation


% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]
afe = 0.286; % latrice parameter [nm]
Vat = afe^3/2; %Atomic volume of bcc iron [nm^3]
rat = (3*Vat/4/pi)^(1/3); % atomic radius [nm]

% time grid - sqrt(t)
nt = 51;
t = linspace(0,sqrt(dt),nt);
t = t.^2;

% N diffusion
D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
D =D0*exp(-Qd./(kb*Ta))*1e+18; %diffusion coefficient [nm^2/s]
gam = D/rat^2; 

% alloy thermodynamics
Xeq = 10.^(2.43 - 1840./Ta) *1e-2; % solubility
Xp = 1/9; % precipitate
gs *= 6.24150913; % convert to eV/nm2
R0 = (2*gs*Vat)./(kb*Ta)/rat; 
dG0 = (4/3)*pi*R0.^2*rat^2*gs/kb./Ta;
Z = 1/20;
S0 = Xp*log(X0./Xeq)+(1-Xp)*log((1-X0)./(1-Xeq));
Rc = R0./S0;
b0 = 4*pi*R0.^2*Z*(rat/afe)^4;

% R grid
% find maximum radius
zfunc = @(r) time2R(r,X0,Xp,Xeq,R0)-dt*gam;
Rmax = fsolve(zfunc,2*Rc)
% define log grid
R = logspace(log10(0.8*Rc),log10(2*Rmax),51);
dR = diff([0 R]);
m = length(R);
Rmid = zeros(size(R));
R3 = zeros(size(R));
Rmid(2:m) = (R(1:m-1)+R(2:m))/2;
R3(2:m) = Rmid(2:m).*(R(1:m-1).^2+R(2:m).^2)/2;

cutoff=1;
dbg=0;
% integrate the PDF
tic
[f, X, dfdt, cutoff] = ...
  pdf_integ(zeros(1,m),X0,t*gam,Xp,Xeq,R,R0,dG0,b0,incub,cutoff,dbg);
toc

% post-process
S = Xp*log(X./Xeq)+(1-Xp)*log((1-X)./(1-Xeq));
R1 = R0/Xp./log(X/Xeq);
b = b0*X;
tau = S.^2/2./b;
Js = b./S.^2 .* exp(-dG0./S.^2);
if incub,
  Js .*= exp(-tau/gam./t);
end

  
Nt = (f*dR')';
Rm = (f*(dR.*Rmid)')'./Nt;
F = (f*(dR.*R3)')';

% plot

figure 2
clf

z = sqrt(t);

subplot(3,2,1)
plot(z,Rm*rat,'.-',z,R1*rat,'.-')
ylabel('R (nm) ');

subplot(3,2,2)
plot(z,[0 diff(Rm)./diff(t)*rat],'.-')
ylabel('dR/dt (nm/s) ');

subplot(3,2,3)
plot(z,X,'.-')
ylabel('Solute mole fraction');

subplot(3,2,4)
plot(z,F,'.-')
ylabel('Transformed volume fraction ');

subplot(3,2,5)
plot(z,Nt,'.-')
xlabel('t^{1/2} (s^{1/2}) ');
ylabel('Clusters per atom');

subplot(3,2,6)
plot(z,Js*gam,'.-')
xlabel('t^{1/2} (s^{1/2}) ');
ylabel('Jn (s^{-1})');

figure 3
clf
subplot(2,1,1)
semilogx(R*rat,f(end,:).*dR,'.-')
title('PDF')
ylabel('f\cdot \Delta R')
hold on
yy=get(gca,'ylim');
semilogx([1 1]*R(cutoff)*rat,[yy(1) 0.9*yy(2)],'k--')
hold off
subplot(2,1,2)
semilogx(R*rat,dfdt(end,:).*dR,'.-')
hold on
yy=get(gca,'ylim');
semilogx([1 1]*R(cutoff)*rat,[yy(1)*1.1 0.9*yy(2)],'k--')
hold off
xlabel('R (nm)')
ylabel('df/dt \Delta R')

