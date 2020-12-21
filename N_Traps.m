clear
clf

Ta = 243; %ageing temperature [K]
dt = 1e9; % ageing time in s

Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]
Kb1 = 8.617e-5; %boltzmann constant [eV/K]


D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
D =D0*exp(-Qd/(Kb1*Ta))*1e+18; %diffusion coefficient [nm^2/s]

P.a = 2.86e-10*1e9; %latrice parameter [nm]
Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
Xc0 = 466e-6/Vat; %initial N concentration
V2_0 = 40e-6/Vat; %intial V2 concentration

Eb = 0.49; % ev binding energy of VN2, Barouh 2015 table V
K1 = 4*pi*1.7*P.a*D;
K2 = K1./Vat.*exp(-Eb/(Kb1*Ta));

% log time grid - 10 pts per decade
% t in s
nt = (log10(dt)+1)*10 + 1;
t = logspace(-1,log10(dt),nt);

ifunc = @(x,t) traps(x,t,Xc0,V2_0,K1,K2);

tic
x = lsode (ifunc, [Xc0], t); 
toc

[xdot, V2, V2N] = traps(x',t,Xc0,V2_0,K1,K2);

subplot(3,1,1)
semilogx(t,x,'.-')
ylabel('Nitrogen (N/Vat)')
xlabel('time (sec)')

subplot(3,1,2)
semilogx(t,V2,'.-')
ylabel('Traps (T/Vat)')
xlabel('time (sec)')

subplot(3,1,3)
semilogx(t,V2N,'.-')
ylabel('Clusters NT (NT/Vat)')
xlabel('time (sec)')