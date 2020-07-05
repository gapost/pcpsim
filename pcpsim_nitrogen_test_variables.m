clear
clf

Ta = [300 307.80 320.55 333.85 347.70 362.10 377.10 392.75 409.05 426.00 ...
      443.65 462.05 481.25 501.20 521.95 543.60 566.15 589.65 614.10 639.55 ...
      666.10, 693.70]'; %ageing temperature [K]
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]
Kb1 = 8.617e-5; %boltzmann constant [eV/K]

% equilibrium solute mole fraction of m in the matrix
P.Xeq1 = 10.^(3.12-2160./Ta) *1e-2;
P.Xeq2 = 10.^(2.48-1770./Ta) *1e-2;
P.Xeqs = 10.^(2.43 - 1840./Ta) *1e-2;


D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
D =D0*exp(-Qd/(Kb1.*Ta))*1e+18; %diffusion coefficient [nm^2/s]

P.a = 2.86e-10*1e9; %latrice parameter [nm]
Z = 1/20; %Zeldovich factor
gs = 0.062*1e-18; %surfuce tension [J/m^2] 0.174

P.Xc0 = 466e-6; %initial value of carbon 
%P.Xc0 = 0.00088; %initial value of carbon, abiko
P.Xp = 1/9;  %mole fraction of precipitation

Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
N0 = 1/Vat; %Number of sites per unit volume [1/nm^3]

P.R0s = (2*gs*Vat)/(Kb.*Ta); %[nm] 

P.b0 = (4*pi.*P.R0s.^2*P.Xc0*Z)./P.a^2;
%P.dG0 = (4/3)*pi.*P.R0s.^2*gs/Kb./Ta;
P.dG0 = (4/3)*pi.*P.R0s.^2.*gs;
P.h = (4/3)*pi*N0.*P.R0s.^3;
P.S0 = P.Xp*log(P.Xc0./P.Xeqs) +(1-P.Xp).*log((1-P.Xc0)./(1-P.Xeqs));

Rc = P.R0s ./ P.S0'; 
dGc = P.dG0 ./ P.S0'.^2;
bc = (4*pi.*Rc.^2*P.Xc0.*D*Z)./P.a^4; 
tc = 1 ./ (2.*bc);

A = N0.*bc.*exp(-dGc./Ta'/Kb);

%u = D*t'./P.a^2;
u = logspace(0, 15, 51);
%t = u*P.a^2./D/60;
t = logspace(-4, 4, 22);
B = exp(-tc./t);

figure 1
subplot(3,2,1)
semilogy(Ta,P.R0s,'.-')
hold on
semilogy(Ta,abs(Rc),'.-')
hold off
xlabel('T (K) ');
ylabel('R (nm) ');

subplot(3,2,2)
semilogy(Ta,tc,'.-')
xlabel('T(K) ');
ylabel('Incubation time ');

subplot(3,2,3)
semilogy(Ta,dGc,'.-')
xlabel('T (K)');
ylabel('Critical Gibbs energy');

subplot(3,2,4)
semilogy(Ta,A,'.-')
xlabel('T (K)');
ylabel(' N0\cdotb*exp(-dGc./KTa)');

subplot(3,2,5)
semilogy(Ta,P.Xeqs,'.-')
hold on
semilogy([300 700],[P.Xc0 P.Xc0],'.-')
xlabel('T (K)');
ylabel('C_{eq} ');

subplot(3,2,6)
plot(Ta,B,'.-')
xlabel('T (K)');
ylabel('exp(-\tau / t)');

