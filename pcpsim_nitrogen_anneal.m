clear
clf

##Ta = [300 307.80 320.55 333.85 347.70 362.10 377.10 392.75 409.05 426.00 ...
##      443.65 462.05 481.25 501.20 521.95 543.60 566.15 589.65 614.10 639.55 ...
##      666.10, 693.70]';

Ta = [300 307.80 320.55 333.85 347.70 362.10 377.10 392.75 409.05 426.00 443.65 462.05]';

Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]
Kb1 = 8.617e-5; %boltzmann constant [eV/K]
D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
P.a = 2.86e-10*1e9; %latrice parameter [nm]
Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
N0 = 1/Vat; %Number of sites per unit volume [1/nm^3]
Z = 1/20; %Zeldovich factor
gs = 0.062*1e-18; %surfuce tension [J/m^2] 0.174
P.Xc0 = 466e-6; %initial value of carbon
P.Xp = 1/9;  %mole fraction of precipitation

N = length(Ta);
P.Xeqs = zeros(N,1); % equilibrium solute mole fraction of m in the matrix
D = zeros(N,1); %diffusion coefficient [nm^2/s]
P.R0s = zeros(N,1); %[nm] 
P.b0 = zeros(N,1);
P.dG0 = zeros(N,1);
P.F0 = zeros(N,1);
P.S0 = zeros(N,1);
N1 = zeros(N+1,1);
C1 = zeros(N+1,1);
R1 = zeros(N+1,1);
S1 = zeros(N,1);

uf = zeros(N,1);
tf = zeros(N+1,1);
for k=1:N; 
C1(1) = P.Xc0; 
 
P.Xeqs(k) = 10.^(2.43 - 1840./Ta(k)) *1e-2;
D(k) =D0*exp(-Qd/(Kb1*Ta(k)))*1e+18; 
P.R0s(k) = (2*gs*Vat)/(Kb*Ta(k)); 

P.b0(k) = (4*pi*P.R0s(k)^2*C1(k)*Z)./P.a^2;
P.dG0(k) = (4/3)*pi*P.R0s(k)^2*gs/Kb/Ta(k);
P.F0(k) = (4/3)*pi*N0*P.R0s(k)^3;
P.S0(k) = P.Xp*log(C1(k)./P.Xeqs(k)) +(1-P.Xp).*log((1-C1(k))./(1-P.Xeqs(k)));
##
##u = logspace(0, 8, 51);
##t = u*P.a^2/D(k)/60;
##i = find (t<=10);
##u = u(i);
##t = t(i);

tf(1) = 0.01;
%t = logspace(log10(tf(k)), log10(k*600), 50);
t = linspace(tf(k), k*600, 50);
u = t.*D(k)/P.a^2;


uf(k) = u(end);
tf(k+1) = t(end);

%ifunc = @(x,u) nuclea_coars(x,u,P);
%ifunc = @(x,u) nuclea_anneal(x,u,P,k);
ifunc = @(x,u) nuclea_coars_calder_anneal(x,u,P,k);

B = P.S0(1).^2 / 2 / P.b0(1);
A = (P.b0(1)./P.S0(1).^2) .*exp(-P.dG0(1)./P.S0(1).^2);

N1(1) = -A*B*expint(B/u(1))+ A* exp(-B/u(1))*u(1);
R1(1) = 1.05/P.S0(1);

tic
x = lsode (ifunc, [N1(k) 1.05/P.S0(k) C1(k)], u); %
toc

%[xdot, fcoars, S] = nuclea_coars(x',u,P);
%[xdot, F, S] = nuclea_anneal(x',u,P,k);
[xdot, fcoars, S] = nuclea_coars_calder_anneal(x',u,P,k);

N1(k+1) = x(end,1);
C1(k+1) = x(end,3);
R1(k+1) = x(end,2);
S1(k) = S(end);

figure 1
subplot(3,1,1)
semilogy(t/60,x(:,2).*P.R0s(k),'o-')
hold on
xlabel('t (min) ');
ylabel('R (nm) ');

subplot(3,1,2)
plot(t/60,x(:,3),'.-')
hold on
xlabel('t (min)');
ylabel('Solute mole fraction');

subplot(3,1,3)
plot(t/60,x(:,1)*N0*1e9,'.-')
hold on
xlabel('t (min)');
ylabel('Density (\mu m^3)');
end
hold off
hold off
hold off

figure 2
subplot(4,1,1)
semilogy(Ta,R1(2:end).*P.R0s,'.-')
hold on
semilogy(Ta,P.R0s./S1,'.-')
hold off
xlabel('T (K) ');
ylabel('R (nm) ');

subplot(4,1,2)
plot(Ta,C1(2:end)*1e6,'.-')
xlabel('T (K) ');
ylabel('Solute mole fraction (ppm)');

subplot(4,1,3)
plot(Ta,N1(2:end)*N0*1e9,'.-')
xlabel('T (K)');
ylabel('Density (\mu m^3)');

subplot(4,1,4)
semilogy(Ta,C1(2:end)./P.Xeqs,'.-')
xlabel('T (K)');
ylabel('C / C_{eq}');
