clear
clf

##Ta = [300 307.80 320.55 333.85 347.70 362.10 377.10 392.75 409.05 426.00 ...
##      443.65 462.05 481.25 501.20 521.95 543.60 566.15 589.65 614.10 639.55 ...
##      666.10, 693.70]';

Ta = [300 307.80 320.55 333.85 347.70 362.10 377.10 392.75 409.05 426.00 443.65 462.05 ...
      465 470 475 480 487.485]';

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
Dt = 600;

N = length(Ta);
D =D0.*exp(-Qd./(Kb1.*Ta))*1e+18; %diffusion coefficient [nm^2/s]
P.Xeqs = 10.^(2.43 - 1840./Ta) *1e-2; % equilibrium solute mole fraction of m in the matrix
P.R0s = (2*gs*Vat)./(Kb.*Ta); %[nm] 
P.dG0 = (4/3)*pi.*P.R0s.^2*gs/Kb./Ta;
P.F0 = (4/3)*pi*N0.*P.R0s.^3;
P.S0 = P.Xp*log(P.Xc0./P.Xeqs) +(1-P.Xp).*log((1-P.Xc0)./(1-P.Xeqs));
P.b0 = (4*pi.*P.R0s.^2*P.Xc0*Z)./P.a^2; 
P.a_R0 = P.a^2./P.R0s.^2;

B = P.S0(1).^2 / 2 / P.b0(1);
A = (P.b0(1)./P.S0(1).^2) .*exp(-P.dG0(1)./P.S0(1).^2);

N1 = zeros(N+1,1);
C1 = zeros(N+1,1);
R1 = zeros(N+1,1);
S1 = zeros(N,1);

uf = zeros(N,1);
tf = zeros(N+1,1);

C1(1) = P.Xc0;

R1(1) = 1.05/P.S0(1);
for k=13:N; 

tf(13) = 0.01;
t = logspace(log10(tf(k)), log10((k-12)*480), 100);
%u = sqrt(t.*D(k)/P.a^2);
u = t.*D(k)/P.a^2;

uf(k) = u(end);
tf(k+1) = t(end);

R1(13) = 1.7026e+02; %163.9614;
C1(13) = 3.4146e-04; %0.00030616; 
N1(13) = 5.4044e-11;%7.7648e-11;

V = [P.Xeqs(k) P.F0(k) P.S0(k) P.b0(k) P.Xc0 P.a_R0(k) P.Xp N1(13)];

ifunc = @(x2,u) dissolve_anneal_2(x2,u,V);

frqT(1) = 1;
for w=2:N; 
frqT(w) = Ta(w)/Ta(w-1);
end

tic
x2 = lsode (ifunc, [R1(k)*frqT(k)], u); 
toc

[xdot2, Xc, F] = dissolve_anneal_2(x2',u,V);

R1(k+1) = x2(end,1);
C1(k+1) = Xc(end);
F1(k) = F(end);

figure 1
subplot(2,1,1)
semilogy(t/60,x2(:,1).*P.R0s(k),'o-')
hold on
xlabel('t (min) ');
ylabel('R (nm) ');

subplot(2,1,2)
%plot(t/60,x2(:,2),'.-')
plot(t/60,Xc,'.-')
hold on
xlabel('t (min)');
ylabel('Solute mole fraction');
end
figure 2
subplot(2,1,1)
plot(Ta,R1(2:end).*P.R0s,'.-')
xlabel('T (K) ');
ylabel('R (nm) ');

subplot(2,1,2)
plot(Ta,C1(2:end)*1e6,'.-')
##hold on
##plot(Ta,P.Xeqs*1e6,'.-')
##hold off
xlabel('T (K) ');
ylabel('Solute mole fraction (ppm)');
