clear
clf

##Ta = [300 307.80 320.55 333.85 347.70 362.10 377.10 392.75 409.05 426.00 ...
##      443.65 462.05 481.25 501.20 521.95 543.60 566.15 589.65 614.10 639.55 ...
##      666.10, 693.70]';

Ta = [300 307.80 320.55 333.85 347.70 362.10 377.10 392.75 409.05 426.00 443.65 462.05 ...
      481.25 501.2 521.95 543.60]';
N = length(Ta);
Rgas = 8.314; % gas constant [J/K*mol]
Kb = 1.38e-23; %boltzmann constant [J/K]
Kb1 = 8.617e-5; %boltzmann constant [eV/K]
D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
P.a = 2.86e-10*1e9; %latrice parameter [nm]
Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
N0 = 1/Vat; %Number of sites per unit volume [1/nm^3]
Z = 1/20; %Zeldovich factor
%gs = 0.062*1e-18; %surfuce tension [J/m^2] 0.174  0.062
gs1 = [0.057 0.059 0.060 0.061 0.062]*1e-18; %surfuce tension [J/m^2] 0.174  0.062
P.Xc0 = 466e-6; %initial value of carbon
P.Xp = 1/9;  %mole fraction of precipitation
%Xeqs = [10.^(3.12-2160./Ta) 10.^(2.43 - 1840./Ta) 10.^(2.48-1770./Ta) ] *1e-2; % equilibrium solute mole fraction of m in the matrix
Xeqs = [0.8 .*10.^(2.43 - 1840./Ta) 10.^(2.43 - 1840./Ta) 1.2 .*10.^(2.43 - 1840./Ta)] *1e-2; % equilibrium solute mole fraction of m in the matrix
D1 = D0.*exp(-Qd./(Kb1.*Ta))*1e+18; %diffusion coefficient [nm^2/s]
Dt1 = [420 480 540];

Ctot = zeros(N+1,length(gs1));
for m = 1:length(gs1),
  
gs = gs1(m);
P.Xeqs = Xeqs(:,2); 
D = (0.7 + 0.15*2) .* D1;
Dt = Dt1(2)
  

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


   
for k=1:N; 

tf(1) = 0.01;
t = logspace(log10(tf(k)), log10(k*Dt), 50);
u = t.*D(k)/P.a^2;

uf(k) = u(end);
tf(k+1) = t(end);

V = [P.Xeqs(k) P.F0(k) P.S0(k) P.b0(k) P.dG0(k) P.a_R0(k) P.Xp];

ifunc = @(x,u) nuclea_anneal(x,u,V);


N1(1) = -A*B*expint(B/u(1))+ A* exp(-B/u(1))*u(1);

frqT(1) = 1;
for w=2:N; 
frqT(w) = Ta(w)/Ta(w-1);
end

tic
x = lsode (ifunc, [N1(k) R1(k)*frqT(k) C1(k)], u); 
toc


[xdot, F, S] = nuclea_anneal(x',u,V);


N1(k+1) = x(end,1);
C1(k+1) = x(end,3);
R1(k+1) = x(end,2);
S1(k) = S(end);
F1(k) = F(end);



##figure 1
##
##plot(t/60,x(:,3),'.-')
##hold on
##xlabel('t (min)');
##ylabel('Solute mole fraction');
end
##hold off
Ctot(:,m) = C1;
##figure 2
plot(Ta,C1(2:end)*1e6,'o-')
xlabel('T (K) ');
ylabel('Solute mole fraction (ppm)');
hold on
end
hold off

title('Step annealing, \Delta t = 8 mins, \gamma = 0.060 J/m^2');
title('Step annealing, \gamma = 0.060 J/m^2')
legend('\gamma = 0.057 J/m^2', '\gamma = 0.059 J/m^2', '\gamma = 0.060 J/m^2',...
'\gamma = 0.061 J/m^2', '\gamma = 0.062 J/m^2', 'location', 'southwest' )
%legend('0.8 * X_{eq}','X_{eq}', '1.2 X_{eq}', 'location', 'southwest' )
%legend('0.85 * D','D', '1.15 * D', 'location', 'southwest' )
%legend('\Deltat = 7 mins','\Deltat = 8 mins', '\Deltat = 9 mins', 'location', 'southwest' )

A5 = [Ta Ctot(2:end,:)];
save -ascii CvsT.dat A5