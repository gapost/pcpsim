clear
clf

% Model parameters/options
dt = 8*60; % annealing time in s

% annealing conditions
% Annealing temperatures (K)
Ta = [231.6 241.2 251.2 261.6 272.45 283.75 295 307.80 320.55 ...
      333.85 347.70 362.10 377.10 392.75 409.05 426.00 ...
      443.65 462.05 481 500 520];
Ta = [231.6 241.2 251.2 261.6 272.45 283.75 295 307.80 320.55 333.85 ];

nTa = length(Ta);

Rgas = 8.314; % gas constant [J/K*mol]
kb = 8.617e-5; %boltzmann constant [eV/K]


D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
D =D0*exp(-Qd./(kb*Ta))*1e+18; %diffusion coefficient [nm^2/s]

P.a = 2.86e-10*1e9; %latrice parameter [nm]
Vat = P.a^3/2; %Atomic volume of bcc iron [nm^3]
rat = (3*Vat/4/pi)^(1/3); % atomic radius [nm]

Xc0 = 466e-6/Vat; %initial N concentration
def = [6 31 64]*1e-6/Vat; %intial V2 concentration

Eb = 0.49; % ev binding energy of VN2, Barouh 2015 table V

K1 = 4*pi*1.73*P.a.*D;
K2 = K1./Vat.*exp(-Eb./(kb*Ta));

% time grid for each annealing T (step = 20 s) 
t = 1:20:dt; 
nt = length(t);

sol = zeros(nTa,2,3);

for j=1:3,
V2_0 = def(j);  

for i=1:nTa

ifunc = @(x,t) traps(x,t,Xc0,V2_0,K1(i),K2(i));

tic
x = lsode (ifunc, [Xc0], t); 
toc

sol(i,1,j) = x(end);

[ sol(i,2,j)] = ifunc(x(end,:)',t(end));
end
end

N = sol(:,1,:);
dNdt = sol(:,2,:);
V2 = V2_0 - (Xc0 - N); % traps V2
V2N = Xc0 - N;        %cluster of V2N

figure 1
subplot(4,1,1)
plot(Ta,N(:,:,1),'.-',Ta,N(:,:,2),'.-',Ta,N(:,:,3),'.-')
ylabel('Nitrogen (N/Vat)')
xlabel('Ta (K)')

subplot(4,1,2)
plot(Ta,-dNdt(:,:,1),'.-',Ta,-dNdt(:,:,2),'.-',Ta,-dNdt(:,:,3),'.-')
ylabel('dN/dt (N/Vat)')
xlabel('Ta (K)')

subplot(4,1,3)
plot(Ta,V2(:,:,1),'.-',Ta,V2(:,:,2),'.-',Ta,V2(:,:,3),'.-')
ylabel('Traps V2(T/Vat)')
xlabel('Ta (K)')

subplot(4,1,4)
plot(Ta,V2N(:,:,1),'.-',Ta,V2N(:,:,2),'.-',Ta,V2N(:,:,3),'.-')
ylabel('Clusters V2N (NT/Vat)')
xlabel('Ta (K)')

figure 2
%plot(Ta,-dNdt(:,:,1)/Xc0,'.-',Ta,-dNdt(:,:,2)/Xc0,'.-',Ta,-dNdt(:,:,3)/Xc0,'.-')
plot(Ta,-dNdt(:,:,1)/def(1),'.-',Ta,-dNdt(:,:,2)/def(2),'.-',Ta,-dNdt(:,:,3)/def(3),'.-')
ylabel('-1/traps \cdot dN/dt (N/Vat)')
xlabel('Ta (K)')