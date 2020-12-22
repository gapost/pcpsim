%% Calculate annealing of Fe - 467ppm N  with mean radius equations
clear

% Model parameters/options
fname = 'data/FeN58_8min_incub.dat';
dt = 8*60; % annealing time in s
gs = 0.058; %surface tension [J/m^2]
incub=1; % Calculate incubation time for nucleation
dbg=0; % debug level 

% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]

% annealing conditions
% Annealing temperatures (K)
Ta = [231.6 241.2 251.2 261.6 272.45 283.75 295 307.80 320.55 333.85 347.70  ...
      362.10 377.10 392.75 409.05 426.00 443.65 462.05 ...
     481 500 520];
%Ta = [231.6 241.2 251.2 261.6 272.45 283.75 295 307.80 320.55 333.85 347.70];
nTa = length(Ta);


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
X0 = 4.67e-4; % Nominal N concentration
Xeq = 10.^(2.43 - 1840./Ta) *1e-2; % solubility
Xp = 1/9; % precipitate

% Traps
V2_0 = 0e-6; %intial V2 traps concentration
Eb = 0.49; % ev binding energy of VN2, Barouh 2015 table V

K1 =1.73; % Trapping coefficient,T Jourdan 2011 table 3 
K2 = K1.*exp(-Eb./(kb*Ta)); % De-trapping coefficient

% derived quantities
gs *= 6.24150913; % convert to eV/nm2
R0 = (2*gs*Vat)./(kb*Ta)/rat; 
dG0 = 0.5*R0.^3;
Z = 1/20;
S = Xp*log(X0./Xeq)+(1-Xp)*log((1-X0)./(1-Xeq));
Rc = R0./S;
b0 = 4*pi*R0.^2*Z*(rat/afe)^4;

% time grid for each annealing T (step = 20 s) 
t = 0:20:dt; 

nt = length(t);

sol = zeros(nTa,8);
V2t = zeros(1,nTa);

figure 1
clf

for i=1:nTa
  
  if i==1,
    x = [0 1.05*R0(1)/S(1) X0 0];
  else
    x = sol(i-1,[1:3 7]);
    
    % check the following:
    %   1) we have nucleation (dN/dt>0) 
    %   2) current R is below the new Rc   
    % if yes, delete all nuclei (they are unstable)
    s = Xp*log(x(3)./Xeq(i))+(1-Xp)*log((1-x(3))./(1-Xeq(i)));
    jn = x(3)*b0(i)/s^2*exp(-dG0(i)/s^2);
    if x(2) < 1.05*R0(i)/s && jn>0,
     % x = [0 1.05*R0(i)/S(i) X0]; 
      x = [0 1.05*R0(i)/S(i) sol(i-1,3) sol(i-1,7)];
      % disp(['R, R*, S, Jn = ' num2str([x(2) 1.05*R0(i)/s s jn])])
    end
  end
  
  disp(Ta(i))
  
  odeopt = odeset('InitialStep', 0.1*gam(i),...
    'AbsTol',[1e-32, 1e-6, 1e-6 1e-6]',...
    'NonNegative',[1 1 1 1]');
  
  ifunc = @(t,x) mean_radius_ng_traps(t,x,Xp,Xeq(i),b0(i),dG0(i),R0(i),incub,dbg,V2_0,K1,K2(i));
  tic
  [ttt,x] = ode23(ifunc,t*gam(i),x,odeopt);
  toc
  sol(i,[1:3 7]) = x(end,:);
  sol(i,[4:6 8]) = ifunc(t(end)*gam(i),x(end,:)');
  
    if sol(i,7)>V2_0,
    V2t(i) = 0;,
    else 
    V2t(i) = V2_0 - sol(i,7);
    end

end

Nt = sol(:,1)';
R = sol(:,2)';
X = sol(:,3)';
NTr = sol(:,7)';
dNdt = sol(:,4)';
dRdt = sol(:,5)';
dXdt = sol(:,6)';
dNTrdt = sol(:,8)';

res = zeros(nTa,1);
res(:,1) = X*0.7 + V2t*1.5 + NTr*0.397*2.2;

F = R.^3.*Nt;
S = Xp*log(X./Xeq)+(1-Xp)*log((1-X)./(1-Xeq));
Rc = R0./S;
i=find(S<0);
Rc(i) = NaN;
i=find(R<0);
R(i) = NaN;

figure 1
subplot(3,2,1)
semilogy(Ta,R*rat,'.-',Ta,Rc*rat,'.-')
ylabel('R, R* (nm) ');

subplot(3,2,2)
plot(Ta,Nt,'.-')
ylabel('Clusters per atom');

subplot(3,2,3)
plot(Ta,X/X0,'.-')
ylabel('X / X_0');
xlabel('Ta (K)');

subplot(3,2,4)
plot(Ta,F,'.-')
ylabel('Transformed volume fraction ');
xlabel('Ta (K)');

subplot(3,2,5)
plot(Ta,NTr,'.-')
ylabel('Trapped N per atom');

subplot(3,2,6)
%plot(Ta,dNTrdt,'.-')
plot(Ta(2:end),diff(NTr)./diff(Ta),'.-')
ylabel('Trapping rate');
xlabel('Ta (K)');

figure 2
plot(Ta,res,'.-')
ylabel('Resistivity (10^{-11} \Omega m/ppm )');
xlabel('Ta (K)');


A = [Ta' X' F' Nt' Rc'*rat];
A = [Ta' res];
%save('-ascii',fname,'A'); 
save('-ascii','res0','A'); 

##print2pdf(gcf,[20 20],'FeN_meanR_anneal')