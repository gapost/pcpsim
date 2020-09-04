%% Calculate Fe-N annealing with PSD evolution
clear

% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]

% annealing conditions
% Annealing temperatures (K)
Ta = [295 307.80 320.55 333.85 347.70  ...
      362.10 377.10 392.75 409.05 426.00 443.65 462.05 ...
     481 500 520];
nTa = length(Ta);
dt = 5*60; % annealing time in s

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
gs = 0.0555; %surface tension [J/m^2]
X0 = 4.67e-4; % Nominal N concentration
Xeq = 10.^(2.43 - 1840./Ta) *1e-2; % solubility
Xp = 1/9; % precipitate

% Options
incub=1; % Calc. incubation time for nucleation
dbg=0; % debug level (0: off, 1: messages, 2: msg+plots, 3: msg+more plots)

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

% log R grid 0.5 - 100 nm, in units of rat
R = logspace(log10(0.5/rat),log10(120/rat),101); 


m = length(R);
dR = diff([0 R]);
dR2 = 0.5*diff([0 R.^2]);
dR4 = 0.25*diff([0 R.^4]);

f = zeros(nTa,m);
dfdt = zeros(nTa,m);
cutoff=1;

Nt_t = ones(1,nTa*nt)*NaN;
F_t = ones(1,nTa*nt)*NaN;
X_t = ones(1,nTa*nt)*NaN;
Rm_t = ones(1,nTa*nt)*NaN;
R1_t = ones(1,nTa*nt)*NaN;

k = 1;

for i=1:nTa

  clear fi dfdti Xi
  
  if i==1,
    fi = zeros(1,m); 
    Xi = X0;
  else
    fi = f(i-1,:);
    Xi = X(i-1);
  end
  
  
  tic
  [fi,Xi,dfdti,cutoff,iter] = ...
  psd_integ(fi,Xi,t*gam(i),Xp,Xeq(i),R,R0(i),dG0(i),b0(i),incub,cutoff,dbg);
  cput = toc();
  
  
  disp(['Ta = ' num2str(Ta(i),'%.1f') ...
  ', cpu = '  num2str(cput,3) ...
  ', iter = '  num2str(iter) ...
  ', cutoff = ' num2str(cutoff)]);
  
  f(i,:) = fi(end,:);
  dfdt(i,:) = dfdti(end,:);
  X(i) = Xi(end);
  
  
  for j=2:nt, 
    k += 1; 
    Nt_t(k) = (fi(j,:)*dR')';
    F_t(k) = (fi(j,:)*dR4')';
    Rm_t(k) = (fi(j,:)*dR2')'./Nt_t(k);
    X_t(k) = Xi(j);
    R1_t(k) = R0(i)/Xp./log(Xi(j)./Xeq(i));    
  end
  
end

Nt = (f*dR')';
Rm = (f*dR2')'./Nt;
F = (f*dR4')';
R1 = R0/Xp./log(X./Xeq);

figure 1
clf

subplot(2,2,1)
plot(Ta,Rm*rat,'.-',Ta,R1*rat,'.-')
ylabel('R, R* (nm) ');

subplot(2,2,2)
plot(Ta,Nt,'.-')
ylabel('Clusters per atom');

subplot(2,2,3)
plot(Ta,X/X0,'.-')
ylabel('Solute mole fraction');
xlabel('Ta (K)');

subplot(2,2,4)
plot(Ta,F,'.-')
ylabel('Transformed volume fraction ');
xlabel('Ta (K)');

figure 2
clf
idx=3:length(Ta);
plot_psd

figure 3
clf
x = (1:nt*nTa)*(t(2)-t(1))/60; % continous annealing time
subplot(2,2,1)
plot(x,Rm_t*rat,'.-')
title('Rm, R* (nm)')
subplot(2,2,2)
plot(x,X_t,'.-')
title('X')
subplot(2,2,3)
plot(x,Nt_t,'.-')
title('Ntot')
xlabel('t (min)')
subplot(2,2,4)
plot(x,F_t,'.-')
title('F')
xlabel('t (min)')

A = [Ta' X' F' Nt' Rm'*rat];
save -ascii FeN555.dat A 





