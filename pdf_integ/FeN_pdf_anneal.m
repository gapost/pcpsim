%% Calculate Fe-N annealing with PDF evolution
clear

kb = 8.617e-5; %boltzmann constant [eV/K]

afe = 0.286; % latrice parameter [nm]
Vat = afe^3/2; %Atomic volume of bcc iron [nm^3]
rat = (3*Vat/4/pi)^(1/3); % atomic radius [nm]

gs = 0.062; %surface tension [J/m^2] 0.062
gs *= 6.24150913; % convert to eV/nm2

% Annealing temperatures (K)
Ta = [295 307.80 320.55 333.85 347.70  ...
      362.10 377.10 392.75 409.05 426.00 443.65 462.05 ...
     481 500 520];

nTa = length(Ta);

% annealing time and time grid (step = 20 s)   
dt = 8*60; % s
nt = 48/2+1;
t = linspace(0,dt,nt);

D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; %energy for diffusion [ev]
D =D0*exp(-Qd./(kb*Ta))*1e+18; %diffusion coefficient [nm^2/s]

gam = D/rat^2; 

X0 = 4.67e-4;

Xeq = 10.^(2.43 - 1840./Ta) *1e-2;

Xp = 1/9;

R0 = (2*gs*Vat)./(kb*Ta)/rat; 

S = Xp*log(X0./Xeq)+(1-Xp)*log((1-X0)./(1-Xeq));

Rc = R0./S;

dG0 = (4/3)*pi*R0.^2*rat^2*gs/kb./Ta;

Z = 1/20;

b0 = 4*pi*R0.^2*Z*(rat/afe)^4;

% log R grid 0.5 - 80 nm
R = logspace(log10(0.5/rat),log10(80/rat),101); 
dR = diff([0 R]);
m = length(R);
Rmid = zeros(size(R));
R3 = zeros(size(R));
Rmid(2:m) = (R(1:m-1)+R(2:m))/2;
R3(2:m) = Rmid(2:m).*(R(1:m-1).^2+R(2:m).^2)/2;

f = zeros(nTa,m);

Nt_t = ones(1,nTa*(nt))*NaN;
F_t = ones(1,nTa*(nt))*NaN;
X_t = ones(1,nTa*(nt))*NaN;
Rm_t = ones(1,nTa*(nt))*NaN;
R1_t = ones(1,nTa*(nt))*NaN;

k = 1;
dk = (t(2)-t(1))/60;
cutoff=1;
dbg=0;
incub=0; % incubation time off



for i=1:nTa
  
  if i==1,
    fi = zeros(1,m); 
    Xi = X0;
  else
    fi = f(i-1,:);
    Xi = X(i-1);
  end
  
  
  tic
  [f_t,X_,dfdt,cutoff] = ...
  pdf_integ(fi,Xi,t*gam(i),Xp,Xeq(i),R,R0(i),dG0(i),b0(i),incub,cutoff,dbg);
  cput = toc();
  
  disp(['Ta = ' num2str(Ta(i),'%.1f') ...
  ', cpu = '  num2str(cput,3) ...
  ', cutoff = ' num2str(cutoff)]);
  
  f(i,:) = f_t(end,:);
  X(i) = X_(end);
  
  
  for j=2:nt,  
    Nt_t(k) = (f_t(j,:)*dR')';
    F_t(k) = (f_t(j,:)*(dR.*R3)')';
    Rm_t(k) = (f_t(j,:)*(dR.*Rmid)')'./Nt_t(k);
    X_t(k) = X_(j);
    R1_t(k) = R0(i)/Xp./log(X_t(k)./Xeq(i));
    k += 1;   
  end
  
  figure 4
  x = (1:k)*dk;
  subplot(2,2,1)
  plot(x,Rm_t(1:k)*rat,'.-')
  title('Rm, R* (nm)')
  subplot(2,2,2)
  plot(x,X_t(1:k),'.-')
  title('X')
  subplot(2,2,3)
  plot(x,Nt_t(1:k),'.-')
  title('Ntot')
  xlabel('t (min)')
  subplot(2,2,4)
  plot(x,F_t(1:k),'.-')
  title('F')
  xlabel('t (min)')
  
  figure 3
  subplot(2,1,1)
  semilogx(R,f(i,:),'.-')
  hold on
  semilogx([1 1]*R(cutoff),get(gca,'ylim'),'k--')
  hold off
  title(num2str([Ta(i) X(i) cutoff],3))
  subplot(2,1,2)
  semilogx(R,dfdt(end,:),'.-')
  
  drawnow  
   
end

Nt = (f*dR')';
Rm = (f*(dR.*Rmid)')'./Nt;
F = (f*(dR.*R3)')';
R1 = R0/Xp./log(X./Xeq);

figure 2
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

figure 3
clf
semilogx(R*rat,sqrt(f.*dR),'.-')
xlabel('R (nm)')
ylabel('(f \DeltaR)^{1/2}')



