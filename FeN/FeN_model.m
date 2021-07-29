function [x,F] = FeN_model(t,Ta,X0,gam,incub)
  
% add path to pcpsim functions
addpath("../m");

% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]
 
% Fe lattice data
afe = 0.286e-9; % lattice parameter [m]

% N diffusion in Fe
D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; % migration energy for diffusion [eV]
D =D0*exp(-Qd./(kb*Ta)); %diffusion coefficient [m^2/s]

% alloy data
Xeq = 10.^(2.43 - 1840./Ta) *1e-2; % N conc. in equil. with a" from ???
Xp = 1/9; % precipitate N conc.

% get calculation parameters
[rat,tau,b0,dG0,R0] = ngparam(afe,2,gam,D,Ta);


nt = length(t);
nTa = length(Ta);

if nTa==1, % isothermal
    % initial condition
    S = Xp*log(X0/Xeq) +(1-Xp).*log((1-X0)/(1-Xeq));
    x0 = [0; 1.05*R0/S; X0];

    % call dae solver
    x = ngdae(t/tau,x0,...
      Xp,Xeq,b0,dG0,R0,incub);

else % nTa>1 -> isochronal
  if nt>1 && nt~=nTa,
    error('For isochronal give a single time or a t(1:n) where n=length(Ta)')
  end
  if nt==1,
    t = ones(1,nTa)*t;
  end

  x = zeros(nTa,3);
    
      
    for i=1:nTa
  
      if i==1,
        s = Xp*log(X0./Xeq(1))+(1-Xp)*log((1-X0)./(1-Xeq(1)));
        y = [0 1.05*R0(1)/s X0]';    
      else
        y = x(i-1,:)';
        
        % check the following:
        %   1) we have nucleation (dN/dt>0) 
        %   2) current R is below the new Rc   
        % if yes, delete all nuclei (they are unstable)
        s = Xp*log(y(3)./Xeq(i))+(1-Xp)*log((1-y(3))./(1-Xeq(i)));
        jn = y(3)*b0(i)/s^2*exp(-dG0(i)/s^2);
        if y(2) < 1.05*R0(i)/s && jn>0,
          y = [0 1.05*R0(i)/s X0]'; 
        end
      end

      y = ngdae(linspace(0,t(i)/tau(i),21),y,...
        Xp,Xeq(i),b0(i),dG0(i),R0(i),incub);
      
      x(i,:) = y(end,:);
      
      
    end
  end

  F = x(:,1).*x(:,2).^3;