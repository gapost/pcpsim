function [x,F,R] = FeNirrad_model(t,Ta,X0,XVN0,gam,incub,Eb)
  
% add path to pcpsim functions
addpath("../m");

% turn debugging off
dbg = 0;

% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]
 
% Fe lattice data
afe = 0.286e-9; % lattice parameter [m]

% N diffusion in Fe
D0 = 1.26e-7; %Pre-exp. diffusion [m^2/s]
Qd = 0.76; % migration energy for diffusion [eV]
Dn = D0*exp(-Qd./(kb*Ta)); %diffusion coefficient [m^2/s]

% V diffusion in Fe
D0v = 1e-6; %Pre-exp. diffusion [m^2/s] Domain 2004 pg5
Em = 0.65; % migration energy for diffusion [eV] Domain 2004 pg5
Dv = D0v*exp(-Em./(kb*Ta)); %diffusion coefficient [m^2/s]

% alloy data
Xeq = 10.^(2.43 - 1840./Ta) *1e-2; % N conc. in equil. with a" from ???
Xp = 1/9; % precipitate N conc.

% get nucleation & growth parameters
[rat,tau,b0,dG0,R0] = ngparam(afe,2,gam,Dn,Ta);

% calculate reaction rate coefficients
rd1 = 1.2*1.7*afe*1e9; % reaction distance of V + N = VN [nm]
rd2 = 1.2*1.9*afe*1e9; % reaction distance of VN + N = VN2 [nm]
rd3 = 3.3*afe*1e9;     % reaction distance of V + In = In [nm]

k = zeros(5,length(Ta));
k(1,:) = 3*rd1/rat*(Dv./Dn+1); % rate coef of V + N -> VN
k(2,:) = k(1,:).*exp(-Eb(1)/kb./Ta); % rate coef of VN -> V + N
k(3,:) = 3*rd2/rat; % rate coef of VN + N -> VN2
k(4,:) = k(3,:).*exp(-Eb(2)/kb./Ta); % rate coef of VN2 -> VN + N
k(5,:) = 3*rd3/rat*Dv./Dn;  % rate coef of V + In -> In

nt = length(t);
nTa = length(Ta);

%x(1) = precipitates, x(2) = Radius, 
%x(3) = N, x(4)=V, x(5)=VN, x(6)=VN2, x(7) = In

% initial condition
s = Xp*log(X0/Xeq(1)) +(1-Xp).*log((1-X0)/(1-Xeq(1)));
x0 = [0 1.05*R0(1)/s X0 0 XVN0 0 XVN0/4]';

if nTa==1, % isothermal
    
    % call dae solver
    x = ngdae_irrad(t/tau,x0,...
      Xp,Xeq,b0,dG0,R0,incub,k,dbg);

else % nTa>1 -> isochronal
  if nt>1 && nt~=nTa,
    error('For isochronal give a single time or a t(1:n) where n=length(Ta)')
  end
  if nt==1,
    t = ones(1,nTa)*t;
  end

  x = zeros(nTa,7);
  R = zeros(nTa,3);  
      
  for i=1:nTa

    if i==1,
      y = x0;  
    else
      y = x(i-1,:)';
      
      % check the following:
      %   1) we have nucleation (dN/dt>0) 
      %   2) current R is below the new Rc   
      % if yes, delete all nuclei (they are unstable)
      s = Xp*log(y(3)./Xeq(i))+(1-Xp)*log((1-y(3))./(1-Xeq(i)));
      jn = y(3)*b0(i)/s^2*exp(-dG0(i)/s^2);
      if y(2) < 1.05*R0(i)/s && jn>0,
        y(3) = y(3) + y(1)*y(2)^3*Xp; % dissolve all N from precipitates
        y(2) = 1.05*R0(i)/s;
        y(1) = 0; 
      end
    end

    if dbg,
      disp(['Ta = ' num2str(Ta(i))])
    end
    
    y = ngdae_irrad(linspace(0,t(i)/tau(i),41),y,...
      Xp,Xeq(i),b0(i),dG0(i),R0(i),incub,k(:,i),dbg);
      
    x(i,:) = y(end,:);

  end
end

F = x(:,1).*x(:,2).^3;
R = zeros(size(x,1),3);
R(1,:) = reaction_rates(x0,k(1,:))/tau(1);
for i=2:nTa,
  R(i,:) = reaction_rates(x(i-1,:),k(:,i))/tau(i);
end
  
end % function FeNirrad_model

function x = ngdae_irrad(t,x0,Xp,Xeq,b0,dG0,R0,incub,k,dbg)

  % x(1) = precipitates, x(2) = Radius, 
  % x(3) = N, x(4)=V, x(5)=VN, x(6)=VN2, x(7) = In
    
  % check initial conditions
  F0 = x0(1)*x0(2)^3;
  X0 = (x0(3) + x0(5) + 2*x0(6))*(1-F0) + Xp*F0;
  if F0==0 && X0==0,
      error(['No solutes in initial condition x0=' num2str(x0)]);
  end

  % calc initial xdot & rates
  xdot0 = daesystem_irrad(t(1),x0,zeros(size(x0)),...
    X0,Xp,Xeq,b0,dG0,R0,incub,k,dbg);
  
  options = odeset('RelTol',1e-6,...
    'AbsTol',[1e-23; 1e-3; 1e-9; 1e-9; 1e-9; 1e-9; 1e-9]);
    
  if dbg,
    odeset(options,'Stats','on');
    disp( '===> Starting ode15i,')
    disp(['  t = ' num2str(t(1)) '..' num2str(t(end))])
    disp(['  x0 = ' num2str(x0')])
    disp(['  xdot0 = ' num2str(xdot0')])
  end
  
  dae_func = @ (t, x, xdot) daesystem_irrad(t,x,xdot,...
    X0,Xp,Xeq,b0,dG0,R0,incub,k,dbg);
    
        
  [t,x] = ode15i(dae_func,t,x0,xdot0,options);
 
end % ngdae_irrad

function ret = daesystem_irrad(t,x,xdot,X0,Xp,Xeq,b0,dG0,R0,incub,k,dbg)
    %x(1) = precipitates, x(2) = Radius, 
    %x(3) = N, x(4)=V, x(5)=VN, x(6)=VN2, x(7) = In  
  
  [dNdt, dRdt] = ng(t,x(1),x(2),xdot(1),x(3),Xp,Xeq,b0,dG0,R0,incub,dbg);
  
  F = x(1)*x(2)^3; 
 
  R = reaction_rates(x,k);
    
  ret = [-xdot(1) + dNdt;        % Np
         -xdot(2) + dRdt;        % Rp
         -X0 + (x(3) + x(5) + 2*x(6))*(1-F) + Xp*F; % N
         -xdot(4) - R(1) - R(3); % V
         -xdot(5) + R(1) - R(2); % VN
         -xdot(6) + R(2);        % VN2
         -xdot(7) - R(3)/4];     % In
         
  if dbg, disp(num2str([t x' xdot' ret'],3)); end
  
end

function R = reaction_rates(x,k)
  R = [k(1)*x(3)*x(4) - k(2)*x(5), ... % N + V  -> VN 
       k(3)*x(3)*x(5) - k(4)*x(6), ... % N + VN -> VN2
       k(5)*x(4)*x(7)];                % V + In -> In-1
end
