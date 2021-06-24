function [dNdt, dRdt] = ng(t,N,R,dNdt,X,Xp,Xeq,b0,dG0,R0,incub,dbg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [dNdt, dRdt, S] = ng(t,N,R,X,Xp,Xeq,b0,dG0,R0,incub,dbg)
%
% Define the ODEs describing the time evolution of precipitate concentration
% and mean radius during homogeneous nucleation and growth 
%
% Input:
%  t       : time (in units D*t/rat^2)
%  N       : precipitate concentration (atomic) 
%  R       : mean radius R (in units of rat)
%  dNdt    : time derivative of N 
%  X, Xp, Xeq : solute concentration in the matrix, precipitate and 
%               matrix at equilibrium
%  R0,dG0,b0 : nucleation & growth physical parameters
%  incub   : if 1 then incubation time is calculated
%  dbg     : if 1 turn on debugging
%
% Output
%  dNdt    : time derivative of N
%  dRdt    : time derivative of R
%  S       : nucleation entropy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  
  % check input 
  if N<0 && abs(N)>1e-23, 
    error(['Negative N value, N=' num2str(N)]);
  end
  if X<0,
    error(['Negative X value, X=' num2str(X)]);
  end
  if R<0 && abs(R)>1e-3,
    warning(['Negative R value, R=' num2str(R)]);
  end

  % nucleation entropy 
  S = Xp*log(X/Xeq) +(1-Xp).*log((1-X)/(1-Xeq)); 
  
  % solute conc. at prec./matrix interface
  Xr = Xeq*exp(R0./R*(1-Xeq)/(Xp-Xeq));
  
  % precipite radius growth rate
  if Xr<0.5*Xp,
    dRdt = (X-Xr) / (Xp-Xr) / R * (R>0);
  else
    Rp = R0*(1-Xeq)/(Xp-Xeq)/log(Xp/2/Xeq);
    dRdt = (X-0.5*Xp)/0.5/Xp/Rp*(R/Rp)*(R>0);
  end

  % reduction of R due to new nuclei
  Rs = 1.05*R0/S;
  if N>1e-23 && R>Rs, % only if appreciable N > 1 per cm^3
    dRdt = dRdt - dNdt/N*(R-Rs);
  end
 
  % compute nucleation rate
  dNdt = 0;  
  if S>0, % nucleation 
    S2 = S^2;
    b = X * b0 / S2;
    dNdt = b .*exp(-dG0./S2); % nucleation rate, no incubation
    if incub, % with incubation
      if t>0,
        dNdt = dNdt * exp(-1/2/b/t);
      else % for t=0, nucl rate = 0
        dNdt = 0;
      end
    end
  end

end