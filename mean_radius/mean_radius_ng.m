function [xdot, F, S] = mean_radius_ng(t,x,Xp,Xeq,b0,dG0,R0,incub,dbg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xdot, F, S] = mean_radius_ng(t,x,Xp,Xeq,b0,dG0,R0,incub,dbg)
%
% Define the ODEs describing mean precipitate radius during nucleation and
% growth 
%
% Input:
%  t       : time (D*t/rat^2) 
%  x(1,3)  : ODE variables, 
%            x(1): precipitate atomic concentration N, 
%            x(2): mean radius R (in units of rat), 
%            x(3): atomic concentration of solute in matrix X
%  Xp, Xeq : solute conc. in the precipitate and in the matrix at equilibrium
%  R0,dG0,b0 : nucleation & growth physical parameters
%  incub   : if 1 then incubation time is calculated
%  dbg     : if 1 turn on debugging
%
% Output
%   xdot   : dx/dt
%   F      : precipitate volume fraction
%   S      : nucleation entropy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  N = x(1); R = x(2); X = x(3);
  if N<0, 
    error(['Negative N value, N=' num2str(N)]);
  endif
  if X<0,
    error(['Negative X value, X=' num2str(X)]);
  endif
  
  xdot = zeros(size(x));
  y=0; % (1/N)dN/dt = xdot(1)/x(1)
  
  S = Xp*log(X/Xeq) +(1-Xp).*log((1-X)/(1-Xeq)); % nucleation entropy
  F = R^3 * N; % precipitate volume fraction
  
  % solute conc. at prec./matrix interface
  Xr = Xeq*exp(R0./R*(1-Xeq)/(Xp-Xeq));
  
  % precipite radius growth rate
  if Xr<0.5*Xp,
    Vr = (X-Xr) / (Xp-Xr) / R * (R>0);
  else
    Xr = 0.5*Xp;
    Rp = R0*(1-Xeq)/(Xp-Xeq)/log(Xp/2/Xeq);
    Vr = (X-Xr) / (Xp-Xr) / Rp * (R/Rp);
  end
  
  if S>0, % nucleation 
    S2 = S^2;
    b = X * b0 / S2;
    xdot(1) = b .*exp(-dG0./S2); % nucleation rate, no incubation
    if incub, % with incubation
      if t>0,
        xdot(1) = xdot(1) .*exp(-1./2/b/t);
        if N>1e-23, % only if appreciable N > 1 per cm^3
          y = xdot(1,:)/N;
        else
          y = 0;
        endif
      else % for t=0, nucl rate = 0
        xdot(1) = 0;
      end
    else % no incubation
      if N>0,
        y=xdot(1)/N;
      endif
    endif
  end
  
  % calc dR/dt & dX/dt only for R>0
  % R<0 means precipitates have dissolved 
  if R>0, 
    xdot(2) = (X-Xr) / (Xp-Xr) / R;
    xdot(3) = (X - Xp) * F / (1-F) * 3 * xdot(2) / R;
  endif 
  
  if dbg, disp(num2str([t x' xdot'])); end
    
endfunction 