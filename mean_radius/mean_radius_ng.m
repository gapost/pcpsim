function [xdot, F, S] = mean_radius_ng(x,t,Xp,Xeq,b0,dG0,R0,incub,dbg)

  N = x(1); R = x(2); X = x(3);
  if N<0,
    error(['Negative N value, N=' num2str(N)]);
  endif
  if X<0,
    error(['Negative X value, X=' num2str(X)]);
  endif

  xdot = zeros(size(x));

  S = Xp*log(X/Xeq) +(1-Xp).*log((1-X)/(1-Xeq));
  F = R^3 * N;
  
  if S>0, % nucleation & growth
    if R<0,
      error(['Negative R in nucleation, R=' num2str(R)]);
    endif
    S2 = S^2;
    b = X * b0 / S2;
    y=0;
    xdot(1) = b .*exp(-dG0./S2);
    if incub,
      if t>0,
        xdot(1) = xdot(1) .*exp(-1./2/b/t);
        if N>1e-23,
          y = xdot(1,:)/N;
        else
          y = 0;
        endif
      end
    else
      if N>0,
        y=xdot(1)/N;
      endif
    endif
    Xr = Xeq*exp(R0./R/Xp);
    Xr=min(0.9*Xp,Xr);
    xdot(2) = (X-Xr) / (Xp-Xr) / R + y*(1.05*R0./S - R);
    if xdot(2)<0, xdot(2)=0; end
    xdot(3) = (X - Xp) * F / (1-F) * (3*xdot(2)/R + y );    
  else % dissolution
    if R>0, 
      Xr = Xeq*exp(R0./R/Xp);
      Xr=min(0.9*Xp,Xr);
      xdot(2) = (X-Xr) / (Xp-Xr) / R;
      xdot(3) = (X - Xp) * F / (1-F) * 3 * xdot(2) / R;
    endif
  endif


  if dbg, disp(num2str([t x' xdot'])); end

endfunction 