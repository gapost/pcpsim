function [f, X, dfdt, cutoff] = pdf_integ(f0,X0,t,Xp,Xeq,R,R0,dG0,b0,incub,...
  cutoff0,dbg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f, X, dfdt, cutoff] = pdf_integ(f0,X0,t,Xp,Xeq,R,R0,dG0,b0,incub,...
%                                      cutoff0,dbg)
%
% Integrate the precipitate distribution function (PDF) time evolution during 
% nucleation & growth or dissolution. 
%
% Input:
%  f0(1,M) : initial distribution on M grid points of precipitate radius
%  X0      : initial solute concentration in the matrix_type
%  t(1,N)  : vector of N time values where the solution is to be returned. t(1) 
%            is the starting point t0 where f=f0 
%  Xp, Xeq : solute conc. in the precipitate and in the matrix at equilibrium
%  R(1,M)  : the grid of precipitate radii
%  R0,dG0,b0 : nucleation & growth physical parameters
%  incub   : if 1 then incubation time is calculated
%  cutoff0 : initial lower cuttoff index
%  dbg     : if 1 turn on debugging
%
% Output
%   f(N,M) : each row corresponds to the PDF at the corresponding t(i)
%   X(N)   : colute concentration at each time point
%   dfdt(N,M) : time difference of the PDF at t(i)
%   cutoff : the new cutoff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  nt = length(t);  
  m = size(R,2);
  cutoff = cutoff0;
  
  dR = diff([0 R]);
  R3 = zeros(size(R)); 
  R3(2:m) = dR(2:m).*(R(1:m-1)+R(2:m)).*(R(1:m-1).^2+R(2:m).^2)/4;
  
  f = zeros(nt,m);
  dfdt = zeros(nt,m);
  X = zeros(1,nt);
  
  f(1,:) = f0;
  X(1) = X0;
  
  for i=2:nt,
    [f(i,:) X(i) dfdt(i,:) cutoff] = pdf_step(f(i-1,:),X(i-1),[t(i-1) t(i)],...
                                 Xp,Xeq,R0,dG0,R,dR,R3,b0,incub,cutoff,dbg);
  endfor

endfunction


% time step for the PDF integration
% Integrates from t(1) to t(2)
function [f, X, dfdt, cutoff] = pdf_step(f0,X0,t,Xp,Xeq,...
  R0,dG0,R,dR,R3,b0,incub,cutoff0,dbg)
    
  m = size(R,2);
  cutoff = cutoff0;
  f = f0;
  X = X0;
  F = R3*f0';
  ti = t(1); % ti is the current time
  
  dfdt = zeros(size(f0));
  vr = zeros(size(R));
  
  do 
    % calculate S and critical radius (simplified)
    S = Xp*log(X./Xeq)+(1-Xp)*log((1-X)./(1-Xeq));
    R1 = R0/Xp/log(X/Xeq);
    
    % calculate nucleation current
    Js = 0;
    is = 1;
    if S>0,
      b = b0*X;
      Js = b/S^2 * exp(-dG0/S^2);
      is = lookup(R,1.06*R1,'rl')+1; % the nuclei are added at f(is)
      if incub,
        tau = S^2/2/b;
        if ti>0, Js *= exp(-tau/ti); else Js=0; end
      end
    end

    % call the PDF difference equations        
    [dfdt, vr, J] = pdf_diff(f,X,Xp,Xeq,R0,R,dR,Js,is,cutoff);

    % estimate dt    
    i = find(f>0);
    if isempty(i), i=2:m; end
    dt = min(0.5*abs(dR(i)./vr(i))); % Courant stability criterion
    dt = min([dt t(2)-ti]); % do not extend beyond t(2)
    
    % check for cutoff changes
    if J(cutoff)<0 && abs(J(cutoff)*dt)>(1e-23 + 1e-3*(f*dR')),
      % decrease the cutoff - dissolution
      dfdt(cutoff) = -J(cutoff)/dR(cutoff);
      if cutoff>1, cutoff -= 1; end
    else
      if (f(cutoff+1)*dR(cutoff+1))<1e-23 && cutoff<m-2,
        cutoff += 1;
      endif
    endif

    % zero out
    f(1:cutoff)=0;
    dfdt(1:cutoff)=0;

    % re-estimate dt
    i = find(f>0);
    if isempty(i), i=2:m; end
    dt = min(0.5*abs(dR(i)./vr(i)));
    dt = min([dt t(2)-ti]);

    % plotting for debugging    
    if dbg 
      clf
      subplot(2,1,1)
      semilogx(R,f,'.-')
      title(num2str([dt/(t(2)-t(1)) (ti+dt-t(1))/(t(2)-t(1)) X Js cutoff],3))
      subplot(2,1,2)
      semilogx(R,dfdt,'.-')
      drawnow
    end
    
    % advance the solution
    dfdt *= dt;
    f = f + dfdt;
    dF = R3*dfdt';
    X += (X-Xp)*dF/(1-F);
    F += dF;
    ti += dt;
    
  % check if we are finished
  until ti>=t;
    
endfunction

% calculate the PDF finite difference equations
function [dfdt, vr, J] = pdf_diff(f,X,Xp,Xeq,R0,R,dR,Js,is,cutoff)
  
  dfdt = zeros(size(f)); 
  vr = zeros(size(R));
  J = zeros(size(R));
  m = size(R,2);

  % concentration at precipitate radius  
  Xr = Xeq.*exp(R0/Xp./R);
  % impose upper bound Xr < 0.5 (Xp + X)
  Xmax = 0.5*(Xp+X);
  i = find(Xr>Xmax);
  Xr(i) = Xmax;

  % growth rate  
  vr = (X-Xr)./(Xp-Xr)./R;
  
  % current
  J = vr.*f;
 
  if any(vr<0),
    % in case there is dissolution (vr<0) at some points
    % the current is calculated differently    
    i = find(vr<0);
    if i(end)==m, i(end)=[]; end   
    J(i) = f(i+1).*vr(i); % i < m
  end

  % observe the boundary condition
  if vr(cutoff)>0, J(cutoff)=0; end
 
  % only the points above cutoff are considered 
  dfdt(cutoff+1:end) = -J(cutoff+1:end)+J(cutoff:end-1);
 
  % add the nucleation  
  dfdt(is) += Js;

  % scale all to dR 
  dfdt = dfdt./dR;
  
endfunction