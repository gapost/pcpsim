function [f, X, dfdt, cutoff, total_iter] = psd_integ(f0,X0,t,Xp,Xeq,R,R0,dG0,b0,incub,...
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
  
  % concentration at precipitate radius  
  Xr = Xeq.*exp(R0/Xp./R);
  
  if dbg>0,
    printf('step  \titer  \tkc     \tX/Xeq-1  \tF      \tN      \n');
  end
  
  total_iter = 0;
  
  for i=2:nt,
    [f(i,:) X(i) dfdt(i,:) cutoff iter] = psd_step(f(i-1,:),X(i-1),[t(i-1) t(i)],...
                                 Xr,Xp,Xeq,R0,dG0,R,dR,R3,b0,incub,cutoff,dbg);
                                 
    total_iter += iter;
                                 
    if dbg>0,
      printf('%d/%d\t%d\t%d\t%f\t%f\t%e\n',...
        [i nt iter cutoff X(i)/Xeq-1 f(i,:)*R3' f(i,:)*dR']);
    end
    
    if dbg==2,
      R1 = R0/Xp./log(X(i)/Xeq);
      clf
      subplot(2,1,1)
      semilogx(R,f(i,:),'.-')
      hold on
      yy=get(gca,'ylim');
      semilogx([1 1]*R1,[yy(1) 0.9*yy(2)],'k--')
      hold off
      title(num2str([i X(i) cutoff],3))
      subplot(2,1,2)
      semilogx(R,dfdt(i,:),'.-')
      drawnow
    end
    
  endfor

endfunction


% time step for the PSD integration
% Integrates from t(1) to t(2) with adaptive step
% iter = # of iterations/steps
function [f, X, dfdt, cutoff, iter] = psd_step(f0,X0,t,Xr,Xp,Xeq,...
  R0,dG0,R,dR,R3,b0,incub,cutoff0,dbg)
    
  m = size(R,2);
  cutoff = cutoff0;
  f = f0;
  X = X0;
  F = R3*f0';
  X00 = X*(1-F)+Xp*F; % total solute concentration
  ti = t(1); % ti is the current time
  
  dfdt = zeros(size(f0));
  vr = zeros(size(R));
  vrq = 1./(Xp-Xr)./R;
  
  iter = 0;
  
  do 
    % calculate S, critical radius (simplified), vr
    S = Xp*log(X./Xeq)+(1-Xp)*log((1-X)./(1-Xeq));
    R1 = R0/Xp/log(X/Xeq);
    vr = (X-Xr).*vrq;
    
    % calculate nucleation current
    Js = 0;
    is = 1;
    if S>0,
      b = b0*X;
      Js = b/S^2 * exp(-dG0/S^2);
      is = lookup(R,R1,'rl')+2; % the nuclei are added at f(is)
      if is>m, is=m; end
      if incub,
        tau = S^2/2/b;
        if ti>0, Js *= exp(-tau/ti); else Js=0; end
      end
    end

    % call the PSD difference equations        
    [dfdt, J] = psd_diff(f,vr,dR,Js,is,cutoff);

    % estimate adaptive dt    
    i = find(f.*dR>1e-23);
    if isempty(i), 
      i=cutoff+1:m; 
      dt = min(0.5*abs(dR(i)./vr(i))); % Courant stability criterion
      dt = min([dt t(2)-ti]); % do not extend beyond t(2)
    else 
      dt = min(0.5*abs(dR(i)./vr(i))); % Courant stability criterion
      dt = min([dt t(2)-ti]); % do not extend beyond t(2)
      
      % check for cutoff changes
      cutoff_changed = 0;
      if J(cutoff)<0 && abs(J(cutoff)*dt)>(1e-23 + 1e-3*(f*dR')),
        % decrease the cutoff - dissolution
        if cutoff>1, 
          dfdt(cutoff) = -J(cutoff)/dR(cutoff);
          cutoff -= 1; 
          cutoff_changed = 1;         
        end
      else
        if (f(cutoff+1)*dR(cutoff+1))<1e-23 && cutoff<m-2,
          cutoff += 1;
          cutoff_changed = 1;
        endif
      endif
      
      % if cutoff changed, re-estimate dt
      if cutoff_changed,
        % zero out
        f(1:cutoff)=0;
        dfdt(1:cutoff)=0;
        
        % re-estimate dt
        i = find(f>0);
        dt = min(0.5*abs(dR(i)./vr(i)));
        dt = min([dt t(2)-ti]);
      end
    end

    % check is particles flow past the last grid point    
    if J(end)*dt > 1e-23,
      error('particles flow past Rmax')
    end

    % plotting for debugging    
    if dbg==3, 
      clf
      subplot(2,1,1)
      semilogx(R,f,'.-')
      title(num2str([dt/(t(2)-t(1)) (ti+dt-t(1))/(t(2)-t(1)) X Js cutoff],3))
      subplot(2,1,2)
      semilogx(R,dfdt,'.-')
      drawnow
    end

    % 2nd criterion
    % check if R* changes more than 1%   
    F = R3*(f + dfdt*dt)';
    X = (X00-Xp*F)/(1-F);    
    if R1>0, % only if R* exists - i.e. supersaturation
      dRrel = abs(R0/Xp/log(X/Xeq)/R1-1);     
      while dRrel>0.01,
        dt *= 0.5;
        F = R3*(f + dfdt*dt)';
        X = (X00-Xp*F)/(1-F);
        dRrel = abs(R0/Xp/log(X/Xeq)/R1-1);
      endwhile
    endif
       
    % advance the solution
    f = f + dfdt*dt;
    ti += dt;
    iter++;
    
  % check if we are finished
  until ti>=t || f*dR'<1e-23;
    
endfunction

% calculate the PDF finite difference equations
function [dfdt, J] = psd_diff(f,vr,dR,Js,is,cutoff)
  
  dfdt = zeros(size(f)); 
  J = zeros(size(vr));
  m = size(vr,2);

  % current
  J = vr.*f;
 
  if any(vr<0),
    % in case there is dissolution (vr<0) at some points
    % the current is calculated differently    
    i = find(vr<0);
    if i(end)==m, i(end)=[]; end   
    J(i) = f(i+1).*vr(i); % i < m
  end

  % observe the low boundary condition
  if vr(cutoff)>0, J(cutoff)=0; end
 
  % only the points above cutoff are considered 
  dfdt(cutoff+1:end) = -J(cutoff+1:end)+J(cutoff:end-1);
 
  % add the nucleation  
  if is>cutoff,
    dfdt(is) += Js;
  end

  % scale all to dR 
  dfdt = dfdt./dR;
  
endfunction