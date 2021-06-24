function [x,F,S] = ngdae(t,x0,Xp,Xeq,b0,dG0,R0,incub,solver,dbg)

  if nargin<8,
    error("No enough input arguments");
  end

  if nargin<9 || isempty(solver), 
    solver = "ode15i";
  end

  if nargin<10 || isempty(dbg), 
    dbg=0;
  end
  
  % check initial conditions
  F0 = x0(1)*x0(2)^3;
  X0 = x0(3)*(1-F0) + Xp*F0;
  if F0==0 && X0==0,
    error(['No solutes in initial condition x0=' num2str(x0)]);
  end

  % calc initial xdot
  xdot0 = zeros(size(x0));  
  [xdot0(1), xdot0(2)] = ng(t(1),x0(1),x0(2),xdot0(1),x0(3),Xp,Xeq,b0,dG0,R0,...
    incub,dbg);

  if strcmp(solver,'daspk'),

    dae_func = @ (x, xdot, t) daesystem(t,x,xdot,...
      X0,Xp,Xeq,b0,dG0,R0,incub,dbg);
      
    daspk_options("absolute tolerance",[1e-23; 1e-3; 1e-9]);
    daspk_options("relative tolerance",[1; 1; 1]*1e-9);
    daspk_options("enforce inequality constraints",0);
    daspk_options("inequality constraint types",[1; 1; 1]);
      
    [x, xdot, istate, msg] = daspk (dae_func, x0, xdot0, t);

  elseif strcmp(solver,'ode15i')
  
    options = odeset('RelTol',1e-9,'AbsTol',[1e-23; 1e-3; 1e-9]);
    
    dae_func = @ (t, x, xdot) daesystem(t,x,xdot,...
      X0,Xp,Xeq,b0,dG0,R0,incub,dbg);
    
    [t,x] = ode15i(dae_func,t,x0,xdot0,options);

  else
    error(["Unknown solver:" solver]);
  end
  
  F = x(:,1).*x(:,2).^3;
  S = Xp*log(x(:,3)./Xeq)+(1-Xp)*log((1-x(:,3))./(1-Xeq));

end % ngdae




