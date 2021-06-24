function ret = daesystem(t,x,xdot,X0,Xp,Xeq,b0,dG0,R0,incub,dbg)
  [dNdt, dRdt] = ng(t,x(1),x(2),xdot(1),x(3),Xp,Xeq,b0,dG0,R0,incub,dbg);
  F = x(1)*x(2)^3;
  ret = [xdot(1) - dNdt; 
         xdot(2) - dRdt; 
         X0 - x(3)*(1-F) - Xp*F];
  if dbg, disp(num2str([t x' xdot' ret'],3)); end
end