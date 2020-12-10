function [xdot, S, F] = ng(x,u,P)

xdot = zeros(size(x));



S = P.Xp*log(x(3,:)./P.Xeqs) +(1-P.Xp).*log((1-x(3,:))./(1-P.Xeqs));

F = P.h * x(2,:).^3 .* exp(x(1,:));

xdot(1,:) = (P.b0./S.^2) .*exp(- P.dG0./S.^2 - (S.^2)./(2*P.b0*u) - x(1,:)); 

xdot(2,:) = (P.a^2/P.R0s^2./x(2,:)) .*((x(3,:) - P.Xeqs*exp(1/P.Xp./x(2,:))) ...
  ./ (P.Xp - P.Xeqs*exp(1/P.Xp./x(2,:)))) ...
  - xdot(1,:).*(1.05 ./S - x(2,:));
  

% xdot(3,:) = P.h*(x(3,:) - P.Xp) .* (xdot(1,:).*x(2,:).^3 + 3*x(1,:).*x(2,:).^2 .*xdot(2,:)) ./ (1 - P.h* x(1,:).*x(2,:).^3);

xdot(3,:) = (x(3,:) - P.Xp) .* F ./ (1-F) .* (3*xdot(2,:)./x(2,:) + xdot(1,:) );

endfunction