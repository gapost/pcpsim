function [xdot, F, S] = nuclea_v2(x,u,P)

xdot = zeros(size(x));
S = P.Xp*log(x(3,:)./P.Xeqs) +(1-P.Xp).*log((1-x(3,:))./(1-P.Xeqs));

F = P.F0 * x(2,:).^3 .* x(1,:);

A.a = S.^2;
A.b = A.a./P.b0;
A.c = P.a^2/P.R0s^2;
A.d = exp(1/P.Xp./x(2,:));
A.e = (x(3,:) - P.Xeqs.*A.d) ./ (P.Xp - P.Xeqs*A.d);


y = xdot(1,:)./ x(1,:);

xdot(1,:) = (1./A.b) .*exp(-P.dG0./A.a) .*exp(-A.b./(2*u)); 

xdot(2,:) = (A.c./x(2,:)) .* A.e - y.*(1.05 ./S - x(2,:));

xdot(3,:) = (x(3,:) - P.Xp) .* F ./ (1-F) .* (3*xdot(2,:)./x(2,:) + y );


endfunction