function [xdot, F, S] = nuclea_anneal(x,u,P,k)

xdot = zeros(size(x));

S = P.Xp*log(x(3,:)./P.Xeqs(k)) +(1-P.Xp).*log((1-x(3,:))./(1-P.Xeqs(k)));

F = P.F0(k) * x(2,:).^3 .* x(1,:);

B = P.S0(k).^2 / 2 / P.b0(k);

y = xdot(1,:)./ x(1,:);

xdot(1,:) = (P.b0(k)./S.^2) .*exp(-P.dG0(k)./S.^2) .*exp(-(S.^2)./(2*P.b0(k)*u)); 

xdot(2,:) = (P.a^2/P.R0s(k)^2./x(2,:)) .*((x(3,:) - P.Xeqs(k)*exp(1/P.Xp./x(2,:))) ./ (P.Xp - P.Xeqs(k)*exp(1/P.Xp./x(2,:)))) - y.*(1.05 ./S - x(2,:));

xdot(3,:) = (x(3,:) - P.Xp) .* F ./ (1-F) .* (3*xdot(2,:)./x(2,:) + y );


endfunction