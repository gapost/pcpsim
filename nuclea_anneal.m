function [xdot, fsat, S] = nuclea_anneal(x,u,P,k)

xdot = zeros(size(x));
S = P.Xp*log(x(3,:)./P.Xeqs(k)) +(1-P.Xp).*log((1-x(3,:))./(1-P.Xeqs(k)));
F = P.F0(k) .* x(2,:).^3 .* x(1,:);
B = P.S0(k).^2 ./ 2 ./ P.b0(k);

if k < 12,
fsat = 1;
else
fsat = 0;
endif   

y = u/B;
if y<0.005,
  y = B^2/u.^2;
else
  y = xdot(1,:)./ x(1,:);
endif

Rgrowth = (P.a^2/P.R0s(k)^2./x(2,:)) .*((x(3,:) - P.Xeqs(k)*exp(1/P.Xp./x(2,:))) ./ (P.Xp - P.Xeqs(k)*exp(1/P.Xp./x(2,:)))) - fsat*y.*(1.05 ./S - x(2,:));

Ngrowth = (P.b0(k)./S.^2) .*exp(-P.dG0(k)./S.^2) .*exp(-(S.^2)./(2*P.b0(k)*u)); 


xdot(1,:) = fsat*Ngrowth;

xdot(2,:) = Rgrowth;

%Cgrowth = P.F0(k)*(x(3,:) - P.Xp) .* (xdot(1,:).*x(2,:).^3 + 3*x(1,:).*x(2,:).^2 .*xdot(2,:)) ./ (1 - P.F0(k)* x(1,:).*x(2,:).^3);
xdot(3,:) = (x(3,:) - P.Xp) .* F ./ (1-F) .* (3.*xdot(2,:)./x(2,:) + fsat*y );
%xdot(3,:) = Cgrowth;


endfunction