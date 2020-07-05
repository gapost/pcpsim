function [xdot, fcoars, S] = nuclea_coars_calder_anneal(x,u,P,k)

xdot = zeros(size(x));
S = P.Xp*log(x(3,:)./P.Xeqs(k)) +(1-P.Xp).*log((1-x(3,:))./(1-P.Xeqs(k)));

%fcoars  = 1 - erf(4.*(x(2,:).*AS - 1)) ;
fcoars  = 1 - erf(4.*(x(3,:)./P.Xeqs(k).*x(2,:).*S - 1)) ;

B = P.S0(k).^2 / 2 / P.b0(k);

y = u/B;
if y<0.005,
  y = B^2/u.^2;
else
  y = xdot(1,:)./ x(1,:);
endif

Rcoars = (4./27) .* (P.Xeqs(k)*(1-P.Xeqs(k))./(P.Xp-P.Xeqs(k)).^2) .* (P.a^2./(P.R0s(k).^2 .* x(2,:).^2));
Rgrowth = (P.a^2/P.R0s(k)^2./x(2,:)) .*((x(3,:) - P.Xeqs(k)*exp(1/P.Xp./x(2,:))) ./ (P.Xp - P.Xeqs(k)*exp(1/P.Xp./x(2,:)))) - y.*(1.05 ./S - x(2,:));

Ngrowth = (P.b0(k)./S.^2) .*exp(-P.dG0(k)./S.^2) .*exp(-(S.^2)./(2*P.b0(k)*u)); 
Ncoars = (3/4)*P.F0(k)./(4*pi*P.R0s(k).^3.*x(2,:).^4) .*xdot(2,:);

xdot(1,:) = (1 - fcoars).*Ngrowth + fcoars.*Ncoars;
xdot(2,:) = (1 - fcoars).*Rgrowth + fcoars.*Rcoars;

Ccoars = P.Xeqs(k)*(1-P.Xeqs(k))./(P.Xeqs(k) - P.Xp) ./ x(2,:).^2 .*xdot(2,:);
Cgrowth = P.F0(k)*(x(3,:) - P.Xp) .* (xdot(1,:).*x(2,:).^3 + 3*x(1,:).*x(2,:).^2 .*xdot(2,:)) ./ (1 - P.F0(k)* x(1,:).*x(2,:).^3);

xdot(3,:) = (1 - fcoars).*Cgrowth + fcoars.*Ccoars;
endfunction