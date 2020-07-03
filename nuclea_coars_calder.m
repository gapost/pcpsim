function [xdot, fcoars, S] = nuclea_coars_calder(x,u,P)

xdot = zeros(size(x));
S = P.Xp*log(x(3,:)./P.Xeqs) +(1-P.Xp).*log((1-x(3,:))./(1-P.Xeqs));

%fcoars  = 1 - erf(4.*(x(2,:).*AS - 1)) ;
fcoars  = 1 - erf(4.*(x(3,:)./P.Xeqs.*x(2,:).*S - 1)) ;

B = P.S0.^2 / 2 / P.b0;

y = u/B;
if y<0.005,
  y = B^2/u.^2;
else
  y = xdot(1,:)./ x(1,:);
endif

Rcoars = (4./27) .* (P.Xeqs*(1-P.Xeqs)./(P.Xp-P.Xeqs).^2) .* (P.a^2./(P.R0s.^2 .* x(2,:).^2));
Rgrowth = (P.a^2/P.R0s^2./x(2,:)) .*((x(3,:) - P.Xeqs*exp(1/P.Xp./x(2,:))) ./ (P.Xp - P.Xeqs*exp(1/P.Xp./x(2,:)))) - y.*(1.05 ./S - x(2,:));
%Rgrowth = (P.a^2/P.R0s^2./x(2,:)) .*((x(3,:) - P.Xeqs*exp(1/P.Xp./x(2,:))) ./ (P.Xp - P.Xeqs*exp(1/P.Xp./x(2,:))));

Ngrowth = (P.b0./S.^2) .*exp(-P.dG0./S.^2) .*exp(-(S.^2)./(2*P.b0*u)); 
##Ncoars = Rcoars./x(2,:).*( (x(3,:).*(1 - x(3,:))./(x(3,:)-P.Xp).^2 ./x(2,:)).* ...
##(1/P.h./x(2,:).^3 - x(1,:)) -(3.*x(1,:)) );
Ncoars = (3/4)*P.h./(4*pi*P.R0s.^3.*x(2,:).^4) .*xdot(2,:);

xdot(1,:) = (1 - fcoars).*Ngrowth + fcoars.*Ncoars;

xdot(2,:) = (1 - fcoars).*Rgrowth + fcoars.*Rcoars;

##Ccoars = (x(3,:).*(1-x(3,:))./(x(3,:)-P.Xp)) .* (1./x(2,:).^2) .* xdot(2,:);
Ccoars = P.Xeqs*(1-P.Xeqs)./(P.Xeqs - P.Xp) ./ x(2,:).^2 .*xdot(2,:);
Cgrowth = P.h*(x(3,:) - P.Xp) .* (xdot(1,:).*x(2,:).^3 + 3*x(1,:).*x(2,:).^2 .*xdot(2,:)) ./ (1 - P.h* x(1,:).*x(2,:).^3);

xdot(3,:) = (1 - fcoars).*Cgrowth + fcoars.*Ccoars;
endfunction