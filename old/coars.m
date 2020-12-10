function [ydot, fcoars, AS] = coars(y,u,P)

ydot = zeros(size(y));
AS = P.Xp*log(y(3,:)./P.Xeqs) +(1-P.Xp).*log((1-y(3,:))./(1-P.Xeqs));

fcoars  = 1 - erf(4.*(y(2,:).*AS - 1)) ;

##ydot(1,:) = fcoars.*(ydot(2,:)./y(2,:).*( (y(3,:).*(1 - y(3,:))./(y(3,:)-P.Xp).^2 ./y(2,:)).* ...
##(1/P.h./y(2,:).^3 - y(1,:)) -(3.*y(1,:)) ));



Rcoars = (4./27) .* (P.Xeqs./(P.Xp-P.Xeqs)) .* (P.a^2./(P.R0s.^2 .* y(2,:).^2));
Rgrowth = (P.a^2/P.R0s^2./y(2,:)) .*((y(3,:) - P.Xeqs*exp(1/P.Xp./y(2,:))) ./ (P.Xp - P.Xeqs*exp(1/P.Xp./y(2,:))));

ydot(1,:) = fcoars.*(Rcoars./y(2,:).*( (y(3,:).*(1 - y(3,:))./(y(3,:)-P.Xp).^2 ./y(2,:)).* ...
(1/P.h./y(2,:).^3 - y(1,:)) -(3.*y(1,:)) ));

ydot(2,:) = (1 - fcoars).*Rgrowth + fcoars.*Rcoars;

Ccoars = (y(3,:).*(1-y(3,:))./(y(3,:)-P.Xp)) .* (1./y(2,:).^2) .* ydot(2,:);

Cgrowth = P.h*(y(3,:) - P.Xp) .* (ydot(1,:).*y(2,:).^3 + 3*y(1,:).*y(2,:).^2 .*ydot(2,:)) ./ (1 - P.h* y(1,:).*y(2,:).^3);

ydot(3,:) = (1 - fcoars).*Cgrowth + fcoars.*Ccoars;
%ydot(3,:) = Cgrowth;
%ydot(3,:) = Ccoars;
endfunction