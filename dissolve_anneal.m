function [xdot2, Xc, F] = dissolve_anneal(x2,u,P,N1,k)

xdot2 = zeros(size(x2));
F = P.F0(k) .* x2(1,:).^3 .* N1(13);
Xc = (P.Xc0 - P.Xp.*F) ./ (1 - F);
A.d = exp(1/P.Xp./x2(1,:));

%xdot2(1,:) = (2.*u).*(P.a^2/P.R0s(k)^2./x2(1,:)) .*(x2(2,:) - P.Xeqs(k)*A.d) ./ (P.Xp - P.Xeqs(k)*A.d) ;
%xdot2(2,:) = (x2(2,:) - P.Xp) .* F ./ (1-F) .* 3.*xdot2(1,:)./x2(1,:);

xdot2(1,:) = (2.*u).*(P.a^2/P.R0s(k)^2./x2(1,:)) .*(Xc - P.Xeqs(k)*A.d) ./ (P.Xp - P.Xeqs(k)*A.d) ;

endfunction