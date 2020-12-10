function [xdot2, Xc, F] = dissolve_anneal_2(x2,u,V)
%V = [P.Xeqs(k) P.F0(k) P.S0(k) P.b0(k) P.Xc0 P.a_R0(k) P.Xp N1(13)];

xdot2 = zeros(size(x2));
F = V(2) .* x2.^3 .* V(8);
Xc = (V(5) - V(7).*F) ./ (1 - F);

e_R = exp(1/V(7)./x2);
Xq = V(1) .* e_R;

%xdot2 = (2.*u).*(V(6)./x2) .*(Xc - Xq) ./ (V(7) - Xq) ;
xdot2 = (V(6)./x2) .*(Xc - Xq) ./ (V(7) - Xq) ;

endfunction