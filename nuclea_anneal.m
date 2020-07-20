function [xdot, F, S] = nuclea_anneal(x,u,V)
%V = [P.Xeqs(k) P.F0(k) P.S0(k) P.b0(k) P.dG0(k) P.a_R0(k) P.Xp];

xdot = zeros(size(x));
S = V(7)*log(x(3,:)./V(1)) +(1-V(7)).*log((1-x(3,:))./(1-V(1)));
F = V(2) .* x(2,:).^3 .* x(1,:);
B = V(3).^2 ./ 2 ./ V(4);

S2 = S.^2;
b0_S2 = V(4)./S2;
e_R = exp(1/V(7)./x(2,:));
Xq = V(1) .* e_R;
e_G_S = exp(-V(5)./S.^2);


y = u/B;
if y<0.005,
  y = B^2/u.^2;
else
  y = xdot(1,:)./ x(1,:);
endif

Rgrowth = (V(6)./ x(2,:)) .*((x(3,:) - Xq) ./ (V(7) - Xq)) + y.*(1.05 ./S - x(2,:));

Ngrowth = b0_S2 .* e_G_S .*exp(-1 ./ (2*b0_S2.*u)); 


xdot(1,:) = Ngrowth;

xdot(2,:) = Rgrowth;

xdot(3,:) = (x(3,:) - V(7)) .* F ./ (1-F) .* (3.*xdot(2,:)./x(2,:) + y );



endfunction