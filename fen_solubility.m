% Calculate solubility of N in equilibrium with gamma-FeN and a"-FeN
% from Wriedt 1987 "The Fe-N system"
% equations in page 360

clear
clf

T=50:10:500;
T=T+273;

Cn_g = 10.^(1.69 - 1810./T);

Cn_a1 = 10.^(3.12-2160./T);
Cn_a2 = 10.^(2.48-1770./T);
Cn_a3 = 10.^(2.43 - 1840./T);

plot([Cn_g; Cn_a1; Cn_a2; Cn_a3],T,'.-')
legend('\gamma','\alpha 1','\alpha 2','\alpha 3')