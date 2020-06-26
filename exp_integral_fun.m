clear
clf

x = logspace(-3, 0, 100);

y = - expint(1./x) + exp(-1./x).* x;
ydot = exp(-1./x);

dy = ydot ./ y;

ya = 1 ./ x.^2;

dif = dy - ya  ;
relerr = dif ./ dy *100;

loglog(x,dif,'.-')