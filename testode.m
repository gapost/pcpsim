clear

P.Ta = 423;
P.km = 0.2;

function xdot = func(x,t,P)
  xdot = P.km*t;
endfunction


func(1,1,P)

ifunc = @(x,t) func(x,t,P);

t = 0:0.1:1;

x = lsode (ifunc, 0, t);

plot(t,x)

