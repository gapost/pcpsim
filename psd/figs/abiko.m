clear
clf

A = load('abiko.csv');
t_ab = A(:,1);
RR_ab = A(:,2);
R_ab = RR_ab*1.3/100 + 1.3;

A = load('tetsu.csv');
t_te = A(:,1);
R_te = A(:,2);

R77_fe = 0.612;
semilogx(t_ab,R_ab-R77_fe,'o-',t_te,R_te-R77_fe,'^-')

return

A = load('tvsX.dat');
t_sim = A(:,1);
X_sim = A(:,2);

R77_fe = 0.69; %??-cm
X0 = 8.8e-4; % Initial N concentration, Abiko

%ax = plotyy (t_ab,R_ab, t_te,R_te, @semilogx, @semilogx);
semilogx(t_ab,(R_ab-R77_fe)/(R_ab(1)-R77_fe),...
         t_te,(R_te-R77_fe)/(R_te(1)-R77_fe),...
         t_sim/60,(X_sim/X0));
xlim([0.1, 1e4])        
xlabel ("t (min) ");
ylabel ("DR/R0 -or- X/X0");
legend('abiko','tetsu','sim');
title("Aging at T=100 C")