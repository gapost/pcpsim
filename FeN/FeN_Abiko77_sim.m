% Calculate Fe - 880 ppm N nucleation & growth at 373 K
% Exp: Abiko & Imai 1977
% model by mean R equations
clear

% annealing conditions
Ta = 373; % ageing T (K)
dt = 1e6; % ageing time in s

% alloy exp data 
% [N] = 880 ppm
% [C] = 140 ppm
% rho(77 K) = 1.3 micro-Ohm-cm
%
% rho(Fe,77K) = 0.612 micro-Ohm-cm
% Defect resistivity = 1.3 - 0.612 = 0.688 = 688 nOhm-cm
%
% We assume that variation of resistivity is solely due to
% [N] precipitation
X0 = 8.8e-4;  % [N]

% Options
incub=1; % Calc. incubation time for nucleation

% log time grid - 10 pts per decade
% t in s
nt = (log10(dt)-1)*10 + 1;
t = logspace(1,log10(dt),nt);

% model calculation for 3 values of gamma (surface tension)
solver = 'daspk'; % 'ode15i', 'daspk' (octave only)
x = zeros(nt,3,3);
F = zeros(nt,3);
[x(:,:,1),F(:,1)] = FeN_model(t,Ta,X0,0.058,incub,solver);
[x(:,:,2),F(:,2)] = FeN_model(t,Ta,X0,0.060,incub,solver);
[x(:,:,3),F(:,3)] = FeN_model(t,Ta,X0,0.062,incub,solver);

% model resistivity calculation
Xp = 1/9;
rhoNm = 0.7; % nOhm-cm / ppm of [N], Wagenblast 1968
rhoNp = rhoNm*0.24; % to fit the data

rho=zeros(nt,3);
for i=1:3,
  rho(:,i) = 1e6*rhoNm*(x(:,3,i).*(1-F(:,i)) + rhoNp/rhoNm*Xp*F(:,i) - X0);
end

clf

A = load('abiko.csv');
t_ab = A(:,1);
RR_ab = A(:,2);
R_ab = (RR_ab-RR_ab(1))/100*1300;

semilogx(t_ab,R_ab,'o-',t/60,rho)

xlabel('Annealing time (min)')
ylabel('\rho - \rho_0 (n\Omega-cm)')
legend('Abiko & Imai 1977', ...
    '\gamma = 0.058 J/m^2', ...
    '\gamma = 0.060 J/m^2', ...
    '\gamma = 0.062 J/m^2',...
    'location','southwest')
    
text(0.15,-200,'Fe - 880 ppm N','fontsize',20)
text(0.15,-240,'Isothermal, T=373 K','fontsize',14)
text(0.15,-280,'\rho_N^p / \rho_N^m = 0.24','fontsize',14)

print -dpdfcrop FeN_Abiko77_sim