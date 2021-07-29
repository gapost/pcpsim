% Simulate isochronal annealing of irradiated Fe-N
% Model:
%   - Fe16N2 nucleation by classical MLS model
%   - V - N and V - SIA reactions
clear

% experimental annealing conditions
Ta = [256 266 276	285	295	305	314	325	335	344	354	364	375	...
      385	395	406	416	432	446	460	474	488	504	518 532 548 ...
      564 570 586 602 616 632 650 700];
dt = 8*60; % ageing time in s

% alloy data
X0 = 4.8e-4; % Initial N concentration,
Xp = 1/9; % precipitate N concentration

% nucleation parameters
gam = 0.058; % Surface tension [J/m^2]
incub=0; % Calc. incubation time for nucleation

% initial VN cluster concentration
% scales with irradiation dose levels 1:5:10
XVN0 = [0 1 5 10]*8e-6; 

% V-N Binding energies - Barouh 2015
Eb(1) = 0.7; % eV binding energy of VN
Eb(2) = 0.8; % eV binding energy of VN2

% model resistivity calculation
rhoNm = 0.6; % nOhm-cm / ppm of [N], Wagenblast 1968
rhoNp = rhoNm*0.20; % validated from annealing exprmnts
rhoD = 1.5; % V & SIA resistivity
rhoVN = 1.5; % Adjusted to fit experiment
rhoVN2 = 0.2; % Adjusted to fit experiment 

rho = zeros(length(Ta),4);

  % x(1) = precipitates, x(2) = Radius, 
  % x(3) = N, x(4)=V, x(5)=VN, x(6)=VN2, x(7) = In

figure 1
clf

for i=1:4,
  
  [x,F,R] = FeNirrad_model(dt,Ta,X0-XVN0(i),XVN0(i),gam,incub,Eb);

  rho(:,i) = (x(:,3:7)*[rhoNm rhoD rhoVN rhoVN2 rhoD*4]').*(1-F) + ...
    rhoNp*Xp*F - rhoNm*X0;
  
  subplot(2,2,i)
  semilogy(Ta,[x(:,[3 5]).*(1-F) 2*x(:,6).*(1-F)  Xp*F])
  xlabel(' Ta (K)')
  ylabel('N atomic concentration')
  title(['X_{VN}^0 = ' num2str(XVN0(i))])
  ylim([1e-6 1e-3])
  legend('Fe matrix','VN','VN2','Precipitates',...
    'location','southeast')
  drawnow
  
end

rho = 1e6*rho;

figure 2
clf 
plot(Ta,rho,'o-');
hold on
plot([250 650], [0 0],'--k')
hold off
title('Resistivity')
xlabel(' Ta (K)')
ylabel ('\rho - \rho_0 (n\Omega-cm)');
xlim([250 650])
lbls = {};
for i=1:4, lbls(i) = ['X_{VN}^0 = ' num2str(XVN0(i))]; end
legend(lbls,'location','southeast')