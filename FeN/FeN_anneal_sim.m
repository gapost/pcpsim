% Simulate isochronal annealing of Fe-N
% Exp: Our 5 min annealing
% model by mean R equations
clear

% annealing conditions
Ta = [276	285	295	305	314	325	335	344	354	364	375	...
      385	395	406	416	432	446	460	474	488	504	518];
dt = 5*60; % ageing time in s

% Options
incub=0; % Calc. incubation time for nucleation

% alloy data
X0 = 4.8e-4; % Initial N concentration, 

% model calculation for 3 values of gamma (surface tension)
solver = 'ode15i'; % 'ode15i', 'daspk' (octave only)
nTa = length(Ta);
x = zeros(nTa,3,3);
F = zeros(nTa,3);
[x(:,:,1),F(:,1)] = FeN_model(dt,Ta,X0,0.056,incub,solver);
[x(:,:,2),F(:,2)] = FeN_model(dt,Ta,X0,0.058,incub,solver);
[x(:,:,3),F(:,3)] = FeN_model(dt,Ta,X0,0.060,incub,solver);

% model resistivity calculation
Xp = 1/9;
rhoNm = 0.7; % nOhm-cm / ppm of [N], Wagenblast 1968
rhoNp = rhoNm*0.20; % to fit the data

rho=zeros(nTa,3);
for i=1:3,
  rho(:,i) = 1e6*rhoNm*(x(:,3,i).*(1-F(:,i)) + rhoNp/rhoNm*Xp*F(:,i) - X0);
end

clf

A = load('FeN_I6_1.ascii');
TaI6 = A(:,5);
RI6 = A(:,3); 
i = find(TaI6<540);
plot(TaI6(i),RI6(i)-RI6(1),'o-',Ta,rho)

xlabel('Annealing Temperature (K)')
ylabel('\rho - \rho_0 (n\Omega-cm)')
legend('Isochronal, \Deltat = 5min','\gamma = 0.056 J/m^2', ...
    '\gamma = 0.058 J/m^2', ...
    '\gamma = 0.060 J/m^2',...
    'location','southwest')
    
text(257,-100,'Fe - 480 ppm N','fontsize',20)
text(257,-120,'Isochronal, \Deltat=5 min','fontsize',14)
text(257,-140,'\rho_N^p / \rho_N^m = 0.20','fontsize',14)
    
print -dpng FeN_anneal_sim

