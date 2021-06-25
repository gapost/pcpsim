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
gam = 0.058; % Surface tension [J/m^2]

% model calculation for 3 values of gamma (surface tension)
solver = 'ode15i'; % 'ode15i', 'daspk' (octave only)
[x,F] = FeN_model(dt,Ta,X0,gam,incub,solver);

clf

subplot(2,2,1)
plot(Ta,x(:,1),'.-')
title('Clusters per atom');

subplot(2,2,2)
plot(Ta,x(:,2),'.-')
ylabel('nm');
title('radius R')

subplot(2,2,3)
plot(Ta,x(:,3),'.-')
title('matrix X_N');
xlabel('Annealing T (K)');

subplot(2,2,4)
plot(Ta,F,'.-')
title('Transformed volume fraction ');
xlabel('Annealing T (K)');
    
print -dpng FeN_isochronal_nucl