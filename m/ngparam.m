function [rat,tau,b0,dG0,R0] = ngparam(a,na,gam,D,T)

% Constants
kb = 8.617e-5; %boltzmann constant [eV/K]

% lattice data
a *= 1e9; % convert m to nm
Vat = a.^3./na; % Atomic volume [nm^3]
rat = (3*Vat/4/pi).^(1/3); % atomic radius [nm]

% Diffusion time const
tau = rat.^2./(D*1e18); % sec

% derived quantities
gam *= 6.24150913; % convert to eV/nm2
R0 = (2*gam.*Vat)./(kb*T)./rat; 
dG0 = 0.5*R0.^3;
Z = 1/20;
b0 = (4*pi*Z)*R0.^2.*(rat/a).^4;