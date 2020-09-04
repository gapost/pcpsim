clear

figure 4
clf

A = load('FeN_I6_1.ascii');
TaI6 = A(:,5);
RI6 = A(:,3); 
RRI6 = (RI6)./(RI6(1))*100;

A = load('FeN555.dat');
Ta = A(:,1);
X = A(:,2); X0 = A(1,2);
F = A(:,3); Xp = 1/9;
C = (X.*(1-F) + Xp*F*0.4)/X0*100; 



% data from Irrad8 RvsTa 2nd take
A = load('resistancevsTa2.dat'); 
Ta8 = A(:,1);
R8_2 = A(:,2);
i = [2:2:size(R8_2)];
R8_2 = R8_2(i)*1594;
Ta8 = Ta8(i);
RR8 = R8_2 ./ R8_2(1)*100;

plot(TaI6(1:26),RRI6(1:26),'o-g')
hold on
plot(Ta8(1:20),RR8(1:20),'o-')
plot(Ta,C,'o-')
hold off

title('Experimental vs Model')
xlabel('Ta (K)')
ylabel('R/R_0 & C/C_0(%)')
legend('Un - Irrad (exp)','Irra8-2nd (exp) ',  
  '\gamma = 0.0555 J/m^2, \rho_p / \rho_m = 0.4',...
  'location', 'southwest' )
