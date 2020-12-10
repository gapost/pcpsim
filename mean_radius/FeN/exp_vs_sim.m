clear
clf

A = load('FeN_I6_1.ascii');
TaI6 = A(:,5);
RI6 = A(:,3); 
RRI6 = (RI6)./(RI6(1))*100;

B = load('../FeN_XvsTa.dat');
Ta = B(:,1);
C = B(:,2:4); 
CC = (C)./(C(1))*100;

% data from Irrad8 RvsTa 2nd take
C = load('resistancevsTa2.dat'); 
Ta8 = C(:,1);
R8_2 = C(:,2);
i = [2:2:size(R8_2)];
R8_2 = R8_2(i)*1594;
Ta8 = Ta8(i);
RR8 = R8_2 ./ R8_2(1)*100;

plot(TaI6(1:26),RRI6(1:26),'o-g')
hold on
plot(Ta8(1:20),RR8(1:20),'o-')
plot(Ta,CC,'o-')
hold off

title('Experimental vs Model')
xlabel('Ta (K)')
ylabel('R/R_0 & C/C_0(%)')
legend('Un - Irrad (exp)','Irra8-2nd (exp) ',  '\gamma = 0.059 J/m^2 (sim)',...
 '\gamma = 0.060 J/m^2 (sim)','\gamma = 0.061 J/m^2 (sim)', 'location', 'southwest' )
