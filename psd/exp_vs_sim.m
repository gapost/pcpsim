clear

figure 4
clf

A = load('FeN_I6_1.ascii');
TaI6 = A(:,5);
RI6 = A(:,3); 
RRI6 = (RI6)./(RI6(1))*100;

A = load('FeN555.dat');
Ta = A(:,1);
C = A(:,2)/A(1,2)*100; 

A = load('FeN57.dat');
C = [C A(:,2)/A(1,2)*100]; 
A = load('FeN58.dat');
C = [C A(:,2)/A(1,2)*100];
##A = load('FeN59noincub.dat');
##C = [C A(:,2)/A(1,2)*100];
##A = load('FeN60incub.dat');
##C = [C A(:,2)/A(1,2)*100];
##A = load('FeN60noincub.dat');
##C = [C A(:,2)/A(1,2)*100];


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
  '\gamma = 0.0555 J/m^2 (sim)',...
  '\gamma = 0.057 J/m^2 ',...
  '\gamma = 0.058 J/m^2 ',...
##  '\gamma = 0.059 J/m^2 (no incub.)',...
##  '\gamma = 0.060 J/m^2 ',...
##  '\gamma = 0.060 J/m^2 (no incub.)',...
  'location', 'southwest' )
