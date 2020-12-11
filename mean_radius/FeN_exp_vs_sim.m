clear
clf

A = load('FeN_I6_1.ascii');
TaI6 = A(:,5);
RI6 = A(:,3); 
RRI6 = (RI6)./(RI6(1))*100;

% data from Irrad8 RvsTa 2nd take
C = load('resistancevsTa2.dat'); 
Ta8 = C(:,1);
R8_2 = C(:,2);
i = [2:2:size(R8_2)];
R8_2 = R8_2(i)*1594;
Ta8 = Ta8(i);
RR8 = R8_2 ./ R8_2(1)*100;

C = [];
A = load('FeN55_8min_incub.dat');
Ta = A(:,1); C = [C A(:,2)/A(1,2)*100];
A = load('FeN56_8min_incub.dat');
Ta = A(:,1); C = [C A(:,2)/A(1,2)*100];
A = load('FeN57_8min_incub.dat');
Ta = A(:,1); C = [C A(:,2)/A(1,2)*100]; 


plot(TaI6(1:26),RRI6(1:26),'o-g')
hold on
plot(Ta8(1:20),RR8(1:20),'o-')
plot(Ta,C,'o-')
hold off

title('Experimental vs Model (Incub. On, 8 min anneal)')
xlabel('Ta (K)')
ylabel('R/R_0 & C/C_0(%)')
legend('Un - Irrad (exp)','Irra8-2nd (exp) ',  ...
  '\gamma = 0.055 J/m^2 (sim)',...
  '\gamma = 0.056 J/m^2 (sim)',...
  '\gamma = 0.057 J/m^2 (sim)',...
  'location', 'southwest' )
  
print2pdf(gcf,[20 15],'FeN_exp_vs_sim')
