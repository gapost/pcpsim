clear
clf


A = load('res0');
Ta = A(:,1);
res = zeros(size(Ta), 4);
res(:,1) = A(:,2);

A = load('res6');
res(:,2) = A(:,2);

A = load('res31');
res(:,3) = A(:,2);

A = load('res64');
res(:,4) = A(:,2);


plot(Ta,res(:,1)*1e6-467*0.7,'.-')
hold on
plot(Ta,res(:,2)*1e6-467*0.7,'.-')
plot(Ta,res(:,3)*1e6-467*0.7,'.-')
plot(Ta,res(:,4)*1e6-467*0.7,'.-')
hold off

ylabel('\rho - \rho_0 (n\Omega cm)')
xlabel('Ta (K)')
