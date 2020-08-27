
subplot(2,1,1)

function [x,y,z] = select_dat(R,dR,f,dfdt)
  j = find(f.*dR>1e-23);
  x = R(j);
  y = f(j);
  z = dfdt(j);
endfunction

[xx,yy] = select_dat(R,dR,f(idx(1),:),dfdt(idx(1),:));
loglog(xx*rat,yy/rat,'.-')

hold on
for kk=2:length(idx),
  [xx,yy] = select_dat(R,dR,f(idx(kk),:),dfdt(idx(kk),:));
  loglog(xx*rat,yy/rat,'.-')  
endfor
hold off
title('PSD')
ylabel('f (nm^{-1})')
if exist("lbls"), legend(lbls); end

subplot(2,1,2)

[xx,yy,zz] = select_dat(R,dR,f(idx(1),:),dfdt(idx(1),:));
semilogx(xx*rat,zz/rat,'.-')

hold on
for kk=2:length(idx),
  [xx,yy,zz] = select_dat(R,dR,f(idx(kk),:),dfdt(idx(kk),:));
  semilogx(xx*rat,zz/rat,'.-')  
endfor
hold off

xlabel('R (nm)')
ylabel('df/dt (nm^{-1}s^{-1})')
