function t = time2R(R,X,Xp,Xeq,R0)
% find the time t for a precipitate to grow from R* to R
  
  Rc = R0/Xp/log(X/Xeq);
  
  fxr = @(r) Xeq*exp(R0/Xp./r);
  fivr = @(r) r./(X-fxr(r)).*(Xp-fxr(r));

  t = zeros(size(R));
  for i=1:length(R),  
    t(i) = quad(fivr,1.05*Rc,R(i));
  end
  
  
endfunction