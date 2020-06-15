function sol = coupled_ode(E, dfuns, steps, a, b, ini, method)
  range = b-a;
  h=range/steps;  
  rows = (range/h)+1;
  columns = size(dfuns)(2)+1;
  sol= zeros(abs(rows),columns);
  heun=zeros(1,columns-1);
  for i=1:abs(rows)
    if i==1
      sol(i,1)=a;
    else
      sol(i,1)=sol(i-1,1)+h;      
    end  
    for j=2:columns
      if i==1
        sol(i,j)=ini(j-1);
      else
        if strcmp("euler",method)
          sol(i,j)=sol(i-1,j)+h*dfuns{j-1}(E, sol(i-1,1:end));      
        elseif strcmp("heun",method)
          heun(j-1)=sol(i-1,j)+h*dfuns{j-1}(E, sol(i-1,1:end));          
        elseif strcmp("rk4",method)
          k1=h*dfuns{j-1}(E, [sol(i-1,1), sol(i-1,2:end)]);
          k2=h*dfuns{j-1}(E, [sol(i-1,1)+(0.5*h), sol(i-1,2:end)+(0.5*h*k1)]);
          k3=h*dfuns{j-1}(E, [sol(i-1,1)+(0.5*h), sol(i-1,2:end)+(0.5*h*k2)]);
          k4=h*dfuns{j-1}(E, [sol(i-1,1)+h, sol(i-1,2:end)+(h*k3)]); 
          sol(i,j)=sol(i-1,j)+((1/6)*(k1+(2*k2)+(2*k3)+k4));       
        end  
      end
    end
    if strcmp("heun",method)
      if i~=1
        for k=2:columns
          sol(i,k)=sol(i-1,k)+(h/2)*((dfuns{k-1}(E, sol(i-1,1:end)))+(dfuns{k-1}(E, [sol(i,1),heun])));
        end 
      end  
    end     
  end
end