function s = intinvarrho(t,r,alphas)
  s=2*((1 - r^2)^(2*alphas))*(1 - cos(t))./(abs(1 -r^2*exp(i*t)).^(2*alphas)-...
      ((1-r^2)^(2*alphas)))./abs(1-r^2*exp(i*t)).^2;
  
endfunction
