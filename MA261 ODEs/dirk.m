function y = dirk (f,Df, t0,y0, h, alpha,beta,gamma)
  m = length(y0);
  s = length(gamma);
  y = y0;
  k = zeros(s,m);
  for i=1:s
    ti = t0+alpha(i)*h;
    yi = y0;
    for j=1:i-1
        yi = yi + (h*beta(i,j)).*k(j,:);
    end
    if beta(i,i) == 0
        k(i,:) = f(ti,yi); 
    else
        [k(i,:), n] = newton (@(k) k-f(ti,yi+(h*beta(i,i)).*k),...
                              @(k) eye(m)-(h*beta(i,i)).*Df(ti,yi+(h*beta(i,i)).*k),...
                              f(t0,y0), 1e-15,1000);
    end
    y = y + (h*gamma(i)).*k(i,:);
  end
end
