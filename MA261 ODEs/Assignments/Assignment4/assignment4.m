function assignment3
  compute(0,true)
  compute(0,false)
  compute(1,true)
  compute(1,false)
  exit(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compute(problem,explicit)
  if explicit
    fprintf("problem %d using explicit methods\n",problem)
  else
    fprintf("problem %d using implicit methods\n",problem)
  end
  format longE
  maxh = 0.05; % maximal step size for adaptive methods

  % setup problem - defines
  % f,Df,T,y0, the tolerance for an adaptive method, and N0 which is the number of steps
  % per unit time, i.e., in [0,1].
  if problem == 0
    % the problem used in the week 9 recording - but also in the singular
    % perturbation recording in week 8 (with constant eps)
    a    = 2;      % choose a=0 for a constant eps
    eps0 = 0.0001;
    eps = @(t) a*(2*t-floor(2*t)) + eps0;
    T  = 10;
    f =  @(t,y) [ 1./eps(t)*(sin(t)-y(1)) ];
    Df = @(t,y) [ -1./eps(t) ];
    y0 = [1.];
    if explicit
      N0 = 500;
    else
      N0 = 50;
    end
    tol = 1e-5;
    finished = @(ydiff,t) true; % don't do any loop
    name = "p0";
  elseif problem == 1
    % van der Pol equation with quite forgiving parameter choice
    mu = 10;
    f =  @(t,y) [ y(2) mu*(1-y(1)*y(1))*y(2)-y(1) ];
    Df = @(t,y) [ 0 1 ; -2*mu*y(1)*y(2)+1 mu*(1-y(1)*y(1)) ];
    y0 = [2. 0.];
    T  = 20;
    freq = 18.863; % approximate frequency for this set of parameters
    if explicit
      N0 = 200;
    else
      N0 = 20;
    end
    tol = 1e-1;
    finished = @(ydiff,t) true; % should be replaced by a function handle to compute the frequency
    name = "p1";
  end

  % setup ode solver
  if explicit
    alpha = [0.    0.5   1.];
    gamma = [1./6. 2./3. 1./6.];
    beta  = [0.   0. 0. ;...
             0.5  0. 0. ;...
            -1.   2. 0. ];
    gammaStar = [0     1.    0.];
    name = strcat(name,"_ex");
    bi_solver = @ode45;
  else
    alpha = [1./2.  2./3.  1./2.  1.    ];
    gamma = [3./2. -3./2.  1./2.  1./2. ];
    beta  = [1./2.  0.     0.     0.    ;...
             1./6.  1./2.  0.     0.    ;...
            -1./2.  1./2.  1./2.  0.    ;...
             3./2. -3./2.  1./2.  1./2. ];
    gammaStar = [1.     0.     0.     0.    ];
    name = strcat(name,"_im");
    bi_solver = @ode23s;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fprintf("fixed: ");
  while 1 % loop added to approximate frequency - needs to define 'finished'
    stepper = @(f,Df,t0,y0,h) dirk(f,Df,t0,y0,h,alpha,beta,gamma);
    y = evolve(stepper,f,Df,0,y0, T,N0*T); 
    t = linspace(0,T,length(y));
    if finished(y-y0,t)
      break
    end
    N0 = N0*2;
  end
  if length(y0)==2  % phase portrait
    plot(y(:,1),y(:,2),'g.-');
  else              % show solution as function of t
    plot(t,y);
  end
  print(strcat(name,'_fixed'),'-depsc');
  fprintf("N=%d h=%d\n",length(y),1/N0);
  % store fixed time step solution to compare with adaptive solvers
  tfixed = t;
  yfixed = y;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fprintf("adaptive: TODO\n");

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fprintf("build-in: ");
  while 1
    [t,y] = bi_solver(@(t,y) f(t,y)',[0 T],y0',...
            odeset('MaxStep',maxh,'NormControl','on','RelTol',tol));
    if finished(y-y0,t)
      break
    end
    tol = tol/10.;
  end
  h = diff(t);
  if length(y0)==2
    plot(y(:,1),y(:,2),'g.-');
    print(strcat(name,'buildin_pp'),'-depsc');

    plot(t,y(:,1),'g.-',t,y(:,2),'b.-');
  else
    plot(t,y,'g.-');
  end
  print(strcat(name,'buildin'),'-depsc');

  % compare on part of the time interval - of interest for problem 0 but
  % you could do something similar around a point in time of interest for
  % other problems....
  if length(y0)==1
    idx = find(t<0.1);
    idxfixed = find(tfixed<0.1);
    plot(t(idx),y(idx),'g.-',tfixed(idxfixed),yfixed(idxfixed),'r.-');
    print(strcat(name,'buildin_cmp'),'-depsc');
  end
  fprintf("N=%d hmin=%d\n",length(y),min(h));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,h] = embeddirk( f,Df, t0,y0, h0, alpha,beta,gamma,gammaStar, p,tol,maxh)
  % implemented embeded rk method here - start with the dirk method given
  % below (copy it here since the idea is to only compute stages once, so
  % you can't call the 'dirk' function directly.
end
function [t,y] = adaptEvolve(phi, f, Df, t0,y0, T,hStart)
  y = zeros([ 1 length(y0) ]);
  t = zeros([1]);
  y(1,:) = y0;
  t(1) = 0;
  h = hStart;
  while 1
    [ynew,h] = phi(f,Df, t(end),y(end,:), h);
    tnew = t(end)+h;
    y = [y ; ynew];
    t = [t ; tnew];
    if tnew > T
      break;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = evolve(phi, f,Df, t0,y0, T,N)
  % compute y_{n+1} = phi(t_n,y_n;h) for n=1,...,N and t_i=ih
  % with h=T/N; so t_1=t0 and t_{N+1}=T
  h = T/N;
  y = zeros([ N+1 length(y0) ]);
  y(1,:) = y0;
  t = 0;
  for i=1:N
    y(i+1,:) = phi(f,Df, t,y(i,:), h);
    t = t+h;
  end
end
function y = dirk( f,Df, t0,y0, h, alpha,beta,gamma)
  m = length(y0);
  s = length(alpha);
  k = zeros(s,m);
  y = y0;
  for i=1:s
    ti = t0+alpha(i)*h;
    yi = y0;
    for j=1:s-1
      yi = yi + (h*beta(i,j)).*k(j,:);
    end
    if beta(i,i) == 0
      k(i,:) = f(ti,yi); 
    else
      [k(i,:), n] = newton (@(ki) ki-f(ti,yi+(h*beta(i,i)).*ki),...
                            @(ki) eye(m)-(h*beta(i,i)).*Df(ti,yi+(h*beta(i,i)).*ki),...
                            f(t0,y0), 1e-15,1000);
    end
    y = y + (h*gamma(i)).*k(i,:);
  end
  assert( any(isnan(y))==0 );
end
function [x,k] = newton(F,DF,x0,eps,K)
  k  = 0;
  x  = x0;
  Fx = F(x);
  while dot(Fx,Fx) > eps*eps && k<K
    x  = x - (DF(x)\Fx')';
    Fx = F(x);
    k  = k+1;
  end
end
