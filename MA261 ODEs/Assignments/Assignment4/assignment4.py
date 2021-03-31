import numpy as np
from numpy.linalg import solve
from numpy import array, diff, exp, sin, log10, linspace, floor
from scipy.integrate import solve_ivp
from matplotlib import pyplot

def newton(F,DF,x0,eps,K):
    k  = 0
    x  = x0.copy().astype(float)  # note: with x=x0 changes to x also changes to x0 with numpy arrays
    delta = np.zeros([len(x0)])
    Fx = F(x)
    while Fx.dot(Fx) > eps*eps and k<K:
        delta[:] = solve(DF(x), Fx)
        x[:] -= delta[:]  # don't construct a new vector in each step - they could be large
        Fx = F(x)
        k += 1
    return [x,k]
  
def dirk(f, Df, t0, y0, h, alpha, beta, gamma):
    s  = gamma.size
    m  = y0.size
    k  = np.zeros((s, m))
    f0 = f(t0,y0)
    y  = y0.copy().astype(float)
    for i in range(s):
        ti = t0+alpha[i]*h
        yi = y0+h*sum( [beta[i,l]*k[l, :] for l in range(i)] )
        if beta[i,i]==0:
          k[i,:] = f(ti,yi) 
        else:
          k[i,:] = newton( lambda d: d - f(ti, yi + h*beta[i,i]*d),
                           lambda d: np.eye(m) - h*beta[i,i]*Df(ti,yi + h*beta[i,i]*d),
                           f0, 1e-15, 1000)[0]
        y[:] += h*gamma[i]*k[i, :]
    assert( not np.any(np.isnan(y)) )
    return y
  
def evolve(phi, f,Df, t0,y0, T,N):
    # compute y_{n+1} = phi(t_n,y_n;h) for n=1,...,N and t_i=ih
    # with h=T/N so t_1=t0 and t_{N+1}=T
    h = T/N
    y = np.zeros( [N+1,len(y0)] )
    y[0] = y0
    t = 0
    for i in range(N):
        y[i+1] = phi(f,Df, t,y[i], h)
        t = t+h
    return y

####################################################################
####################################################################

def embeddirk( f,Df, t0,y0, h0, alpha,beta,gamma,gammaStar, p,tol,maxh):
    # implemented embeded rk method here - start with the dirk method given
    # below (copy it here since the idea is to only compute stages once, so
    # you can't call the 'dirk' function directly.
    # return values are the new value y_1 and time step h used
    # return y,h 
    
def adaptEvolve(phi, f, Df, t0,y0, T,hStart):
    y = np.zeros([ 1,len(y0) ])
    t = np.zeros([1])
    y[0,:] = y0
    t[0] = 0
    h = hStart
    while 1:
        ynew,h = phi(f,Df, t[-1],y[-1,:], h)
        tnew = t[-1]+h
        y = np.append(y, [ynew], 0)
        t = np.append(t, [tnew], 0)
        if tnew > T:
            break
    return t,y

####################################################################
####################################################################

def compute(problem,explicit):
    if explicit:
        print("problem ",problem," using explicit methods")
    else:
        print("problem ",problem," using implicit methods")
    maxh = 0.05 # maximal step size for adaptive methods

    # setup problem - defines
    # f,Df,T,y0, the tolerance for an adaptive method, and N0 which is the number of steps
    # per unit time, i.e., in [0,1].
    if problem == 0:
        # the problem used in the week 9 recording - but also in the singular
        # perturbation recording in week 8 (with constant eps)
        a    = 2      # choose a=0 for a constant eps
        eps0 = 0.0001
        eps = lambda t: a*(2*t-floor(2*t)) + eps0
        f =  lambda t,y: array([ 1./eps(t)*(sin(t)-y[0]) ])
        Df = lambda t,y: array([ -1./eps(t) ])
        T  = 10
        y0 = array([1.])
        if explicit:
            N0 = 500
        else:
            N0 = 50
        tol = 1e-5
        finished = lambda ydiff,t: True
        name = "p0"
    elif problem == 1:
        mu = 10
        T  = 20
        freq = 18.863
        f  = lambda t,y: array([ y[1], mu*(1-y[0]*y[0])*y[1]-y[0] ])
        Df = lambda t,y: array([ [0, 1], [-2*mu*y[0]*y[1]+1, mu*(1-y[0]*y[0])] ])
        y0 = array([2., 0.])
        if explicit:
            N0 = 200
        else:
            N0 = 20
        tol = 1e-1
        finished = lambda ydiff,t: True
        name = "p1"

    # setup ode solver
    if explicit:
        alpha = array([0.,    0.5,   1.])
        gamma = array([1./6., 2./3., 1./6.])
        beta  = array([ [0.,  0., 0.],
                    [0.5, 0., 0.],
                    [-1., 2., 0.] ])
        gammaStar = array([0,     1.,    0.])
        name = name+"_ex"
        bi_solver = "RK45"
    else:
        alpha = array([1./2.,  2./3.,  1./2.,  1.    ])
        gamma = array([3./2., -3./2.,  1./2.,  1./2. ])
        beta  = array([ [1./2.,  0.,     0.,     0.   ],
                    [1./6.,  1./2.,  0.,     0.   ],
                    [-1./2., 1./2.,  1./2.,  0.   ],
                    [3./2., -3./2.,  1./2.,  1./2.] ])
        gammaStar = array([1.,     0.,     0.,     0.    ])
        name = name+"_im"
        bi_solver = "Radau"

    ######################################################

    print("fixed: ",end="")
    while True:
        stepper = lambda f,Df,t0,y0,h: dirk(f,Df,t0,y0,h,alpha,beta,gamma)
        y = evolve(stepper,f,Df,0,y0, T,N0*T)
        t = linspace(0,T,len(y))
        if finished(y-y0,t):
            break
        N0 = N0*2
    if len(y0)==2:    # phase portrait
        pyplot.plot(y[:,0],y[:,1],'g.-')
    else:             # show solution as function of t
        pyplot.plot(t,y)
    pyplot.savefig(name+'_fixed.eps')
    pyplot.clf()
    print("N=",len(y),"h=",1./N0)
    # store fixed time step solution to compare with adaptive solvers
    tfixed = t
    yfixed = y

    ######################################################

    print("adaptive: TODO\n",end="")

    ######################################################

    print("build-in: ",end="")
    while 1:
        res = solve_ivp(f, [0,T], y0, method=bi_solver, max_step=maxh, rtol=tol)
        y = res.y.transpose()
        t = res.t
        if finished(y-y0,t):
            break
        tol = tol/10.
    h = np.diff(t)
    if len(y0)==2:
        pyplot.plot(y[:,0],y[:,1],'g.-')
        pyplot.savefig(name+'buildin_pp.eps')
        pyplot.clf()

        pyplot.plot(t,y[:,0],t,y[:,1])
    else:
        pyplot.plot(t,y[:,0],tfixed,yfixed)
    pyplot.savefig(name+'buildin.eps')
    pyplot.clf()
    print("N=",len(y),"h=",min(h))

##########################################################

if __name__ == "__main__":
    compute(0,True)
    compute(0,False)
    compute(1,True)
    compute(1,False)
