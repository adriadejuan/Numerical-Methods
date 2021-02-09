import numpy as np
from numpy import array
from numpy.linalg import norm
from numpy import sin,cos,pi,exp

def newton (F, DF, x0, eps, K):
    x = x0 - np.matmul(np.linalg.inv(DF(x0)),F(x0))
    k=1
    while ((norm(F(x)) > eps) and (k<K)):
        x = x - np.matmul(np.linalg.inv(DF(x)),F(x))
        k += 1
    return x,k
    
def backwardEuler(f, Df, t0, y0, h):
    F = lambda delta: delta - f(t0+h, y0+h*delta)
    DF = lambda delta: np.identity(len(y0)) - h*Df(t0+h, y0+h*delta)
    d,k = newton(F, DF, y0, 1e-8, 100)
    y = y0 + h*d
    return y
    
# ---------- FORWARDEULER TEST 1 ----------
y = evolve( backwardEuler, lambda t,y: -10*y,
                           lambda t,y: array([-10]),
                            0,array([1.]), 1,20)
print("%d "     % len(y), end="")
print("%1.3e  " % y[0],   end="")
print("%1.3e  " % y[-1],  end="")
# SHOULD RETURN
# 21 1.000e+00  3.007e-04

# ---------- FORWARDEULER TEST 2 ----------
y = evolve( backwardEuler, lambda t,y: -10*y,
                           lambda t,y: array([-10]),
                            0,array([1.]), 1,40)
print("%d "     % len(y), end="")
print("%1.3e  " % y[0],   end="")
print("%1.3e  " % y[-1],  end="")
# SHOULD RETURN
# 41 1.000e+00  1.329e-04

# ---------- FORWARDEULER TEST 3 ----------
y = backwardEuler( lambda t,y: array([-y[0],y[0]-exp(-t) ]),
                   lambda t,y: array([[-1,0],[1,0]]),
                   2,array([exp(-2),1.]), 0.5)
print("%d " % len(y), end="")
print("%1.3e  %1.3e  " % tuple(y),   end="")
# SHOULD RETURN
# 2 9.022e-02  1.004e+00

# ---------- FORWARDEULER TEST 4 ----------
y = backwardEuler( lambda t,y: array([ -y[0],
                                       y[0]-exp(-1*y[2]),
                                       1 ]),
                   lambda t,y: array([[-1,0,0],[1,0,-exp(-y[2])],[0,0,0]]),
                   0,array([exp(-2),1.,2.]), 0.5)
print("%d " % len(y), end="")
print("%1.3e  %1.3e  %1.3e  " % tuple(y),   end="")
# SHOULD RETURN
# 3 9.022e-02  1.004e+00  2.500e+00
