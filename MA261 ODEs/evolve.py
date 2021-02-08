# Write a function 'evolve' to compute the approximation y_n with n=0,...,N
# to the solution Y of an initial value problem

import numpy as np

def forwardEuler(f, Df, t0, y0, h):
      y = y0 + h*f(t0,y0)
      return y
      
def evolve(phi, f,Df, t0,y0, T,N):
    h = T/N
    y = np.zeros(N+1)
    y[0]=y0
    for i in range(1,N+1):
        y[i] = phi(f,Df,h*i,y[i-1],h)
    return y
    
# ---------- EVOLVE TEST 1 ----------
y = evolve( forwardEuler, lambda t,y: -10*y,
                          lambda t,y: array([-10]),
                          0,array([1]), 1,20)
print("%d "     % len(y), end="")
print("%1.3e  " % y[0],   end="")
print("%1.3e  " % y[-1],  end="")
# SHOULD RETURN
# 21 1.000e+00  9.537e-07
