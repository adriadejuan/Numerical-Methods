from math import exp
import numpy as np
from numpy import array
from numpy.linalg import norm
from numpy import sin,cos,pi

def forwardEuler(f, Df, t0, y0, h):
      y = y0 + h*f(t0,y0)
      return y

def evolve(phi, f,Df, t0,y0, T,N):
    h = T/N
    y = [0]*(N+1)
    y[0]=y0
    for i in range(1,N+1):
        y[i] = phi(f,Df,h*i,y[i-1],h)
    return y
      
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
