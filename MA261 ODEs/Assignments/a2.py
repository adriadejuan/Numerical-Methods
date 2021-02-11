from math import exp,sqrt
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt

def evolve(phi, f, Df, t0, y0, T, N):
    h = T/N
    y = [0]*(N+1)
    y[0]=y0
    for i in range(1,N+1):
        y[i] = phi(f,Df,h*i,y[i-1],h)
    return y

def forwardEuler(f, Df, t0, y0, h):
      y = y0 + h*f(t0,y0)
      return y

def heunMethod(f, Df, t0, y0, h):
    x = f(t0,y0)
    y = y0 + h/2*(x + f(t0 + h, y0 +h*x))
    return y

def computeEocs(herr):
    eocs = np.zeros(len(herr)-1)
    for k in range(len(herr)-1):
        eocs[k] = np.log((herr[k+1][1])/(herr[k][1]))/np.log((herr[k+1][0])/(herr[k][0]))
    return eocs

def npPrint(a, format_string ='{0:.4e}'):
    print(*[format_string.format(v,i) for i,v in enumerate(a)], end="")

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

def CrankNicholson (f, Df, t0, y0, h):
    F = lambda delta: delta - f(t0+h, y0+h*delta)
    DF = lambda delta: np.identity(len(y0)) - h*Df(t0+h, y0+h*delta)
    d,k = newton(F, DF, y0, 1e-8, 100)
    y = y0 + h/2*(f(t0,y0) + d)
    return y


c=1.5
T=10
N0 = 20

X = np.linspace(0,T,T+1)
Y1 = np.empty_like(X)
Y2 = np.empty_like(X)

#Q2.1
print("Q2.1: Backward Euler method")
print("Errors at T depending on h (time step):")
for i in range(T+1):
    X[i] = T/(N0*(2**i))
    h = T/(N0*(2**i))
    y = evolve(backwardEuler, lambda t,y: (c-y)**2, lambda t,y: 2*(y-c), 0, np.array([1.]), T, N0*(2**i))
    Y1[i] = abs((1 + T*c*(c-1))/(1+T*(c-1)) - y[-1])
    print("The error at T with h (time step) =", end=" ")
    print("%1.2e" % h, end=" ")
    print("is", end=" ")
    print("%1.4e " % Y1[i], end="")
    print(" and f has been evaluated", N0*(2**i), "times.")
print("Experimental order of convergence (EOCs):")
herr1 = np.zeros((T,2))
for i in range(T):
    herr1[i][0] = T/(N0*2**i)
    herr1[i][1] = Y1[i]
npPrint(computeEocs(herr1))
print("")
print("")
#
# #Q2.2
print("Q2.2: Crank-Nicholson method")
print("Errors at T depending on h (time step):")
for i in range(T+1):
    h = T/(N0*(2**i))
    yother = evolve(CrankNicholson, lambda t,y: (c-y)**2, lambda t,y: 2*(y-c), 0, np.array([1.]), T, N0*(2**i))
    Y2[i] = abs((1 + T*c*(c-1))/(1+T*(c-1)) - yother[-1])
    print("The error at T with h (time step) =", end=" ")
    print("%1.2e" % h, end=" ")
    print("is", end=" ")
    print("%1.4e " % Y2[i], end="")
    print(" and f has been evaluated", 2*N0*(2**i), "times.")
print("Experimental order of convergence (EOCs):")
herr2 = np.zeros((T,2))
for i in range(T):
    herr2[i][0] = T/(N0*2**i)
    herr2[i][1] = Y2[i]
npPrint(computeEocs(herr2))
print("")
print("")

#Q2.3
print("Q2.3")
def fun(t,y):
    p = 1/(sqrt(2))
    if (t < p):
        return -y
    else:
        return y

def dfun(t,y):
    p = 1/(sqrt(2))
    if (t < p):
        return -1
    else:
        return 1

def exact(t):
    p = 1/(sqrt(2))
    if(t <  p):
        return exp(-t)
    else:
        return exp(t-sqrt(2))

N=15
err1 = np.zeros((N,2))
err2 = np.zeros((N,2))
print("Forward Euler:")
for i in range(N):
    print("Error for h = %.6e" %(1/2**i), end=": ")
    y = evolve(forwardEuler, fun, dfun, 0, np.array([1.]), 1, 2**i)
    err1[i,:] = [1/2**i, np.linalg.norm(np.abs(y[-1] - exact(1)))]
    print("%.4e" % err1[i,1])
print("Forward Euler EOC:")
x1 = computeEocs(err1)
npPrint(x1)
print("\nEOC average is %.4e" %np.mean(x1))

print("\nHeun's Method:")

for i in range(N):
    print("Error for h = %.6e" %(1/2**i), end=": ")
    y1 = evolve(heunMethod, fun, dfun, 0, np.array([1.]), 1, 2**i)
    err2[i,:] = [1/2**i, np.linalg.norm(np.abs(y1[-1] - exact(1)))]
    print("%.4e" % err2[i,1])


print("\n Heun's Method EOC:")
x2 = computeEocs(err2)
npPrint(x2)
print("\nEOC average is %.4e" %np.mean(x2))

h = 2**N
X = np.linspace(0,1,h + 1)
Y1 = evolve(forwardEuler, fun, dfun, 0, np.array([1.]), 1,h)
Y2 = np.empty_like(X)
for i in range(h + 1):
    Y2[i] = exact(i/h)
plt.plot(X,Y1, label="Forward Euler")
plt.plot(X,Y2, label="Exact Y")
plt.legend()
plt.show()
