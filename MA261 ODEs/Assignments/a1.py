from math import exp
import numpy as np
import matplotlib.pyplot as plt

def forwardEuler(f, Df, t0, y0, h):
      y = y0 + h*f(t0,y0)
      return y

def otherMethod(f, Df, t0, y0, h):
    x = f(t0,y0)
    y = y0 + h/2*(x + f(t0 + h, y0 +h*x))
    return y

def npPrint(a, format_string ='{0:.4e}'):
    print(*[format_string.format(v,i) for i,v in enumerate(a)], end="")

#Q2.0
def computeEocs(herr):
    eocs = np.zeros(len(herr)-1)
    for k in range(len(herr)-1):
        eocs[k] = np.log((herr[k+1][1])/(herr[k][1]))/np.log((herr[k+1][0])/(herr[k][0]))
    return eocs

def evolve(phi, f, Df, t0,y0, T,N):
    h = T/N
    y = [0]*(N+1)
    y[0]=y0
    for i in range(1,N+1):
        y[i] = phi(f,Df,h*i,y[i-1],h)
    return y

c=1.5
T=10
N0 = 20

X = np.linspace(0,T,T+1)
Y1 = np.empty_like(X)
Y2 = np.empty_like(X)

#Q2.1
print("Q2.1: Forward Eurler method")
print("Errors at T depending on h (time step):")
for i in range(T+1):
    X[i] = T/(N0*(2**i))
    h = T/(N0*(2**i))
    y = evolve(forwardEuler, lambda t,y: (c-y)**2, lambda t,y: 1, 0, np.array([1.]), T, N0*(2**i))
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

#Q2.2
print("Q2.2: Other method")
print("Errors at T depending on h (time step):")
for i in range(T+1):
    h = T/(N0*(2**i))
    yother = evolve(otherMethod, lambda t,y: (c-y)**2, lambda t,y: 1, 0, np.array([1.]), T, N0*(2**i))
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
print("Q2.3: Comparison between methods")
print()
plt.plot(X, Y1, label="Forward Euler method")
plt.plot(X, Y2, label="Heun's method")
plt.xlabel('h (Time Step)')
plt.yscale("log")
plt.ylabel('Error')
plt.title('Error on T vs h (time step)')
plt.legend()
plt.show()

# COMPARING BOTH METHODS:
# From our tested region (h in (5e-4, 0,5)), we can observe that the
# alternative Q2.2 method (Heun's method) converges faster and returns,
# for same time steps, more precise approximations  (less error), at the
# cost of doubling function f evaluations. So we can say the alternative
# method is more precise and equally efficient.
