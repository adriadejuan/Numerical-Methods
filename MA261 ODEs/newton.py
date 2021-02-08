import numpy as np
from numpy.linalg import norm

def newton (F, DF, x0, eps, K):
    x = x0 - np.matmul(np.linalg.inv(DF(x0)),F(x0))
    k=1
    while ((norm(F(x)) > eps) and (k<K)):
        x = x - np.matmul(np.linalg.inv(DF(x)),F(x))
        k += 1
    return x,k
    
# ---------- NEWTON TEST 1 ----------
x,k = newton( lambda y: cos(y/2),
              lambda y: array([-sin(y/2)/2]),
              array([3.]), 1e-8, 100)
relErr = norm(x-pi) / pi
if relErr<1e-10:
  print("CORRECT:\n relative error < 1e-10")
else:
  print("value %g is incorrect\n relative error is %g" % (x,relErr) )
# SHOULD RETURN
# CORRECT:
# relative error < 1e-10

# ---------- NEWTON TEST 2 ----------
A = array([ [1,2], [3,1] ]);
root = array([1,1]);
b = array([3,4]);
x,k = newton( lambda y: A.dot(y)-b,
              lambda y: A,
              array([0.,0.]), 1e-8, 100)
absErr = norm(x-root)
if absErr<1e-8 and k==1:
  print("CORRECT:\n absolute error < 1e-8 computed in single step")
else:
  print("Computed root:",x)
  print("Absolute error is",absErr," Steps needed",k)
# SHOULD RETURN
# CORRECT:
# absolute error < 1e-8 computed in single step

# ---------- NEWTON TEST X ----------
x,k = newton( lambda y: y*y+1,
              lambda y: array([2*y]),
              array([2.]), 1e-12, 100)
if k == 100:
  print("CORRECT:\n no convergence!")
else:
  print("there is no root so k should be K=100 but is %d" % k)
# SHOULD RETURN
# CORRECT:
# no convergence!
