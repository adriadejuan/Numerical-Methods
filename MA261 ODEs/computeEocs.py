# Write a function to compute the experimental order of convergence (EOC) for a given sequence
# of step sizes and errors given in the form of a mx2 matrix ((h_1,e_1),...,(h_m,e_m)).
# Return a vector of EOCs values.

import numpy as np

def npPrint(a, format_string ='{0:.4e}'):
    print(*[format_string.format(v,i) for i,v in enumerate(a)], end="")

def computeEocs(herr):
    eocs = np.zeros(len(herr)-1)
    for k in range(len(herr)-1):
        eocs[k] = np.log((herr[k+1][1])/(herr[k][1]))/np.log((herr[k+1][0])/(herr[k][0]))
    return eocs

# def computeEocs( herr ):
#     # herr = ( (h1,e1),(h2,e2),...,(hm,em) )
#     m = len(herr)
#     eocs = numpy.zeros( m-1 )
#     for i in range(m-1):
#         eocs[i] = numpy.log(herr[i+1,1]/herr[i,1] ) /\
#                   numpy.log(herr[i+1,0]/herr[i,0] );
#     return eocs

# ---------- COMPUTEEOCS TEST 1 ----------
Y = lambda t: exp(-10*t)
herr = np.zeros( (5,2) )
for m in range(5):
    h = 1/(20*2**m)
    y = forwardEuler( lambda t,y: -10*y,
                      lambda t,y: np.array([[-10]]),
                      0,np.array([1.]), h)
    herr[m,:] = [h, np.linalg.norm(y-Y(h))]
npPrint( computeEocs(herr) )
# SHOULD RETURN
# 1.8871e+00 1.9417e+00 1.9704e+00 1.9851e+00

# ---------- COMPUTEEOCS TEST 2 ----------
herr = np.zeros( (5,2) )
for p in range(1,5):
  for m in range(5):
      h = 1/(20*2**m)
      herr[m,:] = [h, h**(p/2)]
  npPrint( computeEocs(herr) )
  print()
# SHOULD RETURN
# 5.0000e-01 5.0000e-01 5.0000e-01 5.0000e-01
# 1.0000e+00 1.0000e+00 1.0000e+00 1.0000e+00
# 1.5000e+00 1.5000e+00 1.5000e+00 1.5000e+00
# 2.0000e+00 2.0000e+00 2.0000e+00 2.0000e+00
  
# ---------- COMPUTEEOCS TEST 3 ----------
Y = lambda t: exp(-10*t)
herr = np.zeros( (5,2) )
for m in range(5):
    h = 1/(20*2**m)
    y = forwardEuler( lambda t,y: -10*y,
                      lambda t,y: np.array([[-10]]),
                      0,np.array([1.]), h)
    herr[m,:] = [h, np.linalg.norm(y-Y(h))]
npPrint( computeEocs(herr) )
# SHOULD RETURN
# 1.8871e+00 1.9417e+00 1.9704e+00 1.9851e+00
