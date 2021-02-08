def forwardEuler(f, Df, t0, y0, h):
      y = y0 + h*f(t0,y0)
      return y
      
# ---------- FORWARDEULER TEST 1 ----------
fprintf("%1.4e",forwardEuler(@(t,y) [-y], @(t,y) [-1], 0,1,0.1) )
# SHOULD RETURN
# 9.0000e-01

# ---------- FORWARDEULER TEST 2 ----------
fprintf("%1.4e",forwardEuler(@(t,y) [-100*y], @(t,y) [-100], 0,1,0.1) )
# SHOULD RETURN
# -9.0000e+00

# ---------- FORWARDEULER TEST 3 ----------
fprintf("%1.4e",forwardEuler(@(t,y) [-100*y],@(t,y) [-100], 0,1,0.009) )
# SHOULD RETURN
# 1.0000e-01

# ---------- FORWARDEULER TEST 4 ----------
y = forwardEuler(@(t,y) [t*t+1], @(t,y) [0], 3,12, 0.1);
T = 3 + 0.1;
Y = @(t) [t+t^3/3];
fprintf("%1.4e", norm(y-Y(T))/norm(Y(T)) );
# SHOULD RETURN
# 2.3279e-03

# ---------- FORWARDEULER TEST 5 ----------
y = forwardEuler(@(t,y) [-y(1), t*t+1], @(t,y) [-1 0 ; 0 0], 3,[exp(-3),12], 0.1);
T = 3 + 0.05;
Y = @(t) [exp(-t) t+t^3/3];
fprintf("%1.4e", norm(y-Y(T))/norm(Y(T)) );
# SHOULD RETURN
# 3.9373e-02
