def CrankNicholson (f, Df, t0, y0, h):
    F = lambda delta: delta - f(t0+h, y0+h*delta)
    DF = lambda delta: np.identity(len(y0)) - h*Df(t0+h, y0+h*delta)
    d,k = newton(F, DF, y0, 1e-8, 100)
    y = y0 + h/2*(f(t0,y0) + d)
    return y
