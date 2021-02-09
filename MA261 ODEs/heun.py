def heunMethod(f, Df, t0, y0, h):
    x = f(t0,y0)
    y = y0 + h/2*(x + f(t0 + h, y0 +h*x))
    return y
