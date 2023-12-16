import math
import numpy as np
import scipy.optimize

def f(x):
    return ((x - 1) ** 4) * np.exp(-(x**2)/2) if x >= -1 or x <= 2 else 50

def s_quad_interpolation(a, b, c):
    e = 1e-7
    s0 = a * f(b) * f(c) / (e + (f(a) - f(b)) * (f(a) - f(c)))
    s1 = b * f(a) * f(c) / (e + (f(b) - f(a)) * (f(b) - f(c)))
    s2 = c * f(a) * f(b) / (e + (f(c) - f(a)) * (f(c) - f(b)))
    return s0 + s1 + s2

def brent(f, a, b, tol=1e-7):
    if abs(f(a)) < abs(f(b)):
        a, b = b, a

    c = a
    u = True

    while abs(b - a) > tol:
        if f(a) != f(c) or f(b) != f(c):
            s = s_quad_interpolation(a, b, c)
        elif f(a) != f(b):
            s = (b * f(a) - a * f(b)) / (f(a) - f(b))

        if (s - (3 * a + b) / 4) * (s - b) > 0 or (u and abs(s - b) >= abs(b - c) / 2) or (
                not u and abs(s - b) >= abs(c - d) / 2) or (u and abs(b - c) < tol) or (
                not u and abs(c - d) < tol):
            s = (a + b) / 2
            u = True
        else:
            u = False

        d = c
        c = b

        if f(a) * f(s) < 0:
            b = s
        else:
            a = s

        if abs(f(a)) < abs(f(b)):
            a, b = b, a

    return b

def golden_mean_search(f, a, b, tol=1e-7):
    golden_ratio = (1 + math.sqrt(5)) / 2

    c = b - (b - a) / golden_ratio
    d = a + (b - a) / golden_ratio

    while abs(b - a) > tol:
        if f(c) < f(d):
            b = d
        else:
            a = c
        c = b - (b - a) / golden_ratio
        d = a + (b - a) / golden_ratio

    return (a + b) / 2

print("my implementation: {}".format(brent(f, -1, 2)))

print('scipy: {}'.format(scipy.optimize.brent(f)))
