import math
import matplotlib.pyplot as plt
import numpy as np 
from numpy import ones,copy,cos,tan,pi,linspace
import scipy 
from scipy.integrate import quad

def hermite_polynomial(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return 2*x
    else:
        return 2*x*hermite_polynomial(n-1, x) - 2*(n-1)*hermite_polynomial(n-2, x)

def compute_factorial(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n * compute_factorial(n-1)

def wavefunction(n, x):
    return (1/math.sqrt((2**n)*compute_factorial(n)*math.sqrt(math.pi))) * math.exp(-(x**2)/2) * hermite_polynomial(n, x)

def integrand_function(x, n):
    return x**2 * wavefunction(n, x)**2

def gaussian_quadrature_points_weights(N):
    a = linspace(3, 4*N-1, N) / (4*N+2)
    x = cos(pi*a + 1/(8*N*N*tan(a)))
    
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = ones(N, float)
        p1 = copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2*k+1)*x*p1 - k*p0) / (k+1)
        dp = (N+1)*(p0 - x*p1) / (1 - x*x)
        dx = p1 / dp
        x -= dx
        delta = max(abs(dx))
    
    w = 2*(N+1)*(N+1) / (N*N*(1 - x*x)*dp*dp)
    
    return x, w

a = -1
b = 1
p = 100

x_gaussian, w_gaussian = gaussian_quadrature_points_weights(p)
xp_gaussian = 0.5*(b-a)*x_gaussian + 0.5*(b+a)
wp_gaussian = 0.5*(b-a)*w_gaussian

integration_gaussian = 0
for i in range(p):
    z = xp_gaussian[i]
    integration_gaussian += wp_gaussian[i] * integrand_function(z/(1-(z**2)), 5) * (1+z**2)/(1-z**2)**2

uncertainty_gaussian = math.sqrt(integration_gaussian)

print()
print("Uncertainty with Gaussian-quadrature: ", uncertainty_gaussian)

x_hermgauss, w_hermgauss = np.polynomial.hermite.hermgauss(p)
xp_hermgauss = 0.5*(b-a)*x_hermgauss + 0.5*(b+a)
wp_hermgauss = 0.5*(b-a)*w_hermgauss

integration_hermgauss = 0
for i in range(p):
    z = xp_hermgauss[i]
    integration_hermgauss += w_hermgauss[i] * integrand_function(z/(1-(z**2)), 5) * (1+z**2)/(1-z**2)**2

integration_hermgauss *= (math.pi)

print("Uncertainty with Gauss-Hermite quadrature: ", math.sqrt(integration_hermgauss))

quad_result, _ = quad(lambda x: integrand_function(x, 5), -np.inf, np.inf)
print("Uncertainty with SciPy's quad function:", math.sqrt(quad_result))
print()
