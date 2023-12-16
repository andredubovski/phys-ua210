import numpy as np
from numpy import ones,copy,cos,tan,pi,linspace
import math
import matplotlib.pyplot as plt

def V(x):
	return x**4

#function for finding integration points and weights for Gaussian quadrature
#written by Mark Newman, June 4, 2011
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

#period function with L = 2
def T(x,m,a):
    return math.sqrt((16*m)/(V(a) - V(x)))

m = 1	#mass
a = 2	#amplitude
n = 20	#steps

amplitudes = np.arange(0, a, 0.01) 
integrations = np.empty(shape = (0, ))

#Integrate
for i in range(amplitudes.size):
	x, w = gaussxw(n)
	xp = 0.5*(amplitudes[i])*x + 0.5*(amplitudes[i])
	wp = 0.5*(amplitudes[i])*w
	
	integ = 0

	for j in range(n):
		integ += wp[j]*T(xp[j], m, amplitudes[i])

	integrations = np.append(integrations, integ)

plt.plot(amplitudes, integrations)
plt.xlabel('a (m)')
plt.ylabel('T (s)')
plt.title('Amplitude vs. Period for Anharmonic Oscillator')
plt.savefig('exercise2.png')
plt.show()

