from math import sqrt, exp, pi, factorial
import matplotlib.pyplot as plt
import numpy as np

def H(x, n):
    if n == 0:
        return 1
    elif n == 1:
        return 2 * x
    else:
        return 2 * x * H(x, n - 1) - 2 * (n - 1) * H(x, n - 2)

#Non-recursive H polynomial finder:
def H_nr(x,n):
    if n == 0:
        return 1
    elif n == 2:
        return 2*x
    else:
        for i in range (2, n+1):
            h_0 = 1; h_1 = 2*x
            h_2 = 2*x*h_1 - 2*(i-1)*h_0
            h_0 = h_1
            h_1 = h_2
        return h_1

def psi(x, n):
	return (1/sqrt((2**n)*factorial(n)*sqrt(pi))) * exp(-(x**2)/2) * H_nr(x,n)

#x_range = np.arange(-4, 4, 0.01)

#psi_ranges = [np.empty(shape = (0, ))]*4

#for x in x_range:
#	for i in range(len(psi_ranges)):
#		psi_ranges[i] = np.append(psi_ranges[i], psi(x, i))



#Comment this out to do part B:

#for i, psi_range in enumerate(psi_ranges):
	#plt.plot(x_range, psi_range, label='n='+str(i))

#plt.legend()
#plt.title('$\Psi_n(x)$ for First Four Levels in 1D Quantum Oscillator')
#plt.xlabel('x (m)')
#plt.ylabel("$\Psi_n(x)$ ($\\frac{1}{\sqrt{m}})$")
#plt.show()


#Part B:
x_range = np.arange(-10,10, 0.01)
psi24_range = np.empty(shape = (0, ))

for i in range(x_range.size):
    psi24_range = np.append(psi24_range, psi(x_range[i], 30))

plt.plot(x_range, psi24_range)
plt.legend()
plt.title('Wavefunction of 30th Energy Level for 1D Quantum Oscillator')
plt.xlabel('x (m)')
plt.ylabel("$\Psi_n(x)$ ($\\frac{1}{\sqrt{m}})$")
plt.show()
