import numpy as np

def a_quadsolve(a, b, c):
	return [(-b + np.sqrt(b**2-(4*a*c)))/(2*a), (-b - np.sqrt(b**2-(4*a*c)))/(2*a)]

print(a_quadsolve(0.001, 1000, 0.001))

def b_quadsolve(a, b, c):
	return [2*c/(-b - np.sqrt(b**2 - 4*a*c)), 2*c/(-b - np.sqrt(b**2 - 4*a*c))]

print(b_quadsolve(0.001, 1000, 0.001))

def c_quadsolve(a, b, c):
	a, b, c  = np.float64(a), np.float64(b), np.float64(c)
	return a_quadsolve(a, b, c)

print(c_quadsolve(0.001, 1000, 0.001))
