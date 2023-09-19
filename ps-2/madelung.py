import numpy as np
import time


half_L = int(input("What should half of L be (integer)?"))

M = 0
uses_for_loop = False

def V_without_constants(i, j, k):
	return sign(i, j, k) / ((i**2+j**2+k**2)**0.5) if (i**2+j**2+k**2) != 0 else 0

def sign(i, j, k):
        return -1 if (i+j+k) % 2 == 0 else 1

t_i = time.time()


if(uses_for_loop):

	def V_without_constants(i, j, k):
		return sign(i, j, k) / ((i**2+j**2+k**2)**0.5) if (i**2+j**2+k**2) != 0 else 0

	for i in range(-half_L, half_L):
		print(i)
		for j in range(-half_L, half_L):
			for k in range(-half_L, half_L):
				M += V_without_constants(i, j, k)
else:
	i = np.arange(-half_L, half_L)
	j = np.arange(-half_L, half_L)
	k = np.arange(-half_L, half_L)
	
	i, j, k = np.meshgrid(i, j, k)
	
	mask = (i != 0) | (j != 0) | (k != 0)	
	d = (i**2 + j**2 + k**2)**0.5
	sign = (-1.0) ** (i + j + k)

	M = -np.sum(sign[mask] / d[mask])	

t_f = time.time()

print("M:", M)

print("Time:", t_f-t_i)
