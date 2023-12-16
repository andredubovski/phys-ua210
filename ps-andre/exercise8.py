import scipy.optimize as opt
import numpy as np

# Define the function
def f(x):
    return (x[0]-1)**2 + (x[1]-2)**2 + (x[2]+3)**2

# Define the gradient of the function
def grad(x):
    return np.array([2*(x[0]-1), 2*(x[1]-2), 2*(x[2]+3)])

# Initial guess
x0 = np.array([0, 0, 0])

# Use CG method
res_cg = opt.minimize(f, x0, method='CG')
print("CG method:")
print("Solution:", res_cg.x)
print("Number of iterations:", res_cg.nit)

# Use L-BFGS-B method
res_lbfgsb = opt.minimize(f, x0, method='L-BFGS-B')
print("\nL-BFGS-B method:")
print("Solution:", res_lbfgsb.x)
print("Number of iterations:", res_lbfgsb.nit)

# Use Newton-CG method
res_ncg = opt.minimize(f, x0, method='Newton-CG', jac=grad)
print("\nNewton-CG method:")
print("Solution:", res_ncg.x)
print("Number of iterations:", res_ncg.nit)
