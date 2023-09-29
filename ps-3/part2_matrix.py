import numpy as np
import time
import matplotlib.pyplot as plt

# Function to multiply two NxN matrices A and B
def matrix_multiply(A=None, B=None, N=None):
    C = np.zeros([N, N], dtype=float)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i, j] += A[i, k] * B[k, j]
    return C

# Timing the matrix multiplication function for N = 10 to 200
times_lst = []

for N in range(10, 201, 10):
    A = np.ones([N, N])
    B = np.random.randint(10, size=(N, N))
    
    start_time = time.time()
    ans = matrix_multiply(A, B, N)
    end_time = time.time()
    elapsed_time = end_time - start_time
    
    if N != 0:
        times_lst.append(elapsed_time)

# Plotting N vs the cube root of the time
plt.plot(range(10, 201, 10), np.cbrt(times_lst))
plt.title("Times for Matrix Multiplication")
plt.ylabel("Cube Root of Time (s^(1/3))")
plt.xlabel("N for NxN Matrix")
plt.savefig("part2_matrix_mult.png")
plt.show()

# Now timing the dot function of numpy
times_lst = []
dot_lst = []

for N in range(10, 201, 10):
    start_time = time.time()
    dot_lst.append(np.dot(A, B))
    end_time = time.time()
    elapsed_time = end_time - start_time
    times_lst.append(elapsed_time)

plt.plot(range(10, 201, 10), times_lst)
plt.title("Times for Matrix Multiplication Using dot()")
plt.ylabel("Time (s)")
plt.xlabel("N for NxN Matrix")
plt.savefig("part2_np.dot.png")
plt.show()
