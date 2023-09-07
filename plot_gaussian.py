import matplotlib.pyplot as plt
import numpy as np

x_axis = np.arange(-10, 10, 0.01)
mean = 0
sd = 3

def norm_pdf(x, mean, std_dev):
    exponent = -0.5 * ((x - mean) / std_dev) ** 2
    coefficient = 1 / (std_dev * np.sqrt(2 * np.pi))
    return coefficient * np.exp(exponent)

plt.plot(x_axis, norm_pdf(x_axis, mean, sd))
plt.xlabel('x', fontsize=20)
plt.ylabel('y', fontsize=20)


plt.savefig("gaussian.png")
plt.show()

