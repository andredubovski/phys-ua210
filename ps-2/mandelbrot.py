import numpy as np
import matplotlib.pyplot as plt

n = 2222

x = np.linspace(-2, 2, n)
y = np.linspace(-2, 2, n)

x, y = np.meshgrid(x, y)
c = y*1j + x

z = np.zeros_like(c)
mandelbrot = np.zeros_like(c, dtype=int)

for i in range(100):
	mask = np.abs(z) < 2.0
	z[mask] = z[mask] ** 2 + c[mask]
	mandelbrot += mask

plt.imshow(mandelbrot, extent=(-2, 2, -2, 2), interpolation='bilinear')
plt.title("Mandelbrot Set")
plt.xlabel("Real Component")
plt.ylabel("Imaginary Component")
plt.colorbar()
plt.show()

