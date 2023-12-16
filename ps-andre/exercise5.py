import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as linalg

# Function to retrieve data from input file and store them in arrays
def load_data(file_path):
    data = np.genfromtxt(file_path, delimiter="|", skip_header=1)
    time, signal = data[:, 1], data[:, 2]
    return time, signal

# Function to plot data
def plot_data(x, y, color, label, title):
    plt.title(title)
    plt.scatter(x, y, color=color, label=label, marker='.')
    plt.xlabel("Time (s)")
    plt.ylabel("Signal (V)")
    plt.legend()

# Function to perform polynomial fit
def polynomial_fit(order, x, y):
    A = np.vander(x, order + 1)
    u, w, vt = linalg.svd(A, full_matrices=False)
    w_inv = np.where(w < 1e-10, w, 1. / w)
    ainv = vt.T.dot(np.diag(w_inv)).dot(u.T)
    c = ainv.dot(y)
    y_model = A.dot(c)

    print("Condition number for polynomial fit of order", order, ":", np.max(w_inv) / np.min(w_inv))

    return y_model

# Function to plot polynomial fit and residuals
def plot_polynomial_fit_and_residuals(order, x, y, y_model):
    residuals = abs(y - y_model)
    plot_data(x, y, color='green', label='data', title=f'Data Fitting using Polynomial Function of Order {order}')
    plot_data(x, y_model, color='blue', label=f'{order}-order poly fit', title=f'Data Fitting using Polynomial Function of Order {order}')
    plot_data(x, residuals, color='yellow', label='residuals', title=f'Data Fitting using Polynomial Function of Order {order}')

# Function to perform trigonometric fit
def trig_fit(n, q, x, y):
    f = q * (np.max(x))
    A = np.ones((len(x), 2 * n + 1))
    for i in range(1, n + 1):
        A[:, 2 * i - 1] = np.cos(2. * np.pi * i * x / f)
        A[:, 2 * i] = np.sin(2. * np.pi * i * x / f)

    u, w, vt = linalg.svd(A, full_matrices=False)
    w_inv = np.where(w == 0., 0., 1. / w)
    ainv = vt.T.dot(np.diag(w_inv)).dot(u.T)
    c = ainv.dot(y)
    y_model = A.dot(c)

    return y_model

# Function to calculate residuals and r^2 for trigonometric fits
def calculate_trig_fit_statistics(y, y_model):
    residuals = abs(y - y_model)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return residuals, r_squared

# Main part of the code
file_path = 'signal.dat'
time, signal = load_data(file_path)

# Polynomial fits
order_3_ymodel = polynomial_fit(3, time, signal)
plot_polynomial_fit_and_residuals(3, time, signal, order_3_ymodel)
plt.show()

order_9_ymodel = polynomial_fit(9, time, signal)
plot_polynomial_fit_and_residuals(9, time, signal, order_9_ymodel)
plt.show()

# Trigonometric fits
trig_ymodel_half = trig_fit(10, 0.5, time, signal)
plot_data(time, signal, color='green', label='data', title='Data Fitting using Trigonometric Functions')
plot_data(time, trig_ymodel_half, color='red', label='n = 10, q = 0.5', title='Data Fitting using Trigonometric Functions')
plt.show()

trig_ymodel_2 = trig_fit(10, 2, time, signal)
plot_data(time, signal, color='green', label='data', title='Data Fitting using Trigonometric Functions')
plot_data(time, trig_ymodel_2, color='blue', label='n = 10, q = 2', title='Data Fitting using Trigonometric Functions')
plt.show()

# Calculate and print statistics for trigonometric fits
residuals_half, r_squared_half = calculate_trig_fit_statistics(signal, trig_ymodel_half)
residuals_2, r_squared_2 = calculate_trig_fit_statistics(signal, trig_ymodel_2)

print()
print("R^2 value for fit = 1/2 period:", r_squared_half)
print("R^2 value for fit = 2 period:", r_squared_2)
print("Mean residual for 1/2 period:", np.mean(residuals_half))
print("Mean residual value for fit = 2 period:", np.mean(residuals_2))
print()
