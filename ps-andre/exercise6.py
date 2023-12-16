import astropy.io.fits    
import matplotlib.pyplot as plt
import numpy as np
import os
import timeit

# Function to read data from the FITS file
def read_data(filename):
    curr_path = os.path.join(os.getcwd(), filename)
    hdu_list = astropy.io.fits.open(curr_path)
    logwave = hdu_list['LOGWAVE'].data
    flux = hdu_list['FLUX'].data
    hdu_list.close()
    return logwave, flux

# Function to plot the spectrum for the first five galaxies
def plot_spectrum(logwave, flux, labels):
    plt.figure()
    colors = ['blue', 'orange', 'green', 'red', 'purple']
    for i in range(len(labels)):
        plt.plot(logwave, flux[i, :], label=f'{labels[i]}', color=colors[i])
    plt.legend()
    plt.title('Spectrum of the First Five Galaxies')
    plt.xlabel('Log Wavelength (log(Angstrom))')
    plt.ylabel('Flux (10^-17 erg/(s*cm^2*Angstrom))')
    plt.show()

# Function to normalize flux
def normalize_flux(flux):
    flux_sum = flux.sum(axis=1)
    flux_norm = flux / flux_sum[:, None]
    return flux_norm

# Function to subtract the mean from the normalized flux
def subtract_mean(flux_norm):
    flux_mean = np.mean(flux_norm, axis=1)
    flux_without_offset = flux_norm - flux_mean[:, None]
    return flux_without_offset

# Function to calculate eigenvectors using the covariance matrix
def calculate_eigenvectors(flux_without_offset):
    r = flux_without_offset
    c = np.dot(r.T, r)
    eig_val_cov, eig_vec_cov = np.linalg.eig(c)
    return eig_vec_cov

# Function to plot eigenvectors obtained using linalg.eig
def plot_eigenvectors(logwave, eigenvectors):
    plt.figure()
    colors = ['blue', 'orange', 'green', 'red', 'purple']
    for i in range(5):
        plt.plot(logwave, eigenvectors[:, i], label=f'Eigenvector {i + 1}', color=colors[i])
    plt.title("Eigenvectors using Covariance Matrix")
    plt.xlabel('Log Wavelength (log(Angstrom))')
    plt.ylabel('Eigenvector Amplitude (no units)')
    plt.legend()
    plt.show()

# Function to calculate eigenvectors using Singular Value Decomposition (SVD)
def calculate_eigenvectors_svd(flux_without_offset):
    u, w, vt = np.linalg.svd(flux_without_offset, full_matrices=False)
    eig_vec_svd = vt.transpose()
    return eig_vec_svd, w  # Return both eigenvectors and singular values

# Function to plot eigenvectors obtained using SVD
def plot_eigenvectors_svd(logwave, eig_vec_svd):
    plt.figure()
    colors = ['blue', 'orange', 'green', 'red', 'purple']
    for i in range(5):
        plt.plot(logwave, eig_vec_svd[:, i], label=f'Eigenvector {i + 1}', color=colors[i])
    plt.title("Eigenvectors using Singular Value Decomposition (SVD)")
    plt.xlabel('Log Wavelength (log(Angstrom))')
    plt.ylabel('Eigenvector Amplitude (no units)')
    plt.legend()
    plt.show()

# Function to compare runtimes of eig and svd methods
def compare_runtimes(start1, stop1, start2, stop2):
    print("\nRuntime Comparison:")
    print("Runtime for eigenvectors using Covariance Matrix: ", stop1 - start1)
    print("Runtime for eigenvectors using SVD: ", stop2 - start2)
    print()

# Function to compare condition numbers of eig and svd methods
def compare_condition_numbers(eig_val_cov, w, r, c):
    condition_number_eig = np.max(eig_val_cov) / np.min(eig_val_cov)
    condition_number_svd = np.max(w) / np.min(w)
    print("Condition Number Comparison:")
    print("Condition number of R (data matrix) = ", np.linalg.cond(r))
    print("Condition number of C (covariance matrix) = ", np.linalg.cond(c))
    print()

# Function to perform Principal Component Analysis (PCA)
def perform_pca_analysis(flux_without_offset, eig_vec_svd, flux_mean, flux_sum):
    c0 = np.dot(flux_without_offset, eig_vec_svd[:, 0])
    c1 = np.dot(flux_without_offset, eig_vec_svd[:, 1])
    c2 = np.dot(flux_without_offset, eig_vec_svd[:, 2])
    c3 = np.dot(flux_without_offset, eig_vec_svd[:, 3])
    c4 = np.dot(flux_without_offset, eig_vec_svd[:, 4])

    weights = np.vstack((c0, c1, c2, c3, c4))
    eigenvectors = np.vstack((eig_vec_svd[:, 0], eig_vec_svd[:, 1], eig_vec_svd[:, 2], eig_vec_svd[:, 3], eig_vec_svd[:, 4]))

    approx_spectra = np.dot(weights.T, eigenvectors)
    approx_spectra = approx_spectra + flux_mean[:, None]
    approx_spectra = approx_spectra * flux_sum[:, None]
    return c0, c1, c2, c3, c4, approx_spectra

# Function to plot PCA results
def plot_pca_results(c1, c0, c2, logwave):
    fig, axes = plt.subplots(2, 1, figsize=(6, 8))

    axes[0].scatter(c1, c0, color='skyblue')
    axes[0].set_xlabel('Principal Component 1 (no units)')
    axes[0].set_ylabel('Principal Component 0 (no units)')
    axes[0].set_title('PCA: Principal Component 0 vs Principal Component 1')

    axes[1].scatter(c2, c0, color='lightcoral')
    axes[1].set_xlabel('Principal Component 2 (no units)')
    axes[1].set_ylabel('Principal Component 0 (no units)')
    axes[1].set_xlim(-0.003, 0.003)
    axes[1].set_title('PCA: Principal Component 0 vs Principal Component 2')

    plt.tight_layout()
    plt.show()

# Function to perform PCA analysis with a varying number of components
def pca_with_varied_components(flux_without_offset, eig_vec_svd, logwave):
    for i in range(1, 21):
        eigen_spectra = PCA(flux_without_offset, eig_vec_svd, i)
        fractional_residuals = ((np.abs(flux_without_offset) - np.abs(eigen_spectra)) / np.abs(flux_without_offset)) ** 2
        residual_mean = np.mean(fractional_residuals, axis=0)
        plt.scatter(logwave, residual_mean, label=f'Components = {i}', color='mediumseagreen')
    plt.legend()
    plt.xlabel('Log Wavelength (log(Angstrom))')
    plt.ylabel('Fractional Mean Residuals (no units)')
    plt.title('PCA: Fractional Mean Residuals vs Log Wavelength')
    plt.show()

def PCA(flux_without_offset, eig_vec_svd, n):
    weights = np.zeros((flux.shape[0], n))
    eigenvectors = np.zeros((flux.shape[1], n))

    for i in range(0, n):
        c = np.dot(flux_without_offset, eig_vec_svd[:, i])
        weights[:, i] = c
        eigenvectors[:, i] = eig_vec_svd[:, i]

    eigen_spectra = np.dot(weights, eigenvectors.T)
    return eigen_spectra

if __name__ == "__main__":
    filename = 'specgrid.fits'
    logwave, flux = read_data(filename)
    labels = [f'Galaxy {i + 1}' for i in range(5)]

    plot_spectrum(logwave, flux, labels)

    flux_norm = normalize_flux(flux)
    flux_without_offset = subtract_mean(flux_norm)

    start1 = timeit.default_timer()
    eig_vec_cov = calculate_eigenvectors(flux_without_offset)
    stop1 = timeit.default_timer()

    plot_eigenvectors(logwave, eig_vec_cov)

    start2 = timeit.default_timer()
    eig_vec_svd, w = calculate_eigenvectors_svd(flux_without_offset)  # Capture singular values
    stop2 = timeit.default_timer()

    plot_eigenvectors_svd(logwave, eig_vec_svd)

    compare_runtimes(start1, stop1, start2, stop2)
    compare_condition_numbers(eig_vec_cov, w, flux_without_offset, np.dot(flux_without_offset.T, flux_without_offset))

    perform_pca_analysis(flux_without_offset, eig_vec_svd, np.mean(flux_norm, axis=1), flux.sum(axis=1))
    pca_with_varied_components(flux_without_offset, eig_vec_svd, logwave)
