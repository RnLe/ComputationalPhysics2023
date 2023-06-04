import typing
import numpy as np

# Gauss-Seidel algorithm

# Perform one iteration of the Gauss-Seidel method. Changes phi.
def _gauss_seidel(phi: np.ndarray, rho: np.ndarray, N: int, delta: float) -> None:
    # Start at 1 and end at N-1 to neglect the boundaries
    for i in range(1, N-1):
        for j in range(1, N-1):
            phi[i, j] = 0.25 * (phi[i-1, j] + phi[i+1, j] + phi[i, j-1] + phi[i, j+1] + delta**2 * rho[i, j])


# Iterate until convergence. Changes phi.
def iterate(phi: np.ndarray, rho: np.ndarray, kappa: float, N: int, delta: float) -> None:
    while True:
        phi_old = phi.copy()
        _gauss_seidel(phi, rho, N, delta)
        if np.max(np.abs(phi - phi_old)) < kappa:
            break
        
        
# Comparation of numerical with the analytical result
# Analytical function
def _analytical_phi(x: float, y: float) -> float:
    steps = 10000
    result = 0
    eps = 1e-10  # Define a threshold for the minimum term size
    for n in range(1, steps + 1):
        npi = n * np.pi
        term = 2 * (1 - np.cos(npi)) / (npi * np.sinh(npi)) * np.sin(npi * x) * np.sinh(npi * y)
        if abs(term) < eps:  # If the term is smaller than the threshold, break the loop
            break
        result += term
    return result



# Calculate the complete grid
def calc_ana_array(N: int, delta: float) -> np.ndarray:
    result = np.zeros((N, N))
    for i in range(N):
        for k in range(N):
            result[N - k - 1, i] = _analytical_phi(i * delta, k * delta)
    return result  