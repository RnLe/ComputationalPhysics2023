import numpy as np
import matplotlib.pyplot as plt

# Parameters
delta = 0.05
kappa = 1e-5
N = int(1/delta) + 1

# Initialize arrays
phi = np.ones((N, N))  # potential
rho = np.zeros((N, N))  # charge density

# Boundary conditions
phi[0, :] = 0
phi[-1, :] = 0
phi[:, 0] = 0
phi[:, -1] = 0

def gauss_seidel(phi, rho):
    """Perform one iteration of the Gauss-Seidel method."""
    for i in range(1, N-1):
        for j in range(1, N-1):
            phi[i, j] = 0.25 * (phi[i-1, j] + phi[i+1, j] + phi[i, j-1] + phi[i, j+1] + delta**2 * rho[i, j])

def iterate(phi, rho):
    """Iterate until convergence."""
    while True:
        phi_old = phi.copy()
        gauss_seidel(phi, rho)
        if np.max(np.abs(phi - phi_old)) < kappa:
            break

# Test the Gauss-Seidel method
iterate(phi, rho)

# Plot the potential
plt.imshow(phi, cmap='viridis', extent=(0, 1, 0, 1))
plt.colorbar(label='Potential')
plt.show()

# Reset potential
phi = np.zeros((N, N)) 

# Boundary conditions
phi[0, :] = 0
phi[-1, :] = 1
phi[:, 0] = 0
phi[:, -1] = 0

# Iterate with Gauss-Seidel method
iterate(phi, rho)

# Plot the potential
plt.imshow(phi, cmap='viridis', extent=(0, 1, 0, 1))
plt.colorbar(label='Potential')
plt.show()

# Reset potential and charge density
phi = np.zeros((N, N)) 
rho = np.zeros((N, N)) 

# Place charge at the center
rho[N//2, N//2] = 1

# Boundary conditions
phi[0, :] = 0
phi[-1, :] = 0
phi[:, 0] = 0
phi[:, -1] = 0

# Iterate with Gauss-Seidel method
iterate(phi, rho)

# Calculate electric field
E = np.gradient(-phi)

# Plot the potential
plt.imshow(phi, cmap='viridis', extent=(0, 1, 0, 1))
plt.colorbar(label='Potential')
plt.show()

# Plot the electric field
plt.quiver(E[1], E[0])
plt.show()

# Reset potential and charge density
phi = np.zeros((N, N)) 
rho = np.zeros((N, N)) 

# Place charges
rho[int(0.25/delta), int(0.25/delta)] = 1
rho[int(0.75/delta), int(0.75/delta)] = 1
rho[int(0.25/delta), int(0.75/delta)] = -1
rho[int(0.75/delta), int(0.25/delta)] = -1

# Boundary conditions
phi[0, :] = 0
phi[-1, :] = 0
phi[:, 0] = 0
phi[:, -1] = 0

# Iterate with Gauss-Seidel method
iterate(phi, rho)

# Calculate electric field
E = np.gradient(-phi)

# Plot the potential
plt.imshow(phi, cmap='viridis', extent=(0, 1, 0, 1))
plt.colorbar(label='Potential')
plt.show()

# Plot the electric field
plt.quiver(E[1], E[0])
plt.show()
