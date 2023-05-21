import numpy as np
import matplotlib.pyplot as plt

# ::::::::::::::::::::::::::::::::::::::::::::
# a)

# Parameters
delta = 0.05
kappa = 1e-5
N = int(1 / delta) + 1

print(f"Step size Delta: {delta}")
print(f"Error bound Kappa: {kappa}")
print(f"Discretization N: {N}")

# ::::::::::::::::::::::::::::::::::::::::::::
# Gauss-Seidel algorithm

# Perform one iteration of the Gauss-Seidel method.
def gauss_seidel(phi, rho):
    # Start at 1 and end at N-1 to neglect the boundaries
    for i in range(1, N-1):
        for j in range(1, N-1):
            phi[i, j] = 0.25 * (phi[i-1, j] + phi[i+1, j] + phi[i, j-1] + phi[i, j+1] + delta**2 * rho[i, j])

# Iterate until convergence.
def iterate(phi, rho):
    while True:
        phi_old = phi.copy()
        gauss_seidel(phi, rho)
        if np.max(np.abs(phi - phi_old)) < kappa:
            break
        
# ::::::::::::::::::::::::::::::::::::::::::::
# b)
        
# Initialize arrays
phi = np.ones((N, N))  # potential
rho = np.zeros((N, N))  # charge density

# Boundary conditions
phi[0, :] = 0
phi[-1, :] = 0
phi[:, 0] = 0
phi[:, -1] = 0

# Perform the Gauss-Seidel method
iterate(phi, rho)

# Plot the potential
plt.imshow(phi, cmap='viridis', extent=(0, 1, 0, 1))
plt.colorbar(label='Potential')
plt.show()

# Calculate electric field
E = np.gradient(-phi)

# Plot the electric field
plt.quiver(E[1], E[0])
plt.show()

# ::::::::::::::::::::::::::::::::::::::::::::
# c)

# Reset potential
phi = np.zeros((N, N)) 

# Boundary conditions
phi[0, :] = 1
phi[-1, :] = 0
phi[:, 0] = 0
phi[:, -1] = 0

# Iterate with Gauss-Seidel method
iterate(phi, rho)

# Plot the potential
plt.imshow(phi, cmap='viridis', extent=(0, 1, 0, 1))
plt.colorbar(label='Potential')
plt.show()

# Calculate electric field
E = np.gradient(-phi)

# Plot the electric field
plt.quiver(E[1], -E[0])
# Invert y-axis
# Because quiver interprets the starting point at the bottom left corner, while imshow starts at the top left corner
plt.gca().invert_yaxis()
plt.show()

# Comparation of numerical with the analytical result
# Analytical function
def analytical_phi(x, y):
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
def calc_ana_array(N):
    result = np.zeros((N, N))
    for i in range(N):
        for k in range(N):
            result[N - k - 1, i] = analytical_phi(i * delta, k * delta)
    return result  

ana = calc_ana_array(N)

# Plot the analytical solution
plt.imshow(ana, cmap='viridis', extent=(0, 1, 0, 1))
plt.colorbar(label='Analytical solution')
plt.show()

# Calculate electric field
E = np.gradient(-ana)

# Plot the electric field
plt.quiver(E[1], -E[0])
# Invert y-axis
# Because quiver interprets the starting point at the bottom left corner, while imshow starts at the top left corner
plt.gca().invert_yaxis()
plt.show()

# Plot the difference
plt.imshow(ana - phi, cmap='viridis', extent=(0, 1, 0, 1))
plt.colorbar(label='Analytical - Numerical solution')
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
#plt.show()

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

# Calculate electric field
E = np.gradient(-phi)

# Plot the electric field
plt.quiver(E[1], -E[0])
# Invert y-axis
# Because quiver interprets the starting point at the bottom left corner, while imshow starts at the top left corner
plt.gca().invert_yaxis()
plt.show()
