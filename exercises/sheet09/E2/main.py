import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Function to minimize
def f(r):
    return (r[..., 0] - 1.9)**2 + (r[..., 1] - 2.1)**2 + 2*np.cos(4*r[..., 0] + 2.8) + 3*np.sin(2*r[..., 1] + 0.6)


# Constants
N = 10
n_iterations = 100
omega = 0.9
c1 = 0.1
c2 = 0.05

# Initialize positions and velocities
positions = np.random.uniform(-5, 5, (N, 2))
velocities = np.random.uniform(-1, 1, (N, 2))

# Initialize personal best positions and global best position
personal_best_positions = np.copy(positions)
personal_best_values = np.apply_along_axis(f, 1, personal_best_positions)
global_best_position = personal_best_positions[np.argmin(personal_best_values)]
global_best_values = [f(global_best_position)]

personal_best_values_per_iteration = []  # list to store personal best values for each iteration

snapshots = []  # list to store positions and velocities for snapshots
snapshot_iterations = [0, 20, 40, 60, 80, 99]  # iterations at which to take snapshots

for i in range(n_iterations):
    # Update velocities
    r1, r2 = np.random.rand(2, N, 1)
    velocities = omega * velocities \
        + c1 * r1 * (personal_best_positions - positions) \
        + c2 * r2 * (global_best_position - positions)

    # Update positions
    positions += velocities

    # Update personal bests
    values = np.apply_along_axis(f, 1, positions)
    improved = values < personal_best_values    # Mask
    personal_best_positions[improved] = positions[improved]
    personal_best_values[improved] = values[improved]

    # Save personal best values for this iteration
    personal_best_values_per_iteration.append(personal_best_values.tolist())

    # Update global best
    if np.min(values) < f(global_best_position):
        global_best_position = positions[np.argmin(values)]
        
    global_best_values.append(f(global_best_position))

    # Take snapshot if this iteration is in snapshot_iterations
    if i in snapshot_iterations:
        snapshots.append((np.copy(positions), np.copy(velocities)))
        
print(f"Global minimum at {global_best_position} with f(r) = {f(global_best_position)}")

# Plotting the best function values
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(global_best_values, label="Global Best")
ax.plot([np.min(pers_best) for pers_best in personal_best_values_per_iteration], label="Personal Best")
ax.set_xlabel('Iteration')
ax.set_ylabel('Function Value')
ax.set_title('Best Function Values Over Iterations')
ax.legend()
plt.show()

# Plotting the course of the swarm
fig, axs = plt.subplots(2, 3, figsize=(18, 12))

# Generate a grid for the function
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)
Z = f(np.stack((X, Y), axis=2))

for ax, (positions, velocities), iteration in zip(axs.flatten(), snapshots, snapshot_iterations):
    ax.contourf(X, Y, Z, 50, cmap='viridis')
    ax.quiver(*positions.T, *velocities.T, color='r')
    ax.plot(*positions.T, 'bo')
    ax.plot(*personal_best_positions.T, 'go', label="Personal bests")
    ax.plot(*global_best_position, 'ro', label="Global best")
    ax.set_title(f'Iteration {iteration}')

plt.legend()
plt.tight_layout()
plt.show()
