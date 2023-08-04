import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm

# Function to minimize
def f(r):
    return (r[..., 0] - 1.9)**2 + (r[..., 1] - 2.1)**2 + 2*np.cos(4*r[..., 0] + 2.8) + 3*np.sin(2*r[..., 1] + 0.6)

N = 100  # number of particles
omega = 0.95
c1 = 0.1
c2 = 0.05

n_iterations = 100

# Initialize particles
particles = np.random.uniform(low=-5, high=5, size=(N, 2))
velocities = np.random.uniform(low=-1, high=1, size=(N, 2))

# Initialize personal bests
personal_bests = particles.copy()
personal_best_values = f(personal_bests)

# Initialize global best
global_best = personal_bests[np.argmin(personal_best_values)]
global_best_value = np.min(personal_best_values)

# Prepare for animation
fig, ax = plt.subplots(figsize=(10, 10))
X, Y = np.meshgrid(np.linspace(-5, 5, 100), np.linspace(-5, 5, 100))
Z = f(np.stack((X, Y), axis=2))
contour = ax.contourf(X, Y, Z, levels=100, cmap='viridis') # removed norm=LogNorm()
particles_plot, = ax.plot(particles[:, 0], particles[:, 1], 'ro')
global_best_plot, = ax.plot(global_best[0], global_best[1], 'bo')
quivers = ax.quiver(particles[:, 0], particles[:, 1], velocities[:, 0], velocities[:, 1], color='r')
ax.set_aspect("equal")

def update(i):
    global particles, velocities, personal_bests, personal_best_values, global_best, global_best_value
    r1 = np.random.uniform(size=(N, 1))
    r2 = np.random.uniform(size=(N, 1))
    
    velocities = omega * velocities + c1 * r1 * (personal_bests - particles) + c2 * r2 * (global_best - particles)
    particles += velocities
    
    new_values = f(particles)
    
    improved_mask = new_values < personal_best_values
    personal_bests[improved_mask] = particles[improved_mask]
    personal_best_values[improved_mask] = new_values[improved_mask]
    
    if np.min(personal_best_values) < global_best_value:
        global_best = personal_bests[np.argmin(personal_best_values)]
        global_best_value = np.min(personal_best_values)
    
    particles_plot.set_data(particles[:, 0], particles[:, 1])
    global_best_plot.set_data(global_best[0], global_best[1])
    
    quivers.set_UVC(velocities[:, 0], velocities[:, 1])
    quivers.set_offsets(particles)
    
    return particles_plot, global_best_plot, quivers,

ani = animation.FuncAnimation(fig, update, frames=n_iterations, interval=20, blit=True)

plt.show()