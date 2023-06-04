import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def lorenz(X, Y, Z, s=10, r=28, b=8/3):
    X_dot = -s*X + s*Y
    Y_dot = -X*Z + r*X - Y
    Z_dot = X*Y - b*Z
    return X_dot, Y_dot, Z_dot

def runge_kutta_step(X, Y, Z, dt, s, r, b):
    dx1, dy1, dz1 = lorenz(X, Y, Z, s, r, b)
    dx2, dy2, dz2 = lorenz(X + dx1*dt/2, Y + dy1*dt/2, Z + dz1*dt/2, s, r, b)
    dx3, dy3, dz3 = lorenz(X + dx2*dt/2, Y + dy2*dt/2, Z + dz2*dt/2, s, r, b)
    dx4, dy4, dz4 = lorenz(X + dx3*dt, Y + dy3*dt, Z + dz3*dt, s, r, b)

    X += dt*(dx1 + 2*dx2 + 2*dx3 + dx4)/6
    Y += dt*(dy1 + 2*dy2 + 2*dy3 + dy4)/6
    Z += dt*(dz1 + 2*dz2 + 2*dz3 + dz4)/6
    return X, Y, Z

def generate_trajectories(r_values, X_init, Y_init, Z_init, dt, T):
    trajectories = []
    for r in r_values:
        X, Y, Z = X_init, Y_init, Z_init
        X_values, Y_values, Z_values = [], [], []
        for _ in np.arange(0, T, dt):
            X_values.append(X)
            Y_values.append(Y)
            Z_values.append(Z)
            X, Y, Z = runge_kutta_step(X, Y, Z, dt, s=10, r=r, b=8/3)
        trajectories.append((X_values, Y_values, Z_values))
    return trajectories

r_values = [20, 28]
X_init, Y_init, Z_init = 1, 1, 1  # You can adjust the initial conditions here
dt = 0.01  # Time step
T = 50  # Total time

trajectories = generate_trajectories(r_values, X_init, Y_init, Z_init, dt, T)

fig, axs = plt.subplots(2, 1, figsize=(7, 10))
for ax, (r, trajectory) in zip(axs, zip(r_values, trajectories)):
    X_values, Y_values, Z_values = trajectory
    ax.plot(X_values, Y_values)
    ax.set_title(f'Projection onto the X-Y plane for r={r}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
plt.tight_layout()
plt.show()

fig, axs = plt.subplots(2, 1, figsize=(7, 10))
for ax, (r, trajectory) in zip(axs, zip(r_values, trajectories)):
    X_values, Y_values, Z_values = trajectory
    X_cut, Y_cut = [], []
    for i in range(len(Z_values) - 1):
        if Z_values[i] > 20 >= Z_values[i+1]:
            X_cut.append((X_values[i] + X_values[i+1])/2)
            Y_cut.append((Y_values[i] + Y_values[i+1])/2)
    ax.plot(X_cut, Y_cut, '.')
    ax.set_title(f'Poincaré cut for r={r}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
plt.tight_layout()
plt.show()

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
for r, trajectory in zip(r_values, trajectories):
    X_values, Y_values, Z_values = trajectory
    ax.plot(X_values, Y_values, Z_values, label=f'r={r}')
ax.legend()
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
