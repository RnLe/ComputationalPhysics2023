import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import signal
from tqdm import tqdm

# Constants
L = 1
dx = 0.01
D = 1
dt = 0.0005  # needs to be small enough to satisfy the stability criterion
N = int(L / dx) + 1
steps = 10000

# Function for initialization
def init_u(condition: str) -> np.ndarray:
    x = np.linspace(0, L, N)
    if condition == 'const':
        return np.ones(N)
    elif condition == 'delta':
        return (np.abs(x - 0.5) < dx).astype(int)
    elif condition == 'heaviside':
        return (x > 0.5).astype(int)
    elif condition == 'dirac':
        return sum([signal.unit_impulse(N, int(0.1 * n / dx)) for n in range(1, 10)]) / 9
    else:
        return np.zeros(N)

# The main function for solving the task
def solve_task(condition: str):
    global u, ax, fig

    # Initialize
    u = init_u(condition)

    # Create a new figure for each animation
    fig, ax = plt.subplots()

    # Create the animation
    ani = animation.FuncAnimation(fig, update, frames=range(steps), interval=100)

    # Save the animation
    print(f"Saving {condition}_animation.mp4")
    pbar = tqdm(total=100)  # 100 corresponds to 100% completion.

    def update_progress(current_frame_number: int, total_frames: int):
        pbar.update(np.round((current_frame_number / total_frames) * 100, 2) - pbar.n)  # pbar.n is the current value

    ani.save(f'{condition}_animation.mp4', writer='ffmpeg', fps=25, progress_callback=update_progress)
    pbar.close()
    print("Saved!\n")

# Update function for the animation
def update(i: int):
    global u

    # Apply the FTCS method
    u_new = np.zeros_like(u)
    u_new[1:-1] = u[1:-1] + D * dt / dx**2 * (u[:-2] - 2*u[1:-1] + u[2:])

    # Boundary conditions (isolating)
    u_new[0] = u_new[1]
    u_new[-1] = u_new[-2]

    # Update the plot
    ax.clear()
    ax.plot(np.linspace(0, L, N), u_new, 'b')
    ax.set_ylim(0, 2)

    # Update u
    u = u_new

# Run the task for all initial conditions
for condition in ['const', 'delta', 'heaviside', 'dirac']:
    solve_task(condition)
