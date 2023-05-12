import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from tqdm import tqdm

# List of position data files
files_r = ["c)r0.01.dat", "c)r1.dat", "c)r100.dat"]

# List of set data files
files_set = ["c)set0.01.dat", "c)set1.dat", "c)set100.dat"]

# Create a figure and a set of subplots
fig, axs = plt.subplots(2, 3, figsize=(15, 10))

# Set up the scatter plots and text boxes for each subplot
scatters = []
texts = []

# Set up the line plots for each subplot in the second row
lines = []
# Define twin axes globally here
axs_twin1 = []
axs_twin2 = []

for ax in axs[0, :]:
    scatter = ax.scatter([], [], animated=True)
    text = ax.text(0.05, 0.95, '', transform=ax.transAxes, va='top', ha='left')
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 16)    
    ax.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    scatters.append(scatter)
    texts.append(text)

# Plot the lines on the twin axes, not just ax
for ax in axs[1, :]:
    ax2 = ax.twinx()
    ax3 = ax.twinx()
    axs_twin1.append(ax2)
    axs_twin2.append(ax3)
    line1, = ax.plot([], [], 'x', label='T', color='tab:blue', markersize=3)
    line2, = ax2.plot([], [], label='Ekin', color='tab:orange')
    line3, = ax3.plot([], [], label='Epot', color='tab:green')
    lines.append((line1, line2, line3))

# Initialize the data arrays
data_r = []
data_set = []

# Read the data from the files
for file_r, file_set in zip(files_r, files_set):
    data_r.append(np.genfromtxt(file_r, skip_header=True, delimiter=','))
    data_set.append(np.genfromtxt(file_set, skip_header=True, delimiter=',', unpack=True))

# Rewrite this init function to use the global twin axes
def init():
    # Set up the labels and legends for the lower plots
    for ax, single_data_set, (line1, line2, line3), ax2, ax3 in zip(axs[1, :], data_set, lines, axs_twin1, axs_twin2):
        time, T, Ekin, Epot, vx_, vy_ = single_data_set
        ax.set_xlim(time[0], time[-1])
        ax.set_ylim(np.min(T), np.max(T))
        ax.set_xlabel('Time')
        ax.set_ylabel('T', color='tab:blue')
        ax2.set_ylabel('Ekin', color='tab:orange')
        ax2.set_ylim(np.min(Ekin), np.max(Ekin))
        ax3.set_ylabel('Epot', color='tab:green')
        ax3.spines['right'].set_position(('outward', 60))
        ax3.set_ylim(np.min(Epot), np.max(Epot))
        ax.legend(loc='upper left')

    for scatter in scatters:
        scatter.set_offsets(np.empty((0, 2)))  # Provide an empty 2D array
    for text in texts:
        text.set_text('')
    for line1, line2, line3 in lines:
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
    return scatters + texts + [item for sublist in lines for item in sublist]  # Include the new text objects in the return statement

def update(frame):
    for scatter, data, text, single_data_set, (line1, line2, line3) in zip(scatters, data_r, texts, data_set, lines):
        xpos = data[frame, ::2]
        ypos = data[frame, 1::2]
        scatter.set_offsets(np.c_[xpos, ypos])
        time, T, Ekin, Epot, vx_, vy_ = single_data_set
        text.set_text(f't={time[frame]:.2f} T={T[frame]:.2f}\nEkin={Ekin[frame]:.2f} Epot={Epot[frame]:.2f}')
        line1.set_data(time[2:frame], T[2:frame])
        line2.set_data(time[2:frame], Ekin[2:frame])
        line3.set_data(time[2:frame], Epot[2:frame])

    return scatters + texts + [item for sublist in lines for item in sublist]  # Include the new text objects in the return statement

ani = animation.FuncAnimation(fig, update, frames=min(len(data) for data in data_r), init_func=init, blit=True, interval=1)

# Adjust layout to avoid overlapping
fig.tight_layout()

plt.show()