import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle, ConnectionPatch

# Read the data from the file
x, y = np.genfromtxt('a.txt', delimiter=', ').T
x2, y2 = np.genfromtxt('a_stable.txt', delimiter=', ').T

# Create a plot with log scales
plt.plot(x, y, label="Unstable")
plt.plot(x2, y2, label="Stable")
plt.xscale('log')
plt.yscale('log')

# Set axis labels
plt.xlabel('x')
plt.ylabel('y')

plt.legend()

#plt.savefig("a.png", dpi=300)
plt.savefig("a.pdf")

plt.close()

# Next plot b)

# Read the data from the file
x, y = np.genfromtxt('b.txt', delimiter=', ').T
x2, y2 = np.genfromtxt('b_stable.txt', delimiter=', ').T

# Create a plot with log scales
fig, ax = plt.subplots()
ax.plot(x, y, label="Unstable")
ax.plot(x2, y2, label="Stable")
ax.set_xscale('log')
ax.set_yscale('log')

# Set axis labels
ax.set_xlabel('x')
ax.set_ylabel('y')

ax.legend()

######################################
# This is solely for the inset plot (the magnification)
# Create inset plot
ax_inset = inset_axes(ax, width="30%", height="30%", loc="upper center")
ax_inset.plot(x, y, label="Unstable")
ax_inset.plot(x2, y2, label="Stable")

# Set the limits for the inset plot
ax_inset.set_xlim(0.002, 0.0003)
ax_inset.set_ylim(0.00101328, 0.000198682)

# Connect the inset plot with a magnification square
# Draw a square indicating the origin of the inset plot
x1_range = (0.002, 0.0003)
y1_range = (0.00101328, 0.000198682)
rect = Rectangle((x1_range[0], y1_range[0]), x1_range[1] - x1_range[0], y1_range[1] - y1_range[0],
                 edgecolor='black', facecolor='none', linestyle='dashed', linewidth=1)
ax.add_patch(rect)

# Add lines connecting the square with the inset plot
xyA = (x1_range[1], y1_range[0])
xyB = (0.0, 0.0)
coordsA = "data"
coordsB = "axes points"
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA=coordsA, coordsB=coordsB, axesA=ax, axesB=ax_inset, color="black", linestyle='dashed', linewidth=1)
ax.add_artist(con)

xyA = (x1_range[1], y1_range[1])
xyB = (0.0, 1.0)
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA=coordsA, coordsB=coordsB, axesA=ax, axesB=ax_inset, color="black", linestyle='dashed', linewidth=1)
ax.add_artist(con)

######################################

# Save the figure
#plt.savefig("b.png", dpi=300)
plt.savefig("b.pdf")

plt.close()

# Next plot c)

# Read the data from the file
x, y = np.genfromtxt('c.txt', delimiter=', ').T
x2, y2 = np.genfromtxt('c_stable.txt', delimiter=', ').T

# Create a plot with log scales
fig, ax = plt.subplots()
ax.plot(x, y, label="Unstable")
ax.plot(x2, y2, label="Stable")
ax.set_xscale('log')
ax.set_yscale('log')

# Set axis labels
ax.set_xlabel('$\delta$')
ax.set_ylabel('y')

ax.legend()

######################################
# This is solely for the inset plot (the magnification)
# Create inset plot
ax_inset = inset_axes(ax, width="30%", height="30%", loc="upper center")
ax_inset.plot(x, y, label="Unstable")
ax_inset.plot(x2, y2, label="Stable")

# Set the limits for the inset plot
ax_inset.set_xlim(5e-08, 2.05e-07)
ax_inset.set_ylim(1.1e-07, 1.9e-07)

# Connect the inset plot with a magnification square
# Draw a square indicating the origin of the inset plot
x1_range = (6e-08, 2e-07)
y1_range = (1.1e-07, 1.9e-07)
rect = Rectangle((x1_range[0], y1_range[0]), x1_range[1] - x1_range[0], y1_range[1] - y1_range[0],
                 edgecolor='black', facecolor='none', linestyle='dashed', linewidth=1)
ax.add_patch(rect)

# Add lines connecting the square with the inset plot
xyA = (x1_range[1], y1_range[0])
xyB = (0.0, 0.0)
coordsA = "data"
coordsB = "axes points"
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA=coordsA, coordsB=coordsB, axesA=ax, axesB=ax_inset, color="black", linestyle='dashed', linewidth=1)
ax.add_artist(con)

xyA = (x1_range[1], y1_range[1])
xyB = (0.0, 1.0)
con = ConnectionPatch(xyA=xyA, xyB=xyB, coordsA=coordsA, coordsB=coordsB, axesA=ax, axesB=ax_inset, color="black", linestyle='dashed', linewidth=1)
ax.add_artist(con)

######################################

plt.savefig("c.png", dpi=300)
plt.savefig("c.pdf")