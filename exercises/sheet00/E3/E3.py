import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# LaTeX in axis-labels
matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

# Analytical solution for the given differential equation
def y(t):
    return np.exp(-t)

# Setup
interval = (0, 10)
steps = 100
step = interval[-1] / steps

# a)
# Initial conditions
y_e_0 = 1

# Euler Method
y_e_a = [y_e_0]
for i in range(steps):
    y_e_a.append(y_e_a[i] * (1 - step))

# Symmetric Euler Method
y_se_a = [y_e_0, np.exp(-step)]
for i in range(steps - 1):
    y_se_a.append(-2 * step * y_se_a[i + 1] + y_se_a[i])

# Plot results for part a)
t = np.linspace(interval[0], interval[-1], 200)
t_d = np.array(range(steps + 1)) / interval[-1]

plt.plot(t_d, y_se_a, 'k.-', label="Symmetric Euler Method")
plt.plot(t_d, y_e_a, 'x', color='tab:orange', label="Euler Method")
plt.plot(t, y(t), label="Analytical solution")
plt.legend()
plt.xlabel(r'$\Delta t$')
plt.ylabel('arbitrary units')
plt.savefig('3a.pdf')
plt.close()

# b)
# Initial conditions
y_e_0 = 1 - step

# Euler Method with alternative initial condition
y_e_b = [y_e_0]
for i in range(steps):
    y_e_b.append(y_e_b[i] * (1 - step))

# Symmetric Euler Method with alternative initial conditions
y_se_b = [1, 1 - step]
for i in range(steps - 1):
    y_se_b.append(-2 * step * y_se_b[i + 1] + y_se_b[i])

fig = plt.figure()

# Plot Symmetric Euler Method
ax1 = fig.add_subplot(2, 1, 1)
ax1.plot(t, y(t), label="Analytical solution")
ax1.plot(t_d, y_se_a, 'k.-', label="Symmetric Euler Method")
ax1.plot(t_d, y_se_b, '.-', color='tab:grey', label="Symmetric Euler Method, b)")
ax1.legend()
ax1.set_xlabel(r'$\Delta t$')
ax1.set_ylabel('arbitrary units')

# Plot Euler Method
ax2 = fig.add_subplot(2, 1, 2)
ax2.plot(t, y(t), label="Analytical solution")
ax2.plot(t_d, y_e_a, 'x', color='tab:orange', label="Euler Method")
ax2.plot(t_d, y_e_b, 'x', color='#dea768', label="Euler Method b)")
ax2.legend()
ax2.set_xlabel(r'$\Delta t$')
ax2.set_ylabel('arbitrary units')

plt.tight_layout()
plt.savefig('3b.pdf')
plt.close()