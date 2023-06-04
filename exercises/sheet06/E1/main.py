import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def generate_bifurcation_data(r_values, map_func, x_init=0.5, N_warmup=10000, N_iter=64) -> pd.DataFrame:
    df = pd.DataFrame()
    for r in r_values:
        x = iterate_map(r, x_init, map_func, N_warmup)
        
        # Count orbits by storing unique values in a set
        orbits = []
        for _ in range(N_iter):
            x = map_func(r, x)
            orbits.append(np.round(x, 5))  # Floats need to be rounded to make them comparable

        # Save all data related to current r value
        # print(f"r: {r}, orbits: {orbits}\n")
        count = len(orbits)
        df_temp = pd.DataFrame(
            {"r": np.full_like(orbits, r),
             "y": orbits,
             "orbits": np.full_like(orbits, count)})
        df_temp = df.drop_duplicates()
        df = pd.concat([df, df_temp])
    return df

def clean_bifurcation_data(r, x_init, orbits):
    """Cleans the data, deleting duplicates and appending relevant parameters."""
    df = pd.DataFrame()
    return df


def logistic_map(r, x):
    return r * x * (1 - x)

def cubic_map(r, x):
    return r * x - x ** 3

def iterate_map(r, x, map_func, N):
    """Perform N iterations of the given map function."""
    for _ in range(N):
        x = map_func(r, x)
    return x

r_values = np.linspace(0, 4, num=20)

# Generate data for the logistic map
orbit_counts_logistic = generate_bifurcation_data(r_values, logistic_map)
orbit_counts_logistic.plot.scatter(x='r', y='y')

# # Generate data for the cubic map
# r_values_cubic = np.linspace(0, (2/3)**0.5, num=1000)
# orbit_counts_cubic = generate_bifurcation_data(r_values_cubic, cubic_map)

# # Plot bifurcation diagrams
# plt.figure(figsize=(10, 7))
# plt.subplot(2, 1, 1)
# plt.plot(r_values, orbit_counts_logistic, '.')
# plt.title('Bifurcation diagram of the logistic map')
# plt.xlabel('r')
# plt.ylabel('Number of Orbits')

# plt.subplot(2, 1, 2)
# plt.plot(r_values_cubic, orbit_counts_cubic, '.')
# plt.title('Bifurcation diagram of the cubic map')
# plt.xlabel('r')
# plt.ylabel('Number of Orbits')

# plt.tight_layout()
# plt.show()

# # Plot bifurcation diagrams
# plt.figure(figsize=(10, 7))
# plt.subplot(2, 1, 1)
# plt.plot(r_values, orbit_counts_logistic, '.')
# plt.title('Bifurcation diagram of the logistic map')
# plt.xlabel('r')
# plt.ylabel('Number of Orbits')

# plt.subplot(2, 1, 2)
# plt.plot(r_values_cubic, orbit_counts_cubic, '.')
# plt.title('Bifurcation diagram of the cubic map')
# plt.xlabel('r')
# plt.ylabel('Number of Orbits')

# plt.tight_layout()
# plt.show()