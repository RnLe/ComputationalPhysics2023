import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

def generate_bifurcation_data(r_values, map_func, x_init=0.5, N_warmup=1000, N_iter=64) -> pd.DataFrame:
    df = pd.DataFrame()

    # Define ANSI escape codes for bold and green text
    bold_green_start = "\033[1;32m"
    bold_end = "\033[0m"

    for r in tqdm(r_values, bar_format=f'{bold_green_start}{{l_bar}}{{bar:20}}{{r_bar}}{{bar:-20b}}{bold_end}'):

        x = iterate_map(r, x_init, map_func, N_warmup)
        
        # Count orbits by storing unique values in a set
        orbits = set()
        for _ in range(N_iter):
            x = map_func(r, x)
            orbits.add(np.round(x, 6))  # Floats need to be rounded to make them comparable

        # Save all data related to current r value
        count = len(orbits)
        orbits = np.array(list(orbits))
        df_temp = pd.DataFrame(
            {"r": np.full_like(orbits, r),
             "y": orbits,
             "orbits": np.full_like(orbits, count, dtype=int)})
        df_temp = df_temp.drop_duplicates()
        df = pd.concat([df, df_temp], ignore_index=True)
    return df

def iterate_map(r, x, map_func, N):
    """Perform N iterations of the given map function."""
    for _ in range(N):
        x = map_func(r, x)
    return x