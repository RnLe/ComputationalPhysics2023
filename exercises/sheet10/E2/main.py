# C++ implementation: exercises/sheet10/E1/metropolis.h
# Remove the next line if no pybind11 is installed
import Solver

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal    # For convolution

# Create a 2D grid of spins
def create_lattice(n, random=True):
    if not random:
        return -np.ones((n, n), dtype=np.int8)
    else:
        return np.random.choice([1, -1], size=(n, n))

# Calculate energy of a configuration
def energy(lattice):
    # We use the convolution theorem to calculate the energy
    # Define kernel
    kernel = np.ones((3, 3))
    # Set corners to 0
    kernel[0, 0] = 0
    kernel[0, 2] = 0
    kernel[2, 0] = 0
    kernel[2, 2] = 0
    # Set center to 0
    kernel[1, 1] = 0    # Set center to 0
    
    # Convolve with lattice
    neighbors_sum = signal.convolve2d(lattice, kernel, mode='same', boundary='wrap')
    
    return -np.sum(lattice * neighbors_sum)

# Metropolis algorithm
def metropolis(lattice, T):
    n = lattice.shape[0]
    
    # Perform sweeps N = n**2 times (in one sweep)
    for _ in range(n**2):
        # Choose a random spin
        # (Possible optimization: loop over all spins instead of choosing a random one; random number generation is expensive)
        x, y = np.random.randint(0, n, 2)

        # Calculate energy change if this spin is flipped
        # Here, we use periodic boundary conditions by using the modulo operator to wrap around the lattice indices if they are out of bounds (e.g. if x = 0, then x-1 = -1, but we want it to be n-1)
        E = lattice[x, y] * (lattice[(x+1)%n, y] + lattice[(x-1)%n, y]                      # Right and left neighbor
                                + lattice[x, (y+1)%n] + lattice[x, (y-1)%n]                 # Top and bottom neighbor
                                + lattice[(x+1)%n, (y+1)%n] + lattice[(x-1)%n, (y-1)%n]     # Top right and bottom left neighbor
                                + lattice[(x+1)%n, (y-1)%n] + lattice[(x-1)%n, (y+1)%n])    # Top left and bottom right neighbor
        
        delta_E = -2 * E
        # Metropolis condition
        if delta_E < 0 or np.random.rand() < np.exp(-delta_E/T):
            lattice[x, y] = -lattice[x, y]

# Simulation
def simulation_a(T, sweeps):
    # Initialize lattice
    lattice = create_lattice(100)
    energies = np.empty(sweeps)
    
    # First snapshot
    plt.imshow(lattice, cmap='gray')
    plt.title(f"Initial configuration for T={T}")
    plt.tight_layout()
    plt.savefig(f"Initial_configuration_T={T}.pdf")
    plt.show()
    
    # Equilibration
    # Perform 100 sweeps to equilibrate the system
    for _ in range(0):
        lattice = np.array(Solver.metropolis(lattice, T))
        
        # Print progress (in a single line)
        print(f'\r\033[1;32mEquilibration progress: {int(100 * (_ + 1) / sweeps)}%\033[0m        ', end='')

    # Perform sweeps and save energy
    # Save lattice for 5 snapshots, equidistant in time (sweeps)
    for i in range(sweeps):
        lattice = np.array(Solver.metropolis(lattice, T))
        energies[i] = energy(lattice)
            
        if i % (sweeps // 5) == 0:
            plt.imshow(lattice, cmap='gray')
            plt.title(f"Configuration for T={T} at sweep {i}\nNo equilibration sweeps were performed")
            plt.tight_layout()
            plt.savefig(f"snapshot{i}, T={T}.pdf")
            plt.show()
            
        # Print progress (in a single line)
        print(f'\r\033[1;32mSimulation progress: {int(100 * (i + 1) / sweeps)}%\033[0m          ', end='')

    return energies / lattice.size, lattice

# Examining equilibration phase
def equilibration_only(T, sweeps=10000):
    lattice = create_lattice(100)
    energies = np.empty(sweeps)
    print()
    
    # Equilibration
    # Use two criteria for equilibration: Tolerance or time
    # If energy changes by less than 0.1% in 500 sweeps, we assume that the system has equilibrated
    # Otherwise, we abort equilibration after 30 seconds
    tolerance = 0.001
    start_time = time.time()
    end_sweep = sweeps
    found = False
    for i in range(sweeps):
        lattice = np.array(Solver.metropolis(lattice, T))
        energies[i] = energy(lattice)
        # Normalize energy by number of spins
        energies[i] /= lattice.size
        
        # Initialize energy change
        energy_change = 100
        
        # Check if energy has changed by less than 0.1% in the last 500 sweeps
        if i > 500:
            energy_change = np.abs(energies[i] - energies[i-500]) / np.abs(energies[i])
        
        # Print progress (dual output: a progress for tolerance and time (shown as 'energy_change in %' and 'time in seconds / 60'))
        # Also, reset line after printing
        print(f'\r\033[1;32mEquilibration progress for T = {T}; Energy change: {energy_change*100:.2f}%, must fall below {tolerance*100}%; Time limit: {(time.time() - start_time):.2f}s / 60 s\033[0m                       ', end='')

                 
        if not found and energy_change < tolerance:
            found = True
            end_sweep = i
            # break
        
        # # Check if 30 seconds have passed
        # if time.time() - start_time > 60:
        #     end_sweep = i
        #     break
        

    print()
    print(f'\033[1;32mEquilibration finished after {end_sweep} sweeps\033[0m')
    return energies, end_sweep

# Calculate magnetization of a configuration
def magnetization(lattice):
    return np.sum(lattice)

# Simulation
# In this method, we calculate the average energy, average magnetization and absoulte average magnetization for a given temperature as function of time
# This means we need to calculate the energy and magnetization for each sweep
def simulation_c(T, sweeps=10000):
    # Initialize lattice
    n = 100
    lattice = create_lattice(n)
    energies = np.empty(sweeps)
    energy_squares = np.empty(sweeps)
    magnetizations = np.empty(sweeps)
    
    print()
    
    # Equilibration
    # We somewhat observed equilibration after 500 sweeps for T = 3.0 and 1500 for T = 1.0
    # Linearly interpolate between these values (changing 500 to 1000 and 1500 to 2500 for more precision)
    equilibration_sweeps = int(1000 + (T - 1.0) * (2500 - 1000) / (3.0 - 1.0))
    for _ in range(1):
        lattice = np.array(Solver.metropolis(lattice, T))
        
        # Print progress (in a single line)
        print(f'\r\033[1;32mEquilibration progress for T={T}: {int(100 * (_ + 1) / equilibration_sweeps)}%\033[0m', end='')

    print()
    
    # Perform sweeps and save energy
    for i in range(sweeps):
        lattice = np.array(Solver.metropolis(lattice, T))
        E = energy(lattice)
        energies[i] = E
        energy_squares[i] = E**2
        magnetizations[i] = magnetization(lattice)
            
        # Print progress (in a single line)
        print(f'\r\033[1;32mSimulation progress for T={T}: {int(100 * (i + 1) / sweeps)}%\033[0m', end='')

    return energies, energy_squares, magnetizations, equilibration_sweeps

if __name__ == '__main__':
    # Get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', action='store_true', help='Task 1')
    parser.add_argument('-b', action='store_true', help='Task 2')
    parser.add_argument('-c', action='store_true', help='Task 3')
    parser.add_argument('-d', action='store_true', help='Task 4')
    parser.add_argument('--all', action='store_true', help='All tasks')
    # Argument to either simulate or load results from file
    parser.add_argument('--load', action='store_true', help='Plot results from file')
    args = parser.parse_args()
    
    # Show command line arguments if no task is specified
    if not args.a and not args.b and not args.c and not args.d:
        print("Please specify a task with -a, -b, -c, -d or --all")
    
    if args.a or args.all:
        # 1)
        # Perform simulations for T = 1, 3
        T = [1, 3]
        sweeps = 5000

        energies, lattice = simulation_a(T[0], sweeps=sweeps)

        # Plot the lattice
        plt.imshow(lattice, cmap='gray')
        plt.title(f"Final configuration of the lattice for T={T[0]}")
        plt.tight_layout()
        plt.savefig(f"final_configuration_T={T[0]}.pdf")
        plt.close()

        energies, lattice = simulation_a(T[1], sweeps=sweeps)

        # Plot the lattice
        plt.imshow(lattice, cmap='gray')
        plt.title(f"Final configuration of the lattice for T={T[1]}")
        plt.tight_layout()
        plt.savefig(f"final_configuration_T={T[1]}.pdf")
        plt.close()
        
    if args.b or args.all:
        # 2)
        # Define boundaries for the temperatures
        T = [1.5, 3]
        sweeps = 10000

        # Perform simulations for different temperatures (only 2)
        for temperatures in T:
            energies, end_sweep = equilibration_only(temperatures, sweeps=sweeps)
            
            # Plot the average energy per spin as a function of time
            plt.plot(energies, label=f"T={temperatures}")
            plt.plot(end_sweep, energies[end_sweep], 'o', label=f"T={temperatures}, equilibrated with 0.1% tolerance")
            
        plt.title("Average energy per spin during equilibration phase")
        plt.xlabel("Sweep")
        plt.ylabel("Average energy per spin")
        plt.legend()
        plt.tight_layout()
        plt.savefig("equilibration_b.pdf")
        plt.close()

    if args.c or args.all:
        # 3)
        # Define boundaries for the temperatures
        # (T_critial - 0.3, T_critial + 0.3) is a high resolution range around the critical temperature
        T = [1.5, 3]
        T_critial = 2.269
        T_range = np.linspace(T[0], T[1], 10)
        T_range = np.append(T_range, T_critial)
        # Sort the temperatures
        T_range.sort()
        
        # Create arrays for averages to hold the averages for each temperature
        E_avgs = []
        m_avgs = []
        m_abs_avgs = []
        c_avgs = []
        equi_sweeps = []
        
        sweeps = 200000
        n = 100
        # Perform simulations for different temperatures
        # Respect flag to either perform simulations or load data
        if not args.load:   
            for temperatures in T_range:
                energies, energy_squares, magnetizations, equilibration_sweeps = simulation_c(temperatures, sweeps=sweeps)
                equi_sweeps.append(equilibration_sweeps)
                # Print the current sweep
                print("\nCalculating averages for T =", temperatures, "\n")
                # Calculate the average energy, average magnetization and absoulte average magnetization as a function of time
                E_avg = np.cumsum(energies) / (np.arange(len(energies)) + 1)
                E_avgs.append(E_avg / n**2)
                m_avg = np.cumsum(magnetizations) / (np.arange(len(magnetizations)) + 1)
                m_avgs.append(m_avg / n**2)
                m_abs_avg = np.cumsum(np.abs(magnetizations)) / (np.arange(len(magnetizations)) + 1)
                m_abs_avgs.append(m_abs_avg / n**2)
                # Calculate the specific heat
                # $c(T)=\frac{\left\langle\mathcal{H}^2\right\rangle-\langle\mathcal{H}\rangle^2}{k_{\mathrm{B}} T^2 N}$, where $N$ is the number of spins
                c = (np.cumsum(energy_squares) / (np.arange(len(energy_squares)) + 1))
                c -= E_avg**2
                c /= temperatures**2
                c_avgs.append(c / n**2)
                
            # Write data to file for each average
            # The columns are the temperatures and the rows are the sweeps
            with open("data_c_E.txt", "w") as f:
                for t in T_range:
                    f.write(f"{t} ")
                f.write("\n")
                for i in range(len(E_avgs[0])):
                    for j in range(len(E_avgs)):
                        f.write(f"{E_avgs[j][i]} ")
                    f.write("\n")
            with open("data_c_m.txt", "w") as f:
                for t in T_range:
                    f.write(f"{t} ")
                f.write("\n")
                for i in range(len(m_avgs[0])):
                    for j in range(len(m_avgs)):
                        f.write(f"{m_avgs[j][i]} ")
                    f.write("\n")
            with open("data_c_m_abs.txt", "w") as f:
                for t in T_range:
                    f.write(f"{t} ")
                f.write("\n")
                for i in range(len(m_abs_avgs[0])):
                    for j in range(len(m_abs_avgs)):
                        f.write(f"{m_abs_avgs[j][i]} ")
                    f.write("\n")
            with open("data_c_c.txt", "w") as f:
                for t in T_range:
                    f.write(f"{t} ")
                f.write("\n")
                for i in range(len(c_avgs[0])):
                    for j in range(len(c_avgs)):
                        f.write(f"{c_avgs[j][i]} ")
                    f.write("\n")
            
        # Create plot for each average
        # T_range holds multiple temperatures. Plot all averages of one type in one plot.
        # Format the lines until equilibrium as dashed. The time sweep where equilibrium is reached is stored in equi_sweeps.
        # In total, there are 4 plots. Each .._avgs array holds the averages for one type.
        
        # Define colors for the plots
        colors = plt.cm.coolwarm(np.linspace(0, 1, len(T_range)))
        
        equilibration_sweeps = np.array(1000 + (T_range - 1.0) * (2500 - 1000) / (3.0 - 1.0), dtype=int)
        
        # Flag to either plot from simutlation or from file
        if args.load:
            # Load data from file
            E_avgs = np.loadtxt("data_c_E.txt", skiprows=1, unpack=True)
            m_avgs = np.loadtxt("data_c_m.txt", skiprows=1, unpack=True)
            m_abs_avgs = np.loadtxt("data_c_m_abs.txt", skiprows=1, unpack=True)
            c_avgs = np.loadtxt("data_c_c.txt", skiprows=1, unpack=True)
            equi_sweeps = equilibration_sweeps
        
        for avgs, name in zip([E_avgs, m_avgs, m_abs_avgs, c_avgs], ["E", "m", "|m|", "c"]):
            for i, temperature in enumerate(T_range):
                # Plot the average until equilibrium as dashed
                # Plot the critical temperature differently
                if temperature == T_critial:
                    plt.plot(avgs[i], label=f"T={round(temperature, 3)}", color='k')
                else:
                    plt.plot(np.arange(equi_sweeps[i]), avgs[i][:equi_sweeps[i]], '--', color=colors[i])
                    plt.plot(np.arange(equi_sweeps[i], len(avgs[i]) - 1), avgs[i][equi_sweeps[i]:-1], label=f"T={round(temperature, 3)}", color=colors[i])
                    
            # For the specific heat, show only the x range 10e4 to end
            if name == "c":
                plt.title(f"Average ${name}$ per spin")
                plt.xlabel("Time sweep")
                plt.xscale("log")
                plt.xlim(10e4, len(avgs[0]) - 1)
                plt.ylim(0, 12)
                plt.ylabel(f"Average ${name}$ per spin")
                plt.legend()
                plt.tight_layout()
                plt.savefig(f"equilibration_c_d_{name}.pdf")
                plt.show()
            else:
                plt.title(f"Average ${name}$ per spin")
                plt.xlabel("Time sweep")
                plt.xscale("log")
                plt.ylabel(f"Average ${name}$ per spin")
                plt.legend()
                plt.tight_layout()
                plt.savefig(f"equilibration_c_d_{name}.pdf")
                plt.show()