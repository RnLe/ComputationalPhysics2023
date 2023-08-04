# C++ implementation: exercises/sheet10/E1/metropolis.h
# Remove the next line if no pybind11 is installed
import Solver
# Python implementation: exercises/sheet10/E1/metropolis.py
from metropolis import metropolis_range

import numpy as np
import matplotlib.pyplot as plt
import argparse
import time
import sys

# This script has three modes:
# 1. Run the Monte Carlo simulation in C++ (and plot the results)
# 2. Run the Monte Carlo simulation in Python (and plot the results)
# 3. Run the Monte Carlo simulation in Python and C++ and compare the runtimes (and plot the results)

# Main function (with command line arguments)
if __name__ == '__main__':
    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--python', action='store_true', help='Run the Python version of the Monte Carlo simulation')
    parser.add_argument('--cpp', action='store_true', help='Run the C++ version of the Monte Carlo simulation')
    parser.add_argument('--compare', action='store_true', help='Compare the Python and C++ versions of the Monte Carlo simulation')
    
    # List options of parameters when no command line arguments are given
    if len(sys.argv) == 1:
        print('No command line arguments given. Use one of the following options:')
        print('python main.py --python')
        print('python main.py --cpp')
        print('python main.py --compare')
        sys.exit()
    
    # Set parameters
    # HINT: For demonstation purposes, the step size is set to 10**2. For the actual exercise, use 10**4.
    H_values = np.linspace(-5, 5, 10**4)
    T = 1.0  # temperature
    steps = 10**5  # number of Monte Carlo steps

    # Initialize array to store magnetization values
    m_values = np.zeros(H_values.shape)
    
    # Case differentiation
    args = parser.parse_args()
    if args.python:
        # Run the Python version of the Monte Carlo simulation
        m_values = metropolis_range(T, steps, H_values)
        
    # Remove this if no pybind11 is installed
    # :::::::::::::::::::::::::::::::::::::::
    elif args.compare:
        # Compare the Python and C++ versions of the Monte Carlo simulation
        print()
        start = time.time()
        m_values = metropolis_range(T, steps, H_values)
        end = time.time()
        # Print both times in bold yellow
        print('\033[93m' + 'Python runtime: {:.3f} s'.format(end - start) + '\033[0m')
        start = time.time()
        m_values = Solver.metropolis_range(T, steps, H_values)
        end = time.time()
        print('\033[93m' + 'C++ runtime: {:.3f} s'.format(end - start) + '\033[0m')
        print()
    elif args.cpp:
        # Run the C++ version of the Monte Carlo simulation
        m_values = Solver.metropolis_range(T, steps, H_values)

    # :::::::::::::::::::::::::::::::::::::::
    # Remove until here if no pybind11 is installed
    
    
    
    # Calculate analytical result
    beta = 1 / T
    m_analytical = np.tanh(beta * H_values)

    # Plot numerical and analytical results
    plt.figure()
    plt.plot(H_values, m_values, label='Numerical')
    plt.plot(H_values, m_analytical, '--', label='Analytical')
    plt.xlabel('$H$')
    plt.ylabel('$m$')
    plt.grid()
    plt.legend()
    plt.savefig('magnetization.pdf')
    plt.show()