import numpy as np

# Metropolis algorithm for a single spin
def metropolis(T, steps, H):
    kb = 1.0  # Boltzmann constant
    beta = 1.0 / (kb * T)

    # Initialize spin
    sigma = 1

    # Initialize magnetization
    m = 0.0

    for i in range(steps):
        # Proposed flip
        sigma_new = -sigma

        # Calculate energy change
        delta_E = -sigma_new*H - (-sigma*H)

        # Metropolis condition
        if delta_E < 0:
            sigma = sigma_new
        elif np.random.rand() < np.exp(-beta*delta_E):
            sigma = sigma_new

        # Update magnetization
        m += sigma

    return m / steps  # return average magnetization

# Perform Monte Carlo simulations
def metropolis_range(T, steps, H_values):
    m_values = np.zeros(H_values.shape)

    for i, H in enumerate(H_values):
        m_values[i] = metropolis(T, steps, H)
        
        # Print progress (in a single line)
        print(f'\r\033[1;32mProgress: {int(100 * (i + 1) / len(H_values))}%\033[0m', end='')
    
    # Print newline
    print()
    return m_values