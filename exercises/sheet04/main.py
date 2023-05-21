import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
from calculation import iterate, calc_ana_array
from plotting import plot_potential, plot_field

rc('text', usetex=True)

# LaTeX in axis-labels
matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

def plot_results(phi: np.ndarray, E: np.ndarray, N: int, title: str, filename: str, label_potential: str, label_field: str):
    # Create a mosaic subplot
    fig, axs = plt.subplot_mosaic("AB", figsize=(10, 6))

    # Plot the potential
    im = plot_potential(phi, axs['A'], fig, label_potential)

    # Plot the electric field
    plot_field(E, axs['B'], fig, im, N, label_field)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(filename)
    plt.show()


def main():  
    # ::::::::::::::::::::::::::::::::::::::::::::
    # a)

    # Parameters
    delta = 0.05
    kappa = 1e-5
    N = int(1 / delta) + 1

    print(f"Step size Delta: {delta}")
    print(f"Error bound Kappa: {kappa}")
    print(f"Discretization N: {N}")
        
    # ::::::::::::::::::::::::::::::::::::::::::::
    # b)
            
    # Initialize arrays
    phi = np.ones((N, N))  # potential
    rho = np.zeros((N, N))  # charge density

    # Boundary conditions
    phi[0, :] = 0
    phi[-1, :] = 0
    phi[:, 0] = 0
    phi[:, -1] = 0

    # Perform the Gauss-Seidel method
    iterate(phi, rho, kappa, N, delta)

    # Calculate the electric field
    E = np.gradient(-phi)
    
    # Plotting
    # ::::::::::::::::::::::::::::::::::::::::::::
    plot_results(phi, E,  N,"Uniform, positive potential with uniform boundary conditions\n" + r"\textbf{Gauss-Seidel algorithm}", "b.pdf", "Potential", "Electric Field")

    # ::::::::::::::::::::::::::::::::::::::::::::
    # c)

    # Reset potential
    phi = np.zeros((N, N)) 

    # Boundary conditions
    phi[0, :] = 1
    phi[-1, :] = 0
    phi[:, 0] = 0
    phi[:, -1] = 0

    # Iterate with Gauss-Seidel method
    iterate(phi, rho, kappa, N, delta)

    # Calculate the electric field
    E = np.gradient(-phi)
    
    # Plotting
    # ::::::::::::::::::::::::::::::::::::::::::::
    plot_results(phi, E, N, "Upper boundary potential\n" + r"\textbf{Gauss-Seidel algorithm}", "c_numerical.pdf", "Potential - Numerical", "Electric Field - Numerical")

    # :::::::::::::::::::::::::::::::::::::::::::
    
    # Calculate the analytical solution
    ana = calc_ana_array(N, delta)

    # Calculate the electric field
    E = np.gradient(-ana)
    
    # Plotting
    # ::::::::::::::::::::::::::::::::::::::::::::
    plot_results(ana, E, N, "Analytical Solution to the upper boundary potential", "c_analytical.pdf", "Potential - Analytical", "Electric Field - Analytical")
    
    # Plot the difference in a single plot
    plt.imshow(ana - phi, cmap='viridis', extent=(0, 1, 0, 1))
    plt.colorbar(label='Analytical - Numerical solution')
    plt.title("Difference: Analytical - Numerical\nUpper boundary potential")
    plt.savefig("c_difference.pdf")
    plt.show()
    
    # :::::::::::::::::::::::::::::::::::::::::::
    # d)

    # Reset potential and charge density
    phi = np.zeros((N, N)) 
    rho = np.zeros((N, N)) 

    # Place charge at the center
    rho[N//2, N//2] = 1

    # Boundary conditions
    phi[0, :] = 0
    phi[-1, :] = 0
    phi[:, 0] = 0
    phi[:, -1] = 0

    # Iterate with Gauss-Seidel method
    iterate(phi, rho, kappa, N, delta)

    # Calculate the electric field
    E = np.gradient(-phi)
    
    # Plotting
    # ::::::::::::::::::::::::::::::::::::::::::::
    plot_results(phi, E, N, "Single Charge, center\n" + r"\textbf{Gauss-Seidel algorithm}", "d_single.pdf", "Potential", "Electric Field")

    # ::::::::::::::::::::::::::::::::::::::::::::
    # e)
    
    # Reset potential and charge density
    phi = np.zeros((N, N)) 
    rho = np.zeros((N, N)) 

    # Place charges
    rho[int(0.25/delta), int(0.25/delta)] = 1
    rho[int(0.75/delta), int(0.75/delta)] = 1
    rho[int(0.25/delta), int(0.75/delta)] = -1
    rho[int(0.75/delta), int(0.25/delta)] = -1

    # Boundary conditions
    phi[0, :] = 0
    phi[-1, :] = 0
    phi[:, 0] = 0
    phi[:, -1] = 0

    # Iterate with Gauss-Seidel method
    iterate(phi, rho, kappa, N, delta)

    # Calculate the electric field
    E = np.gradient(-phi)
    
    # Plotting
    # ::::::::::::::::::::::::::::::::::::::::::::
    plot_results(phi, E, N, "Four Charges\n" + r"\textbf{Gauss-Seidel algorithm}", "e_four_charges.pdf", "Potential - Four Charges", "Electric Field - Four Charges")

    # :::::::::::::::::::::::::::::::::::::::::::


if __name__ == "__main__":
    main()