import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
from calculation import iterate, calc_ana_array
from plotting import plot_potential, plot_field

rc('text', usetex=True)

# LaTeX in axis-labels
matplotlib.rcParams.update({'font.size': 12, 'text.usetex': True})

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

    # Plotting
    # :::::::::::::::::::::::::::::::::::::::::::
    # Create a mosaic subplot
    fig, axs = plt.subplot_mosaic("AB", figsize=(10, 6))

    # Plot the potential
    im = plot_potential(phi, axs['A'], fig, "Potential")

    # Calculate electric field
    E = np.gradient(-phi)

    # Plot the electric field
    plot_field(E, axs['B'], fig, im, N, "Electric Field")

    fig.suptitle("Uniform, positive potential with uniform boundary conditions\n" + r"\textbf{Gauss-Seidel algorithm}")
    fig.tight_layout()
    fig.savefig("b.pdf")
    plt.show()

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

    # Plotting
    # :::::::::::::::::::::::::::::::::::::::::::
    # Create a mosaic subplot
    fig, axs = plt.subplot_mosaic("AB", figsize=(10, 6))

    # Plot the potential
    plot_potential(phi, axs['A'], fig, "Potential - Numerical")

    # Calculate electric field
    E = np.gradient(-phi)

    # Plot the electric field
    plot_field(E, axs['B'], fig, im, N, "Electric Field - Numerical")

    fig.suptitle("Upper boundary potential\n" + r"\textbf{Gauss-Seidel algorithm}")
    fig.tight_layout()
    fig.savefig("c_numerical.pdf")

    plt.show()
    # :::::::::::::::::::::::::::::::::::::::::::
    
    # Calculate the analytical solution
    ana = calc_ana_array(N, delta)

    # Plotting
    # :::::::::::::::::::::::::::::::::::::::::::
    # Plot the analytical solution
    fig, axs = plt.subplot_mosaic("AB", figsize=(10, 6))

    plot_potential(ana, axs['A'], fig, "Potential - Analytical")

    # Calculate electric field
    E = np.gradient(-ana)

    # Plot the electric field
    plot_field(E, axs['B'], fig, im, N, "Electric Field - Analytical")

    fig.suptitle("Analytical Solution to the upper boundary potential")
    fig.tight_layout()
    fig.savefig("c_analytical.pdf")

    plt.show()
    # :::::::::::::::::::::::::::::::::::::::::::
    
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

    # Plotting
    # :::::::::::::::::::::::::::::::::::::::::::
    # Create a mosaic subplot
    fig, axs = plt.subplot_mosaic("AB", figsize=(10, 6))

    # Plot the potential
    plot_potential(phi, axs['A'], fig, "Potential")

    # Calculate electric field
    E = np.gradient(-phi)

    # Plot the electric field
    plot_field(E, axs['B'], fig, im, N, "Electric Field")

    fig.suptitle("Single Charge, center\n" + r"\textbf{Gauss-Seidel algorithm}")
    fig.tight_layout()
    fig.savefig("d_single.pdf")

    plt.show()
    # :::::::::::::::::::::::::::::::::::::::::::
    
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

    # Plotting
    # :::::::::::::::::::::::::::::::::::::::::::
    # Create a mosaic subplot
    fig, axs = plt.subplot_mosaic("AB", figsize=(10, 6))

    # Plot the potential
    plot_potential(phi, axs['A'], fig, "Potential - Four Charges")

    # Calculate electric field
    E = np.gradient(-phi)

    # Plot the electric field
    plot_field(E, axs['B'], fig, im, N, "Electric Field - Four Charges")

    fig.suptitle("Four Charges\n" + r"\textbf{Gauss-Seidel algorithm}")
    fig.tight_layout()
    fig.savefig("e_four_charges.pdf")

    plt.show()
    # :::::::::::::::::::::::::::::::::::::::::::
    
    
if __name__ == "__main__":
    main()