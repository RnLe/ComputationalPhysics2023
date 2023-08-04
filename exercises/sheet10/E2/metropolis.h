#include <iostream>
#include <random>
#include <cmath>
#include <vector>

// Metropolis algorithm for the 2D Ising model
// lattice: 2D array of spins
// T: temperature

int mod(int a, int b) {
    int ret = a % b;
    if (ret < 0) {
        ret += b;
    }
    return ret;
}

// This function is made available in Python via the pybind11 library

std::vector<std::vector<int>> metropolis(std::vector<std::vector<int>> lattice, double T) {
    double beta = 1.0 / T;
    int n = lattice.size();
    int steps = n * n;

    // Random number generation (Mersenne Twister; slow but good)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::random_device rd2;
    std::mt19937 gen2(rd2());
    std::uniform_int_distribution<> dis2(0, n-1);

    for (int i = 0; i < steps; ++i) {
        // Choose a random spin
        // (Possible optimization: loop over all spins instead of choosing a random one; random number generation is expensive)
        int x = dis2(gen2);
        int y = dis2(gen2);

        // Calculate energy change if this spin is flipped
        // Here, we use periodic boundary conditions by using the modulo operator to wrap around the lattice indices if they are out of bounds (e.g. if x = 0, then x-1 = -1, but we want it to be n-1)
        // Custom modulo function to handle negative numbers
        
        // int E = -lattice[x][y] * (lattice[mod(x+1, n)][y] + lattice[mod(x-1, n)][y]                             // Right and left neighbor
        //                         + lattice[x][mod(y+1, n)] + lattice[x][mod(y-1, n)]                             // Top and bottom neighbor
        //                         + lattice[mod(x+1, n)][mod(y+1, n)] + lattice[mod(x-1, n)][mod(y-1, n)]         // Top right and bottom left neighbor
        //                         + lattice[mod(x+1, n)][mod(y-1, n)] + lattice[mod(x-1, n)][mod(y+1, n)]);       // Top left and bottom right neighbor

        int E = -lattice[x][y] * (lattice[mod(x+1, n)][y] + lattice[mod(x-1, n)][y]                             // Right and left neighbor
                                 + lattice[x][mod(y+1, n)] + lattice[x][mod(y-1, n)]);                             // Top and bottom neighbor
        
        float delta_E = -2. * E;
        // Metropolis condition
        if (delta_E < 0 || dis(gen) < exp(-beta*delta_E)) {
            lattice[x][y] = -lattice[x][y];
        }
        // std::cout << "Random number: " << dis(gen) << "delta_E: " << delta_E << "exp(-beta*delta_E): " << exp(-beta*delta_E) << "Flip accepted: " << (delta_E < 0 || dis(gen) < exp(-beta*delta_E)) << std::endl;

        

    }
    return lattice;
}