#include <iostream>
#include <random>
#include <cmath>
#include <vector>

// Metropolis algorithm for the 2D Ising model
// H: external magnetic field
// T: temperature
// steps: number of steps
// returns: average magnetization

// These functions are made available in Python via the pybind11 library

double metropolis(double H, double T, int steps) {
    double kb = 1.0;  // Boltzmann constant
    double beta = 1.0 / (kb * T);

    // Initialize spin (up or down)
    int sigma = 1;

    // Initialize magnetization
    double m = 0.0;

    // Random number generation (Mersenne Twister; slow but good)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < steps; ++i) {
        // Proposed flip
        int sigma_new = -sigma;

        // Calculate energy change
        double delta_E = -sigma_new*H - (-sigma*H);

        // Metropolis condition
        if (delta_E < 0) {
            sigma = sigma_new;
        } else if (dis(gen) < exp(-beta*delta_E)) {
            sigma = sigma_new;
        }

        // Update magnetization
        m += sigma;
    }

    return m / steps;  // return average magnetization
}

// Iterate over metropolis() for a given range of magnetic fields (come in as an array)
std::vector<double> metropolis_range(double T, int steps, std::vector<double> H_range) {
    // Initialize vector and reserve memory
    // (reserve vs. resize: https://stackoverflow.com/questions/7397768/what-is-the-difference-between-vector-resize-and-vector-reserve)
    std::vector<double> m_range;
    m_range.reserve(H_range.size());

    for (int i = 0; i < H_range.size(); ++i) {
        m_range.push_back(metropolis(H_range[i], T, steps));

        // Progress bar (every 1%)
        if (i % (H_range.size() / 100) == 0) {
            int progress = round((double)i / H_range.size() * 100);
            std::cout << "\033[1;32mProgress: " << progress << "%\033[0m\r";
            std::cout.flush();
        }
    }

    std::cout << "\033[1;32mProgress: 100%\033[0m\n";
    return m_range;
}