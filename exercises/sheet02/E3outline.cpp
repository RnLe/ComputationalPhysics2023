// Basic headers
#include <iostream>
#include <cmath>
// Containers
#include <vector>
#include <tuple>
// File handling
#include <fstream>
// For the time measurement
#include <chrono>
// For mean and std calculation
#include <algorithm>
#include <numeric>

using namespace std;

const double G = 1.0;
const double m1 = 1.0;
const double m2 = 2.0;

tuple<double, double> force(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    double distance = sqrt(dx * dx + dy * dy);
    double distance_cubed = distance * distance * distance;
    double f = -G * m1 * m2 / distance_cubed;
    return {f * dx, f * dy};
}

vector<tuple<double, double, double, double, double, double, double, double>> simulate(string method, double dt, double T, double t = 0.,
                                                                                        double x1 = 0., double y1 = 1., double x2 = 0., double y2 = -0.5,
                                                                                        double vx1 = 0.8, double vy1 = 0., double vx2 = -0.4, double vy2 = 0.) {
    // Calculate initial acceleration and previous positions for Verlet method
    // ...

    vector<tuple<double, double, double, double, double, double, double, double>> results;

    while (dt > 0 ? t <= T : t >= T) {      // Checks whether dt is positive or negative (ternary conditional operator)
        // This loop contains a case distinction between the two methods
        
        // Add current positions and velocities to results vector
        // ...

        // Update positions and velocities based on chosen method (Euler or Verlet)
        // ...
    }

    return results;
}

void write_to_file(string filename, auto results) {
    // Open output file, check for errors, write data, and close the file
    // ...
}

template <typename T>
double mean(const std::vector<T> &data) {
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

template <typename T>
double standard_deviation(const std::vector<T> &data) {
    // Calculate standard deviation using mean() function
    // ...
}

int main() {
    double dt = 0.1;
    double T = 100;
    int N = 1000;
    vector<int64_t> elapsed_times_euler;
    vector<int64_t> elapsed_times_verlet;

    // a) Calculate positions, coarse and adequate step size, and write results to files
    // ...

    // b) Measure computation time for Euler and Verlet methods
    // ...

    // c) Reverse the simulation and write results to files
    // ...

    return 0;
}
