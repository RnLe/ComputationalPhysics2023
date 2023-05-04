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

// Gravitational constant
const double G = 1.0;

// Masses
const double m1 = 1.0;
const double m2 = 2.0;

// Force function
// The task is reduced to a 2-dimensional problem
tuple<double, double> force(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    double distance = sqrt(dx * dx + dy * dy);
    double distance_cubed = distance * distance * distance;
    double f = -G * m1 * m2 / distance_cubed;
    return {f * dx, f * dy};
}

// Simulate function
// "If you've written the same code twice or more, you can probably wrap it in a function." - Felix
// Insofar a nested apporoach is used
vector<tuple<double, double, double, double, double, double, double, double>> simulate(string method, double dt, double T, double t = 0.,
                                                                                        double x1 = 0., double y1 = 1., double x2 = 0., double y2 = -0.5,
                                                                                        double vx1 = 0.8, double vy1 = 0., double vx2 = -0.4, double vy2 = 0.) {

    // Calculate initial acceleration
    auto [ax1, ay1] = force(x1, y1, x2, y2);
    ax1 /= m1;
    ay1 /= m1;
    double ax2 = -ax1 * (m1 / m2);
    double ay2 = -ay1 * (m1 / m2);

    // Previous positions for Verlet method
    double prev_x1 = x1 - vx1 * dt + 0.5 * ax1 * dt * dt;
    double prev_y1 = y1 - vy1 * dt + 0.5 * ay1 * dt * dt;
    double prev_x2 = x2 - vx2 * dt + 0.5 * ax2 * dt * dt;
    double prev_y2 = y2 - vy2 * dt + 0.5 * ay2 * dt * dt;

    vector<tuple<double, double, double, double, double, double, double, double>> results;

    while (dt > 0 ? t <= T : t >= T) {    // Checks whether dt is positive or negative (ternary conditional operator)
        results.push_back({x1, y1, x2, y2, vx1, vy1, vx2, vy2});

        // The formulas are provided in the submission
        // For the acceleration a, Newtons Law F=ma is used
        if (method == "Euler") {            // Easy and fast, but not accurate
            auto [fx, fy] = force(x1, y1, x2, y2);
            vx1 += (fx / m1) * dt;
            vy1 += (fy / m1) * dt;
            vx2 -= (fx / m2) * dt;
            vy2 -= (fy / m2) * dt;

            x1 += vx1 * dt;
            y1 += vy1 * dt;
            x2 += vx2 * dt;
            y2 += vy2 * dt;
        } else if (method == "Verlet") {    // More accurate and stable, but slightly slower than euler
            double next_x1 = 2. * x1 - prev_x1;
            double next_y1 = 2. * y1 - prev_y1;
            double next_x2 = 2. * x2 - prev_x2;
            double next_y2 = 2. * y2 - prev_y2;

            auto [fx, fy] = force(x1, y1, x2, y2);
            next_x1 += (fx / m1) * dt * dt;
            next_y1 += (fy / m1) * dt * dt;
            next_x2 -= (fx / m2) * dt * dt;
            next_y2 -= (fy / m2) * dt * dt;

            prev_x1 = x1;
            prev_y1 = y1;
            prev_x2 = x2;
            prev_y2 = y2;

            x1 = next_x1;
            y1 = next_y1;
            x2 = next_x2;
            y2 = next_y2;

            // Only calculated to be used for the reverse calculation
            vx1 += (fx / m1) * dt;
            vy1 += (fy / m1) * dt;
            vx2 -= (fx / m2) * dt;
            vy2 -= (fy / m2) * dt;
        }

        t += dt;
    }

    return results;
}

// Write an array of values to a file (used for plotting in Python)
void write_to_file(string filename, auto results) {
    ofstream outfile(filename);

    if (!outfile.is_open()) {
        cerr << "Error: Unable to open the output file." << endl;
        return;
    }
    outfile << "x1" << ',' << "y1" << ',' << "x2" << ',' << "y2" << "vx1" << ',' << "vy1" << ',' << "vx2" << ',' << "vy2" << endl;
    for (const auto &[x1, y1, x2, y2, vx1, vy1, vx2, vy2] : results) {
        outfile << x1 << ", " << y1 << ", " << x2 << ", " << y2 << ", " << vx1 << ", " << vy1 << ", " << vx2 << ", " << vy2 << endl;
    }
    cout << "File \"" << filename <<  "\" written." << endl;
    cout << results.size() << " lines written." << endl;

    outfile.close();
}

// Used to calculate the mean and standard deviation. I miss the ease of np.mean() and np.std() here..
// These functions are introduced to make a reliable statement about the computation time of the two integration methods
template <typename T>
double mean(const std::vector<T> &data) {
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

template <typename T>
double standard_deviation(const std::vector<T> &data) {
    double mean_value = mean(data);
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0,
                                       std::plus<>(),
                                       [mean_value](const T &x, const T &y) { return (x - mean_value) * (y - mean_value); });
    return std::sqrt(sq_sum / data.size());
}


int main() {
    double dt = 0.1;      // Time step
    double T = 100;       // Total simulation time

    // Run the simulations for both methods
    int N = 1000;
    vector<int64_t> elapsed_times_euler;        // int64_t is the underlying data type of the count() function
    vector<int64_t> elapsed_times_verlet;

    cout << "a)" << endl;
    // a)
    // Calculate the positions, coarse step size
    vector<tuple<double, double, double, double, double, double, double, double>> euler_positions  = simulate("Euler",  dt, T);
    vector<tuple<double, double, double, double, double, double, double, double>> verlet_positions = simulate("Verlet", dt, T);

    // Write the results
    write_to_file("Euler_dt00_1.csv", euler_positions);
    write_to_file("Verlet_dt00_1.csv", verlet_positions);

    // Adequate step size
    dt = 0.01;

    // Calculate the positions, coarse step size
    euler_positions  = simulate("Euler",  dt, T);
    verlet_positions = simulate("Verlet", dt, T);

    // Write the results
    write_to_file("Euler_dt00_01.csv", euler_positions);
    write_to_file("Verlet_dt00_01.csv", verlet_positions);


    cout << endl << "b)" << endl;
    // b)
    // We calculate the two methods serially to mitigate fluctuations in processing speed
    for (int i = 0; i < N; i++) {
        // Measure time for Euler method
        auto start_euler = chrono::steady_clock::now();
        vector<tuple<double, double, double, double, double, double, double, double>> euler_positions = simulate("Euler", dt, T);
        auto end_euler = chrono::steady_clock::now();
        auto elapsed_euler = chrono::duration_cast<chrono::microseconds>(end_euler - start_euler).count();
        elapsed_times_euler.push_back(elapsed_euler);

        // Measure time for Verlet method
        auto start_verlet = chrono::steady_clock::now();
        vector<tuple<double, double, double, double, double, double, double, double>> verlet_positions = simulate("Verlet", dt, T);
        auto end_verlet = chrono::steady_clock::now();
        auto elapsed_verlet = chrono::duration_cast<chrono::microseconds>(end_verlet - start_verlet).count();
        elapsed_times_verlet.push_back(elapsed_verlet);
    }

    double mean_euler = mean(elapsed_times_euler), mean_verlet = mean(elapsed_times_verlet);
    double std_euler = standard_deviation(elapsed_times_euler), std_verlet = standard_deviation(elapsed_times_verlet);

    cout << "Time needed for Euler in μs: " << int(mean_euler) << " +-" << int(std_euler) << endl;
    cout << "Time needed for Verlet in μs: " << int(mean_verlet) << " +-" << int(std_verlet) << endl;
    cout << "Iterations, N = " << N << endl;
    cout << "Euler's method is " << (1-mean_euler/mean_verlet)*100 << "\% faster." << endl;


    cout << endl << "c)" << endl;
    // For c), pass the latest calculated points as starting values and define a negative time step, as well as T=0
    // New initial conditions, reverse everything
    dt = -dt;
    double t = T;
    T = 0;

    // New positions and velocities
    double x1 = get<0>(euler_positions.back()), y1 = get<1>(euler_positions.back()), x2 = get<2>(euler_positions.back()), y2 = get<3>(euler_positions.back());
    double vx1 = get<4>(euler_positions.back()), vy1 = get<5>(euler_positions.back()), vx2 = get<6>(euler_positions.back()), vy2 = get<7>(euler_positions.back());

    vector<tuple<double, double, double, double, double, double, double, double>> euler_positions_reverse  = simulate("Euler",  dt, T, t, x1, y1, x2, y2, vx1, vy1, vx2, vy2);
    vector<tuple<double, double, double, double, double, double, double, double>> verlet_positions_reverse = simulate("Verlet",  dt, T, t, x1, y1, x2, y2, vx1, vy1, vx2, vy2);

    // Merge the vectors
    euler_positions.insert(euler_positions.end(), euler_positions_reverse.begin(), euler_positions_reverse.end());
    verlet_positions.insert(verlet_positions.end(), verlet_positions_reverse.begin(), verlet_positions_reverse.end());

    // Write the results
    write_to_file("Euler_reverse_dt00_01.csv", euler_positions);
    write_to_file("Verlet_reverse_dt00_01.csv", verlet_positions);

    return 0;
}