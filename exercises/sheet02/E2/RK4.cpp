#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <fstream>
#include <functional>

using namespace std;

// Define the system of ODEs as a function
tuple<double, double, double> shape_equations(double S, double R, double Z, double Psi, double p_La_gamma, double rho_ga2_gamma) {
    double dR_dS = cos(Psi);
    double dZ_dS = sin(Psi);
    double dPsi_dS;
    if (R == 0) {
        dPsi_dS = (p_La_gamma - rho_ga2_gamma * Z) / 2;
    } else {
        dPsi_dS = p_La_gamma - rho_ga2_gamma * Z - sin(Psi) / R;
    }
    
    return {dR_dS, dZ_dS, dPsi_dS};
}

// Implement the RK4 algorithm
vector<tuple<double, double, double, double>> rk4(const function<tuple<double, double, double>(double, double, double, double, double, double)> &f,
                                                  double S0, double R0, double Z0, double Psi0,
                                                  double p_La_gamma, double rho_ga2_gamma,
                                                  double h, double S_max) {

    vector<tuple<double, double, double, double>> results;
    double S = S0, R = R0, Z = Z0, Psi = Psi0;
    
    results.push_back({S, R, Z, Psi});
    
    while (S < S_max && R < 0.5) {
        // RK4 algorithm
        // Note that the variable S is only passed to f to resemble the RK4 algorithm from the literature
        auto [k1_R, k1_Z, k1_Psi] = f(S, R, Z, Psi, p_La_gamma, rho_ga2_gamma);
        auto [k2_R, k2_Z, k2_Psi] = f(S + h/2, R + h*k1_R/2, Z + h*k1_Z/2, Psi + h*k1_Psi/2, p_La_gamma, rho_ga2_gamma);
        auto [k3_R, k3_Z, k3_Psi] = f(S + h/2, R + h*k2_R/2, Z + h*k2_Z/2, Psi + h*k2_Psi/2, p_La_gamma, rho_ga2_gamma);
        auto [k4_R, k4_Z, k4_Psi] = f(S + h, R + h*k3_R, Z + h*k3_Z, Psi + h*k3_Psi, p_La_gamma, rho_ga2_gamma);

        double dR = (k1_R + 2*k2_R + 2*k3_R + k4_R) * h / 6;
        double dZ = (k1_Z + 2*k2_Z + 2*k3_Z + k4_Z) * h / 6;
        double dPsi = (k1_Psi + 2*k2_Psi + 2*k3_Psi + k4_Psi) * h / 6;

        R += h*dR;
        Z += h*dZ;
        Psi += h*dPsi;
        S += h;
        
        results.push_back({S, R, Z, Psi});
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
    outfile << 'S' << ',' << 'R' << ',' << 'Z' << ',' << "Psi" << endl;
    for (const auto &[S, R, Z, Psi] : results) {
        outfile << S << ',' << R << ',' << Z << ',' << Psi << endl;
    }

    outfile.close();
}

// Write an array of values to a file (used for plotting in Python)
void write_to_file_d(string filename, auto results) {
    ofstream outfile(filename);

    if (!outfile.is_open()) {
        cerr << "Error: Unable to open the output file." << endl;
        return;
    }
    outfile << "p_La_gamma" << ',' << "volumes" << endl;
    for (const auto &[p_La_gamma, volumes] : results) {
        outfile << p_La_gamma << ',' << volumes << endl;
    }

    outfile.close();
}

double volume(vector<tuple<double, double, double, double>> SRZPsi, double h) {
    double volume_ = 0;
    for (int i = 0; 0.5 > get<1>(SRZPsi[i]); i++) {
        double r1 = get<1>(SRZPsi[i]);
        double r2 = get<1>(SRZPsi[i+1]);
        double z1 = get<2>(SRZPsi[i]);
        double z2 = get<2>(SRZPsi[i+1]);

        // Calculate the volume of the infinitesimal disk formed by rotating each segment of the curve around the Z-axis
        double dV = M_PI * (r1 * r1 + r1 * r2 + r2 * r2) * (z2 - z1) / 3;
        volume_ += dV;
    }
    return volume_;
}


int main() {
    // Define parameters
    double p_La_gamma = 2, rho_ga2_gamma = 0.1;
    double S0 = 0, R0 = 0, Z0 = 0, Psi0 = 0;
    double h =1e-2;
    double S_max = 10000;

    cout << "Parameters: p_La_gamma: " << p_La_gamma << ", rho_ga2_gamma: " << rho_ga2_gamma << ", h: " << h << ", S_max: " << S_max << endl;

    // c) Run the RK4 algorithm
    cout << "c) Running RK4 algorithm..." << endl;
    auto results = rk4(shape_equations, S0, R0, Z0, Psi0, p_La_gamma, rho_ga2_gamma, h, S_max);
    cout << "RK4 algorithm finished.\n";

    // Print the first two and last two lines of data
    cout << "S\tR\tZ\tPsi" << endl;
    for (size_t i = 0; i < results.size(); i++) {
        if (i < 2 || i >= results.size() - 2) {
            auto &[S, R, Z, Psi] = results[i];
            cout << S << '\t' << R << '\t' << Z << '\t' << Psi << endl;
        }
    }

    string filename_c = "rk4_gamma2.txt";
    write_to_file(filename_c, results);
    cout << "Data written to file: " << filename_c << endl;

    // d)
    // Redefine parameters and explore the parameter space
    p_La_gamma = 0, rho_ga2_gamma = 0.5;
    double p_La_gamma_max = 5, N = 200;     // Introducing N as number of volumes calculated
    vector<tuple<double, double>> volumes;  // The volumes vector should also know at what parameter it was evaluated

    cout << "\nd) Exploring the parameter space..." << endl;
    cout << "N = " << N << " iterations of the RK4 algorithm." << endl;

    for (int i = 0; i < N; i++) {
        auto results = rk4(shape_equations, S0, R0, Z0, Psi0, p_La_gamma, rho_ga2_gamma, h, S_max);
        volumes.push_back({p_La_gamma, volume(results, h)});
        p_La_gamma = i / N * p_La_gamma_max;
    }

    string filename_d = "volumes_N200.txt";
    write_to_file_d(filename_d, volumes);
    cout << "Data written to file: " << filename_d << endl;

    return 0;
}