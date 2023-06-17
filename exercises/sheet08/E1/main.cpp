#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

typedef unsigned long long int ulli;

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// CLASSES
class LCG {
    public:
                                                LCG(ulli seed, ulli a, ulli c, ulli m) : seed(seed), a(a), c(c), m(m) {};
        float                                   random_f();          
    private:
        ulli                                    seed;
        ulli                                    a;
        ulli                                    c;
        ulli                                    m;

        ulli                                    iteration_step();

};

float LCG::random_f() {
    ulli result = iteration_step();
    return static_cast<float>(result) / static_cast<float>(m);
}

ulli LCG::iteration_step() {
    seed = (a * seed + c) % m;
    return seed;
}

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// CLASSES END

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// SAMPLINGS

double box_muller_transform(LCG& generator) {
    double U = generator.random_f();
    double V = generator.random_f();

    double X = sqrt(-2.0 * log(U)) * cos(2.0 * M_PI * V);
    // double Y = sqrt(-2.0 * log(U)) * sin(2.0 * M_PI * V); // For the second random number

    return X;
}

double central_limit_transform(LCG& generator, int N) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += generator.random_f();
    }

    // Normalize: Mean 0, Std 1
    return (sum - N/2.0) / sqrt(N/12.0);
}

double rejection_sampling(LCG& generator) {
    double pi = 3.14159265358979323846;
    while (true) {
        double u = generator.random_f() * pi;  // Scaled to [0, pi]
        double v = generator.random_f() * 0.5;  // Scaled to [0, 0.5]
        if (v <= std::sin(u) / 2) {
            return u;
        }
    }
}

double inversion_sampling(LCG& generator) {
    return std::cbrt(generator.random_f());  // cbrt() is the cubic root
}

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// WRITE TO CSV

void write_to_csv_box_muller(std::vector<LCG>& randomGenerators, int N) {
    std::ofstream file;
    file.open("box_muller_numbers.csv");
    file << "Gaussian\n";
    for (int i = 0; i < N; i++) {
        file << box_muller_transform(randomGenerators[3]) << "\n";  // Generator with index 3 was asked for
    }
    file.close();
}

void write_to_csv_central_limit(std::vector<LCG>& randomGenerators, int N, int n_samples) {
    std::ofstream file;
    file.open("central_limit_numbers.csv");
    file << "Gaussian\n";
    for (int i = 0; i < n_samples; i++) {
        file << central_limit_transform(randomGenerators[3], N) << "\n";
    }
    file.close();
}

void write_to_csv(std::vector<LCG>& randomGenerators, int N) {
    std::ofstream file;
    file.open("random_numbers.csv");
    file << "a),b),c),d)\n";
    for (int i = 0; i < N; i++) {
        for (LCG& generator : randomGenerators) {
            file << generator.random_f();
            if (&generator != &randomGenerators.back()) {
                file << ", ";
            }
        }
        file << "\n";
    }
    file.close();
}

void write_to_csv_rejection_sampling(std::vector<LCG>& randomGenerators, int n_samples) {
    std::ofstream file;
    file.open("rejection_sampling_numbers.csv");
    file << "Rejection\n";
    for (int i = 0; i < n_samples; i++) {
        file << rejection_sampling(randomGenerators[3]) << "\n";
    }
    file.close();
}

void write_to_csv_inversion_sampling(std::vector<LCG>& randomGenerators, int n_samples) {
    std::ofstream file;
    file.open("inversion_sampling_numbers.csv");
    file << "Inversion\n";
    for (int i = 0; i < n_samples; i++) {
        file << inversion_sampling(randomGenerators[3]) << "\n";
    }
    file.close();
}

// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




int main() {
    std::vector<ulli> seed = {1234, 1234, 123456789, 1234};
    std::vector<ulli> a = {20, 137, 65539, 16807};
    std::vector<ulli> c = {120, 187, 0, 0};
    std::vector<ulli> m = {6075, 256, (1ULL << 31), (1ULL << 31) - 1};

    std::vector<LCG> randomGenerators = {};
    for (int i = 0; i < seed.size(); i++) {
        randomGenerators.push_back(LCG(seed[i], a[i], c[i], m[i]));
    }
    
    std::cout << "a) 50 Numbers:\n";
    for (int i = 0; i < 50; i++) {
        std::cout << randomGenerators[0].random_f();
        if ((i + 1) % 10 != 0) std::cout << ", ";
        else std::cout << "\n";
    }
    
    write_to_csv(randomGenerators, 100000);
    write_to_csv_box_muller(randomGenerators, 100000);
    write_to_csv_central_limit(randomGenerators, 12, 100000);  // N=12 because it's close to sqrt 100
    write_to_csv_rejection_sampling(randomGenerators, 100000);
    write_to_csv_inversion_sampling(randomGenerators, 100000);

    return 0;
}