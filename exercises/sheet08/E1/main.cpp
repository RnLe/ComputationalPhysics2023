#include <iostream>
#include <vector>
#include <fstream>

typedef unsigned long long int ulli;

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

    return 0;
}