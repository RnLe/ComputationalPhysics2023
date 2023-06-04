#include <iostream>
#include <fstream>
#include <string>
#include "a.hpp"
#include "b.hpp"
#include "c.hpp"

using namespace std;
// Function to write values to a file for plotting (omitted for brevity)
// All limits are chosen by try and error to demonstrate the cut-off point
int main() {
    // Define variables and flags for different calculations
    int steps;
    float limit_a, limit_a_stable, limit_b, limit_b_stable, limit_c, limit_c_stable;
    float* a, *a_stable, *b, *b_stable, *c, *c_stable;

    // a) Calculate unstable and stable results for a)
    steps = 1000;
    limit_a = 2e7;
    limit_a_stable = 2e25;
    a = a_func(limit_a, steps);  
    a_stable = a_func_stable(limit_a_stable, steps);   

    // Print and/or write results for a)

    // Free memory for a)
    delete[] a;
    delete[] a_stable;

    // b) Calculate unstable and stable results for b)
    steps = 100;
    limit_b = 2e-3;
    limit_b_stable = 2e-36;
    b = b_func(limit_b, steps);  
    b_stable = b_func_stable(limit_b_stable, steps);  

    // Print and/or write results for b)

    // Free memory for b)
    delete[] b;
    delete[] b_stable;

    // c) Calculate unstable and stable results for c)
    steps = 100;
    limit_c = 2e-7;
    limit_c_stable = 2e-36;
    c = c_func(limit_c, steps);  
    c_stable = c_func_stable(limit_c_stable, steps);  

    // Print and/or write results for c)

    // Free memory for c)
    delete[] c;
    delete[] c_stable;

    return 0;
}