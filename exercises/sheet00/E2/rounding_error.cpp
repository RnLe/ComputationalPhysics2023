#include <iostream>
#include <fstream>  // Only used to write a file to plot the data in python externally
#include <string>
// .hpps are used to organize the code
#include "a.hpp"
#include "b.hpp"
#include "c.hpp"

using namespace std;

// Only used to plot the data in python
void write_to_file(float* values, int steps, string filename, float limit) {
    // Open a file for writing
    std::ofstream outfile(filename);

    // Check if the file is opened successfully
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open the output file." << std::endl;
        return;
    }

    // Write the array values to the file
    for (int i = 0; i < steps; i++) {
        outfile << (i+1)*(limit/steps) << ", " << values[i] << std::endl;
    }

    // Close the file
    outfile.close();

    return;
}

int main() {
    // Flags for easier code handling
    bool quickView = true;      // Only the last value of an array is shown, instead of the whole array.
    bool writeToFile = true;

    int steps = 1000;
    float limit_a = 2e7;            // 2e7 is a good value, because cancellation is happening right before reaching the limit
    float limit_a_stable = 2e25;    // ~1.175e38 is the upper limit for floats
    float* a = a_func(limit_a, steps);  
    float* a_stable = a_func_stable(limit_a_stable, steps);   


    if (quickView)
    {
        cout << endl << "a)" << endl << "Unstable calculation: " << endl;

        cout << a[steps-1];

        cout << endl << "Stable calculation: " << endl;

        cout << a_stable[steps-1];

        cout << endl;
    } else
    {
        cout << endl << "a)" << endl << "Unstable calculation: " << endl;

        for (int i = 0; i < steps; i++)
        {
            cout << a[i] << ", ";
        }

        cout << endl << "Stable calculation: " << endl;

        for (int i = 0; i < steps; i++)
        {
            cout << a_stable[i] << ", ";
        }

        cout << endl;
    }
    
    
    if (writeToFile)
    {
        write_to_file(a, steps, "a.txt", limit_a);
        write_to_file(a_stable, steps, "a_stable.txt", limit_a_stable);
    }
    
    delete[] a;
    delete[] a_stable;

    // b)
    steps = 100;
    float limit_b = 2e-3;
    float limit_b_stable = 2e-36; 
    float* b = b_func(limit_b, steps);  
    float* b_stable = b_func_stable(limit_b_stable, steps);  

    if (quickView)
    {
        cout << endl << "b)" << endl << "Unstable calculation: " << endl;

        cout << b[0];

        cout << endl << "Stable calculation: " << endl;

        cout << b_stable[0];

        cout << endl;
    } else
    {
        cout << endl << "b)" << endl << "Unstable calculation: " << endl;

        for (int i = 0; i < steps; i++)
        {
            cout << b[i] << ", ";
        }

        cout << endl << "Stable calculation: " << endl;

        for (int i = 0; i < steps; i++)
        {
            cout << b_stable[i] << ", ";
        }

        cout << endl;
    }
    
    
    if (writeToFile)
    {
        write_to_file(b, steps, "b.txt", limit_b);
        write_to_file(b_stable, steps, "b_stable.txt", limit_b_stable);
    }

    delete[] b;
    delete[] b_stable;

    // c)
    steps = 100;
    float limit_c = 2e-7;
    float limit_c_stable = 2e-36; 
    float* c = c_func(limit_c, steps);  
    float* c_stable = c_func_stable(limit_c_stable, steps);  

    if (quickView)
    {
        cout << endl << "c)" << endl << "Unstable calculation: " << endl;

        cout << c[0];

        cout << endl << "Stable calculation: " << endl;

        cout << c_stable[0];

        cout << endl;
    } else
    {
        cout << endl << "c)" << endl << "Unstable calculation: " << endl;

        for (int i = 0; i < steps; i++)
        {
            cout << c[i] << ", ";
        }

        cout << endl << "Stable calculation: " << endl;

        for (int i = 0; i < steps; i++)
        {
            cout << c_stable[i] << ", ";
        }

        cout << endl;
    }
    
    
    if (writeToFile)
    {
        write_to_file(c, steps, "c.txt", limit_c);
        write_to_file(c_stable, steps, "c_stable.txt", limit_c_stable);
    }

    delete[] c;
    delete[] c_stable;

    return 0;
}