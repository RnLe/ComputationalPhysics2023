#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

// Analytical function declarations for f1 and f2
float f_1(float x);
float df_1(float x);
float ddf_1(float x);
float f_2(float x);
float df_2(float x);

// Function declarations for numerical differentiation methods
float two_point(float (*func)(float), float x, float h);
float four_point(float (*func)(float), float x, float h);

// Function declarations for writing data to files
void write_to_file(float* values, int steps, string filename);
void write_to_file_err(float* values, int steps, string filename);

// Function wrappers for numerical differentiation methods
// The wrapper is necessary, because two_point cannot call itself
float two_point_wrapper_x(float x);
float two_point_wrapper_h(float h);
float four_point_wrapper_x(float x);

int main(int argc, char* argv[]) {

    // Flag to write files into txts
    bool writeFiles = false;
    if (argc > 1 && (string(argv[1]) == "--write" || string(argv[1]) == "-w")) writeFiles = true;

    // Defining constants and parameters
    float x_a = M_PIf/4;
    float h = 0.8f;
    int steps_a = 100000;
    int steps_b = 100000;
    int steps_c = 100000;
    int steps_d = 100000;

    // Allocate memory for arrays to store data and errors
    float *a = new float[steps_a];
    float *a_error = new float[steps_a];
    float *b = new float[steps_b];
    float *b_error = new float[steps_b];
    float *c_error = new float[steps_c];
    float *d_error_f1 = new float[steps_d];
    float *d_error_f2 = new float[steps_d];

    // a)
    // Two-point method for f_1, determining best interval h
    // Comparison to analytical solution
    for (int i = 0; i < steps_a; i++) a[i] = two_point(f_1, x_a, (float)i/steps_a);
    for (int i = 0; i < steps_a; i++) a_error[i] = df_1(((float)i*2*M_PIf/steps_a - M_PIf)) - two_point(f_1, ((float)i*2*M_PIf/steps_a - M_PIf), h);

    // b)
    // Two-point method for the first derivative of f_1, determining best interval h
    // Comparison to analytical solution
    for (int i = 0; i < steps_b; i++) b[i] = two_point(two_point_wrapper_h, x_a, (float)i/steps_b);
    for (int i = 0; i < steps_b; i++) b_error[i] = ddf_1(((float)i*2*M_PIf/steps_b - M_PIf)) - two_point(two_point_wrapper_x, ((float)i*2*M_PIf/steps_b - M_PIf), h);

    // c) 
    // Four-point method for f_1
    // Comparison to analytical solution
    for (int i = 0; i < steps_c; i++) c_error[i] = df_1(((float)i*2*M_PIf/steps_c - M_PIf)) - four_point(f_1, ((float)i*2*M_PIf/steps_c - M_PIf), h);

    // d) Four-point method for f_2 and comparison with two-point method
    // Comparison to analytical solution
    for (int i = 0; i < steps_d; i++) d_error_f1[i] = df_2(((float)i*2*M_PIf/steps_d - M_PIf)) - four_point(f_2, ((float)i*2*M_PIf/steps_d - M_PIf), h);
    for (int i = 0; i < steps_d; i++) d_error_f2[i] = df_2(((float)i*2*M_PIf/steps_d - M_PIf)) - two_point(f_2, ((float)i*2*M_PIf/steps_d - M_PIf), h);

    // Write files if flag is set
    // Plots are created in python, matplotlib
    if (writeFiles) {
        write_to_file(a, steps_a, "a.txt");
        write_to_file_err(a_error, steps_a, "a_error.txt");
        write_to_file(b, steps_b, "b.txt");
        write_to_file_err(b_error, steps_b, "b_error.txt");
        write_to_file_err(c_error, steps_c, "c_error.txt");
        write_to_file_err(d_error_f1, steps_d, "d_error2_f2.txt");
        write_to_file_err(d_error_f2, steps_d, "d_error4_f2.txt");
    }

    // Free memory
    delete[] a;
    delete[] a_error;
    delete[] b;
    delete[] b_error;
    delete[] c_error;
    delete[] d_error_f1;
    delete[] d_error_f2;

    return 0;
}