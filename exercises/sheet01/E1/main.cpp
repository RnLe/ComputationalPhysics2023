#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

// Functions f1 and f2
//////////////////////////////////////////

// Function f_1: sin(x)
float f_1(float x) {
    return sinf32(x);
}

// First derivative of f_1: cos(x)
float df_1(float x) {
    return cosf32(x);
}

// Second derivative of f_1: -sin(x)
float ddf_1(float x) {
    return -sinf32(x);
}

// Function f_2: Piecewise function based on x and pi
float f_2(float x) {
    if (x >= 0 && x != M_PIf)
    {
        return (2*floorf32(x/M_PIf) - cosf32(fmod(x, M_PIf)) + 1);
    }
    else if (x < 0 && x != M_PIf)
    {
        return (2*floorf32(x/M_PIf) + cosf32(fmod(x, M_PIf)) + 1);
    }

    return NAN;
}

// First derivative of f_2: Piecewise function based on x and pi
float df_2(float x) {
    if (x >= 0 && x != M_PIf)
    {
        return (sinf32(fmod(x, M_PIf)));
    }
    else if (x < 0 && x != M_PIf)
    {
        return (-sinf32(fmod(x, M_PIf)));
    }

    return NAN;
}

//////////////////////////////////////////
// Functions for numerical differentiation

// Two-point finite difference method
float two_point(float (*func)(float), float x, float h) {
    float result = (func(x+h)-func(x-h))/(2*h);
    return result;
}

// Four-point finite difference method
float four_point(float (*func)(float), float x, float h) {
    float result = (-func(x+2*h)+8*func(x+h)-8*func(x-h)-func(x-2*h))/(12*h);
    return result;
}

/////////////////////////////////////////////////////////////

// Write an array of values to a file (used for plotting in Python)
void write_to_file(float* values, int steps, string filename) {
    ofstream outfile(filename);

    if (!outfile.is_open()) {
        cerr << "Error: Unable to open the output file." << endl;
        return;
    }

    for (int i = 0; i < steps; i++) {
        outfile << 10.0f*i/steps << ", " << values[i] << endl;
    }

    outfile.close();
}

// Write error values to a file
void write_to_file_err(float* values, int steps, string filename) {
    ofstream outfile(filename);

    if (!outfile.is_open()) {
        cerr << "Error: Unable to open the output file." << endl;
        return;
    }

    for (int i = 0; i < steps; i++) {
        outfile << ((float)i*2*M_PIf/steps - M_PIf) << ", " << values[i] << endl;
    }

    outfile.close();
}

// Wrappers for numerical differentiation methods
//////////////////////////////////////////

// Two-point method wrapper for f_1 with a fixed h value
float two_point_wrapper_x(float x) {
    return two_point(f_1, x, 0.01f);
}

// Two-point method wrapper for f_1 with a fixed x value
float two_point_wrapper_h(float h) {
    return two_point(f_1, M_PIf/4, h);
}

// Four-point method wrapper for f_1 with a fixed h value
float four_point_wrapper_x(float x) {
    return four_point(f_1, x, 0.01f);
}

//////////////////////////////////////////////////////////////


int main(int argc, char* argv[]) {

    bool writeFiles = false;
    if (argc > 1 && (string(argv[1]) == "--write" || string(argv[1]) == "-w")) writeFiles = true;

    float x_a = M_PIf/4;
    float h = 0.01f;
    int steps_a = 10000;
    int steps_b = 10000;
    int steps_c = 10000;
    int steps_d = 10000;

    float *a = new float[steps_a];
    float *a_error = new float[steps_a];
    float *b = new float[steps_b];
    float *b_error = new float[steps_b];
    float *c_error = new float[steps_c];
    float *d_error_f1 = new float[steps_d];
    float *d_error_f2 = new float[steps_d];

    // a) Two-point method for f_1
    for (int i = 0; i < steps_a; i++) a[i] = two_point(f_1, x_a, (float)i/steps_a);
    for (int i = 0; i < steps_a; i++) a_error[i] = two_point(f_1, ((float)i*2*M_PIf/steps_a - M_PIf), h) / df_1(((float)i*2*M_PIf/steps_a - M_PIf));

    // b) Two-point method for the first derivative of f_1

    for (int i = 0; i < steps_b; i++) b[i] = two_point(two_point_wrapper_h, x_a, (float)i/steps_b);
    for (int i = 0; i < steps_b; i++) b_error[i] = two_point(two_point_wrapper_x, ((float)i*2*M_PIf/steps_b - M_PIf), h) / ddf_1(((float)i*2*M_PIf/steps_b - M_PIf));

    // c) Four-point method for f_1

    for (int i = 0; i < steps_c; i++) c_error[i] = four_point(f_1, ((float)i*2*M_PIf/steps_c - M_PIf), h) / df_1(((float)i*2*M_PIf/steps_c - M_PIf));

    // d) Four-point method for f_2 and comparison with two-point method

    for (int i = 0; i < steps_d; i++) d_error_f1[i] = four_point(f_2, ((float)i*2*M_PIf/steps_d - M_PIf), h) / df_2(((float)i*2*M_PIf/steps_d - M_PIf));
    for (int i = 0; i < steps_d; i++) d_error_f2[i] = two_point(f_2, ((float)i*2*M_PIf/steps_d - M_PIf), h) / df_2(((float)i*2*M_PIf/steps_d - M_PIf));

    if (writeFiles) {
        write_to_file(a, steps_a, "a.txt");
        write_to_file_err(a_error, steps_a, "a_error.txt");
        write_to_file(b, steps_b, "b.txt");
        write_to_file_err(b_error, steps_b, "b_error.txt");
        write_to_file_err(c_error, steps_c, "c_error.txt");
        write_to_file_err(d_error_f1, steps_d, "d_error2_f2.txt");
        write_to_file_err(d_error_f2, steps_d, "d_error4_f2.txt");
    }

    delete[] a;
    delete[] a_error;
    delete[] b;
    delete[] b_error;
    delete[] c_error;
    delete[] d_error_f1;
    delete[] d_error_f2;


    return 0;
}