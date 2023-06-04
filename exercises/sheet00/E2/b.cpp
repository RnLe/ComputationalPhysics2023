#include "b.hpp"
#include <cmath>

// This function calculates an array of floats with equidistant (input)values in a given range from 0 to 'limit'
float* b_func(float limit, int steps) {

    float* array = new float[steps];        // Dynamically allocate an array of floats

    float interval = limit/steps;           // Interval from 0 to limit in #steps

    for (int i = 0; i < steps; i++)         // Math Magic
    {
        float x = (i+1)*interval;
        float val = (1-cosf32(x))/sinf32(x);// Careful with the choice of the function. The standard sin(x) function returns a double.
        array[i] = val;
    }
    return array;                           // Note that only the pointer is passed
}

// This is the numerically stable function, analogously implemented
float* b_func_stable(float limit, int steps) {

    float* array = new float[steps];

    float interval = limit/steps;

    for (int i = 0; i < steps; i++)
    {
        float x = (i+1)*interval;
        float val = pow(x,2)/(2*sinf32(x)) ;       // Only this line is different
        array[i] = val;
    }
    return array;
}