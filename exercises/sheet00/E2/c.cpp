#include "c.hpp"
#include <cmath>

// Introduce a random, fixed value, which represents the x in c)'s expression. We are only interested in a variation of delta.
float x_ = 1;

// This function calculates an array of floats with equidistant (input)values in a given range from 0 to 'limit'
float* c_func(float limit, int steps) {

    float* array = new float[steps];        // Dynamically allocate an array of floats

    float interval = limit/steps;           // Interval from 0 to limit in #steps

    for (int i = 0; i < steps; i++)         // Math Magic
    {
        float x = (i+1)*interval;
        float val = sinf32(x_ + x) - sinf32(x_);// Careful with the choice of the function. The standard sin(x) function returns a double.
        array[i] = val;
    }
    return array;                           // Note that only the pointer is passed
}

// This is the numerically stable function, analogously implemented
float* c_func_stable(float limit, int steps) {

    float* array = new float[steps];

    float interval = limit/steps;

    for (int i = 0; i < steps; i++)
    {
        float x = (i+1)*interval;
        float val = x * cosf32(x_);       // Only this line is different
        array[i] = val;
    }
    return array;
}