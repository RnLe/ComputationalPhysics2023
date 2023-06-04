#include "a.hpp"
#include <cmath> // Used for basic math

// This function calculates an array of floats with equidistant (input)values in a given range from 0 to 'limit'
float* a_func(float limit, int steps) {

    float* array = new float[steps];        // Dynamically allocate an array of floats

    float interval = limit/steps;           // Interval from 0 to limit in #steps

    for (int i = 0; i < steps; i++)         // Math Magic
    {
        float x = (i+1)*interval;
        float val = 1 / sqrt(x) - 1 / sqrt(x + 1);
        array[i] = val;
    }
    return array;                           // Note that only the pointer is passed
}

// This is the numerically stable function, analogously implemented
float* a_func_stable(float limit, int steps) {

    float* array = new float[steps];

    float interval = limit/steps;

    for (int i = 0; i < steps; i++)
    {
        float x = (i+1)*interval;
        float val = 1 / (sqrt(x) * sqrt(x + 1) * (sqrt(x) + sqrt(x + 1)));       // Only this line is different
        array[i] = val;
    }
    return array;
}