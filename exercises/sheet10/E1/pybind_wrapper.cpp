#include "metropolis.h"

#include <pybind11/pybind11.h>
// Automatic conversion headers
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// Metropolis Header (functions only, no class; metropolis.h)

PYBIND11_MODULE(Solver, m) {

    // Add functions
    // Metropolis algorithm
    m.def("metropolis", &metropolis, "Metropolis algorithm for the 2D Ising model");    // Actually not needed in Python, but for completeness
    m.def("metropolis_range", &metropolis_range, "Iterate over metropolis() for a given range of magnetic fields (come in as an array)");

}
