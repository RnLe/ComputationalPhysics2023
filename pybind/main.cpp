#include <vector>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class Solver {
  public:
  int parameter;
  std::vector<double> solution;

  Solver(int p) : parameter(p) {
    for (int i = 0; i < p; i++) {
      solution.push_back((double) i / (double) p);
    }
  }

  Solver(int p, int c) : parameter(p) {
    for (int i = 0; i < p; i++) {
      solution.push_back((double) i / (double) p + c);
    }
  }

  std::string repr() {
    std::stringstream stream;
    stream << "Solver Object:" << std::endl
           << "parameter = " << parameter;
    return stream.str();
  }
};

PYBIND11_MODULE(Solver, m) {
  py::class_<Solver>(m, "Solver")
    .def(py::init<int>())
    .def(py::init<int, int>())
    .def("__repr__", &::Solver::repr)
    .def_readwrite("parameter", &Solver::parameter)
    .def_readwrite("solution", &Solver::solution);
}
