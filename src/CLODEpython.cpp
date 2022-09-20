//
// Created by Wolf on 18/09/2022.
//

#define PY_SSIZE_T_CLEAN
#define CONFIG_64
#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "clODE_struct_defs.cl"
#include "CLODEfeatures.hpp"

#include "logging/PythonSink.hpp"
#include "spdlog/spdlog.h"

namespace py = pybind11;

PYBIND11_MODULE(clode, m) {

    m.doc() = "CLODE Python interface"; // optional module docstring

    auto python_sink = std::make_shared<PythonSink_mt>();
    auto python_logger = std::make_shared<spdlog::logger>("python", python_sink);
    spdlog::register_logger(python_logger);

    py::class_<ProblemInfo>(m, "problem_info")
            .def(py::init<const std::string &,
                 int,
                 int,
                 int,
                 int,
                 const std::vector<std::string> &,
                 const std::vector<std::string> &,
                 const std::vector<std::string> &>
                 ()
            );

    py::class_<SolverParams<double>>(m, "solver_params")
    .def(py::init<double,
                double,
                double,
                double,
                int,
                int,
                int>
                ()
    );

    py::class_<ObserverParams<double>>(m, "observer_params")
    .def(py::init<int,
            int,
            int,
            double,
            double,
            double,
            double,
            double,
            double,
            double,
            double>
            ()
    );

    py::class_<OpenCLResource>(m, "opencl_resource")
    .def(py::init<>());

    py::class_<CLODEfeatures>(m, "clode_features")
            .def(py::init<ProblemInfo &,
                          std::string &,
                          std::string &,
                          bool,
                          OpenCLResource &,
                          std::string &>())
            .def("initialize", static_cast<void (CLODEfeatures::*)
                                                    (std::vector<double>,
                                                    std::vector<double>,
                                                    std::vector<double>,
                                                    SolverParams<double>,
                                                    ObserverParams<double>)>
                                                    (&CLODEfeatures::initialize), "Initialize CLODEfeatures")
            .def("seed_rng", static_cast<void (CLODEfeatures::*)(int)>(&CLODEfeatures::seedRNG))
            .def("seed_rng", static_cast<void (CLODEfeatures::*)()>(&CLODEfeatures::seedRNG))
            .def("transient", &CLODEfeatures::transient)
            .def("features", static_cast<void (CLODEfeatures::*)()>(&CLODEfeatures::features))
            .def("get_tspan", &CLODEfeatures::getTspan)
            .def("get_f", &CLODEfeatures::getF)
            .def("get_n_features", &CLODEfeatures::getNFeatures);

}
