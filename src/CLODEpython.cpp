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
#include "CLODEtrajectory.hpp"

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
    ).def_readwrite("e_var_ix", &ObserverParams<double>::eVarIx)
    .def_readwrite("f_var_ix", &ObserverParams<double>::fVarIx)
    .def_readwrite("maxEventCount", &ObserverParams<double>::maxEventCount);

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
            .def("build_cl", &CLODEfeatures::buildCL)
            .def("transient", &CLODEfeatures::transient)
            .def("features", static_cast<void (CLODEfeatures::*)(bool)>(&CLODEfeatures::features))
            .def("features", static_cast<void (CLODEfeatures::*)()>(&CLODEfeatures::features))
            .def("get_tspan", &CLODEfeatures::getTspan)
            .def("get_f", &CLODEfeatures::getF)
            .def("get_n_features", &CLODEfeatures::getNFeatures)
            .def("get_feature_names", &CLODEfeatures::getFeatureNames)
            .def("getXf", &CLODEfeatures::getXf)
            .def("shift_tspan", &CLODEfeatures::shiftTspan)
            .def("shift_x0", &CLODEfeatures::shiftX0);

    py::class_<CLODEtrajectory>(m, "clode_trajectory")
            .def(py::init<ProblemInfo &,
                    std::string &,
                    bool,
                    OpenCLResource &,
                    std::string &>())
            .def("initialize", static_cast<void (CLODEtrajectory::*)
                    (std::vector<double>,
                     std::vector<double>,
                     std::vector<double>,
                     SolverParams<double>)>
            (&CLODEtrajectory::initialize), "Initialize CLODEtrajectory")
            .def("seed_rng", static_cast<void (CLODEtrajectory::*)(int)>(&CLODEtrajectory::seedRNG))
            .def("seed_rng", static_cast<void (CLODEtrajectory::*)()>(&CLODEtrajectory::seedRNG))
            .def("build_cl", &CLODEtrajectory::buildCL)
            .def("transient", &CLODEtrajectory::transient)
            .def("trajectory", &CLODEtrajectory::trajectory)
            .def("get_t", &CLODEtrajectory::getT)
            .def("get_x", &CLODEtrajectory::getX)
            .def("get_x0", &CLODEtrajectory::getX0)
            .def("get_dx", &CLODEtrajectory::getDx)
            .def("get_aux", &CLODEtrajectory::getAux)
            .def("get_n_stored", &CLODEtrajectory::getNstored)
            .def("shift_x0", &CLODEtrajectory::shiftX0);
}
