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

template <typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

PYBIND11_MODULE(clode_cpp_wrapper, m) {

    m.doc() = "CLODE C++/Python interface"; // optional module docstring

    py::enum_<spdlog::level::level_enum>(m, "log_level")
        .value("trace", spdlog::level::trace)
        .value("debug", spdlog::level::debug)
        .value("info", spdlog::level::info)
        .value("warn", spdlog::level::warn)
        .value("err", spdlog::level::err)
        .value("critical", spdlog::level::critical)
        .value("off", spdlog::level::off)
        .export_values();

    struct LoggerSingleton {
        std::shared_ptr<PythonSink_mt> sink;
        std::shared_ptr<spdlog::logger> python_logger;
        LoggerSingleton() {
            spdlog::set_level(spdlog::level::info);
            sink = std::make_shared<PythonSink_mt>();
            python_logger = std::make_shared<spdlog::logger>("python", sink);
            spdlog::set_default_logger(python_logger);
        }

        static LoggerSingleton& instance()
        {
            static LoggerSingleton just_one;
            return just_one;
        }

        void set_log_level(spdlog::level::level_enum level){
            python_logger->set_level(level);
        };

        void set_log_pattern(std::string &pattern){
            python_logger->set_pattern(pattern);
        };

        spdlog::level::level_enum get_log_level() {
            return python_logger->level();
        };
    };

    py::class_<LoggerSingleton>(m, "LoggerSingleton")
        .def("set_log_level", &LoggerSingleton::set_log_level)
        .def("set_log_pattern", &LoggerSingleton::set_log_pattern)
        .def("get_log_level", &LoggerSingleton::get_log_level);

    m.def("get_logger",
        &LoggerSingleton::instance,
        py::return_value_policy::reference,
        "Get logger singleton instance");

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
    .def_readwrite("maxEventCount", &ObserverParams<double>::maxEventCount)
    .def("__repr__", [](const ObserverParams<double> &p) {
        return "<observer_params(e_var_ix=" + std::to_string(p.eVarIx) +
               ", f_var_ix=" + std::to_string(p.fVarIx) +
               ", maxEventCount=" + std::to_string(p.maxEventCount) +
               ", minXamp=" + std::to_string(p.minXamp) +
               ", minIMI=" + std::to_string(p.minIMI) +
               ", nHoodRadius=" + std::to_string(p.nHoodRadius) +
               ", xUpThresh=" + std::to_string(p.xUpThresh) +
               ", xDownThresh=" + std::to_string(p.xDownThresh) +
               ", dxUpThresh=" + std::to_string(p.dxUpThresh) +
               ", dxDownThresh=" + std::to_string(p.dxDownThresh) +
               ", eps_dx=" + std::to_string(p.eps_dx) +
               ")>";
    });

    py::enum_<cl_vendor>(m, "cl_vendor")
        .value("VENDOR_ANY", VENDOR_ANY)
        .value("VENDOR_NVIDIA", VENDOR_NVIDIA)
        .value("VENDOR_AMD", VENDOR_AMD)
        .value("VENDOR_INTEL", VENDOR_INTEL)
        .export_values();

    enum cl_device_type_wrapper {
        DEVICE_TYPE_ALL = CL_DEVICE_TYPE_ALL,
        DEVICE_TYPE_CPU = CL_DEVICE_TYPE_CPU,
        DEVICE_TYPE_GPU = CL_DEVICE_TYPE_GPU,
        DEVICE_TYPE_ACCELERATOR = CL_DEVICE_TYPE_ACCELERATOR,
        DEVICE_TYPE_DEFAULT = CL_DEVICE_TYPE_DEFAULT,
        DEVICE_TYPE_CUSTOM = CL_DEVICE_TYPE_CUSTOM
    };

    py::enum_<cl_device_type_wrapper>(m, "cl_device_type")
        .value("DEVICE_TYPE_ALL", DEVICE_TYPE_ALL)
        .value("DEVICE_TYPE_CPU", DEVICE_TYPE_CPU)
        .value("DEVICE_TYPE_GPU", DEVICE_TYPE_GPU)
        .value("DEVICE_TYPE_ACCELERATOR", DEVICE_TYPE_ACCELERATOR)
        .value("DEVICE_TYPE_DEFAULT", DEVICE_TYPE_DEFAULT)
        .value("DEVICE_TYPE_CUSTOM", DEVICE_TYPE_CUSTOM)
        .export_values();

    py::class_<OpenCLResource>(m, "opencl_resource")
        .def(py::init<>())
        .def(py::init<cl_device_type>())
        .def(py::init<cl_vendor>())
        .def(py::init<cl_device_type, cl_vendor>())
        .def(py::init<unsigned int, unsigned int>())
        .def(py::init<unsigned int, std::vector<unsigned int>>())
        .def("get_double_support", &OpenCLResource::getDoubleSupport, "Get double support")
        .def("get_max_memory_alloc_size", &OpenCLResource::getMaxMemAllocSize, "Get max memory alloc size")
        .def("get_device_cl_version", &OpenCLResource::getDeviceCLVersion, "Get device CL version")
        .def("print_devices", &OpenCLResource::print, "Print device info to log");

    py::class_<deviceInfo>(m, "device_info")
        .def_readwrite("name", &deviceInfo::name)
        .def_readwrite("vendor", &deviceInfo::vendor)
        .def_readwrite("version", &deviceInfo::version)
        .def_readwrite("device_type", &deviceInfo::devType)
        .def_readwrite("device_type_str", &deviceInfo::devTypeStr)
        .def_readwrite("compute_units", &deviceInfo::computeUnits)
        .def_readwrite("max_clock", &deviceInfo::maxClock)
        .def_readwrite("max_work_group_size", &deviceInfo::maxWorkGroupSize)
        .def_readwrite("device_memory_size", &deviceInfo::deviceMemSize)
        .def_readwrite("max_memory_alloc_size", &deviceInfo::maxMemAllocSize)
        .def_readwrite("extensions", &deviceInfo::extensions)
        .def_readwrite("double_support", &deviceInfo::doubleSupport)
        .def_readwrite("device_available", &deviceInfo::deviceAvailable)
        .def("__repr__", [](const deviceInfo &d) {
            return "<device_info(name=" + d.name +
                   ", vendor=" + d.vendor +
                   ", version=" + d.version +
                   ", device_type=" + d.devTypeStr +
                   ", compute_units=" + std::to_string(d.computeUnits) +
                   ", max_clock=" + std::to_string(d.maxClock) +
                   ", max_work_group_size=" + std::to_string(d.maxWorkGroupSize) +
                   ", device_memory_size=" + std::to_string(d.deviceMemSize) +
                   ", max_memory_alloc_size=" + std::to_string(d.maxMemAllocSize) +
                   ", extensions=" + d.extensions +
                   ", double_support=" + std::to_string(d.doubleSupport) +
                   ", device_available=" + std::to_string(d.deviceAvailable) +
                   ")>";
        }, "Device info string representation");

    py::class_<platformInfo>(m, "platform_info")
        .def_readwrite("name", &platformInfo::name)
        .def_readwrite("vendor", &platformInfo::vendor)
        .def_readwrite("version", &platformInfo::version)
        .def_readwrite("device_info", &platformInfo::device_info)
        .def_readwrite("device_count", &platformInfo::nDevices)
        .def("__repr__", [](const platformInfo &p) {
            return "<platform_info(name=" + p.name +
                   ", vendor=" + p.vendor +
                   ", version=" + p.version +
                   ", device_count=" + std::to_string(p.nDevices) +
                   ")>";
        }, "Platform info string representation");

    m.def("query_opencl", &queryOpenCL, "Query OpenCL devices");
    m.def("print_opencl", overload_cast_<>()(&printOpenCL), "Print OpenCL devices");

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
            .def("set_tspan", static_cast<void (CLODEfeatures::*)(std::vector<double>)>(&CLODEfeatures::setTspan))
            .def("set_problem_data", static_cast<void (CLODEfeatures::*)(std::vector<double>, std::vector<double>)>(&CLODEfeatures::setProblemData))
            .def("set_x0", static_cast<void (CLODEfeatures::*)(std::vector<double>)>(&CLODEfeatures::setX0))
            .def("set_pars", static_cast<void (CLODEfeatures::*)(std::vector<double>)>(&CLODEfeatures::setPars))
            .def("set_solver_params", static_cast<void (CLODEfeatures::*)(SolverParams<double>)>(&CLODEfeatures::setSolverParams))
            .def("shift_tspan", &CLODEfeatures::shiftTspan)
            .def("shift_x0", &CLODEfeatures::shiftX0)
            .def("get_tspan", &CLODEfeatures::getTspan)
            .def("get_x0", &CLODEfeatures::getX0)
            .def("get_xf", &CLODEfeatures::getXf)
            .def("get_available_steppers", &CLODEfeatures::getAvailableSteppers)
            .def("get_problem_info", &CLODEfeatures::getProblemInfo)
            .def("get_program_string", &CLODEfeatures::getProgramString)
            .def("print_status", &CLODEfeatures::printStatus)                                      //end of CLODE methods
            .def("features", static_cast<void (CLODEfeatures::*)(bool)>(&CLODEfeatures::features)) //CLODEfeatures specializations
            .def("features", static_cast<void (CLODEfeatures::*)()>(&CLODEfeatures::features))
            .def("__repr__", [](const CLODEfeatures &c) {
                return "<CLODEfeatures (observer="
                + c.getObserverName()
                + ", n_features="
                + std::to_string(c.getNFeatures())
                + ")>";
            }, "CLODEfeatures string representation")
            .def("set_observer_params", static_cast<void (CLODEfeatures::*)(ObserverParams<double>)>(&CLODEfeatures::setObserverParams))
            .def("get_observer_name", &CLODEfeatures::getObserverName)
            .def("get_n_features", &CLODEfeatures::getNFeatures)
            .def("get_feature_names", &CLODEfeatures::getFeatureNames)
            .def("get_f", &CLODEfeatures::getF)
            .def("get_available_observers", &CLODEfeatures::getAvailableObservers);

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
            .def("set_tspan", static_cast<void (CLODEtrajectory::*)(std::vector<double>)>(&CLODEtrajectory::setTspan))
            .def("set_problem_data", static_cast<void (CLODEtrajectory::*)(std::vector<double>, std::vector<double>)>(&CLODEtrajectory::setProblemData))
            .def("set_x0", static_cast<void (CLODEtrajectory::*)(std::vector<double>)>(&CLODEtrajectory::setX0))
            .def("set_pars", static_cast<void (CLODEtrajectory::*)(std::vector<double>)>(&CLODEtrajectory::setPars))
            .def("set_solver_params", static_cast<void (CLODEtrajectory::*)(SolverParams<double>)>(&CLODEtrajectory::setSolverParams))
            .def("shift_tspan", &CLODEtrajectory::shiftTspan)
            .def("shift_x0", &CLODEtrajectory::shiftX0)
            .def("get_tspan", &CLODEtrajectory::getTspan)
            .def("get_x0", &CLODEtrajectory::getX0)
            .def("get_xf", &CLODEtrajectory::getXf)
            .def("get_available_steppers", &CLODEfeatures::getAvailableSteppers)
            .def("get_problem_info", &CLODEtrajectory::getProblemInfo)
            .def("get_program_string", &CLODEtrajectory::getProgramString)
            .def("print_status", &CLODEtrajectory::printStatus) //end of CLODE methods
            .def("trajectory", &CLODEtrajectory::trajectory)    //CLODEtrajectory specializations
            .def("get_t", &CLODEtrajectory::getT)
            .def("get_x", &CLODEtrajectory::getX)
            .def("get_dx", &CLODEtrajectory::getDx)
            .def("get_aux", &CLODEtrajectory::getAux)
            .def("get_n_stored", &CLODEtrajectory::getNstored);
}
