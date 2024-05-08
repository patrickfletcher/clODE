//
// Created by Wolf on 18/09/2022.
//

#define PY_SSIZE_T_CLEAN
#define CONFIG_64
#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "clODE_struct_defs.cl"
#include "CLODE.hpp"
#include "CLODEfeatures.hpp"
#include "CLODEtrajectory.hpp"

#include "logging/PythonSink.hpp"
#include "spdlog/spdlog.h"

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = py::detail::overload_cast_impl<Args...>;

std::string vector_to_string(const std::vector<std::string> &vec)
{
    std::string out = "[";
    for (const std::string &s : vec)
    {
        out += s + ", ";
    }
    out += "]";
    return out;
}

PYBIND11_MODULE(clode_cpp_wrapper, m)
{

    m.doc() = "CLODE C++/Python interface"; // optional module docstring

    // logging
    /****************************************************/

    py::enum_<spdlog::level::level_enum>(m, "LogLevel")
        .value("trace", spdlog::level::trace)
        .value("debug", spdlog::level::debug)
        .value("info", spdlog::level::info)
        .value("warn", spdlog::level::warn)
        .value("err", spdlog::level::err)
        .value("critical", spdlog::level::critical)
        .value("off", spdlog::level::off)
        .export_values();

    struct LoggerSingleton
    {
        std::shared_ptr<PythonSink_mt> sink;
        std::shared_ptr<spdlog::logger> python_logger;
        LoggerSingleton()
        {
            spdlog::set_level(spdlog::level::info);
            sink = std::make_shared<PythonSink_mt>();
            python_logger = std::make_shared<spdlog::logger>("python", sink);
            spdlog::set_default_logger(python_logger);
        }

        static LoggerSingleton &instance()
        {
            static LoggerSingleton just_one;
            return just_one;
        }

        void set_log_level(spdlog::level::level_enum level)
        {
            python_logger->set_level(level);
        };

        void set_log_pattern(std::string &pattern)
        {
            python_logger->set_pattern(pattern);
        };

        spdlog::level::level_enum get_log_level()
        {
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

    // OpenCL runtime
    /****************************************************/

    py::enum_<cl_vendor>(m, "CLVendor")
        .value("VENDOR_ANY", VENDOR_ANY)
        .value("VENDOR_NVIDIA", VENDOR_NVIDIA)
        .value("VENDOR_AMD", VENDOR_AMD)
        .value("VENDOR_INTEL", VENDOR_INTEL)
        .export_values();

    py::enum_<e_cl_device_type>(m, "CLDeviceType")
        .value("DEVICE_TYPE_ALL", DEVICE_TYPE_ALL)
        .value("DEVICE_TYPE_CPU", DEVICE_TYPE_CPU)
        .value("DEVICE_TYPE_GPU", DEVICE_TYPE_GPU)
        .value("DEVICE_TYPE_ACCELERATOR", DEVICE_TYPE_ACCELERATOR)
        .value("DEVICE_TYPE_DEFAULT", DEVICE_TYPE_DEFAULT)
        .value("DEVICE_TYPE_CUSTOM", DEVICE_TYPE_CUSTOM)
        .export_values();

    py::class_<OpenCLResource>(m, "OpenCLResource")
        .def(py::init<>())
        .def(py::init<cl_device_type>())
        .def(py::init<cl_vendor>())
        .def(py::init<cl_device_type, cl_vendor>())
        .def(py::init<e_cl_device_type, cl_vendor>())
        .def(py::init<unsigned int, unsigned int>())
        .def(py::init<unsigned int, std::vector<unsigned int>>())
        .def("get_double_support", &OpenCLResource::getDoubleSupport, "Get double support")
        .def("get_max_memory_alloc_size", &OpenCLResource::getMaxMemAllocSize, "Get max memory alloc size")
        .def("get_device_cl_version", &OpenCLResource::getDeviceCLVersion, "Get device CL version")
        .def("print_devices", &OpenCLResource::print, "Print device info to log");

    py::class_<deviceInfo>(m, "DeviceInfo")
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
        .def(
            "__repr__", [](const deviceInfo &d)
            { return "<device_info(name=" + d.name +
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
                     ")>"; },
            "Device info string representation");

    py::class_<platformInfo>(m, "PlatformInfo")
        .def_readwrite("name", &platformInfo::name)
        .def_readwrite("vendor", &platformInfo::vendor)
        .def_readwrite("version", &platformInfo::version)
        .def_readwrite("device_info", &platformInfo::device_info)
        .def_readwrite("device_count", &platformInfo::nDevices)
        .def(
            "__repr__", [](const platformInfo &p)
            { return "<platform_info(name=" + p.name +
                     ", vendor=" + p.vendor +
                     ", version=" + p.version +
                     ", device_count=" + std::to_string(p.nDevices) +
                     ")>"; },
            "Platform info string representation");

    m.def("query_opencl", &queryOpenCL, "Query OpenCL devices");
    m.def("_print_opencl", overload_cast_<>()(&printOpenCL), "Print OpenCL devices");

    // core clODE solver
    /****************************************************/

    py::class_<ProblemInfo>(m, "ProblemInfo")
        .def(py::init<const std::string &,
                      const std::vector<std::string> &,
                      const std::vector<std::string> &,
                      const std::vector<std::string> &,
                      int>(),
             py::arg("src_file"),
             py::arg("vars"),
             py::arg("pars"),
             py::arg("aux") = std::vector<std::string>(),
             py::arg("num_noise") = 1)
        .def(py::init<>())
        .def_readwrite("src_file", &ProblemInfo::clRHSfilename)
        .def_readwrite("num_var", &ProblemInfo::nVar)
        .def_readwrite("num_par", &ProblemInfo::nPar)
        .def_readwrite("num_aux", &ProblemInfo::nAux)
        .def_readwrite("num_noise", &ProblemInfo::nWiener)
        .def_property("vars", &ProblemInfo::getVarNames, &ProblemInfo::setVarNames)
        .def_property("pars", &ProblemInfo::getParNames, &ProblemInfo::setParNames)
        .def_property("aux", &ProblemInfo::getAuxNames, &ProblemInfo::setAuxNames)
        .def(
            "__repr__", [](const ProblemInfo &p)
            { return "<problem_info(src_file=" + p.clRHSfilename +
                     ", vars=" + vector_to_string(p.varNames) +
                     ", pars=" + vector_to_string(p.parNames) +
                     ", aux=" + vector_to_string(p.auxNames) +
                     ", num_noise=" + std::to_string(p.nWiener) +
                     ")>"; },
            "clode Problem info string representation");

    py::class_<SolverParams<double>>(m, "SolverParams")
        .def(py::init<double,
                      double,
                      double,
                      double,
                      unsigned int,
                      unsigned int,
                      unsigned int>(),
             py::arg("dt") = 0.1,
             py::arg("dtmax") = 0.5,
             py::arg("abstol") = 1e-6,
             py::arg("reltol") = 1e-3,
             py::arg("max_steps") = 1000000,
             py::arg("max_store") = 1000000,
             py::arg("nout") = 1)
        .def_readwrite("dt", &SolverParams<double>::dt)
        .def_readwrite("dtmax", &SolverParams<double>::dtmax)
        .def_readwrite("abstol", &SolverParams<double>::abstol)
        .def_readwrite("reltol", &SolverParams<double>::reltol)
        .def_readwrite("max_steps", &SolverParams<double>::max_steps)
        .def_readwrite("max_store", &SolverParams<double>::max_store)
        .def_readwrite("nout", &SolverParams<double>::nout)
        .def(
            "__repr__", [](const SolverParams<double> &s)
            { return "<solver_params(dt=" + std::to_string(s.dt) +
                     ", dtmax=" + std::to_string(s.dtmax) +
                     ", abstol=" + std::to_string(s.abstol) +
                     ", reltol=" + std::to_string(s.reltol) +
                     ", max_steps=" + std::to_string(s.max_steps) +
                     ", max_store=" + std::to_string(s.max_store) +
                     ", nout=" + std::to_string(s.nout) +
                     ")>"; },
            "Solver params string representation");

    py::class_<CLODE>(m, "SimulatorBase")
        .def(py::init<ProblemInfo &,
                      std::string &,
                      bool,
                      OpenCLResource &,
                      std::string &>(),
             py::arg("problem_info"),
             py::arg("stepper"),
             py::arg("cl_single_precision"),
             py::arg("opencl_resource"),
             py::arg("clode_root"))
        .def("set_problem_info", static_cast<void (CLODE::*)(ProblemInfo)>(&CLODE::setProblemInfo))
        .def("set_stepper", static_cast<void (CLODE::*)(std::string)>(&CLODE::setStepper))
        .def("set_precision", static_cast<void (CLODE::*)(bool)>(&CLODE::setPrecision))
        .def("set_opencl", static_cast<void (CLODE::*)(OpenCLResource opencl)>(&CLODE::setOpenCL))
        .def("set_opencl", static_cast<void (CLODE::*)(unsigned int platformID, unsigned int deviceID)>(&CLODE::setOpenCL))
        .def("build_cl", &CLODE::buildCL)
        .def("set_problem_data", static_cast<void (CLODE::*)(std::vector<double>, std::vector<double>)>(&CLODE::setProblemData))
        .def("set_tspan", static_cast<void (CLODE::*)(std::vector<double>)>(&CLODE::setTspan))
        .def("set_x0", static_cast<void (CLODE::*)(std::vector<double>)>(&CLODE::setX0))
        .def("set_pars", static_cast<void (CLODE::*)(std::vector<double>)>(&CLODE::setPars))
        .def("set_solver_params", static_cast<void (CLODE::*)(SolverParams<double>)>(&CLODE::setSolverParams))
        .def("seed_rng", static_cast<void (CLODE::*)()>(&CLODE::seedRNG), "Seed RNG")
        .def("seed_rng", static_cast<void (CLODE::*)(int)>(&CLODE::seedRNG), "Seed RNG", py::arg("seed"))
        .def("transient", &CLODE::transient)
        .def("shift_tspan", &CLODE::shiftTspan)
        .def("shift_x0", &CLODE::shiftX0)
        .def("get_problem_info", &CLODE::getProblemInfo)
        .def("get_tspan", &CLODE::getTspan)
        .def("get_solver_params", &CLODE::getSolverParams)
        .def("get_pars", &CLODE::getPars)
        .def("get_x0", &CLODE::getX0)
        .def("get_xf", &CLODE::getXf)
        .def("get_available_steppers", &CLODE::getAvailableSteppers)
        .def("get_program_string", &CLODE::getProgramString)
        .def("print_status", &CLODE::printStatus);

    // clODE features specialization
    /****************************************************/

    py::class_<ObserverParams<double>>(m, "ObserverParams")
        .def(py::init<unsigned int,
                      unsigned int,
                      unsigned int,
                      unsigned int,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double>(),
             py::arg("e_var_ix") = 0,
             py::arg("f_var_ix") = 0,
             py::arg("max_event_count") = 100,
             py::arg("max_event_timestamps") = 0,
             py::arg("min_amp") = 0.,
             py::arg("min_imi") = 0.,
             py::arg("nhood_radius") = 0.05,
             py::arg("x_up_threshold") = 0.2,
             py::arg("x_down_threshold") = 0.2,
             py::arg("dx_up_threshold") = 0.,
             py::arg("dx_down_threshold") = 0.,
             py::arg("eps_dx") = 0.)
        .def_readwrite("e_var_ix", &ObserverParams<double>::eVarIx)
        .def_readwrite("f_var_ix", &ObserverParams<double>::fVarIx)
        .def_readwrite("max_event_count", &ObserverParams<double>::maxEventCount)
        .def_readwrite("max_event_timestamps", &ObserverParams<double>::maxEventTimestamps)
        .def_readwrite("min_amp", &ObserverParams<double>::minXamp)
        .def_readwrite("min_imi", &ObserverParams<double>::minIMI)
        .def_readwrite("nhood_radius", &ObserverParams<double>::nHoodRadius)
        .def_readwrite("x_up_threshold", &ObserverParams<double>::xUpThresh)
        .def_readwrite("x_down_threshold", &ObserverParams<double>::xDownThresh)
        .def_readwrite("dx_up_threshold", &ObserverParams<double>::dxUpThresh)
        .def_readwrite("dx_down_threshold", &ObserverParams<double>::dxDownThresh)
        .def_readwrite("eps_dx", &ObserverParams<double>::eps_dx)
        .def("__repr__", [](const ObserverParams<double> &p)
             { return "<observer_params(e_var_ix=" + std::to_string(p.eVarIx) +
                      ", f_var_ix=" + std::to_string(p.fVarIx) +
                      ", max_event_count=" + std::to_string(p.maxEventCount) +
                      ", max_event_timestamps=" + std::to_string(p.maxEventTimestamps) +
                      ", min_amp=" + std::to_string(p.minXamp) +
                      ", min_imi=" + std::to_string(p.minIMI) +
                      ", nhood_radius=" + std::to_string(p.nHoodRadius) +
                      ", x_up_threshold=" + std::to_string(p.xUpThresh) +
                      ", x_down_threshold=" + std::to_string(p.xDownThresh) +
                      ", dx_up_threshold=" + std::to_string(p.dxUpThresh) +
                      ", dx_down_threshold=" + std::to_string(p.dxDownThresh) +
                      ", eps_dx=" + std::to_string(p.eps_dx) +
                      ")>"; });

    py::class_<CLODEfeatures, CLODE>(m, "FeatureSimulatorBase")
        .def(py::init<ProblemInfo &,
                      std::string &,
                      std::string &,
                      ObserverParams<double>,
                      bool,
                      OpenCLResource &,
                      std::string &>())
        .def("build_cl", &CLODEfeatures::buildCL)
        .def("set_observer_params", static_cast<void (CLODEfeatures::*)(ObserverParams<double>)>(&CLODEfeatures::setObserverParams))
        .def("set_observer", static_cast<void (CLODEfeatures::*)(std::string)>(&CLODEfeatures::setObserver))
        .def("initialize_observer", &CLODEfeatures::initializeObserver)
        .def("is_observer_initialized", &CLODEfeatures::isObserverInitialized)
        .def("features", static_cast<void (CLODEfeatures::*)(bool)>(&CLODEfeatures::features))
        .def("features", static_cast<void (CLODEfeatures::*)()>(&CLODEfeatures::features))
        .def("get_observer_params", &CLODEfeatures::getObserverParams)
        .def("get_observer_name", &CLODEfeatures::getObserverName)
        .def("get_f", &CLODEfeatures::getF)
        .def("get_n_features", &CLODEfeatures::getNFeatures)
        .def("get_feature_names", &CLODEfeatures::getFeatureNames)
        .def("get_available_observers", &CLODEfeatures::getAvailableObservers)
        .def(
            "__repr__", [](const CLODEfeatures &c)
            { return "<CLODEfeatures (observer=" + c.getObserverName() + ", n_features=" + std::to_string(c.getNFeatures()) + ")>"; },
            "CLODEfeatures string representation");

    // clODE trajectory specialization
    /****************************************************/

    py::class_<CLODEtrajectory, CLODE>(m, "TrajectorySimulatorBase")
        .def(py::init<ProblemInfo &,
                      std::string &,
                      bool,
                      OpenCLResource &,
                      std::string &>())
        .def("build_cl", &CLODEtrajectory::buildCL)
        .def("trajectory", &CLODEtrajectory::trajectory) // CLODEtrajectory specializations
        .def("get_t", &CLODEtrajectory::getT)
        .def("get_x", &CLODEtrajectory::getX)
        .def("get_dx", &CLODEtrajectory::getDx)
        .def("get_aux", &CLODEtrajectory::getAux)
        .def("get_n_stored", &CLODEtrajectory::getNstored);
}
