/*
 * OpenCLResource.cpp
 *
 *  Copyright 2017 Patrick Fletcher <patrick.fletcher@nih.gov>
 *
 */

#include "OpenCLResource.hpp"

#include <stdexcept>
#include <fstream>
#include <filesystem>

#include "spdlog/spdlog.h"

/********************************
 * OpenCLResource Member Functions
 ********************************/

OpenCLResource::OpenCLResource()
{
    getPlatformAndDevices(CL_DEVICE_TYPE_DEFAULT, VENDOR_ANY);
    initializeOpenCL();
}

OpenCLResource::OpenCLResource(cl_deviceType type)
{

    getPlatformAndDevices(type, VENDOR_ANY);
    initializeOpenCL();
}

OpenCLResource::OpenCLResource(cl_vendor vendor)
{

    getPlatformAndDevices(CL_DEVICE_TYPE_DEFAULT, vendor);
    initializeOpenCL();
}

OpenCLResource::OpenCLResource(cl_deviceType type, cl_vendor vendor)
{

    getPlatformAndDevices(type, vendor);
    initializeOpenCL();
}

OpenCLResource::OpenCLResource(e_cl_device_type type, cl_vendor vendor)
{

    getPlatformAndDevices(type, vendor);
    initializeOpenCL();
}

OpenCLResource::OpenCLResource(int argc, char **argv)
{

    // modified from openCLUtilities to include accelerators as a devicetype
    cl_deviceType type = CL_DEVICE_TYPE_ALL;
    cl_vendor vendor = VENDOR_ANY;
    int nValidArgs = 0;

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--device") == 0)
        {
            if (strcmp(argv[i + 1], "cpu") == 0)
            {
                type = CL_DEVICE_TYPE_CPU;
            }
            else if (strcmp(argv[i + 1], "gpu") == 0)
            {
                type = CL_DEVICE_TYPE_GPU;
            }
            else if (strcmp(argv[i + 1], "accel") == 0)
            {
                type = CL_DEVICE_TYPE_ACCELERATOR;
            }
            else
                throw cl::Error(1, "Unkown device type used with --device");
            i++;
            nValidArgs++;
        }
        else if (strcmp(argv[i], "--vendor") == 0)
        {
            if (strcmp(argv[i + 1], "amd") == 0)
            {
                vendor = VENDOR_AMD;
            }
            else if (strcmp(argv[i + 1], "intel") == 0)
            {
                vendor = VENDOR_INTEL;
            }
            else if (strcmp(argv[i + 1], "nvidia") == 0)
            {
                vendor = VENDOR_NVIDIA;
            }
            else
                throw cl::Error(1, "Unkown vendor name used with --vendor");
            i++;
            nValidArgs++;
        }
    }

    if (nValidArgs == 0 && argc > 1)
    {
        spdlog::warn("OpenCLResource didn't recognize the command line arguments. Using default device. ");
    }

    getPlatformAndDevices(type, vendor);
    initializeOpenCL();
}

OpenCLResource::OpenCLResource(unsigned int platformID, unsigned int deviceID)
{

    std::vector<unsigned int> deviceIDs(1, deviceID);
    getPlatformAndDevices(platformID, deviceIDs);
    initializeOpenCL();
}

OpenCLResource::OpenCLResource(unsigned int platformID, std::vector<unsigned int> deviceIDs)
{

    getPlatformAndDevices(platformID, deviceIDs);
    initializeOpenCL();
};

// modified from openCLUtilities
// TODO: check various possibilities on default? eg: any Accel, any non-intel GPU, intel GPU, any CPU
void OpenCLResource::getPlatformAndDevices(cl_deviceType type, cl_vendor vendor)
{

    // query all platform and device info, store as a vector of structs
    //~ std::vector<platformInfo> pinfo=queryOpenCL();

    // Get available platforms
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if (platforms.size() == 0)
        throw cl::Error(1, "No OpenCL platforms were found");

    int tempID = -1;
    if (vendor != VENDOR_ANY)
    {
        std::string vendorStr;
        switch (vendor)
        {
        case VENDOR_NVIDIA:
            vendorStr = "NVIDIA";
            break;
        case VENDOR_AMD:
            vendorStr = "Advanced Micro Devices";
            break;
        case VENDOR_INTEL:
            vendorStr = "Intel";
            break;
        default:
            throw cl::Error(1, "Invalid vendor specified");
            break;
        }

        std::vector<cl::Platform> tempPlatforms;
        for (unsigned int i = 0; i < platforms.size(); i++)
        {
            if (platforms[i].getInfo<CL_PLATFORM_VENDOR>().find(vendorStr) != std::string::npos)
            {
                tempPlatforms.push_back(platforms[i]);
            }
        }

        platforms = tempPlatforms; // keep only the platforms with correct vendor
    }

    std::vector<cl::Device> tempDevices;
    for (unsigned int i = 0; i < platforms.size(); i++)
    {
        try
        {
            platforms[i].getDevices(type, &tempDevices);
            // TODO: apply extra device filters (eg. extensionSupported,enoughMem,...) here?
            tempID = i;
            break;
        }
        catch (cl::Error &e)
        {
            continue;
        }
    }

    if (tempID == -1)
        throw cl::Error(1, "No compatible OpenCL platform found");

    // we found a platform with compatible device(s) to use
    platform = platforms[tempID];
    devices = tempDevices;

    // get info for the selected plaform and device(s)
    platform_info = getPlatformInfo(platform, devices);
}

void OpenCLResource::getPlatformAndDevices(unsigned int platformID, std::vector<unsigned int> deviceIDs)
{

    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    if (platforms.size() == 0)
        throw cl::Error(1, "No OpenCL platforms were found");

    spdlog::debug("Found {} OpenCL platforms", platforms.size());
    spdlog::debug("Using platform {}", platformID);

    if (platformID < platforms.size())
    {
        platform = platforms[platformID];
    }
    else
    {
        throw std::out_of_range("Specified platformID exceeds number of available platforms");
    }

    std::vector<cl::Device> tempDevices;
    platform.getDevices(CL_DEVICE_TYPE_ALL, &tempDevices);
    for (unsigned int i = 0; i < deviceIDs.size(); ++i)
    {

        if (deviceIDs[i] < tempDevices.size())
        {
            devices.push_back(tempDevices[deviceIDs[i]]);
        }
        else
        {
            throw std::out_of_range("Specified deviceID exceeds the number devices on the selected platform");
        }
    }

    // get info for the selected plaform and device(s)
    platform_info = getPlatformInfo(platform, devices);
}

// create the OpenCL context from the device list, and a queue for each device
void OpenCLResource::initializeOpenCL()
{
    try
    {
        context = cl::Context(devices);
        for (unsigned int i = 0; i < devices.size(); ++i)
            queues.push_back(cl::CommandQueue(context, devices[i]));
    }
    catch (cl::Error &er)
    {
        spdlog::error("{}({})", er.what(), CLErrorString(er.err()).c_str());
        throw er;
    }

    spdlog::debug("OpenCLResource Created");
}

// attempt to build OpenCL program given as a string, build options empty if not supplied.
void OpenCLResource::buildProgramFromString(std::string sourceStr, std::string buildOptions)
{
    spdlog::debug("Building program from string");
    spdlog::trace(sourceStr.c_str());
    spdlog::debug(buildOptions.c_str());
    cl::Program::Sources source(1, std::make_pair(sourceStr.c_str(), sourceStr.length()));
    std::string buildLog;
    cl_int builderror;
    try
    {
        program = cl::Program(context, source, &error);
        spdlog::debug("Program Object creation error code: {}", CLErrorString(error).c_str());

        builderror = program.build(devices, buildOptions.c_str());
        spdlog::debug("Program Object build error code: {}", CLErrorString(builderror).c_str());

        std::string kernelnames;
        program.getInfo(CL_PROGRAM_KERNEL_NAMES, &kernelnames);
        spdlog::debug("Kernels built:   {}", kernelnames.c_str());
    }
    catch (cl::Error &er)
    {
        spdlog::error("{}({})", er.what(), CLErrorString(er.err()).c_str());
        if (er.err() == CL_BUILD_PROGRAM_FAILURE)
        {
            // spdlog::info("{}",sourceStr.c_str());
            for (unsigned int i = 0; i < devices.size(); ++i)
            {
                program.getBuildInfo(devices[i], CL_PROGRAM_BUILD_LOG, &buildLog);
                spdlog::error("OpenCL build log, Device {}:", i);
                spdlog::error("{}", buildLog.c_str());
                // spdlog::error(std::__fs::filesystem::current_path().c_str());
            }
        }
        throw er;
    }
}

void OpenCLResource::buildProgramFromSource(std::string filename, std::string buildOptions)
{
    spdlog::debug("Building program from source file");
    spdlog::trace(filename.c_str());
    spdlog::debug(buildOptions.c_str());
    std::string sourceStr = read_file(filename);
    buildProgramFromString(sourceStr, buildOptions);
}

// void OpenCLResource::writeProgramBinary()
// {
//     write_file("cl_program.ptx", program.getInfo<CL_PROGRAM_BINARIES>()[0]); // save binary (ptx file)
// }

// prints the selected platform and devices information (queried on the fly)
void OpenCLResource::print()
{
    std::string tmp;
    spdlog::info("Selected platform and device: ");
    spdlog::info("Platform  --------------------");
    printPlatformInfo(platform_info);
}

/********************************
 * Other functions
 ********************************/

// get info of all devices on all platforms
std::vector<platformInfo> queryOpenCL()
{

    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    if (platforms.size() == 0)
        throw cl::Error(1, "No OpenCL platforms were found");

    std::vector<platformInfo> pinfo;

    for (unsigned int i = 0; i < platforms.size(); ++i)
    {
        pinfo.push_back(getPlatformInfo(platforms[i]));
    }

    return pinfo;
}

// get info of a given cl::Platform and its devices (optionally only devices of a given type, default is CL_DEVICE_TYPE_ALL)
platformInfo getPlatformInfo(cl::Platform platform, std::vector<cl::Device> devices)
{

    platformInfo pinfo;
    platform.getInfo(CL_PLATFORM_NAME, &pinfo.name);
    platform.getInfo(CL_PLATFORM_VENDOR, &pinfo.vendor);
    platform.getInfo(CL_PLATFORM_VERSION, &pinfo.version);

    if (devices.size() == 0)
    { // get all the devices
        platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
    }

    pinfo.nDevices = (unsigned int)devices.size();

    for (unsigned int j = 0; j < pinfo.nDevices; j++)
        pinfo.device_info.push_back(getDeviceInfo(devices[j]));

    return pinfo;
}

// get info for given cl::Device
deviceInfo getDeviceInfo(cl::Device device)
{

    deviceInfo dinfo;

    // Get device info for this compute resource
    device.getInfo(CL_DEVICE_NAME, &dinfo.name);
    device.getInfo(CL_DEVICE_VENDOR, &dinfo.vendor);
    device.getInfo(CL_DEVICE_VERSION, &dinfo.version);
    device.getInfo(CL_DEVICE_TYPE, &dinfo.devType);
    switch (dinfo.devType)
    {
    case CL_DEVICE_TYPE_CPU:
        dinfo.devTypeStr = "CPU";
        break;
    case CL_DEVICE_TYPE_GPU:
        dinfo.devTypeStr = "GPU";
        break;
    case CL_DEVICE_TYPE_ACCELERATOR:
        dinfo.devTypeStr = "Accelerator";
        break;
    default:
        dinfo.devTypeStr = "Unknown";
    }

    device.getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &dinfo.computeUnits);
    device.getInfo(CL_DEVICE_MAX_CLOCK_FREQUENCY, &dinfo.maxClock);
    device.getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &dinfo.maxWorkGroupSize);
    device.getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &dinfo.deviceMemSize);
    device.getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &dinfo.maxMemAllocSize);
    device.getInfo(CL_DEVICE_EXTENSIONS, &dinfo.extensions);

    std::string doubleStr = "fp64";
    dinfo.doubleSupport = dinfo.extensions.find(doubleStr) != std::string::npos;

    device.getInfo(CL_DEVICE_EXTENSIONS, &dinfo.extensions);
    device.getInfo(CL_DEVICE_AVAILABLE, &dinfo.deviceAvailable);

    return dinfo;
}

// TODO: kernel info

// /* cl_kernel_work_group_info */
// #define CL_KERNEL_WORK_GROUP_SIZE                   0x11B0
// #define CL_KERNEL_COMPILE_WORK_GROUP_SIZE           0x11B1
// #define CL_KERNEL_LOCAL_MEM_SIZE                    0x11B2
// #define CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE 0x11B3
// #define CL_KERNEL_PRIVATE_MEM_SIZE                  0x11B4
// #ifdef CL_VERSION_1_2
// #define CL_KERNEL_GLOBAL_WORK_SIZE                  0x11B5
// #endif

// print information about all platforms and devices found
void printOpenCL()
{

    spdlog::info("Querying OpenCL platforms...");
    std::vector<platformInfo> pinfo = queryOpenCL();
    printOpenCL(pinfo);
}

// print information about all platforms and devices found, given pre-queried array of platformInfo structs
void printOpenCL(std::vector<platformInfo> pinfo)
{

    spdlog::info("Number of platforms found: {}", (unsigned int)pinfo.size());
    for (unsigned int i = 0; i < pinfo.size(); ++i)
    {
        spdlog::info("- Platform {} ------------------------------", i);
        printPlatformInfo(pinfo[i]);
    }
    spdlog::info("");
}

// print information about a platform and its devices given pre-queried info in platformInfo struct
void printPlatformInfo(platformInfo pinfo)
{
    spdlog::info("Name:    {}", pinfo.name.c_str());
    spdlog::info("Vendor:  {}", pinfo.vendor.c_str());
    spdlog::info("Version: {}", pinfo.version.c_str());

    for (unsigned int j = 0; j < pinfo.nDevices; j++)
    {
        spdlog::info("- Device {} ------------", j);
        printDeviceInfo(pinfo.device_info[j]);
    }
}

// print info about a specific cl::Device
void printDeviceInfo(cl::Device device)
{
    deviceInfo dinfo = getDeviceInfo(device);
    printDeviceInfo(dinfo);
}

// print info about a specific cl::Device given pre-queried info in deviceInfo struct
void printDeviceInfo(deviceInfo dinfo)
{
    spdlog::info("  Name:   {}", dinfo.name.c_str());
    spdlog::info("  Type:   {}", dinfo.devTypeStr.c_str());
    spdlog::info("  Vendor: {}", dinfo.vendor.c_str());
    spdlog::info("  Version: {}", dinfo.version.c_str());
    spdlog::info("  Compute units (CUs): {}", dinfo.computeUnits);
    spdlog::info("  Clock frequency:     {} MHz", dinfo.maxClock);
    spdlog::info("  Global memory size:  {} MB", (long long unsigned int)(dinfo.deviceMemSize / 1024 / 1024));
    spdlog::info("  Max allocation size: {} MB", (long long unsigned int)(dinfo.maxMemAllocSize / 1024 / 1024));
    spdlog::info("  Max work group/CU:   {}", (int)dinfo.maxWorkGroupSize);
    spdlog::info("  Double support:      {}", (dinfo.doubleSupport ? "true" : "false"));
    spdlog::info("  Device available:    {}", (dinfo.deviceAvailable ? "true" : "false"));
};

std::string CLErrorString(cl_int error)
{
    switch (error)
    {
    case CL_SUCCESS:
        return std::string("Success!");
    case CL_DEVICE_NOT_FOUND:
        return std::string("Device not found.");
    case CL_DEVICE_NOT_AVAILABLE:
        return std::string("Device not available");
    case CL_COMPILER_NOT_AVAILABLE:
        return std::string("Compiler not available");
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:
        return std::string("Memory object allocation failure");
    case CL_OUT_OF_RESOURCES:
        return std::string("Out of resources");
    case CL_OUT_OF_HOST_MEMORY:
        return std::string("Out of host memory");
    case CL_PROFILING_INFO_NOT_AVAILABLE:
        return std::string("Profiling information not available");
    case CL_MEM_COPY_OVERLAP:
        return std::string("Memory copy overlap");
    case CL_IMAGE_FORMAT_MISMATCH:
        return std::string("Image format mismatch");
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:
        return std::string("Image format not supported");
    case CL_BUILD_PROGRAM_FAILURE:
        return std::string("Program build failure");
    case CL_MAP_FAILURE:
        return std::string("Map failure");
    case CL_INVALID_VALUE:
        return std::string("Invalid value");
    case CL_INVALID_DEVICE_TYPE:
        return std::string("Invalid device type");
    case CL_INVALID_PLATFORM:
        return std::string("Invalid platform");
    case CL_INVALID_DEVICE:
        return std::string("Invalid device");
    case CL_INVALID_CONTEXT:
        return std::string("Invalid context");
    case CL_INVALID_QUEUE_PROPERTIES:
        return std::string("Invalid queue properties");
    case CL_INVALID_COMMAND_QUEUE:
        return std::string("Invalid command queue");
    case CL_INVALID_HOST_PTR:
        return std::string("Invalid host pointer");
    case CL_INVALID_MEM_OBJECT:
        return std::string("Invalid memory object");
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
        return std::string("Invalid image format descriptor");
    case CL_INVALID_IMAGE_SIZE:
        return std::string("Invalid image size");
    case CL_INVALID_SAMPLER:
        return std::string("Invalid sampler");
    case CL_INVALID_BINARY:
        return std::string("Invalid binary");
    case CL_INVALID_BUILD_OPTIONS:
        return std::string("Invalid build options");
    case CL_INVALID_PROGRAM:
        return std::string("Invalid program");
    case CL_INVALID_PROGRAM_EXECUTABLE:
        return std::string("Invalid program executable");
    case CL_INVALID_KERNEL_NAME:
        return std::string("Invalid kernel name");
    case CL_INVALID_KERNEL_DEFINITION:
        return std::string("Invalid kernel definition");
    case CL_INVALID_KERNEL:
        return std::string("Invalid kernel");
    case CL_INVALID_ARG_INDEX:
        return std::string("Invalid argument index");
    case CL_INVALID_ARG_VALUE:
        return std::string("Invalid argument value");
    case CL_INVALID_ARG_SIZE:
        return std::string("Invalid argument size");
    case CL_INVALID_KERNEL_ARGS:
        return std::string("Invalid kernel arguments");
    case CL_INVALID_WORK_DIMENSION:
        return std::string("Invalid work dimension");
    case CL_INVALID_WORK_GROUP_SIZE:
        return std::string("Invalid work group size");
    case CL_INVALID_WORK_ITEM_SIZE:
        return std::string("Invalid work item size");
    case CL_INVALID_GLOBAL_OFFSET:
        return std::string("Invalid global offset");
    case CL_INVALID_EVENT_WAIT_LIST:
        return std::string("Invalid event wait list");
    case CL_INVALID_EVENT:
        return std::string("Invalid event");
    case CL_INVALID_OPERATION:
        return std::string("Invalid operation");
    case CL_INVALID_GL_OBJECT:
        return std::string("Invalid OpenGL object");
    case CL_INVALID_BUFFER_SIZE:
        return std::string("Invalid buffer size");
    case CL_INVALID_MIP_LEVEL:
        return std::string("Invalid mip-map level");
    default:
        return std::string("Unknown");
    }
}

// Read source file
std::string read_file(std::string filename)
{
    std::ifstream sourceFile(filename.c_str());
    if (sourceFile.fail())
    {
        spdlog::error("Cannot find file {}", filename.c_str());
        throw cl::Error(1, "Failed to open OpenCL source file");
    }
    std::string sourceStr(std::istreambuf_iterator<char>(sourceFile), (std::istreambuf_iterator<char>()));
    sourceFile.close();
    return sourceStr;
}

// void write_file(const std::string& filename, const std::string& content="") {
// 	std::ofstream file(filename, std::ios::out);
// 	file.write(content.c_str(), content.length());
// 	file.close();
// }