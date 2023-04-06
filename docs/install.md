# Introduction

## C++

To install the C++ library, you will need the following dependencies:

* A C++ compiler (GCC, Clang, MSVC, etc.)
* Bazel (4.0 or later recommended)
* An OpenCL runtime (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)

You can build the C++ libraries using Bazel:

```
bazel build //clode/cpp:cpp
```

There are three libraries that will be built:

* libclode_features.a: The feature extraction library
* libclode_trajectory.a: The trajectory extraction library
* libopencl_resources.a: The OpenCL resources library (to find your OpenCL runtime)

## Python

To install the Python library, you will need the following dependencies:

* A C++ compiler (GCC, Clang, MSVC, etc.)
* Python 3.8 or later
* An OpenCL runtime (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)

You can install the Python library using pip:

```
    pip install clode
```

This will download Bazel to your machine (using Bazelisk)
and build the C++ libraries. It will then install the Python library.

## Windows

On Windows, prior to installing via pip you will need the following dependencies in addition to those listed above:

* The MSVC C++ compiler (e.g., Visual Studio Community installed to default path)
* MSYS2 (add msys64/usr/bin to path)

Bazel will use MSVC to build the C++ libraries.
Further, Bazel will include the OpenCL SDK in the build.
This means that you do not need to install the OpenCL SDK separately.

Should you wish to change this behaviour, you can modify the
library inside bazel/external/opencl_windows.BUILD and
bazel/repository_locations.bzl.

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details

## Verifying the installation

To verify that the installation was successful, you can run the following command:

```
python -c "import clode; print(clode.print_opencl())"
```
