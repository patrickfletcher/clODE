# Installation

## Python

Pre-build binaries are provided via PyPI for Python 3.8-3.12 on MacOS and Windows, which can be installed simply using pip:

```bash
    pip install clode
```

An OpenCL runtime for your device is required. This is often included as part of your
GPU driver (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)

If you are validating or experimenting with the PyOpenCL backend during the current transition period, install the optional dependency as well:

```bash
    pip install clode[pyopencl]
```

The public default backend still remains the current C++ path at the moment. The optional PyOpenCL dependency is for contributor workflows, milestone validation, and the upcoming default-switch work.

For a concrete transition-phase example of selecting the PyOpenCL backend and running a simulation, see `examples/pyopencl_ornstein_uhlenbeck.py`.

### Google Colab

On Google Colab, you need to re-install the nvidia-opencl-dev package
to make the OpenCL runtime work correctly.

```jupyter
!sudo apt-get update
!sudo apt remove nvidia-opencl-dev clinfo -y
!sudo apt install nvidia-opencl-dev clinfo -y
```

### Linux

You might need to install the OpenCL runtime separately.
For example, on Ubuntu, you can install the OpenCL runtime using the following command:

```bash
sudo apt-get update
sudo apt install ocl-icd-libopencl1 ocl-icd-opencl-dev clinfo intel-opencl-icd
```

If you are using an NVIDIA device, install the vendor runtime that exposes OpenCL support through the NVIDIA driver stack instead of the Intel ICD. In either case, verify the available platforms and devices with:

```bash
clinfo -l
```

You can also inspect the OpenCL platforms visible to clODE from Python:

```python
import clode
clode.print_opencl()
```

## Installation from source

To install the Python library from source, you will need the following dependencies:

* A C++ compiler (GCC, Clang, MSVC, etc.)
* Python 3.8 or later
* An OpenCL runtime (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)
* OpenCL development headers on Linux

You can then install the Python library using pip:

```bash
    pip install clode
```

This will download Bazel to your machine (using Bazelisk)
and build the C++ libraries. It will then install the Python library.

If you also want the PyOpenCL backend dependencies in that environment, install:

```bash
    pip install -e .[pyopencl]
```

### Windows

On Windows, prior to installing via pip you will need the following dependencies in addition to those listed above:

* The MSVC C++ compiler (e.g., Visual Studio Community installed to default path)
* MSYS2 (add msys64/usr/bin to path)

Bazel will use MSVC to build the C++ libraries.
Further, Bazel will include the OpenCL SDK in the build.
This means that you do not need to install the OpenCL SDK separately.

Should you wish to change this behaviour, you can modify the
library inside bazel/external/opencl_windows.BUILD and
bazel/repository_locations.bzl.

## C++ stand-alone source installation

To install the C++ library, you will need the following dependencies:

* A C++ compiler (GCC, Clang, MSVC, etc.)
* Bazel (4.0 or later recommended)
* An OpenCL runtime (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)

You can build the C++ libraries using Bazel:

```bash
bazel build //clode/cpp:cpp
```

There are three libraries that will be built:

* libclode_features.a: The feature extraction library
* libclode_trajectory.a: The trajectory extraction library
* libopencl_resources.a: The OpenCL resources library (to find your OpenCL runtime)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details

## Verifying the installation

To verify that the installation was successful, you can run the following command:

```py run
from clode import query_opencl
print(query_opencl())
```

If you are validating the PyOpenCL backend, also verify that `clinfo -l` and `clode.print_opencl()` agree on the runtimes you intend to use before pinning platform and device IDs.
