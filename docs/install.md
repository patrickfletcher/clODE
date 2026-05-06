# Installation

## Python

The default package is a pure-Python distribution for Python 3.8-3.12. Install it with:

```bash
    pip install clode
```

An OpenCL runtime for your device is required. This is often included as part of your
GPU driver (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)

Tagged releases now ship the packaged OpenCL kernel assets needed by the Python-owned runtime path, and PyOpenCL is the default backend. In a source checkout, the legacy C++ backend remains available only as an explicit comparison path after its wrapper extension has been built into `clode/cpp/`.

For a concrete backend-selection example, see `examples/pyopencl_ornstein_uhlenbeck.py`.

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

* Python 3.8 or later
* An OpenCL runtime (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)
* OpenCL development headers on Linux if `pyopencl` needs to build from source in your environment

You can then install the Python library using pip:

```bash
    pip install .
```

This uses the same pure-Python packaging path as the published wheel.

### Legacy C++ compatibility backend

The legacy C++ backend is no longer part of the published wheel or the default Python build path. If you want to compare the legacy runtime against the default PyOpenCL path, use a source checkout.

1. Install the checkout into the active environment:

```bash
pip install -e .
```

1. Build the legacy wrapper with Bazel:

```bash
bazel build //clode/cpp:clode_cpp_wrapper
```

1. Copy the built shared library into `clode/cpp/` using a Python-recognized extension-module suffix for your interpreter:

```bash
python - <<'PY'
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path
import shutil

candidates = []
for root in (Path("bazel-bin"), Path("bazel-out")):
    if root.exists():
        candidates.extend(
            path
            for path in root.rglob("*")
            if path.is_file()
            and path.name
            in {
                "clode_cpp_wrapper",
                "clode_cpp_wrapper.so",
                "libclode_cpp_wrapper.so",
                "libclode_cpp_wrapper.dylib",
                "clode_cpp_wrapper.dll",
                "clode_cpp_wrapper.pyd",
            }
        )

if not candidates:
    raise SystemExit("No Bazel-built clode_cpp_wrapper shared library was found.")

suffix = next(suffix for suffix in EXTENSION_SUFFIXES if suffix.endswith((".so", ".pyd")))
destination = Path("clode/cpp") / f"clode_cpp_wrapper{suffix}"
shutil.copy2(candidates[0], destination)
print(destination)
PY
```

1. Select the legacy backend explicitly when running a comparison script or test:

```bash
_CLODE_BACKEND=cpp python -m pytest test/test_vdp.py -q
```

You can compare the two backends by running the same command twice: once with the default PyOpenCL path or `_CLODE_BACKEND=pyopencl`, and once with `_CLODE_BACKEND=cpp`.

When doing that, choose `platform_id` and `device_id` from the backend you are actually running. PyOpenCL and the legacy wrapper do not necessarily report platforms in the same order.

### Windows

On Windows, prior to building the legacy C++ path from source you will need the following dependencies in addition to those listed above:

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

Compare `clinfo -l` with `clode.print_opencl()` carefully before pinning platform and device IDs. The visible ordering can differ between `clinfo`, the legacy runtime path, and PyOpenCL.
