# Installation

## Python

The default package is a pure-Python distribution for Python 3.10 and newer. Install it with:

```bash
    pip install clode
```

An OpenCL runtime for your device is required. This is often included as part of your
GPU driver (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)

## Verifying the installation

To verify that the installation was successful, you can run the following command:

```py run
from clode import query_opencl
print(query_opencl())
```

Compare `clinfo -l` with `clode.print_opencl()` carefully before pinning platform and device IDs. clODE reports the PyOpenCL-visible ordering, which can differ from `clinfo`.
