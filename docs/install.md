# Installation

## Python

The default package is a pure-Python distribution for Python 3.10 and newer. Install it with:

```bash
    pip install clode
```

An OpenCL runtime for your device is required. This is often included as part of your
GPU driver (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)

Tagged releases ship the packaged OpenCL kernel assets needed by the Python-owned runtime path. PyOpenCL is the only in-tree backend.

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

If you need to distinguish an OpenCL runtime/compiler problem from a clODE kernel problem in a source checkout, run:

```bash
python tools/probe_opencl_runtime.py --platform-id 0 --device-id 0
```

The probe lists visible runtimes, builds a trivial kernel, asks PyOpenCL to match a small diagnostic struct, and then asks it to match clODE's generated `localmax` observer struct. If the first or second probe already fails, the issue is below clODE's kernel layer. The script defaults to `CLODE_TEST_PLATFORM_ID` and `CLODE_TEST_DEVICE_ID` when those environment variables are set.

## Installation from source

To install the Python library from source, you will need the following dependencies:

* Python 3.10 or later
* An OpenCL runtime (AMD APP SDK, Intel OpenCL SDK, NVIDIA CUDA, etc.)
* OpenCL development headers on Linux if `pyopencl` needs to build from source in your environment

You can then install the Python library using pip:

```bash
    pip install .
```

This uses the same pure-Python packaging path as the published wheel.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.md) file for details

## Verifying the installation

To verify that the installation was successful, you can run the following command:

```py run
from clode import query_opencl
print(query_opencl())
```

Compare `clinfo -l` with `clode.print_opencl()` carefully before pinning platform and device IDs. The visible ordering can differ between `clinfo`, the legacy runtime path, and PyOpenCL.
