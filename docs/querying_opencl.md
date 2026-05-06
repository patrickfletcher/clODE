# Querying system OpenCL capabilities

clODE can query the OpenCL capabilities of your machine. This is useful for debugging and for finding the best OpenCL device for your application.

## Print OpenCL capabilities

```python
import clode

clode.print_opencl()
```

`print_opencl()` writes a formatted device report to standard output.

## Query OpenCL capabilities as Python objects

```python
import clode

platforms = clode.query_opencl()
print(platforms)
print(platforms[0].device_info)
```

`query_opencl()` returns a list of `PlatformInfo` objects, each with a `device_info` list of `DeviceInfo` objects.

## Backend-specific ordering note

The platform ordering reported by `clode.print_opencl()` and `clode.query_opencl()` follows the active backend implementation.

- by default, the ordering comes from PyOpenCL
- with `_CLODE_BACKEND=cpp`, the ordering comes from the legacy wrapper runtime

That ordering may differ from `clinfo -l`, so choose platform and device IDs from the same toolchain you will actually use for the simulation.
