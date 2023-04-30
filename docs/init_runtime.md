# (Advanced) Initialise runtimes

There are three ways to initialise the clODE runtime.
These give users fine-grained control over which
OpenCL devices are used by clODE.

Note: Runtimes are never re_initialised. If you wish to
change the runtime, you must first reset the runtime
by calling reset_runtime.

## Automatic initialisation

The simplest way to initialise the clODE runtime is to call
CLODEFeatures and CLODETrajectory without any runtime-specific
arguments. The runtime device and vendor will be selected
automatically.

The default runtime device is cl_device_type.DEVICE_TYPE_DEFAULT.
The default runtime vendor is cl_vendor.ANY.

## Select device and vendor by name

The second way to initialise the clODE runtime is to specify
the device and platform by name. This is done by passing
the `device_type` and `vendor` arguments to initialise_runtime.

You can select from the following devices:
* cl_device_type.DEVICE_TYPE_DEFAULT
* cl_device_type.DEVICE_TYPE_CPU
* cl_device_type.DEVICE_TYPE_GPU
* cl_device_type.DEVICE_TYPE_ACCELERATOR
* cl_device_type.DEVICE_TYPE_CUSTOM
* cl_device_type.DEVICE_TYPE_ALL

You can select from the following vendors:
* cl_vendor.AMD
* cl_vendor.NVIDIA
* cl_vendor.INTEL
* cl_vendor.ANY

```python
import clode

device_type = cl_device_type,DEVICE_TYPE_GPU
vendor = cl_vendor.AMD

clode.initialise_runtime(device_type=device_type, vendor=vendor)
```

## Select platform and device by index

The third way to initialise the clODE runtime is to specify
the platform and device by index. This is done by passing
the `platformID` and `deviceID` arguments
to initialise_runtime_by_platform_id.

```python
import clode

platformID = 0
deviceID = 0

clode.initialise_runtime_by_platform_id(platformID=platformID, deviceID=deviceID)
```

## Selecting platform and multiple devices

You can also select multiple devices on a single platform.
This is done by passing the `deviceIDs` argument
to initialise_runtime_by_device_ids.

```python

import clode

platformID = 0
deviceIDs = [0, 1, 2]

clode.initialise_runtime_by_device_ids(platformID=platformID, deviceIDs=deviceIDs)
```

## Resetting the runtime

You can reset the runtime by calling reset_runtime.

```python

import clode

clode.reset_runtime()
```
