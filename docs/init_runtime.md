# (Advanced) Initialize runtimes

Each CLODEFeature and CLODETrajectory object
maintains its own OpenCL runtime. This is so
that users can initialize multiple runtimes
across different devices.

## Automatic initialization

The simplest way to initialize the clODE runtime is to call
CLODEFeatures and CLODETrajectory without any runtime-specific
arguments. The runtime device and vendor will be selected
automatically.

The default runtime device is cl_device_type.DEVICE_TYPE_DEFAULT.
The default runtime vendor is cl_vendor.ANY.

## Select device and vendor by name

The second way to initialize the clODE runtime is to specify
the device and platform by name. This is done by passing
the `device_type` and `vendor` arguments to CLODEFeatures or CLODETrajectory.

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

input_file: str = "test/van_der_pol_oscillator.cl"
tspan = (0.0, 1000.0)

trajectory = clode.CLODETrajectory(
    src_file=input_file,
    variable_names=["x", "y"],
    parameter_names=["mu"],
    num_noise=0,
    stepper=clode.Stepper.dormand_prince,
    tspan=tspan,
    device_type=device_type,
    vendor=vendor,
)
```

## Select platform and device by index

The third way to initialize the clODE runtime is to specify
the platform and device by index. This is done by passing
the `platformID` and `deviceID` arguments
to CLODEFeatures or CLODETrajectory.

```python
import clode

platform_id = 0
device_id = 0

input_file: str = "test/van_der_pol_oscillator.cl"
tspan = (0.0, 1000.0)

trajectory = clode.CLODETrajectory(
    src_file=input_file,
    variable_names=["x", "y"],
    parameter_names=["mu"],
    num_noise=0,
    stepper=clode.Stepper.dormand_prince,
    tspan=tspan,
    platform_id=platform_id,
    device_id=device_id,
)
```

## Selecting platform and multiple devices

You can also select multiple devices on a single platform.
This is done by passing the `deviceIDs` argument
to CLODEFeatures or CLODETrajectory.

```python

import clode

platformID = 0
deviceIDs = [0, 1, 2]

input_file: str = "test/van_der_pol_oscillator.cl"
tspan = (0.0, 1000.0)

trajectory = clode.CLODETrajectory(
    src_file=input_file,
    variable_names=["x", "y"],
    parameter_names=["mu"],
    num_noise=0,
    stepper=clode.Stepper.dormand_prince,
    tspan=tspan,
    platform_id=platform_id,
    device_ids=device_ids,
)
```

## Printing the devices

You can print the devices that are used
by calling the print_devices method of
CLODEFeature or CLODETrajectory objects.

```python
import clode

features = CLODEFeatures(...)
trajectory = CLODETrajectory(...)

trajectory.print_devices()

features.print_devices()
```
