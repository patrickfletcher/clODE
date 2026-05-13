# (Advanced) Runtime selection

Each simulator instance owns its own OpenCL runtime selection. That means you can construct different simulators for different devices in the same Python process.

Most users do not need to call `initialize_runtime(...)` directly. The common pattern is to pass runtime-selection arguments to `Simulator`, `FeatureSimulator`, or `TrajectorySimulator`.

## Automatic selection

If you do not pass any runtime-selection arguments, clODE chooses a default device automatically.

```python
import clode


simulator = clode.TrajectorySimulator(
    src_file="test/van_der_pol_oscillator.cl",
    variables={"x": 0.0, "y": 1.0},
    parameters={"mu": 1.0},
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, 1000.0),
)
```

## Select by device type and vendor

Use `CLDeviceType` and `CLVendor` when you want a class of device rather than a specific `(platform_id, device_id)` pair.

```python
import clode


simulator = clode.TrajectorySimulator(
    src_file="test/van_der_pol_oscillator.cl",
    variables={"x": 0.0, "y": 1.0},
    parameters={"mu": 1.0},
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, 1000.0),
    device_type=clode.CLDeviceType.DEVICE_TYPE_GPU,
    vendor=clode.CLVendor.VENDOR_NVIDIA,
)
```

Available device types:

- `CLDeviceType.DEVICE_TYPE_DEFAULT`
- `CLDeviceType.DEVICE_TYPE_CPU`
- `CLDeviceType.DEVICE_TYPE_GPU`
- `CLDeviceType.DEVICE_TYPE_ACCELERATOR`
- `CLDeviceType.DEVICE_TYPE_CUSTOM`
- `CLDeviceType.DEVICE_TYPE_ALL`

Available vendors:

- `CLVendor.VENDOR_ANY`
- `CLVendor.VENDOR_NVIDIA`
- `CLVendor.VENDOR_AMD`
- `CLVendor.VENDOR_INTEL`

## Select by platform and device index

Use `platform_id` and `device_id` when you want a specific runtime reported by `clode.query_opencl()` or `clode.print_opencl()`.

```python
import clode


simulator = clode.FeatureSimulator(
    src_file="test/van_der_pol_oscillator.cl",
    variables={"x": 0.0, "y": 1.0},
    parameters={"mu": 1.0},
    observer=clode.Observer.threshold_2,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, 1000.0),
    platform_id=0,
    device_id=0,
)
```

## Notes on `device_ids`

The public constructors still accept `device_ids`, but the high-level simulator API does not implement real multi-device work partitioning. The supported path is explicit single-device selection.

For new work, prefer one of:

- automatic selection
- `device_type` plus `vendor`
- explicit `platform_id` plus `device_id`

## Inspecting the selected runtime

Use the public query helpers when choosing a device tuple:

```python
import clode


platforms = clode.query_opencl()
print(platforms)
clode.print_opencl()
```

You can also print the devices visible to a specific simulator instance:

```python
simulator.print_devices()
```

Like `clode.print_opencl()`, this is an explicit report helper and does not depend on the current logging configuration.

Platform ordering can differ from `clinfo -l`. clODE reports the PyOpenCL-visible ordering, so choose platform and device IDs from `clode.query_opencl()` or `clode.print_opencl()` rather than assuming the `clinfo` order matches.
