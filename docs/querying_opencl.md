# Querying system OpenCL capabilities

clODE can query the OpenCL capabilities of your machine. This is useful for debugging and for finding the best OpenCL device for your application.

## Example - Print OpenCL capabilities

```python
import clode

print(clode.print_opencl())
```

### Output

```
[2023-04-05 17:19:48.614] [info] 
Querying OpenCL platforms...
[2023-04-05 17:19:48.733] [info] Number of platforms found: %u
[2023-04-05 17:19:48.733] [info] 
Platform 0. ------------------------------
[2023-04-05 17:19:48.733] [info] Name:    Apple
[2023-04-05 17:19:48.733] [info] Vendor:  Apple
[2023-04-05 17:19:48.733] [info] Version: OpenCL 1.2 (Nov  4 2022 20:34:31)
[2023-04-05 17:19:48.733] [info] 
Device 0. --------------------
[2023-04-05 17:19:48.733] [info] Name:   Apple M1 Pro
[2023-04-05 17:19:48.733] [info] Type:   GPU
[2023-04-05 17:19:48.733] [info] Vendor: Apple
[2023-04-05 17:19:48.733] [info] Version: OpenCL 1.2 
[2023-04-05 17:19:48.733] [info] Compute units (CUs): 16
[2023-04-05 17:19:48.733] [info] Clock frequency:     1000 MHz
[2023-04-05 17:19:48.733] [info] Global memory size:  10922 MB
[2023-04-05 17:19:48.733] [info] Max allocation size: 2048 MB
[2023-04-05 17:19:48.733] [info] Max work group/CU:   256
[2023-04-05 17:19:48.733] [info] Double support:      false
[2023-04-05 17:19:48.733] [info] Device available:    true
[2023-04-05 17:19:48.733] [info] 
```

## Example - Query OpenCL capabilities

```python
import clode

platforms = clode.query_opencl()
print(platforms)
print(clode.query_open_cl()[0].device_info)
```

### Output

```
[<platform_info(name=Apple, vendor=Apple, version=OpenCL 1.2 (Nov  4 2022 20:34:31), device_count=1)>]

[<device_info(name=Apple M1 Pro, vendor=Apple, version=OpenCL 1.2 , device_type=GPU, compute_units=16, max_clock=1000, max_work_group_size=256, device_memory_size=11453251584, max_memory_alloc_size=2147483648, extensions=cl_APPLE_SetMemObjectDestructor cl_APPLE_ContextLoggingFunctions cl_APPLE_clut cl_APPLE_query_kernel_names cl_APPLE_gl_sharing cl_khr_gl_event cl_khr_byte_addressable_store cl_khr_global_int32_base_atomics cl_khr_global_int32_extended_atomics cl_khr_local_int32_base_atomics cl_khr_local_int32_extended_atomics cl_khr_3d_image_writes cl_khr_image2d_from_buffer cl_khr_depth_images , double_support=0, device_available=1)>]

```
