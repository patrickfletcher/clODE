import clode

# a brief version?

# ocl_info = clode.runtime.query_opencl()
# for i, ocl in enumerate(ocl_info):
#     print(f"Platform {i}: ", ocl.name, f" [{ocl.version}]")
#     for i, dev in enumerate(ocl.device_info):
#         print(f"--Device {i}: ", dev.name, "\n")


# to use the print_opencl (and integrator.print_status) I need to set logger to info
clode.set_log_level(clode.LogLevel.info)
print(clode.runtime.print_opencl()) 