# XPP files

clODE's Python library can read XPP files and convert them to OpenCL source files.

To use this functionality, you must have the Python library installed.
Simply send your .xpp file to the functions
clode_features or clode_trajectory. For a fully-working example, look in `test/test_xpp_converter.py`.

## Example

You can also use the clode functions `read_ode_parameters`
and `format_opencl_rhs`.

### Converting an XPP string to OpenCL

```python
import clode

with open("van_der_pol_oscillator.xpp", "r") as f:
    xpp_str = f.read()

# Read the parameters from the XPP file
parameters, auxiliaries, initial_values, dx, noise, statements = \
    clode.read_ode_parameters(xpp_str)

clode_cl_str = clode.format_opencl_rhs(
    parameters, auxiliaries, initial_values, dx, noise, statements
)
```

### Converting an XPP file to OpenCL

```python
import clode

# The file will be named van_der_pol_oscillator.cl
# and you can find it in the same directory as the .xpp file
clode_cl_filename = clode.convert_xpp_file("van_der_pol_oscillator.xpp")
```

## Limitations

The XPP converter is not perfect. It is not guaranteed to work for all XPP files.
Formatting will not be preserved.

The following features are partially supported:

* exponents
