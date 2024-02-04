# Specifying ODEs

clODE supports many common use cases for specifying systems of ordinary
differential equations (ODEs). This document outlines the various ways to
specify ODEs in clODE.

# Calling helper functions

## Specifying helper functions

It is possible to specify helper functions in Python and call them from the
main ODE function. This can be useful for breaking up complex ODEs into
smaller, more manageable pieces.

```py run
from clode import OpenCLConverter

def helper_function(x: float, y: float) -> float:
    return (x + y) / y

def get_rhs(t: float,
            variables: list[float],
            parameters: list[float],
            derivatives: list[float],
            aux: list[float],
            wiener: list[float]) -> None:
    x: float = variables[0]
    y: float = variables[1]
    
    derivatives[0] = helper_function(x, y)

converter = OpenCLConverter()
converter.convert_to_opencl(helper_function)
print(converter.convert_to_opencl(get_rhs))
```

To load a helper function into the OpenCL program, you must use
the list argument `supplementary_equations` in the Simulator initializers.

Note: The functions must be specified in the order that they are used.
clODE will not automatically resolve dependencies between functions.

```python
simulator = TrajectorySimulator(
    rhs_equation=get_rhs,
    variables={"x": 0.0, "y": 2.0},
    parameters={"a": 1.0},
    supplementary_equations=[helper_function],
    t_span=(0, 10),
    dt=0.01,
)
```

## Using builtins (min, max, etc.)

Certain builtin functions are available in clODE.
In a future version, clODE will automatically determine
which function to use based on the argument types.
Therefore, it is recommended to use the clODE builtins
instead of the Python builtins when possible.

The functions currently available are:
- `OpenCLExp`
- `OpenCLMin`
- `OpenCLMax`

```py run
from clode import OpenCLConverter, OpenCLExp

def get_rhs(
    t: float,
    variables: list[float],
    parameters: list[float],
    derivatives: list[float],
    aux: list[float],
    wiener: list[float],
) -> None:
    x: float = variables[0]
    y: float = variables[1]

    derivatives[0] = (x - y) / OpenCLExp(y)
    derivatives[1] = -x
    

converter = OpenCLConverter()
print(converter.convert_to_opencl(get_rhs))
```